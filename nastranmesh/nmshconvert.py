#!/usr/bin python
import sys
import re
from cell import cell, Face
from sets import Set
from copy import copy
from numpy import zeros, array, sign, cross, dot, ones, arctan, sin, cos, pi, mod, sqrt, inner
from numpy.linalg import norm, solve
from pylab import find
from operator import add
from scipy.optimize import fsolve
from dolfin import File, MeshEditor, Mesh, Cell, facets
from collections import defaultdict
import matplotlib.pyplot as plt

print "Converting from Nastran format to Nek5000, semtex or FEniCS format"
tot_num_curved = 0
# Use regular expressions to identify sections and tokens found in a fluent file
re_dimline  = re.compile(r"\(2\s(\d)\)")
#re_comment  = re.compile(r"\(0\s.*") # Fluent
re_comment  = re.compile(r"^\$\s.*") # Nastran
re_zone0    = re.compile(r"\(10\s\(0\s(\w+)\s(\w+)\s(\d+)\s(\d+)\)\)")
re_zone     = re.compile(r"\(10\s\((\w+)\s(\w+)\s(\w+)\s(\d+)\s(\d)\)(\(|)")
re_face0    = re.compile(r"\(13(\s*)\(0\s+(\w+)\s+(\w+)\s+(0|0 0)\)\)")
re_face     = re.compile(r"\(13(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\)(\s*)(\(|)")
re_periodic = re.compile(r"\(18.*\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\).*\(")
re_pfaces   = re.compile(r"((^\s)|)(\w+)(\s*)(\w+)")
re_cells0   = re.compile(r"\(12(\s*)\(0(\s+)(\w+)(\s+)(\w+)(\s+)(0|0 0)\)\)")
re_cells    = re.compile(r"\(12.*\((\w+)\s+(\w+)\s+(\w+)\s+(\d+)\s+(\d+)\)\)")
re_cells2   = re.compile(r"\(12(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\)(\s*)(\(|)")
re_zones    = re.compile(r"\((45|39)\s+\((\d+)\s+(\S+)\s+(\S+).*\)\((.*|[0-9]+[\.]*[0-9]*)\)\)")
re_parant   = re.compile(r"(^\s*\)(\s*)|^\s*\)\)(\s*)|^\s*\(\s*)")

# Declare some maps that will be built when reading in the lists of nodes and faces:
cell_map = {}               # Maps cell id with nodes
cell_face_map = {}          # Maps cell id with faces
face_cell_map = {}          # Maps face id with two cells
face_list = []              # List of faces [[id, 2-4 nodes, 2 connecting cells and type]]
face_map = {}               # For each cell a dictionary with key=face and val=local face number
nodes = None                # Will be numpy array of nodes
Nodes = []                  # List of nodes (temporary)

# Information about connectivity and boundaries
boundary_map = {}           # For each cell and local face number, get connected cell and local face number
boundary_val_map = {}       # Holds boundary conditions, like V, v, P, etc
boundary_nodes = {}         # List of nodes attached to a boundary. Key is zone id
boundary_nodes_face_map = {}# Map of the faces that belongs to a node on a boundary
boundary_faces = {}         # List of faces attached to a boundary. Key is zone id

# Information about mesh periodicity
periodic_face_map = {}      # Map of face to face periodicity
periodic_cell_face_map = {} # Map (cell, local face number) of periodic face to (cell, local face number) of shadow
periodic_node_map = {}      # Periodic node to node map

# Some global values
num_cells = {}              # Total number of cells in different zones
num_vertices = {}           # Total number of grid points
zones = {}                  # zone information
zone_number_of_faces = {}   # number of faces for each zone

# For Nek5000
temperature_val_map = {}    # Boundary conditions for temperature 
passive_scalar_val_map = [] # Boundary conditions for passive scalars
face_node_map = {}          # Maps face with nodes
mid_point_map = {}          # Store location of midpoint for key face

# For Nek5000 and semtex that can have curved boundaries
curves_map = {}             # Holds curve information
curved_faces = {}           # For each boundary zone store faces that are in curved section
curved_nodes = {}           # For each boundary zone store nodes that are in curved section 
curved_midpoint = {}        # For each cell and curved edge, coordinates for midpoint (empty for straight edge)
# ADDED BY RUD 25.09
curved_east = {}            # For each cell and curved edge, coordinates for the point east of midpoint  
curved_west = {}            # For each cell and curved edge, coordinates for the point west of midpoint  
tot_num_nodes = 0           # Total number of nodes
# RUD PROJ ALG
surf_list = [] # A list containing all elements and the local face that are to be projected on an exact surface from ICEM
# END BY RUD 25.09


######## ADDED BY RUD 25.09.15 #########
# Global dictionaries
#edge2points = {1:[1,2], 2:[2,3], 3:[3,4], 4:[1,4], 5:[5,6], 6:[6,7], 7:[7,8], 8:[8,9], 9:[1,5], 10:[2,6], 11:[3,7], 12:[4,8]}
#face2points = {1:[1,2,6,5], 2:[2,3,7,6], 3:[4,3,7,8], 4:[1,4,8,5], 5:[1,2,3,4], 6:[5,6,7,8]}
#edge2faces = {1:[1,5], 2:[2,5], 3:[3,5], 4:[4,5], 5:[1,6], 6:[2,6], 7:[3,6], 8:[4,6], 9:[4,1], 10:[1,2], 11:[2,3], 12:[3,4]}
#shadow_face = {1:3, 3:1, 2:4, 4:2, 5:6, 6:5}
def find_neighbours():
    global Faces,Cells
    for i in range(len(Faces)):
        #print Faces[i].cell_1,Faces[i].cell_2
        c1 = Faces[i].cell_1 
        c2 = Faces[i].cell_2
        # Add =he neighbours by using the faces
        if c1 != 0 : 
            f1 = Cells[c1-1].face_glloc[i+1] #local face value for cell 1
            Cells[c1-1].nbours[f1] = c2 # zero value indicates no neighbour
        if c2 != 0:
            f2 = Cells[c2-1].face_glloc[i+1] #local face value for cell 2
            Cells[c2-1].nbours[f2] = c1 # zero value indicates no neighbour
    #for i in range(len(Cells)):
        #print Cells[i].nbours

def getreastart(name):
	# A function that read the start of a .rea file and returns it as a string.
	# reads the first 127 lines 
	f = open(name,'r')
	s = ''
	for i in range(127):
		s = s+f.readline()
	return s

def getreaend(name):
	# A function that read the end of a .rea file and returns it as a string.
	# reads the first 127 lines 
	s = ''
	dummy = -1
	for line in open(name).readlines():
#		s = s+full[i]
		test = line.find('PRESOLVE')
		if(test!=-1): dummy = 1
		if(dummy==1): s = s+line
	return s


# This function is made to fix the boundary conditions SYM and v 
def fixbc(name):
    s = ''
    dummy = -1
    for line in open(name).readlines():
           if(dummy==1): #If we are in the BC-section 
               if(line.find('****')!=-1): dummy = -1 # Out of BC'section
               if(dummy==1): 
                   #line = line.replace('M  ','SYM')
                   line = line.replace('V  ','v  ')
           if(line.find('FLUID')!=-1): dummy = 1 # When reaching the BCS
           s = s+line
    newfile = open(name, 'w')
    newfile.write(s)
    return 1

# This function is made to fix the thermal boundary conditions
def fixthermalbc(name,iftemp):
    s = '' # full txt
    thermal = '' # Thermal conditions
    dummy = -1
    for line in open(name).readlines():
           if(dummy==1): #If we are in the FLUID BC-section 
               if(line.find('****')!=-1): 
                   dummy = -1 # Out of BC'section
                   line = '***** THERMAL BOUNDARY CONDITIONS *****\n'
               if(dummy==1): #If we still are in the FLUID BC-section 
                   thermal = thermal+line.replace('v  ','t  ')
           if(line.find('FLUID')!=-1): dummy = 1 # When reaching the BCS
           s = s+line
           if(dummy==-1 and thermal!=''): 
               s = s+thermal # Adding the modified conditions to the termal part 
               thermal = ''
    newfile = open(name, 'w')
    print 'temperature variable = {}'.format(iftemp)
    #if(iftemp=='False'): 
        #print 'temperature variable = {}'.format(iftemp)
        #s = '  ***** NO THERMAL BOUNDARY CONDITIONS ***** \n'
    newfile.write(s) 
    return 1

def points2circ(x1,x2,x3):
    # takes in three 3D points and returns 
    # The radius (r) and center x of the sphere with its center 
    # in the same plane as the three points.
    # INPUT: xi = array([real,real,real])
    # OUTPUT: r = real, x = array([real,real,real])
    n = cross(x2-x1,x3-x1) # 
    phi1 = x2-x1
    psi1 = -sum(x1**2-x2**2)/2
    phi2 = x3-x2
    psi2 = -sum(x2**2-x3**2)/2
    psi3 = dot(n,x1)
    A = array([phi1,phi2,n])
    b = array([psi1,psi2,psi3])
    x = solve(A,b)
    r = sqrt(dot(x-x1,x-x1))
    return r,x

def write_surface_file():

    print 'Writing surface nodes to surf.i !'
    global tot_num_nodes
    ofile  = open('surf.i', "w")
    #ofile.write('{}\n'.format(tot_num_nodes))
    for n in range(tot_num_nodes):
        nn = nodes[:,n] 
        ofile.write(reduce(add, ['{0:.6e} '.format(x) 
        for x in nn]) + '\n')
    ofile.close()

    ### Do necessary changes in SIZE file! 
    ifile = open('SIZE',"r")
    lines = ifile.readlines()
    ifile.close()
    if len(lines) == 0:
        raise IOError("Empty SIZE file")
    ofile = open("SIZE","w")
    # Iterating through the lines
    dummy = 0
    target = '(?<=nsurf\=).+\)'
    repl= '{}'.format(tot_num_nodes)+')'
    for line in lines:
        # want to replace lfdm=XXX with lfdm = 38
        m = re.search(target,line)
        if(m): dummy = 1
        line = re.sub(target,repl,line)
        #re.sub('(?<=lfdm).+)','=38\)',line)
        #print line
        ofile.write(line)
    if (not dummy): ofile.write('      parameter (nsurf={}) ! Number of bdry nodes\n'.format(tot_num_nodes))
    ofile.close()

def write_surf_list():

    print 'Writing surface elements and local faces to bdry.i !'
    ofile = open('bdry.i', "w")
    #ofile.write('{}\n'.format(len(surf_list)))
    for info in surf_list:
        #ofile.write(reduce(add, ['{0:.8e}'.format(x).rjust(16) 
        ofile.write(reduce(add, ['{} '.format(x) 
                                    for x in info]) + '\n')
    ofile.close()

    print 'Appending some essential variables in SIZE as well!! Remeber to call the script from the same folder as your SIZE-file!!'
    ifile = open('SIZE',"r")

    ### Do necessary changes in SIZE file! 
    lines = ifile.readlines()
    ifile.close()
    if len(lines) == 0:
        raise IOError("Empty SIZE file")
    ofile = open("SIZE","w")
    # Iterating through the lines
    dummy = 0
    target = '(?<=nbdry\=).+\)'
    repl= '{}'.format(len(surf_list))+')'
    for line in lines:
        # want to replace lfdm=XXX with lfdm = 38
        m = re.search(target,line)
        if(m): dummy = 1
        line = re.sub(target,repl,line)
        #re.sub('(?<=lfdm).+)','=38\)',line)
        #print line
        ofile.write(line)
    if (not dummy): ofile.write('      parameter (nbdry={}) ! Number of bdry nodes\n'.format(len(surf_list)))
    ofile.close()

######## END RUD 25.09.15 #########

Faces = []
Cells = []

# This part is just a default that goes into the top of the rea-file:
rea_start= """****** PARAMETERS *****
    2.610000     NEKTON VERSION
   {0:2d} DIMENSIONAL RUN
          103 PARAMETERS FOLLOW
   1.00000         p1  DENSITY
  -100.000         p2  VISCOS
   0. 
   0. 
   0. 
   0. 
   1.00000         p7  RHOCP
   1.00000         p8  CONDUCT
   0. 
   0.              p10 FINTIME
   10000           p11 NSTEPS
  -0.10000E-01     p12 DT
   0.              p13 IOCOMM
   0.              p14 IOTIME
   5               p15 IOSTEP
   0.              p16 PSSOLVER
   0. 
  -20.00000        p18 GRID
  -1.00000         p19 INTYPE
   10.0000         p20 NORDER
   0.100000E-05    p21 DIVERGENCE
   0.100000E-06    p22 HELMHOLTZ
   0.              p23 NPSCAL
   0.100000E-01    p24 TOLREL
   0.100000E-01    p25 TOLABS
   1.00000         p26 COURANT/NTAU 
   2.00000         p27 TORDER
   0.              p28 TORDER: mesh velocity (0: p28=p27)
   0.              p29 magnetic visc if > 0, = -1/Rm if < 0
   0.              p30 > 0 ==> properties set in uservp()
   0.              p31 NPERT: #perturbation modes
   0.              p32 #BCs in re2 file, if > 0
   0. 
   0. 
   0. 
   0. 
   0. 
   0. 
   0. 
   0. 
   0.              p41 1-->multiplicative SEMG
   0.              p42 0=gmres/1=pcg
   0.              p43 0=semg/1=schwarz
   0.              p44 0=E-based/1=A-based prec.
   0.              p45 Relaxation factor for DTFS
   0.              p46 reserved
   0.              p47 vnu: mesh matieral prop
   0. 
   0. 
   0. 
   0. 
   0.              p52 IOHIS
   0. 
   0.              p54 1,2,3-->fixed flow rate dir=x,y,z
   0.              p55 vol.flow rate (p54>0) or Ubar (p54<0)
   0. 
   0. 
   0. 
   0.              p59 !=0 --> full Jac. eval. for each el.
   0.              p60 !=0 --> init. velocity to small nonzero
   0. 
   0.              p62 >0 --> force byte_swap for output
   0.              p63 =8 --> force 8-byte output
   0.              p64 =1 --> perturbation restart
   1.00000         p65 #iofiles (eg, 0 or 64); <0 --> sep. dirs
   4.00000         p66 output : <0=ascii, else binary
   4.00000         p67 restart: <0=ascii, else binary
   0.              p68 iastep: freq for avg_all
   0. 
   0. 
   0. 
   0. 
   0. 
   0.              p74 verbose Helmholtz
   0. 
   0. 
   0. 
   0. 
   0. 
   0. 
   0. 
   0. 
   0. 
   0.               p84 !=0 --> sets initial timestep if p12>0
   0.               p85 dt ratio if p84 !=0, for timesteps>0
   0.               p86 reserved
   0. 
   0. 
   0. 
   0. 
   0. 
   0. 
   20.0000          p93 Number of previous pressure solns saved
   3.00000          p94 start projecting velocity after p94 step
   5.00000          p95 start projecting pressure after p95 step
   0. 
   0. 
   0. 
   0.               p99   dealiasing: <0--> off/3--> old/4--> new
   0.              
   0.               p101   No. of additional filter modes
   1.00000          p102   Dump out divergence at each time step
  -1.00000          p103   weight of stabilizing filter (.01) 
      4  Lines of passive scalar data follows2 CONDUCT; 2RHOCP
   1.00000       1.00000       1.00000       1.00000       1.00000
   1.00000       1.00000       1.00000       1.00000
   1.00000       1.00000       1.00000       1.00000       1.00000
   1.00000       1.00000       1.00000       1.00000
           13  LOGICAL SWITCHES FOLLOW
  T     IFFLOW
  F     IFHEAT
  T     IFTRAN
  T F F F F F F F F F F IFNAV & IFADVC (convection in P.S. fields)
  F F T T T T T T T T T T IFTMSH (IF mesh for this field is T mesh)
  F     IFAXIS
  F     IFSTRS
  F     IFSPLIT
  F     IFMGRID
  F     IFMODEL
  F     IFKEPS
  F     IFMVBD
  F     IFCHAR
   5.00000       5.00000      -2.75000      -2.75000     XFAC,YFAC,XZERO,YZERO
"""

# This part goes into the bottom of the rea-file:
rea_end = """            0 PRESOLVE/RESTART OPTIONS  *****
            7         INITIAL CONDITIONS *****
C Default
C Default
C Default
C Default
C Default
C Default
C Default
  ***** DRIVE FORCE DATA ***** BODY FORCE, FLOW, Q
            4                 Lines of Drive force data follow
C
C
C
C
  ***** Variable Property Data ***** Overrrides Parameter data.
            1 Lines follow.
            0 PACKETS OF DATA FOLLOW
  ***** HISTORY AND INTEGRAL DATA *****
            0   POINTS.  Hcode, I,J,H,IEL
  ***** OUTPUT FIELD SPECIFICATION *****
            6 SPECIFICATIONS FOLLOW
  F      COORDINATES
  T      VELOCITY
  T      PRESSURE
  F      TEMPERATURE
  F      TEMPERATURE GRADIENT
            0      PASSIVE SCALARS
  ***** OBJECT SPECIFICATION *****
       0 Surface Objects
       0 Volume  Objects
       0 Edge    Objects
       0 Point   Objects
"""

def read_periodic(lines, periodic_dx):
    """Scan periodic section and create periodic_face_map.
    However, if a dictionary periodic_dx has been supplied then simply skip
    past this section without creating the periodic_face_map."""
    while 1:
        line = lines.pop(0)
        a = re.search(re_pfaces, line)
        if a:
            if not periodic_dx:
                periodic_face_map[int(a.group(3), 16)] = int(a.group(5), 16)
            continue
        break

    if not periodic_dx:
        keys = periodic_face_map.keys()
        vals = periodic_face_map.itervalues()
        for key, val in zip(keys, vals):
            periodic_face_map[val] = key
        
def create_periodic_face_map(periodic_dx):
    """Create a map between periodic zones.
    periodic_dx is a dictionary with key a tuple of the connected zones  
    and value (dx, dy, dz) the actual separation between the zones.
    If periodic_dx is {} the loop is never entered.
    """
    for zone0, zone1 in periodic_dx:
        print 'Generating map for periodic zones ', zone0, zone1, '\n'
        face0_list = []
        face1_list = []
        nodes0 = []
        nodes1 = []
        # Get the faces of mapped zones
        for i, face in enumerate(face_list):
            if face[-1] == zone0:
                face0_list.append((face, i))
                nodes0 += face[1]
                face[-2] = 8 # bc_type is now periodic
            elif face[-1] == zone1:
                face1_list.append((face, i))
                nodes1 += face[1]
                face[-2] = 12 # bc_type is now shadow
        # Get unique lists of nodes
        nodes0 = list(Set(nodes0))
        nodes1 = list(Set(nodes1))
        periodic_node_map[(zone0, zone1)] = {}
        # Get mapping from zone0 to zone1
        dx = array(periodic_dx[(zone0, zone1)])            
        # Go through all nodes in zone0 and find the periodic match
        for node in nodes0:
            original_node = nodes[:, node - 1]
            for shadow_node in nodes1:
                if all(abs(abs(nodes[:, shadow_node - 1] - original_node) 
                           - abs(dx)) < 1.e-7*norm(dx)):
                    periodic_node_map[(zone0, zone1)][node] = shadow_node
                    nodes1.remove(shadow_node)
                    
        # Generate periodic face map
        for face, face_number in face0_list:
            nodes0 = face[1]
            true_nodes_of_shadow = [periodic_node_map[(zone0, zone1)][i] 
                                                        for i in nodes0]            
            for face_of_shadow, shadow_number in face1_list:
                nodes_of_shadow = face_of_shadow[1]
                if len(Set(nodes_of_shadow + true_nodes_of_shadow)) == 4:
                    periodic_face_map[face_number + 1] = shadow_number + 1
                    break
                    
            face1_list.remove((face_of_shadow, shadow_number))
    
def process_cell(cell_no, vertices,tol):
    #print 'processing cell number: {}'.format(cell_no)
    # Update Face dictionary from the nodes in a cell
    # Store cell_map, listing the vertices that describe the cell  
    cell_map[cell_no] = zeros(21, int)
    cell_map[cell_no] = vertices[1:21]
    # A Nek element is numbered different than a Nastran element 
    #cell_map[cell_no][0:8] = vertices[[2, 3, 4, 1, 6, 7, 8, 5]]
    #if (vertices[0] == 20):
    #    cell_map[cell_no][8:20] = vertices[[10, 11, 12, 9, 14, 15, 16, 13, 18, 19, 20, 17]]

    # ADDED BY MAGNUS 06.10 ####### 
    currentCell = cell(cell_no,vertices)
    Cells.append(currentCell)
    #print currentCell.cellID
    #print currentCell.vertices
    # End 

    # Checking edge midpoints for curved edges
    curved_midpoint[cell_no] = {}
    curved_east[cell_no] = {}
    curved_west[cell_no] = {}
    if (vertices[0] == 20):
        check_edge(1, [vertices[1], vertices[2], vertices[9]], cell_no,tol)
        check_edge(2, [vertices[2], vertices[3], vertices[10]], cell_no,tol)
        check_edge(3, [vertices[3], vertices[4], vertices[11]], cell_no,tol)
        check_edge(4, [vertices[4], vertices[1], vertices[12]], cell_no,tol)
        check_edge(5, [vertices[5], vertices[6], vertices[17]], cell_no,tol)
        check_edge(6, [vertices[6], vertices[7], vertices[18]], cell_no,tol)
        check_edge(7, [vertices[7], vertices[8], vertices[19]], cell_no,tol)
        check_edge(8, [vertices[8], vertices[5], vertices[20]], cell_no,tol)
        check_edge(9, [vertices[1], vertices[5], vertices[13]], cell_no,tol)
        check_edge(10, [vertices[2], vertices[6], vertices[14]], cell_no,tol)
        check_edge(11, [vertices[3], vertices[7], vertices[15]], cell_no,tol)
        check_edge(12, [vertices[4], vertices[8], vertices[16]], cell_no,tol)
            

    faceverts = zeros(4)
    face_map[cell_no] = {}

    # Define faces according to Nastran node numbering

    # South face
    faceverts = ([vertices[1], vertices[2], vertices[3], vertices[4]])
    # Check for existing face
    facename = "South"
    face_no = check_face(faceverts, facename, cell_no, len(Faces))
    face_map[cell_no][face_no] = 5  
    (Cells[cell_no-1]).faces[5] = face_no # Added by Magnus 06.10

    # North Face
    faceverts = ([vertices[5], vertices[6], vertices[7], vertices[8]])
    # Check for existing face
    facename = "North"
    face_no = check_face(faceverts, facename, cell_no, len(Faces))
    face_map[cell_no][face_no] = 6
    (Cells[cell_no-1]).faces[6] = face_no # Added by Magnus 06.10

    # West face
    #faceverts = ([vertices[1], vertices[2], vertices[6], vertices[5]])
    faceverts = ([vertices[4], vertices[1], vertices[5], vertices[8]])
    # Check for existing face
    facename = "West"
    face_no = check_face(faceverts, facename, cell_no, len(Faces))
    face_map[cell_no][face_no] = 4
    (Cells[cell_no-1]).faces[4] = face_no # Added by Magnus 06.10

    # East Face
    #faceverts = ([vertices[4], vertices[3], vertices[7], vertices[8]])
    faceverts = ([vertices[2], vertices[6], vertices[7], vertices[3]])
    # Check for existing face
    facename = "East"
    face_no = check_face(faceverts, facename, cell_no, len(Faces))
    face_map[cell_no][face_no] = 2
    (Cells[cell_no-1]).faces[2] = face_no # Added by Magnus 06.10

    # Top face
    #faceverts = ([vertices[1], vertices[4], vertices[8], vertices[5]])
    faceverts = ([vertices[3], vertices[7], vertices[8], vertices[4]])
    # Check for existing face
    facename = "Top"
    face_no = check_face(faceverts, facename, cell_no, len(Faces))
    face_map[cell_no][face_no] = 3
    (Cells[cell_no-1]).faces[3] = face_no # Added by Magnus 06.10

    # Bottom Face
    #faceverts = ([vertices[2], vertices[3], vertices[7], vertices[6]])
    faceverts = ([vertices[1], vertices[2], vertices[6], vertices[5]])
    # Check for existing face
    facename = "Bottom"
    face_no = check_face(faceverts, facename, cell_no, len(Faces))
    face_map[cell_no][face_no] = 1
    (Cells[cell_no-1]).faces[1] = face_no # Added by Magnus 06.10


def check_edge(edge_no, edgeverts, cell_no,tol):
    global tot_num_curved
    
    [xa, ya, za] = [nodes[0, edgeverts[0]-1], nodes[1, edgeverts[0]-1], nodes[2, edgeverts[0]-1]]
    [xb, yb, zb] = [nodes[0, edgeverts[1]-1], nodes[1, edgeverts[1]-1], nodes[2, edgeverts[1]-1]]
    [xc, yc, zc] = [nodes[0, edgeverts[2]-1], nodes[1, edgeverts[2]-1], nodes[2, edgeverts[2]-1]]

    # Test for midpoint
    AC = array([xa-xc, ya-yc, za-zc])
    AB = array([xa-xb, ya-yb, za-zb])
    BC = array([xb-xc, yb-yc, zb-zc])
    AC_abs = norm(AC)
    AB_abs = norm(AB)
    BC_abs = norm(BC)
    dev = abs(AC_abs-BC_abs) / AB_abs
    #if (dev > 1e-03):
        #print "Edge point is not midpoint for edge", edge_no, "cell", cell_no, "(deviation", dev,")"
        #print "A: ", [xa, ya, za]
        #print "B: ", [xb, yb, zb]
        #print "C: ", [xc, yc, zc]
        #print "AC, BC, AB: ", AC_abs, BC_abs, AB_abs
        
    # First test for curvature: Sum of vector lengths 
    if ((AC_abs + BC_abs - AB_abs) > tol * AB_abs):
        #print "Curved line for edge", edge_no, "cell", cell_no
        curved_midpoint[cell_no][edge_no] = edgeverts[2] - 1
        curved_east[cell_no][edge_no] = edgeverts[0] - 1
        curved_west[cell_no][edge_no] = edgeverts[1] - 1
        tot_num_curved += 1
        Cells[cell_no-1].addCurve(edge_no,edgeverts[2]-1,edgeverts[0]-1,edgeverts[1]-1) # ADDED BY RUD 06.10

    # Second test for curvature: Inner product
    #if (abs(inner(AC,BC)) / (norm(AC) * norm(BC)) < 1 - 1e-04):
    #    print "Test 2: Curved line for edge", edge_no, "cell", cell_no

    # Third test for curvature: Vector product
    #if (norm(cross(AC,BC)) > 1e-04):
    #    print "Test 3: Curved line for edge", edge_no, "cell", cell_no

def check_face(faceverts, facename, cell_no, face_no):
    global Faces
    oldFaces = Face.find_by_nodesum(int(sum(faceverts)))
#    print facename, " Length of oldfaces: ", len(oldFaces)
    found = False
#    print facename, " Length of Faces: ", len(Faces)
    if (len(oldFaces) > 0):
        for f in oldFaces:
#            print "Checking face no ", f.getID()
            if f.isEqual(faceverts[0], faceverts[1], faceverts[2], faceverts[3]):
                f.addCell(cell_no)
                found = True
                global_faceno = f.getID() 
    if (not found):
        currentFace = Face()
        currentFace.addData(faceverts[0], faceverts[1], faceverts[2], faceverts[3], cell_no) 
        face_no += 1
        currentFace.addID(face_no)
        Faces.append(currentFace)
        global_faceno = face_no
 #       print facename, " face appended as face no ", face_no

    # Add cell_face_map entry (storing the faces of a given cell)
    if not cell_no in cell_face_map:
        cell_face_map[cell_no] = [global_faceno] 
    else:
        cell_face_map[cell_no].append(global_faceno)

    # Add face_map entry (store ordered faces of a given cell)
#    face_map[cell_no]
        
    return global_faceno

def set_boundary_zone(zone_count, faceverts):
    oldFaces = Face.find_by_nodesum(int(sum(faceverts)))
    found = 0
    if (len(oldFaces) > 0 ):
        for f in oldFaces:
            if f.isEqual(faceverts[0], faceverts[1], faceverts[2], faceverts[3]):
                found += 1
                f.addZone(zone_count)
    else:
        print "Warning: No face found for vertices ", faceverts
    if (found > 1):    
        print "Warning: ", found, " faces found for vertices ", faceverts

def create_face_list():
    for face in Faces:
        zone_id = face.zone_id
        nd = 4
#        nds = [int(x, 16) for x in [face.vert_1, face.vert_2, face.vert_3, face.vert_4]]
        nds = [int(x) for x in [face.vert_1, face.vert_2, face.vert_3, face.vert_4]]
#        cells = [int(x, 16) for x in [face.cell_1, face.cell_2]]
        cells = [int(x) for x in [face.cell_1, face.cell_2]]
        if (min(cells) == 0):  # A boundary zone
            face_number = face.getID()
            if zone_id in boundary_nodes:
                boundary_nodes[zone_id] += nds
                boundary_faces[zone_id] += [face_number - 1]
                for nd in nds:
                    if nd in boundary_nodes_face_map[zone_id]:
                        boundary_nodes_face_map[zone_id][nd] += [face_number - 1]
                    else:
                        boundary_nodes_face_map[zone_id][nd] = [face_number - 1]
            else:
                boundary_nodes[zone_id] = nds
                boundary_faces[zone_id] = [face_number - 1]
                boundary_nodes_face_map[zone_id] = { nd: [face_number - 1]}
        bc_type = 2      # Internal boundary
#        zone_id = 0      # Default interior zone
        face_list.append([nd, copy(nds), copy(cells), bc_type, zone_id])

def create_cell_face_map(dim, mesh_format):
    """Creates maps from cells to faces and gives faces local numbering."""
    if mesh_format == 'fenics':
        create_cell_map(dim)
    else:
        eval('create_cell_face_map_' + str(dim) + 'D()') 

def create_cell_face_map_2D():
    
    for cell, faces in cell_face_map.iteritems():
        
        for face in faces:
            nds = face_list[face - 1][1]
                
            if not cell in cell_map:
                cell_map[cell] = copy(nds)
                
            elif len(cell_map[cell]) == 4:
                # Finished
                pass
            
            elif len(Set(cell_map[cell] + nds)) - len(cell_map[cell] +
                                                    nds) == 0:
                # No common node, no connectivity
                pass
            
            else:
                # Add to list preserving rotational connectivity
                # but not neccessarily right-handed
                if nds[0] in cell_map[cell]:
                    if nds[0] == cell_map[cell][0]:
                        cell_map[cell].insert(0, nds[1])
                    else:
                        cell_map[cell].append(nds[1])
                else:
                    if nds[1] == cell_map[cell][0]:
                        cell_map[cell].insert(0, nds[0])
                    else:
                        cell_map[cell].append(nds[0])

    # Ensure right-handed mesh
    for c, nds in cell_map.iteritems():
        cross_product = crss2d(nds)
        if any([i <= 0 for i in cross_product]):
            cell_map[c].reverse() 
            cross_product = crss2d(nds)
    
    # Generate face_map
    for cell, nds in cell_map.iteritems():
        face_map[cell]={(nds[0], nds[1]) : 1, (nds[1], nds[0]) : 1,
                        (nds[1], nds[2]) : 2, (nds[2], nds[1]) : 2,
                        (nds[2], nds[3]) : 3, (nds[3], nds[2]) : 3,
                        (nds[3], nds[0]) : 4, (nds[0], nds[3]) : 4,
                        1 : (nds[0], nds[1]), 1 : (nds[1], nds[0]),
                        2 : (nds[1], nds[2]), 2 : (nds[2], nds[1]),
                        3 : (nds[2], nds[3]), 3 : (nds[3], nds[2]),
                        4 : (nds[3], nds[0]), 4 : (nds[0], nds[3])}

#def create_cell_face_map_3D():    
#    for cell, faces in cell_face_map.iteritems():


def create_cell_face_map_3D():    
    """From GambiToNek:
    PREFERED NEKTON FORMAT:   - edge/midpt #'s shown in ( )  
         
                       7 ----(6)------- 6      face  corner points                                     
                      / .             / |      1     {1,2,6,5} 
                    (7) .           (5) |      2     {2,3,7,6}     
                    / (11)          / (10)     3     {4,3,7,8}          
                   8 ----(8)------ 5    |      4     {1,4,8,5}                     
                   |    .          |    |      5     {1,2,3,4}                     
                   |    3 .....(2).|... 2      6     {5,6,7,8}                     
                 (12)   .         (9)  /                            
                   | (3)           | (1)
                   | .             | /                                                                   
                   |.              |/                                              
                   4 ----(4)------ 1  
    """
    local_face_number = {'Bottom':1, 1:'Bottom', 
                         'East':2, 2:'East', 
                         'Top':3, 3:'Top', 
                         'West':4, 4:'West', 
                         'South':5, 5:'South', 
                         'North':6, 6:'North'}

    # Local node numbering. 
    # Note that order of points is only relevant for the Bottom face that we 
    # place first because all other faces are placed based on matching
    # nodes with the Bottom face.
    face_dict_numbering = dict(
        Bottom = array([1, 2, 6, 5], int)[::-1] - 1,
        East = array([2, 3, 7, 6], int)[::-1] - 1,
        Top = array([4, 3, 7, 8], int)[::-1] - 1,
        West = array([1, 4, 8, 5], int)[::-1] - 1,
        South = array([1, 2, 3, 4], int)[::-1] - 1,
        North = array([5, 6, 7, 8], int)[::-1] - 1
        )
            
    # Create a map from nodes of the first placed face to the remaining faces.
    # local_node_map maps nodes 1, 2, 6, 5 correctly to nodes 4,3,7,8
    # For new face 'West' node 4 is connected to node 1 of the Bottom face etc.
    # As a (somewhat irrelevant ) rule we always place the first face in Bottom.
    local_node_map = dict(West={1:4, 5:8}, 
                          East={2:3, 6:7}, 
                          Top={1:4, 2:3, 5:8, 6:7}, 
                          South={1:4, 2:3}, 
                          North={5:8, 6:7})

    # Map of opposite faces in the box.
    shadow_face = {1:3, 3:1, 2:4, 4:2, 5:6, 6:5}
    
    for cell, faces in cell_face_map.iteritems():
        print "cell, faces: ", cell, faces
        # Place first face
        new_face_pos = 1               # Put the first face in Bottom
        face = faces[0]                # The first face from the list
        face1 = face_list[face - 1]    # Some more info for the face
        nodes1 = array(face1[1], int)  # The nodes of first face
        face_map[cell] = {}            # cell/face to local face number map
        face_map[cell][face] = new_face_pos        
        cm = cell_map[cell] = zeros(8, int)  
        newface = local_face_number[new_face_pos] # Name of new face        
        # Place the first four nodes of the cell
        cm[face_dict_numbering[newface][:]] = nodes1[:]
        nodes_of_first_face = array(cm, int)
        
        # One face and 4 nodes placed, now to the remaining faces
        remaining_faces = face_dict_numbering.keys()
        remaining_faces.remove(newface)   # already placed     
        faces.remove(face)                
        for face in faces: # Loop over remaining faces to be placed
            face2 = face_list[face - 1]
            nodes2 = array(face2[1], int)  # Nodes of new face
            common_nodes = find(map(lambda x: x in nodes1, nodes2))
            new_nodes = find(map(lambda x: x not in cm, nodes2))
            
            if common_nodes.shape[0] == 0:
                # e.g., opposite Bottom is Top
                opposite_face = shadow_face[new_face_pos]  
                face_map[cell][face] = opposite_face
                remaining_faces.remove(local_face_number[opposite_face])
                # Don't place any nodes because we don't know connectivity
                
            elif common_nodes.shape[0] == 2:
                matching_nodes = nodes2[common_nodes]
                pos = []
                for i, node in enumerate(matching_nodes):
                    pos.append(find(node == cm)[0])
                
                # The two common nodes determine where the face is located
                for rface in remaining_faces:
                    rval = face_dict_numbering[rface]
                    if all(map(lambda x: x in rval, pos)):
                        break # break out of loop when we find the position
                remaining_faces.remove(rface)
                
                # Give the new face its local face number in face_map
                face_map[cell][face] = local_face_number[rface] 
                
                # Now place the new nodes in cell_map
                # The connectivity of nodes is known from face_dict_numbering
                for new_node_orig_pos in new_nodes:
                    if new_node_orig_pos < 3:
                        node_to_right = nodes2[new_node_orig_pos + 1]
                    else:
                        node_to_right = nodes2[0]
                    if new_node_orig_pos > 0:
                        node_to_left = nodes2[new_node_orig_pos - 1]
                    else:
                        node_to_left = nodes2[3]

                        
                    if any(node_to_right == nodes_of_first_face):
                        # if the node to the right is already placed
                        nd = find(node_to_right == nodes_of_first_face)[0]
                        cm[local_node_map[rface][nd + 1] - 1] =  \
                                                      nodes2[new_node_orig_pos]
                        
                    else:
                        nd = find(node_to_left == nodes_of_first_face)[0]
                        cm[local_node_map[rface][nd + 1] - 1] = \
                                                      nodes2[new_node_orig_pos]

                    

        if any(x < 0. for x in crss(cm)):
            # We guessed the wrong direction for first face. Reverse direction.
            cm_copy = copy(cm)
            cm[[0, 1, 5, 4]] = cm_copy[[4, 5, 1, 0]]
            cm[[3, 2, 6, 7]] = cm_copy[[7, 6, 2, 3]]
            # Reversing direction switches North and South faces, face_map must be updated accordingly
            for f, p in face_map[cell].iteritems():
                if (p == 5):
                    f_south = f
                elif (p == 6):
                    f_north = f
            face_map[cell][f_south] = 6
            face_map[cell][f_north] = 5
        
def crss2d(nds):
    """Test that we have the correct orientation for 2D mesh."""
    cross_product = []
    nds = nodes[:, array(nds) - 1]
    for i, j, k in zip([2, 3, 1, 4], [4, 1, 3, 2], [1, 2, 4, 3]):
        xy1 = nds[:, i - 1]
        xy2 = nds[:, j - 1]
        xy0 = nds[:, k - 1]
        v1x = xy1[0] - xy0[0]
        v2x = xy2[0] - xy0[0]
        v1y = xy1[1] - xy0[1]
        v2y = xy2[1] - xy0[1]
        cross_product.append(v1x*v2y - v1y*v2x)
    return cross_product
    
def crss(nds):
    """Test that we have the correct orientation for 3D mesh."""
    cross_product = []
    nds = nodes[:, array(nds) - 1]
    for i, j, k, l in zip([2, 3, 1, 4, 6, 7, 5, 8], [4, 1, 3, 2, 8, 5, 7, 6], 
                          [5, 6, 8, 7, 1, 2, 4, 3], [1, 2, 4, 3, 5, 6, 8, 7]):
        p1 = nds[:, i - 1]
        p2 = nds[:, j - 1]
        p3 = nds[:, k - 1]
        p0 = nds[:, l - 1]        
        u = p1 - p0
        v = p2 - p0
        w = p3 - p0        
        c = array([u[1]*v[2] - u[2]*v[1],
                   u[2]*v[0] - u[0]*v[2],
                   u[0]*v[1] - u[1]*v[0]])
        cross_product.append(dot(w, c))  
        
    for i in range(4, 8):
        cross_product[i] = - cross_product[i]
        
    return cross_product
            
def create_periodic_cell_face_map():
    """Create maps of type local face 1 of cell 2 is periodic to local face 2 
    of cell 10."""
    for f0, f1 in periodic_face_map.iteritems():
        # f0, f1 = periodic face0 - face1
        face0 = face_list[f0 - 1]
        face1 = face_list[f1 - 1] # shadow
        nd, nds, cells, bc_type, zone_id = [0,]*2, [0,]*2, [0,]*2, [0,]*2, [0,]*2
        for i, ff in enumerate([face0, face1]):
            nd[i], nds[i], cells[i], bc_type[i], zone_id[i] = ff
        
        cell_face_pair = []
        for i in range(2):
            c = max(cells[i])
            if len(nds[i]) == 2:
                cell_face_pair.append((c, face_map[c][(nds[i][0], nds[i][1])]))
            else:
                cell_face_pair.append((c, face_map[c][eval('f'+str(i))]))
                
        periodic_cell_face_map[cell_face_pair[0]] = cell_face_pair[1]
        periodic_cell_face_map[cell_face_pair[1]] = cell_face_pair[0]
                   
def create_boundary_section(bcs, temperature, passive_scalars, mesh_format):
    
    if mesh_format == 'fenics':
        return
        
    for i in range(len(passive_scalars)):
        passive_scalar_val_map.append({})
                        
    global bcs_copy
    if not bcs:
        bcs_copy = {}
    else:
        bcs_copy = bcs

    for face_number, face in enumerate(face_list):
        nd, nds, cells, bc_type, zone_id = face

        if min(cells) == 0: # Boundary face
            c = max(cells)
            if len(nds) == 2:
                local_face  = face_map[c][(nds[0], nds[1])]
            else:
                local_face  = face_map[c][face_number + 1]
            if bc_type == 8 or bc_type == 12:
                # Face is periodic.
                boundary_map[(c, local_face)] = \
                                    periodic_cell_face_map[(c, local_face)]                                   
                boundary_val_map[(c, local_face)] = 'P'
                if temperature:
                    temperature_val_map[(c, local_face)] = 'P'
                for ss, scalar in enumerate(passive_scalars):
                    passive_scalar_val_map[ss][(c, local_face)] = 'P' 
                        
            else:
                boundary_map[(c, local_face)] = (0, 0)
                if bcs:
                    boundary_val_map[(c, local_face)] = bcs[zone_id]
                else:
                    if (zones[zone_id][1][-1] == 'S'):
                        # Fix for symmetry boundary condition
                        boundary_val_map[(c, local_face)] = 'SYM'
                    else:
                        boundary_val_map[(c, local_face)] = zones[zone_id][1][-1]
                    bcs_copy[zone_id] = zones[zone_id][1][-1]
                    # RUD PROJ ALG
                    if (zones[zone_id][1][-2] == 'E'):    
                     # A list containing all elements and the local 
                     #face that are to be projected on an exact surface from ICEM
                     surf_list.append([c,local_face]) 
                if temperature:
                    temperature_val_map[(c, local_face)] =  \
                                                    temperature[zone_id]
                for ss, scalar in enumerate(passive_scalars):
                    passive_scalar_val_map[ss][(c, local_face)] = \
                                            passive_scalars[ss][zone_id] 

        else:
            c0 = cells[0]
            c1 = cells[1]
            if len(nds) == 2:
                local_face0  = face_map[c0][(nds[0], nds[1])]
                local_face1  = face_map[c1][(nds[0], nds[1])]
            else: #3D
                local_face0 = face_map[c0][face_number + 1]
                local_face1 = face_map[c1][face_number + 1] 
                
            if bc_type == 2:
                # interior
                boundary_map[(c0, local_face0)] = (c1, local_face1)
                boundary_map[(c1, local_face1)] = (c0, local_face0)
                boundary_val_map[(c0, local_face0)] = 'E'
                boundary_val_map[(c1, local_face1)] = 'E'
                
                if temperature:
                    temperature_val_map[(c0, local_face0)] =  'E'
                    temperature_val_map[(c1, local_face1)] =  'E'
                
                for ss, scalar in enumerate(passive_scalars):
                    passive_scalar_val_map[ss][(c0, local_face0)] = 'E'
                    passive_scalar_val_map[ss][(c1, local_face1)] = 'E'
            
            else:
                raise NotImplementedError
            
def circle_center(x, x0, x1, r):
    return (r**2 - (x[0] - x0[0])**2 - (x[1] - x0[1])**2, 
            r**2 - (x[0] - x1[0])**2 - (x[1] - x1[1])**2)

def barycentric_weight(y):
    N = len(y)
    w = ones(N, float)
    for j in range(1, N):
        for k in range(0, j):
            w[k] = (y[k] - y[j])*w[k]
            w[j] *= (y[j] - y[k])
    w = 1./w
    return w

def fun(x0, x1, y0, y1, xx, yy):
    """w are barycentric weights.
    x is unknown.
    x0, y0 contain x and y in the four nodes used to compute the midpoint."""    

    # Look for point of intersection between interpolated curve between nodes in x, y
    # and the normal to the face between nodes (x0, y0) and (x1, y1)
    # Transform coordinate axes
    # Center of face is xs, ys
    xs = (x0 + x1)/2.
    ys = (y0 + y1)/2.

    if abs(y1 - y0) > abs(x1 - x0):
        theta = arctan((x1 - x0)/(y1 - y0))
        theta2 = arctan((xx - xs)/(yy - ys))
        dy = (yy - ys)/cos(theta2)
        xn = copy(xx)
        yn = copy(yy)
        xn = dy*sin(theta2 - theta)
        yn = dy*cos(theta2 - theta)
        w = barycentric_weight(yn)
        y2 = - yn
        f = zeros(len(y2), float)
        ss = sum(w/y2)
        f[:] = w/y2/ss
        dy = dot(f, xn)
        xny = xs + dy*sin(theta + pi/2.)
        yny = ys + dy*cos(theta + pi/2.)

    else:        
        theta = arctan((y1 - y0)/(x1 - x0))
        theta2 = arctan((yy - ys)/(xx - xs))
        dx = (xx - xs)/cos(theta)
        xn = copy(xx)
        yn = copy(yy)
        xn = dx*cos(theta2 - theta)
        yn = dx*sin(theta2 - theta)
        w = barycentric_weight(xn)
        x2 = - xn
        f = zeros(len(x2), float)
        ss = sum(w/x2)
        f[:] = w/x2/ss
        dy = dot(f, yn)
        xny = xs + dy*cos(theta + pi/2.)
        yny = ys + dy*sin(theta + pi/2.)
    
    return xny, yny

def find_mid_point(curve_zone, curve):
    
    # Get unique list of nodes on boundary
    nodes0 = boundary_nodes[curve_zone]
    face0_list = boundary_faces[curve_zone]
    
    # Loop over faces and compute midpoint using at least two 
    # nodes of two neighbouring faces
    for i in face0_list:
        face = face_list[i] 
        nds = face[1]
        face_node_map[i] = {'left': None, 'right': None}
        for j in face0_list:
            face2 = face_list[j] 
            if len(Set(nds + face2[1])) == 3:
                if nds[0] in face2[1]:
                    face_node_map[i]['left'] = j
                else:
                    face_node_map[i]['right'] = j                   
                    
    for face in face_node_map:
        if face in curved_faces[curve_zone]:
            fl, fr = face_node_map[face]['left'], face_node_map[face]['right']
            if fl is not None and fr is not None:
                common_node_left = map(lambda x: x in face_list[face][1], face_list[fl][1])            
                node0 = face_list[fl][1][common_node_left.index(False)]
                node1 = face_list[fl][1][common_node_left.index(True)]
                common_node_right = map(lambda x: x in face_list[face][1], face_list[fr][1])            
                node2 = face_list[fr][1][common_node_right.index(True)]
                node3 = face_list[fr][1][common_node_right.index(False)]        
                x, y = nodes[:, array([node0, node1, node2, node3]) - 1]            
                x0, y0 = x[1], y[1]
                x1, y1 = x[2], y[2]            
                
            elif fl is None:
                common_node_right = map(lambda x: x in face_list[face][1], face_list[fr][1])            
                node0 = face_list[face][1][common_node_right.index(True)]
                node1 = face_list[fr][1][common_node_right.index(True)]
                node2 = face_list[fr][1][common_node_right.index(False)]        
                x, y = nodes[:, array([node0, node1, node2]) - 1]        
                x0, y0 = nodes[:, array(face_list[face][1]) - 1]            
                x0, x1 = x0[:]
                y0, y1 = y0[:]
                            
            elif fr is None:
                common_node_left = map(lambda x: x in face_list[face][1], face_list[fl][1])            
                node0 = face_list[face][1][common_node_left.index(True)]
                node1 = face_list[fl][1][common_node_left.index(True)]
                node2 = face_list[fl][1][common_node_left.index(False)]        
                x, y = nodes[:, array([node0, node1, node2]) - 1]
                x0, y0 = nodes[:, array(face_list[face][1]) - 1]
                x0, x1 = x0[:]
                y0, y1 = y0[:]
                    
            xn, yn = fun(x0, x1, y0, y1, x, y)                
            mid_point_map[face] = (xn, yn, 0, 0, 0)

def count_lines_startwith(lines, string):
    # Count the occurences of string in the list of lines
    count = 0
    for line in lines:
        if line.startswith(string):
            count += 1
    return(count)

def scan_fluent_mesh(lines):
    """Scan fluent mesh and generate numerous maps."""
    # Warning! Not yet tested for multiple interior zones
    dim = 0
    one = 0
    num_faces = 0
    while 1:
        try:
            line = lines.pop(0)
        except:
            print 'Finished reading file\n'
            break
        if dim == 0: # Dimension usually comes first
            a = re.search(re_dimline, line)
            if a: 
                print 'Reading dimensions\n'
                dim = int(a.group(1))
                print 'Mesh is ' + str(dim) + 'D\n'
                continue
        
        if one == 0: # The total number of nodes
            a = re.search(re_zone0, line)
            if a:
                print 'Reading zone info\n'
                one, num_vertices, dummy1, dummy2 = int(a.group(1)), \
                     int(a.group(2), 16), int(a.group(3), 16), int(a.group(4))
                continue
            
        a = re.search(re_zone, line) # Nodes
        if a:
            zone_id, first_id, last_id, dummy1, dummy2 = int(a.group(1), 16), \
                int(a.group(2), 16), int(a.group(3), 16), int(a.group(4)), \
                int(a.group(5))
            print 'Reading nodes in zone ', zone_id + 1, '\n'
            read_zone_nodes(dim, first_id, last_id, lines)
            lines = lines[(last_id - first_id + 1):]  
            continue
            
        a = re.search(re_zones,line) # Zone info
        if a:
            print 'Reading zone ', line
            dummy, zone_id, zone_type, zone_name, radius =  \
                       int(a.group(1)), int(a.group(2)),  a.group(3), \
                       a.group(4), a.group(5)
            zones[zone_id] = [zone_type, zone_name, radius]
            continue
        
        a = re.search(re_cells0, line) # Get total number of cells/elements
        if a:
            print 'Reading cell info ', line
            first_id, tot_num_cells = int(a.group(3),16), int(a.group(5), 16)
            continue

        a = re.search(re_cells,line) # Get the cell info.
        if a:
            print a.group(1), a.group(2), a.group(3), a.group(4), a.group(5)
            zone_id, first_id, last_id, bc_type, element_type = \
                int(a.group(1),16), int(a.group(2), 16), int(a.group(3), 16), \
                int(a.group(4), 16), int(a.group(5), 16)
            print 'Reading cells in zone ', zone_id, '\n'
            if last_id == 0:
                raise TypeError("Zero elements!")
            num_cells[zone_id] = [first_id, last_id, bc_type, element_type]
            continue

        a = re.search(re_cells2,line) # Get the cell info.
        if a:
            raise TypeError("Wrong cell type. Can only handle one single cell type")

        a = re.search(re_face0, line)
        if a:
            print 'Reading total number of faces\n', line
            num_faces = int(a.group(3),16)
            continue
            
        a = re.search(re_face, line)
        if a:
            print 'Reading faces ', line
            zone_id, first_id, last_id, bc_type, face_type = \
                 int(a.group(2), 16), int(a.group(3),16), int(a.group(4),16), \
                 int(a.group(5), 16), int(a.group(6), 16)
            read_faces(zone_id, first_id, last_id, bc_type, face_type, lines)
            lines = lines[(last_id - first_id + 1):]
            zone_number_of_faces[zone_id] = last_id - first_id + 1
            continue
        
        a = re.search(re_periodic, line)
        if a:
            print 'Reading periodic connectivity\n', line
            read_periodic(lines, periodic_dx)
            continue
        
        print 'Line = ',line
        if any([re.search(st, line) for st in (re_parant, re_comment)]) or \
                                                             not line.strip():   
            continue
                
        # Should not make it here
        print 'Line = ',line
        raise IOError('Something went wrong reading fluent mesh.')
    
def scan_nastran_mesh(lines,tol):
    """Scan nastran mesh and generate numerous maps."""
    # (Warning! Not yet tested for multiple interior zones) (Fluent)
    global tot_num_nodes
    dim = 3    # 3-d assumed
    print 'Mesh is ' + str(dim) + 'D\n'
    num_faces = 0
    num_vertices = count_lines_startwith(lines, 'GRID')
    tot_num_cells = count_lines_startwith(lines, 'CHEXA')   # Only HEXA grids considered

    chexaline = 0
    node_no = 0
    cell_no = 0
    zone_count = -1
    
    while 1:
        try:
            line = lines.pop(0)
        except:
            if ((num_vertices == node_no) and (tot_num_cells == cell_no)):
                tot_num_nodes = node_no
                print "All", tot_num_nodes, "nodes read"
                print "All", cell_no, "cells read"
                print 'Finished reading file\n'
                break
            else:
                raise IOError('Something went wrong reading Nastran mesh, missing nodes or cells')

        if line.startswith('GRID'):
            currentID = int(line[9:23])
            node_no += 1
            # Fixing expontial format (.1-17 --> .1e-17)
            field = line[24:32]
            efield = re.sub('(.+?)(-)(.+)','\\1E-\\3',field)
            currentX = float(efield)
            field = line[32:40]
            efield = re.sub('(.+?)(-)(.+)','\\1E-\\3',field)
            currentY = float(efield)
            field = line[40:48]
            efield = re.sub('(.+?)(-)(.+)','\\1E-\\3',field)
            currentZ = float(efield)
            if (node_no == 1):
                global nodes
                nodes = zeros((dim, num_vertices))
            nodes[:,node_no-1] = [currentX, currentY, currentZ]
            continue
            
        a = re.search(re_zones,line) # Zone info
        if a:
            print 'Reading zone ', line
            dummy, zone_id, zone_type, zone_name, radius =  \
                       int(a.group(1)), int(a.group(2)),  a.group(3), \
                       a.group(4), a.group(5)
            zones[zone_id] = [zone_type, zone_name, radius]
            continue
        
        a = re.search('family', line) # Zone info in Nastran file
        if a:
            zone_count += 1
            ln = line.split()
            zone_name = ln[len(ln)-1]
            print 'zone_name = ', zone_name
            if (zone_count == 0):
                zone_type = 'E'   # Internal boundary
            else:
                zone_type = zone_name[-1]
            radius = 0
            zones[zone_count] = [zone_type, zone_name, radius]
            continue

        if line.startswith('PLOAD'):  # Use this for all outer boundaries
            ln = line.split()
            verts = sorted([int(ln[3]), int(ln[4]), int(ln[5]), int(ln[6])])
            # Update zone number for face
            set_boundary_zone(zone_count, verts)
            continue
            

        a = re.search(re_cells,line) # Get the cell info.
        if a:
            print a.group(1), a.group(2), a.group(3), a.group(4), a.group(5)
            zone_id, first_id, last_id, bc_type, element_type = \
                int(a.group(1),16), int(a.group(2), 16), int(a.group(3), 16), \
                int(a.group(4), 16), int(a.group(5), 16)
            print 'Reading cells in zone ', zone_id, '\n'
            if last_id == 0:
                raise TypeError("Zero elements!")
            num_cells[zone_id] = [first_id, last_id, bc_type, element_type]
            continue

        if line.startswith('CHEXA'):
            vertices = zeros(21, int)
#            cell_no = int(line[8:16])
            cell_no += 1
            vertices[1] = int(line[24:32])
            vertices[2] = int(line[32:40])
            vertices[3] = int(line[40:48])
            vertices[4] = int(line[48:56])
            vertices[5] = int(line[56:64])
            vertices[6] = int(line[64:72])
            chexaline = 1
            continue

        if (chexaline == 1):
            # Continue reading CHEXA node numbers
            vertices[7] = int(line[8:16])
            vertices[8] = int(line[16:24])
            # Check for quadratic element:
            if (len(line) > 24):
                vertices[9] = int(line[24:32])
                vertices[10] = int(line[32:40])
                vertices[11] = int(line[40:48])
                vertices[12] = int(line[48:56])
                vertices[13] = int(line[56:64])
                vertices[14] = int(line[64:72])
                chexaline = 2
            else:
                vertices[0] = 8
                # Process cell
                process_cell(cell_no, vertices,tol)
                chexaline = 0
            continue
            
        if (chexaline == 2):
            # Continue reading CHEXA node numbers
            vertices[15] = int(line[8:16])
            vertices[16] = int(line[16:24])
            vertices[17] = int(line[24:32])
            vertices[18] = int(line[32:40])
            vertices[19] = int(line[40:48])
            vertices[20] = int(line[48:56])
            vertices[0] = 20
            # Process faces
            process_cell(cell_no, vertices,tol)
            chexaline = 0
            continue

        a = re.search(re_cells2,line) # Get the cell info.
        if a:
            raise TypeError("Wrong cell type. Can only handle one single cell type")

        a = re.search(re_face0, line)
        if a:
            print 'Reading total number of faces\n', line
            num_faces = int(a.group(3),16)
            continue
            
        a = re.search(re_face, line)
        if a:
            print 'Reading faces ', line
            zone_id, first_id, last_id, bc_type, face_type = \
                 int(a.group(2), 16), int(a.group(3),16), int(a.group(4),16), \
                 int(a.group(5), 16), int(a.group(6), 16)
            read_faces(zone_id, first_id, last_id, bc_type, face_type, lines)
            lines = lines[(last_id - first_id + 1):]
            zone_number_of_faces[zone_id] = last_id - first_id + 1
            continue
        
        a = re.search(re_periodic, line)
        if a:
            print 'Reading periodic connectivity\n', line
            read_periodic(lines, periodic_dx)
            continue
        
#        print 'Line = ',line.strip()
        if any([re.search(st, line) for st in (re_parant, re_comment)]) or \
                                                             not line.strip():   
#            print '(Ignored)' # CEW
            continue
                
        # Should not make it here
        print 'Line = ',line.strip()
        print 'Unidentified'
        continue
#        raise IOError('Something went wrong reading fluent mesh.')

def write_nek5000_file(dim, ofilename, curves, temperature, passive_scalars,start_of_rea,end_of_rea,curve_type):
    tot_num_cells = len(cell_map)
    ofile  = open(ofilename + '.rea', "w")
    ## Put the mesh in a rea-file
    print 'Create the rea-file: %s\n' %(ofilename+'.rea')
    print 'Writing the elements info\n'
    ofile.write(start_of_rea.format(dim))
    ofile.write(' **MESH DATA**\n')
    ofile.write('       %s       %s       %s      NEL,NDIM,NELV\n' 
                                    %(tot_num_cells, dim, tot_num_cells))
    element_header = '      ELEMENT          %s [    1 ]    GROUP     0\n'    
    for i in range(tot_num_cells):
        ofile.write(element_header %(i + 1))
        if dim == 2:
            n1 = nodes[0, array(cell_map[i + 1]) - 1]
            n2 = nodes[1, array(cell_map[i + 1]) - 1]
            ofile.write(reduce(add, ['{0:.8e}'.format(x).rjust(16) 
                                    for x in n1]) + '\n')
            ofile.write(reduce(add, ['{0:.8e}'.format(x).rjust(16) 
                                    for x in n2]) + '\n')                                                       
        else:
            # Write only the 8 first vertices, not the midside nodes
            nn = nodes[:, array(cell_map[i + 1][0:8]) - 1]
            for i in range(3):
                ofile.write(reduce(add, ['{0:.8e}'.format(x).rjust(16) 
                                        for x in nn[i, :4]]) + '\n')
            for i in range(3):
                ofile.write(reduce(add, ['{0:.8e}'.format(x).rjust(16) 
                                        for x in nn[i, 4:]]) + '\n')
    print 'Writing the curved side info\n'
    ofile.write("  ***** CURVED SIDE DATA *****  \n")
    cc = "{0:6d} Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE \n"
    if tot_num_cells < 1000:
        c1 = "{0:3d}{1:3d}"
    elif tot_num_cells < 1e06:
        c1 = "{0:2d}{1:6d}"
    else:
        c1 = "{0:2d}{1:12d}"
    c2 = "{0:14.6e}{1:14.6e}{2:14.6e}{3:14.6e}{4:14.6e} {5:s}\n"    
    # MIDPOINT NOTATION 
    if(curve_type=='m'):
		ofile.write(cc.format(tot_num_curved))
		for ic in range(tot_num_cells):
                        #for ie in curved_midpoint[ic+1].keys():
                                #xx = nodes[:, curved_midpoint[ic+1][ie]]
                        for ie in Cells[ic].curved_edge:
                                xx = nodes[:,Cells[ic].curved_midpoint[ie]] 
				print 'Curved line for edge {} cell {}'.format(ie,ic)
				ofile.write(c1.format(ie, ic+1))
				# The next 2 lines should provide the identical values
				ofile.write(c2.format(xx[0], xx[1], xx[2], 0.0, 0.0, 'm'))
    elif(curve_type=='c'):
        ofile.write(cc.format(tot_num_curved))
        for ic in range(tot_num_cells):
            for ie in curved_midpoint[ic+1].keys():
                xx = nodes[:, curved_midpoint[ic+1][ie]]
                xxwest= nodes[:, curved_east[ic+1][ie]]
                xxeast= nodes[:, curved_west[ic+1][ie]]
            #for ie in Cells[ic].curved_edge:
                #xx = nodes[:,Cells[ic].curved_midpoint[ie]] 
                #xxwest = nodes[:,Cells[ic].curved_west[ie]] 
                #xxeast = nodes[:,Cells[ic].curved_east[ie]] 
                ofile.write(c1.format(ie, ic+1))
                radius ,center = points2circ(xxwest,xx,xxeast)
                print 'edge number: {}'.format(ie)
                print 'mid node: {}, west nodes: {},east node: {}'.format(Cells[ic].curved_midpoint[ie],Cells[ic].curved_west[ie],Cells[ic].curved_east[ie])
						# SPHERE                                                                               
                #ofile.write(c2.format(center[0],center[1],center[2],radius, 0.0, 's'))
                # CIRCLE 
                ofile.write(c2.format(radius,0.0,0.0,0.0,0.0,'C'))
                #ofile.write(c2.format(-0.05,0.0,0.0,0.0,0.0,'C'))
                ## Printing data
                plt.plot([xxwest[0],xx[0],xxeast[0]],[xxwest[1],xx[1],xxeast[1]],'r')
                #plt.plot(center[0],center[1],'b*')
        plt.show()
    elif(curve_type=='A'):
        ofile.write(cc.format(tot_num_curved))
        for ic in range(tot_num_cells):
            for ie in curved_midpoint[ic+1].keys():
                xx = nodes[:, curved_midpoint[ic+1][ie]]
                xxwest= nodes[:, curved_east[ic+1][ie]]
                xxeast= nodes[:, curved_west[ic+1][ie]]
                ofile.write(c1.format(ie, ic+1))
                radius ,center = points2circ(xxwest,xx,xxeast)
                print 'edge number: {}'.format(ie)
                print 'mid node: {}, west nodes: {},east node: {}'.format(Cells[ic].curved_midpoint[ie],Cells[ic].curved_west[ie],Cells[ic].curved_east[ie])
                ofile.write(c2.format(center[0],center[1],center[2],radius,0.0,'A'))

    # SPHERE NOTATION 
    #####ADDED BY RUD 25.09 ######
    print 'Printing circle information\n'
    #print face_list[0]
    #print face_list[1]
    #print face_list[2]
    #print face_list[3]
   #print face_list[4]
    #print face_list[5]
    #print face_list[6]
    #print face_list[7]
    #print face_list[8]
    #print face_list[9]
    #print face_list[10]
    #print face_list[11]
    #print face_list[12]
    #print face_list[13]
    #print '------------------\n'
    #print parallel_edge_map[1][0]
    #ofile.write(cc.format(tot_num_curved))
    #for ic in range(tot_num_cells):
                  ##print cell_face_map[ic+1]
                ##print curved_midpoint[ic+1] 
                #for ie in curved_midpoint[ic+1].keys():
                        ##print ic+1,ie 
                        ##print face_list[dumdum[0]]
                        ##print ic,'\n'
                        #ofile.write(c1.format(ie, ic+1))
                        #xx = nodes[:, curved_midpoint[ic+1][ie]]
                        #xxwest= nodes[:, curved_east[ic+1][ie]]
                        #xxeast= nodes[:, curved_west[ic+1][ie]]
                        #radius ,center = points2circ(xxwest,xx,xxeast)
                        ## SPHERE 
                        ##ofile.write(c2.format(center[0],center[1],center[2],radius, 0.0, 's'))
                        ## CIRCLE 
                        #ofile.write(c2.format(radius,0.0,0.0,0.0,0.0, 'c'))
                        ## Printing data
                        ##plt.plot([xxwest[0],xx[0],xxeast[0]],[xxwest[1],xx[1],xxeast[1]],'r')
                        ##plt.plot(center[0],center[1],'b*')
    ##plt.show()
                        
    #####END BY RUD 25.09 ######

            
    ## for zone in curves:
    ##     if curves_map[zone][1]['type'] == 'C':
    ##         for i in range(len(curves_map[zone][0])):
    ##             xx = curves_map[zone][1]['x'][i]
    ##             ofile.write(c1.format(curves_map[zone][0][i][1], 
    ##                                 curves_map[zone][0][i][0]))
    ##             ofile.write(c2.format(xx[0], xx[1], xx[2], xx[3], xx[4],
    ##                                 curves_map[zone][1]['type']))
    ##     elif curves_map[zone][1]['type'] == 'm':
    ##         for i in range(zone_number_of_faces[zone]):     
    ##             if curves_map[zone][0][i][2] in curved_faces[zone]:
    ##                 ofile.write(c1.format(curves_map[zone][0][i][1], 
    ##                                     curves_map[zone][0][i][0]))             
    ##                 xx = mid_point_map[curves_map[zone][0][i][2]]
    ##                 ofile.write(c2.format(xx[0], xx[1], xx[2], xx[3], xx[4],
    ##                                     curves_map[zone][1]['type']))

    print 'Writing the boundary and bottom section\n'
    ofile.write('  ***** BOUNDARY CONDITIONS ***** \n')
    ofile.write('  ***** FLUID   BOUNDARY CONDITIONS ***** \n')
    f_str = " {0:s}  "
    f_str3 = " {0:3s}"
    if tot_num_cells < 1000:
        f_str1 = "{0:3d}{1:3d}"
    elif tot_num_cells < 100000:
        f_str1 = "{0:5d}{1:1d}"
    elif tot_num_cells < 1000000:
        f_str1 = "{0:6d}"    
    else:
        f_str1 = "{0:12d}"    
    f_str2 = "{0:14.7e}{1:14.7e}{2:14.7e}{3:14.7e}{4:14.7e}\n"
    for i in range(1, tot_num_cells + 1):
        for j in range(1, 2*dim + 1):
            bm = boundary_map[(i, j)]
            if (len(boundary_val_map[(i, j)]) == 3):
                ofile.write(f_str3.format(boundary_val_map[(i, j)]))
            else:
                ofile.write(f_str.format(boundary_val_map[(i, j)]))
            if tot_num_cells < 1e5:
               ofile.write(f_str1.format(i, j))
            else:
               ofile.write(f_str1.format(i))
            ofile.write(f_str2.format(bm[0], bm[1], 0, 0, 0))
    if not temperature:
        ofile.write('  ***** NO THERMAL BOUNDARY CONDITIONS ***** \n')
    else:
        ofile.write('  ***** THERMAL BOUNDARY CONDITIONS ***** \n')
        for i in range(1, tot_num_cells + 1):
            for j in range(1, 2*dim + 1):
                bm = boundary_map[(i, j)]
                ofile.write(f_str %(temperature_val_map[(i, j)], i, j, bm[0],
                                    bm[1], 0, 0, 0))    
    if len(passive_scalars) == 0:
        pass
    else:
        for si, scalar in enumerate(passive_scalars):
            ofile.write(' ***** PASSIVE SCALAR %s  BOUNDARY CONDITIONS ***** \n' 
                        %(si + 1))
            for i in range(1, tot_num_cells + 1):
                for j in range(1, 2*dim + 1):
                    bm = boundary_map[(i, j)]
                    ofile.write(f_str %(passive_scalar_val_map[si][(i, j)], i, 
                                        j, bm[0], bm[1], 0, 0, 0))
    
    ofile.write(end_of_rea)
    print 'Finished writing the rea-file\n'
    # Close files
    ofile.close()

def convert(nastranmesh, 
        periodic_dx,                            # nek5000 and semtex
        func=None, 
        mesh_format='nek5000',                     # nek5000, semtex or fenics
        curves = {}, bcs = False,                  # nek5000 and semtex
        temperature=False, passive_scalars=[],     # nek5000 only
        cylindrical=1, NZ=1,                       # semtex  only
        reafile='def.rea',outfile='out.rea',curve_type = 'm',tol=1e-02): 
    #Converts a fluent mesh to a mesh format that can be used by Nek5000
       #semtex or FEniCS. 

         #nastranmesh = fluent mesh (*.msh file)

               #func = Optional function of spatial coordinates (x,y) that can 
                      #be used to modify the fluent mesh.
                      #For example, say you have a mesh that is a rectangular 
                      #geometry with -2 < x < 6 and 0 < y 0.5. Now you want to
                      #squeeze this mesh around x = 0 to create a stenosis type
                      #of mesh. This can be accomplished by squeezing the mesh
                      #in the y-direction through:

                      #def func_y(x, y):
                          #if abs(x) < 1.:
                              #return y*(1 - 0.25*(1 + cos(x*pi)))
                          #else:
                              #return y

                      #and supplying this value to the keyword func like:

                      #func={'y': func_y}

                      #Otherwise you might just create this stenosis in your
                      #mesh generation software and supply nothing for func.
                      #Note that in nek5000 you will most likely do this in
                      #subroutines userdat or userdat2.

        #mesh_format = 'nek5000', 'semtex' or 'fenics'

                #bcs = False or dictionary of boundary conditions for 
                      #velocity/pressure (something like {1: 'W', 2: 'v'} for
                      #wall in zone 1 and Dirchlet to be specified in zone 2).
                      #False indicates that dictionary is not used and
                      #in that case we assume the name of the zone ends in
                      #the correct letter, like 'somename_W' for zone 1 and 
                      #'somename_v' for zone 2. Zonenames are easy to modify
                      #at the bottom of the fluent msh files.
                      #Don't include periodic zones here.

        #periodic_dx = Dictionary describing any periodicity in the mesh.
                      #Keys are tuples of the two periodic zones (zone0, zone1)
                      #and values are displacement vectors.
                      #Note that the current program also can read a mesh where
                      #periodicity has been generated in the mesh generator. 
                      #However, this author has still not managed to 
                      #create a periodic 3D mesh correctly in any mesh software.
                      #Hence I prefer to define periodicity through this
                      #dictionary here and do nothing regarding periodicity in
                      #the meshing software. All you need to know are the ids
                      #and displacement of the periodic zones. Connectivity
                      #will then be computed here. Note that each periodic zone
                      #needs to be its own zone. Simply creating a 3D UnitCube
                      #in gambit and not assigning names to the 6 zones won't
                      #work. We need zone identifiers (for now).
                      #Not for FEniCS.

             #curves = Dictionary of curve information. Keys are curve zones 
                      #and value is either 
                               #{'type': 'C', 
                                #'radius': radius, 
                                #'circle_center': (x, y), 
                                #'depth': depth} 
                      #for a circle or 
                               #{'type': 'm'} 
                      #for midpoint with nek5000 or
                               #{'type': 'spline'}
                      #for a curved side with semtex. Here a .geom file containing
                      #the spline information will be created.

                      #The circle may provide the radius or the center of the
                      #circle through 'circle_center'. The curvature may also be 
                      #used in the internal elements inside the surface through
                      #specifying the depth. depth=4 means that the curvature 
                      #is used throughout the first four elements inside that 
                      #surface. This is necessary to get good quality meshes in,
                      #e.g., a cylinder. The radius for an internal face is 
                      #allways computed as the distance to the circle_center.
                      #Not for FEniCS.

        #temperature = False or dictionary of boundary conditions for 
                      #temperature (something like {1: 'W', 2: 'v'}) for
                      #Wall in zone 1 and Dirchlet specified in .usr in zone 2.
                      #False indicates that temperature is not used. 
                      #(nek5000 only)

    #passive_scalars = [] or list of dictionaries of boundary conditions for
                      #passive scalars. Empty means no passive scalars are used.
                      #(nek5000 only)

        #cylindrical = 1 for cylindrical mesh and 0 for regular (semtex only)

                 #NZ = Number of planes in homogeneous direction (semtex only)

# ADDED BY RUD 25.09.15
    global nodes, tot_num_curved
    print tot_num_curved
# This part is just a default that goes into the top of the rea-file:
    if (reafile != 'def.rea'): 
        start_of_rea=getreastart(reafile)
        end_of_rea =getreaend(reafile)
    else:
        print 'REMEMBER TO INPUT A REA-FILE IF YOU HAVE SOME STORED SETTINGS'
        start_of_rea = rea_start
        end_of_rea = rea_end
# END RUD 25.09.15
    #periodic_dx={(3,6):[0,1,0]}
    if(outfile != 'out.rea'): ofilename = outfile[:-4]
    else: ofilename = nastranmesh[:-4]
    ifile  = open(nastranmesh, "r")

    zone_id = 0
    if not nodes:
        # Read all lines of nastran mesh
        lines = ifile.readlines()
        if len(lines) == 0:
            raise IOError("Empty nastran mesh file")

        scan_nastran_mesh(lines,tol)

    dim = nodes.shape[0]
    create_face_list()
    #    create_cell_face_map(dim, mesh_format)
    #print "periodic_dx:", periodic_dx
    #print "periodic_dx keys:", periodic_dx.keys()
    #print "periodic_dx values:", periodic_dx.values()
    create_periodic_face_map(periodic_dx)
    create_periodic_cell_face_map()
    create_boundary_section(bcs, 0, passive_scalars, mesh_format)

# ADDED BY RUD 06.10.15
    for cell_no in range(len(Cells)):
        (Cells[cell_no]).define_face_mappings(Faces)
    find_neighbours()
    # Checking the cell mapping
    #for cell_no in range(len(Cells)):
        #print 'faces for cell number {}:'.format(cell_no+1)
        #print Cells[cell_no].faces
        #print Cells[cell_no].face_glloc
        #for ie in range(6):
            #print Cells[cell_no].faces[ie+1],Cells[cell_no].face_glloc[Cells[cell_no].faces[ie+1]]

    # Method for assigning the neightbours
    if(0):
        print 'propagation is activated \n'
        for cell in Cells:
            cell.create_inf_faces()
            cell.glob_nodes_inf(nodes)
            for i in cell.infliction_faces:
                cell.sendInfo(Cells[cell.nbours[i]-1])
        #Updating info
        for cell in Cells:
                    nodes, tot_num_curved = cell.updateInfo(nodes,tot_num_curved)
        #(Cells[89]).print_info()
        # SINGLE ELEMENT CASE
        #(Cells[357]).print_info()
        #Cells[357].glob_nodes_inf(nodes)
        #for i in Cells[357].infliction_faces:
             #Cells[357].sendInfo(Cells[Cells[357].nbours[i]-1])
        #nodes, tot_num_curved = (Cells[330]).updateInfo(nodes,tot_num_curved)

        #(Cells[330]).print_info()
        # SINGLE ELEM CASE END
        #print 'final totnum {}'.format(tot_num_curved)
        for cell in Cells:
                    if (cell.curved_edge): print 'cell number {} '.format(cell.cellID)
                    for ie in cell.curved_edge:
                            print 'east {}, west {}, mid {}'.format(cell.curved_east[ie],cell.curved_west[ie],cell.curved_midpoint[ie])

        
    #print ('cellID of elem 1 is {}'.format(Cells[0].cellID))
# END RUD 06.10.15
    # Modify the entire mesh using the shape-function func
    if func:
        sz = nodes.shape
        for i in range(sz[1]):
            x, y = nodes[:, i]
            if 'x' in func:
                xnew = func['x'](x, y)
                if abs(xnew - x) > 1e-6:
                    nodes[0, i] = xnew
            if 'y' in func:
                ynew = func['y'](x, y)
                if abs(ynew - y) > 1e-6:
                    nodes[1, i] = ynew
        if mesh_format == 'nek5000':
            print 'Warning!! Consider using userdat/userdat2 instead!'

        #if not mesh_format == 'fenics':
        # read_curved_sides(curves)

    # Generate the mesh files for given mesh format
    if mesh_format == 'nek5000':
        write_nek5000_file(dim, ofilename, curves,0, passive_scalars,start_of_rea,end_of_rea,curve_type)
    elif mesh_format == 'surface':
        write_surface_file();
    elif mesh_format == 'semtex':
        write_semtex_file(dim, ofilename, curves, cylindrical, NZ)
    if mesh_format == 'fenics':
        write_fenics_file(dim, ofilename)

    ifile.close()
	# ADDED BY RUD 25.09.15
    if(mesh_format != 'surface'):
        print 'Fixing symmetry and inflow conditions \n '
        fixbc(ofilename + '.rea')

        print 'Fixing thermal inflow conditions \n '
        fixthermalbc(ofilename + '.rea',temperature)

        if(surf_list): 
            write_surf_list()

	# END ADDED BY RUD 25.09.15
