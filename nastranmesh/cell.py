from numpy import zeros,array,dot,sqrt
from collections import defaultdict
from globals import *
# Global dictionaries
parallel_edge_map = {1:[3,5], 2:[4,6], 3:[1,7], 4:[2,8], 5:[1,7], 6:[2,8], 7:[3,5], 8:[4,6], 9:[10,12], 10:[9,11], 11:[10,12], 12:[9,11]} # Topologically parallel edges by nek5000 edge numbering (not including the "opposite" egde)
edge2points = {1:[1,2], 2:[2,3], 3:[3,4], 4:[1,4], 5:[5,6], 6:[6,7], 7:[7,8], 8:[8,9], 9:[1,5], 10:[2,6], 11:[3,7], 12:[4,8]}
face2points = {1:[1,2,6,5], 2:[2,3,7,6], 3:[4,3,7,8], 4:[1,4,8,5], 5:[1,2,3,4], 6:[5,6,7,8]}
edge2faces = {1:[1,5], 2:[2,5], 3:[3,5], 4:[4,5], 5:[1,6], 6:[2,6], 7:[3,6], 8:[4,6], 9:[4,1], 10:[1,2], 11:[2,3], 12:[3,4]}
shadow_face = {1:3, 3:1, 2:4, 4:2, 5:6, 6:5}
edge2mid = {1:9,2:10,3:11,4:12,5:17,6:18,7:19,8:20,9:13,10:14,11:15,12:16}
class cell:
	def __init__(self,ID,verts):
		self.cellID = ID    # Id number
		self.vertices = verts#List of vertices,[#tot,8 vertices, 12 midpoints] All global
		self.curved_edge = [] # curved edge numbers
		self.curved_midpoint = {} #midpoint info 
		self.curved_west = {} # west corner node info (rel to curved midpoint)
		self.curved_east = {} # east corner node info (rel to curved midpoint)
		self.curved_info = {} # info to be sent to neighbour
		self.curved_edge_info = [] # list of curved edges and their information ready to print in .rea file
		self.nbours= {} # List of actual neighbour cells
		self.info_nbour = {} #each row corresponds to each neighbour, 
		self.rec_info = [] # a list of recieved info
		# info_nbour = [nbour,opposite face or not bool, 2 global verts defining edge on shared face, 
		#................ relative position of midpoint, number of following neighbours to inflict]                   
		self.faces = {}  # A list of Faces
		self.infliction_faces= [] # A list of lists, containing direction of infliction corresponding to each curved edge max 4 directions
		### LOCAL MAPPING FACES ###
		self.face_glloc = {} # Global to local face mapping
		### LOCAL MAPPINGS WITHIN AN ELEMENT DEFINED AS A GLOBAL DICT ### 
		self.n = 2 # The number of cells to inflict
	#-------- Class functions---------#
	#---------------------------------#
	def points_glloc(self,glob):
		'''global to local node mapping'''
		for i in range(20):
			if self.vertices[i+1]== glob:
				return i+1
		return -1
	def points2edge(self,n1,n2):
		''' returns the edge containing the local points n1 and n2'''
		global edge2points
		dummy = [n1,n2]
		dummy.sort()
		for edge in edge2points.keys():
			if (edge2points[edge] == dummy):
				return edge
		return -1
	def addCurve(self,edge,mid,east,west):
            ''' provides the global nodes for the curved edges 
            nodes[:,curved_midpoint[edge]] will give the coordinates to the midpoint'''
            self.curved_edge.append(edge)
            self.curved_midpoint[edge] = mid
            self.curved_east[edge] = east 
            self.curved_west[edge] = west 
	def print_info(self):
		''' Prints the information of a cell'''
		print ('cell nr. {}'.format(self.cellID))
		print ('list of vertices:\n')
		for i in range(21):
			print i , self.vertices[i]
		print 'Curved edges: \n'
		if not self.curved_edge : print 'No registered curved edges'
		else: print self.curved_edge
		print 'Neighbours, corresponding face and global ID: \n'
		for n in range(6):
			print n+1 , self.nbours[n+1]
		print ' neighbours where curvature will be propagated:'
		for n in self.infliction_faces:
			print n
			print '\n'
		print 'info to send to neihgbour:'
		#print self.curved_info
		for key in self.curved_info.keys():
			print 'info : {}'.format(self.curved_info[key])
		print 'global face numbering \n'
		for n in [1,2,3,4,5,6]:
			print n, self.faces[n]
		print 'That was all the information \n'
		print '-------------||----------------'

	def create_inf_faces(self):
		''' initializes the faces where the neighbour that is going to be inflicted lives
		for now it assumes a cylindrical structure of the domain'''
		global edge2faces,shadow_face
		adjFaces = [] # The adjacent faces to the curves
		counter = 0
		for edge in self.curved_edge:
			adjFaces.append(edge2faces[edge][0])
			adjFaces.append(edge2faces[edge][1])
		for i in range(6):
			counter = adjFaces.count(i+1)
			if(counter == 2): # the face to propagate through
                            if(self.nbours[shadow_face[i+1]]>0): # Check that the neightbour is not zero
				self.infliction_faces.append(shadow_face[i+1])

	def glob_nodes_inf(self,nodes):
		''' defines the global nodes which are relevant when inflicting curvature 
		onto another element also assumes a cylindrical structure'''
		inf_edges = []
		inf_offset = []
		for edge in self.curved_edge: # Iterating over the curved egdes
			for i in range(2): # iterating over the parallell edges
				if parallel_edge_map[edge][i] in inf_edges: # if edge already exists the remove it
					inf_edges.remove(parallel_edge_map[edge][i])
					inf_offset.remove(self.calc_midpoint(edge,nodes))
				elif parallel_edge_map[edge][i] not in self.curved_edge: # if edge does not already exist then add it
					inf_edges.append(parallel_edge_map[edge][i])
					inf_offset.append(self.calc_midpoint(edge,nodes))
		for i in range(len(inf_edges)):
			#print 'curved offset for edge {} is equal {}\n'.format(inf_edges[i],inf_offset[i])
			east = self.vertices[edge2points[inf_edges[i]][0]]
			west = self.vertices[edge2points[inf_edges[i]][1]]
			mid  = self.vertices[edge2mid[inf_edges[i]]] 
			self.curved_info[i] = [east,west,mid,inf_offset[i]] # info to be past on
		
	def sendInfo(self,cell2):
		''' sending necessary curved side info to neighobour cell2'''
		for i in self.curved_info.keys():
			cell2.rec_info.append(self.curved_info[i])
                print'sending info from cell {} to cell nr {}'.format(self.cellID,cell2.cellID)

        def updateNodes(self,nodes,info):
            ''' Moving the midside node of a certain edge returning the correct position'''
            x1 = nodes[:,info[0]]
            x2 = nodes[:,info[1]]
            x3 = nodes[:,info[2]]
            a = x2-x1
            k = info[3][0]
            h = info[3][1]
            x3 = k*a+x1 + h*dot(a,a)*array([-a[1],a[0],0])
            return x3
            

        def updateInfo(self,nodes,tot_num_curved):
                 '''updating the info in rec_info so that the curved edges are imposed.'''
                 for info in self.rec_info:
                    #print 'global point {}'.format(info[0])
                         east = self.points_glloc(info[0]) # local east point
                         west = self.points_glloc(info[1]) # local west point
                         mid = self.points_glloc(info[2]) # local mid point
                         #print 'Vertices for this edge : {},{}'.format(east,west)
                         edge = self.points2edge(east,west)
                         if edge == -1: 
                             print 'all global points for element {}:  {}'.format(self.cellID,self.vertices)
                             print 'EDGE NOT FOUND IN points2edge!!'
                             return
                        #print edge
                         self.curved_edge.append(edge)
                         self.curved_midpoint[edge] = info[2]
                         self.curved_east[edge] = info[0]
                         self.curved_west[edge] = info[1]
                         tot_num_curved = tot_num_curved+1
                         nodes[:,info[2]] = self.updateNodes(nodes,info) # Updating the nodes position
                 return nodes,tot_num_curved

	def calc_midpoint(self,edge,nodes):
		'''Calculates the relative position of the given midpoint on edge 
		k is the relative position on edge, and h is the offset from edge.'''
		x1,x2 = self.vertices[edge2points[edge]] # The endpoints of the curved edge
		x1 = x1 -1 # OBS zero and 1-indexing 
		x2 = x2 -1 # OBS zero and 1-indexing 
		x3 = self.curved_midpoint[edge] # The midpoint of the curved edge
		a = nodes[:,x2]-nodes[:,x1]
		b = nodes[:,x3]-nodes[:,x1]
		k = dot(a,b)/dot(a,a) # rel distance in kji direction from corner x1. 
		h = dot(b-k*a,b-k*a)/dot(a,a) # rel distance in eta direction from the line x2-x1
		return k,h/2

	def define_face_mappings(self,Faces):
		'''A function that defines the mappings from global to local faces both ways'''
		if not self.faces: print 'NO FACES REGISTERED! , ERROR: DEFINE_FACE_MAPPING' 
		for count in range(6):
			self.face_glloc[self.faces[count+1]] = count+1
		if count < 5: 
			print ('cell no %d has less than 6 faces' % self.cellID)
			print ('only %d registered faces' % count)

	def opornot(edge,f): 
		'''determining whether an edge belongs to face or the opposite one'''
		count = 0
		for node in edge2points[edge]:
			for node2 in face2points[f]:
				if(node==node2):
					count = count+1
		if(count==2): return 1
		elif(count==0): return 0
		else: print('in element {}, face {} and edge {} has {} common pts'.format(self.cellID,f,edge,count))

	def create_info():
		'''creating the curvature info to be passed on to the next cell'''
		if curved_edge == {}: 
			print 'no curved info to be created'
			return 
		if not infliction_faces: # no faces registered
			for i in curved_edge.keys():
				infliction_faces[i] = edge2faces(curved_edge[i])
		for edge in infliction_faces: #iterating over the edges
			for f in infliction_faces[edge]: # iterating over the faces corresponding to edge
				opposite = opornot(edge,f) 
				pts = self.vertices[edge2points[edge]] #global points
				midpoint = calc_midpoint(curved_edge_info[edge],edge)
				info_nbour.append([f,opposite,pts,midpoint,self.n-1])

	def invoke_midpoint(information):
		'''adding a curvature given the input'''
		# Find the corresponding edge 
			# Find the corresponding face
		face = face_glloc[globface] # The local face
		pts = pts_glloc(points) # The local points
		if(opposite):
			#switch the face and points
			face = opface
			pts = oppts
			# PROBLEM!! faceverts are not consistent !! line 500 nmshconvert.py ask cew
		self.curved_edge.append(points2edge(pts)) # The local edge is now created
			# fix curved_edge_info 
		self.curved_edge_info = transforminfo(information)

	def inflict_midpoint():
		'''inflicting the neighbours with curved edge'''
		for ne in nbours:
			nbours[ne].invoke_midpoint(info_nbour[ne])
			nbours[ne].update_midpoint()
			nbours[ne].inflict_midpoint()

class Face:
    id_index = defaultdict(list)
    nodesum_index=defaultdict(list)
    faceID=0;
    vert_1 = 0
    vert_2 = 0
    vert_3 = 0
    vert_4 = 0
    nodesum = vert_1 + vert_2 + vert_3 + vert_4
    cell_1 = 0
    cell_2 = 0
    zone_id = 0
      
    def addData(self, v1, v2, v3, v4, c1):
	self.vert_1 = v1
	self.vert_2 = v2
	self.vert_3 = v3
	self.vert_4 = v4
	self.cell_1 = c1
        self.nodesum = self.vert_1 + self.vert_2 + self.vert_3 + self.vert_4
        self.verts = sorted([self.vert_1, self.vert_2, self.vert_3, self.vert_4])
#        Face.verts_index[self.verts].append(self)
	Face.nodesum_index[self.nodesum].append(self)
     
    def addCell(self, c2):
	self.cell_2 = c2

    def addZone(self, zid):
	self.zone_id = zid

    def addID(self, ID):
	self.faceID = ID
	Face.id_index[ID].append(self)

    def getID(self):
        return self.faceID

    def isEqual(self,v1,v2,v3,v4):
        return (sorted([v1, v2, v3, v4]) == self.verts)

    @classmethod
    def find_by_ID(cls, ID):
        return Face.id_index[ID]
    @classmethod
    def find_by_nodesum(cls, sum):
        return Face.nodesum_index[sum]   


'''
to be done : 
Making the curved side info, use the global nodes array.
'''
