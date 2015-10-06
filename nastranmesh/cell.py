from numpy import zeros,array,dot, sqrt
from collections import defaultdict

class cell:
	def __init__(self,ID,verts):
		self.cellID = ID    # Id number
		self.vertices = verts#List of vertices,[#tot,8 vertices, 12 midpoints] All global
		self.curved_edge = [] # curved edge numbers
		self.curved_edge_info = [] # list of curved edges and their information ready to print in .rea file
		self.nbours= {} # List of actual neighbour cells
		self.info_nbour = {} #each row corresponds to each neighbour, 
		# info_nbour = [nbour,opposite face or not bool, 2 global verts defining edge on shared face, 
		#................ relative position of midpoint, number of following neighbours to inflict]                   
		self.faces = {}  # A list of Faces
		self.infliction_faces= [] # A list of lists, containing direction of infliction corresponding to each curved edge max 4 directions
		### LOCAL MAPPING FACES ###
		self.face_glloc = {} # Global to local face mapping
		self.points_glloc = {} #
		### LOCAL MAPPINGS WITHIN AN ELEMENT DEFINED AS A GLOBAL DICT ### 
		self.n = 2 # The number of cells to inflict
	#-------- Class functions---------#
	#---------------------------------#
	def addCurve(self,edge):
		self.curved_edge.append(edge)
	def print_info(self):
		''' Prints the information of a cell'''
		print ('cell nr. {}'.format(self.cellID))
		print ('list of vertices:\n')
		for i in range(21):
			print self.vertices[i]
		print 'Curved edges: \n'
		if not self.curved_edge : print 'No registered curved edges'
		else: print self.curved_edge
		print 'Neighbours, corresponding face and global ID: \n'
		for n in range(6):
			print n+1 , self.nbours[n+1]
		print ' neighbours where curvature will be propagated: \n'
		for n in self.infliction_faces:
			print n
		print 'global face numbering \n'
		for n in [1,2,3,4,5,6]:
			print n, self.faces[n]
		print 'That was all the information \n'
		print '-------------||----------------'


	def calc_midpoint(info_array,edge):
		'''Calculates the relative position of the given midpoint on edge 
		k is the relative position on edge, and h is the offset from edge.'''
		global Pts
		x1,x2 = self.vertices[edge2points(edge)]
		x3 = info_array[3:6] 
		a = Pts[x2]-Pts[x1]
		b = Pts[x3]-Pts[x1]
		k = dot(a,a)/dot(a,b)*sqrt(dot(a,a)) # rel distance in kji direction from corner x1. 
		h = dot(b-k*a,b-k*a) # rel distance in eta direction from the line x2-x1
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
