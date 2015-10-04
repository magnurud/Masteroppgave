import numpy as np

# Global dictionaries
edge2points = {1:[1,2], 2:[2,3], 3:[3,4], 4:[1,4], 5:[5,6], 6:[6,7], 7:[7,8], 8:[8,9], 9:[1,5], 10:[2,6], 11:[3,7], 12:[4,8]}
face2points = {1:[1,2,6,5], 2:[2,3,7,6], 3:[4,3,7,8], 4:[1,4,8,5], 5:[1,2,3,4], 6:[5,6,7,8]}
edge2faces = {1:[1,5], 2:[2,5], 3:[3,5], 4:[4,5], 5:[1,6], 6:[2,6], 7:[3,6], 8:[4,6], 9:[4,1], 10:[1,2], 11:[2,3], 12:[3,4]}
shadow_face = {1:3, 3:1, 2:4, 4:2, 5:6, 6:5}
Pts = [] # A list with the position vector of each point

class cell:
    cellID= 0 # Id number
    vertices = np.zeros(21) #List of vertices,[#tot,8 vertices, 12 midpoints] All global
    curved_edge = [] # curved edge numbers
    curved_edge_info = [] # list of curved edges and their information ready to print in .rea file
    nbours = np.zeros(6) # Array of neighbours global numbers corresonding to facenumber
    nbours_list = [] # List of actual neighbour cells
    info_nbour = [] #each row corresponds to each neighbour, 
    # info_nbour = [nbour,opposite face or not bool, 2 global verts defining edge on shared face, 
    #................ relative position of midpoint, number of following neighbours to inflict]                   
    Faces = []  # A list of Faces
    infliction_faces= [] # A list of lists, containing direction of infliction corresponding to each curved edge max 4 directions
    ### LOCAL MAPPING FACES ###
    face_glloc = {} # Global to local face mapping
    face_locgl = {} # Local to global face mapping
    points_glloc = {} #
    ### LOCAL MAPPINGS WITHIN AN ELEMENT DEFINED AS A GLOBAL DICT ### 

    n = 2 # The number of cells to inflict
#----------- Class functions -----------#
#----------- --------------- -----------#
    def print_info():
        ''' Prints the information of a cell'''
        print ('cell nr. {}'.format(self.cellID))
        print ('list of vertices:\n')
        for i in range(21):
            print self.vertices[i]
        print 'Curved edges: \n'
        if not curved_edge: print 'No registered curved edges'
        else: print curved_edge
        print 'Neighbours, corresponding face and global ID: \n'
        for n in range(6):
            print n , nbours[n]
        print ' neighbours where curvature will be propagated: \n'
        for n in infliction_faces:
            print n
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
        k = np.dot(a,a)/np.dot(a,b)*sqrt(np.dot(a,a)) # rel distance in kji direction from corner x1. 
        h = np.dot(b-k*a,b-k*a) # rel distance in eta direction from the line x2-x1
        return k,h/2

    def define_face_mappings():
        '''A function that defines the mappings from global to local faces both ways'''
        count = 0
        for globface in self.Faces:
            count = count+1
            self.face_glloc[globface.faceID] = count
            self.face_locgl[count] = globface.faceID
        if count < 6: 
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
            nbours_list[ne].invoke_midpoint(info_nbour[ne])
            nbours_list[ne].update_midpoint()
            nbours_list[ne].inflict_midpoint()

cells = []

'''
to be done : 

'''
