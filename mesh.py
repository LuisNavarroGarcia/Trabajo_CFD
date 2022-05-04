import numpy as np

class Mesh():
    def __init__(self, num_cells, num_bc):
        
        #load cell data
        self.cells = np.loadtxt(f'mesh/cells_{num_cells}.dat', dtype= np.int32)
        # load node data
        self.Rn = np.loadtxt(f'mesh/nodes_{num_cells}.dat')
        # load boundary conditions

        self.bc = []
        for i in range(num_bc):
            try:
                self.bc.append(np.loadtxt(f'mesh/bc_{i+1}_{num_cells}.dat', dtype= np.int32))
            except Exception:
                print(f'Only found: {i} boundaries')

        
    def preprocess(self):
        
        # initialize arrays
        self.V = np.zeros(len(self.cells))
        self.neighbours = np.zeros((len(self.cells), 3), dtype= np.int32)
        self.areas = np.zeros((len(self.cells), 3))
        self.normals = np.zeros((len(self.cells), 3, 2))
        self.Rc = np.zeros((len(self.cells), 2))
        self.faces = np.zeros((len(self.cells), 3, 2))

        # calculations
        for i in range(len(self.cells)):
          
            # Get indexes of points that define each cell
            indexes = self.cells[i] 

            # Get point that form a cell
            point1 = self.Rn[indexes[0]-1] 
            point2 = self.Rn[indexes[1]-1] 
            point3 = self.Rn[indexes[2]-1]

            # Calculate sides
            delta12 = np.subtract(point2,point1)
            delta23 = np.subtract(point3,point2)
            delta31 = np.subtract(point1,point3)

            # calculate faces
            face1 = np.add(point1, delta12/2)
            face2 = np.add(point2, delta23/2)
            face3 = np.add(point3, delta31/2)

            # construct array of faces
            self.faces[i] = np.array([face1, face2, face3])

            # calculate length of sides (areas)
            area1 = np.linalg.norm(delta12)
            area2 = np.linalg.norm(delta23)
            area3 = np.linalg.norm(delta31)

            # construct array of areas
            self.areas[i] = np.array([area1, area2, area3])

            # construct array of volumes
            self.V[i] = np.array((1/2)*np.linalg.norm(np.cross(delta12, delta23)))

             # construct array of centroids
            self.Rc[i] =  np.array([(point1[0]+point2[0]+point3[0])/3,
                        (point1[1]+point2[1]+point3[1])/3])
            
            # calculate external normal vectors

            #normal for first face
            dx1 = delta12[0]
            dy1 = delta12[1]
            N1 = np.array([-dy1, dx1])/area1
            
            # check if it's not external, and recalculate if it's not
            if np.dot(np.subtract(self.Rc[i],face1), N1) > 0:
                N1 = np.array([dy1, -dx1])/area1
            
            #normal for second face
            dx2 = delta23[0]
            dy2 = delta23[1]
            N2 = np.array([-dy2, dx2])/area2
            
            # check if it's not external, and recalculate if it's not
            if np.dot(np.subtract(self.Rc[i],face2), N2) > 0:
                N2 = np.array([dy2, -dx2])/area2
            
            #normal for third face
            dx3 = delta31[0]
            dy3 = delta31[1]         
            N3 = np.array([-dy3, dx3])/area3
            
            # check if it's not external, and recalculate if it's not
            if np.dot(np.subtract(self.Rc[i],face3), N3) > 0:
                N3 = np.array([dy3, -dx3])/area3

            # construct array of normals
            self.normals[i] = np.array([N1, N2, N3])

            # temporal array for neighbour storage
            neighbours = np.full(3, None)
            
            # loop to check if the other cells are neighbours
            for j in range(len(self.cells)):
                
                # only check for cells that are not the one being preprocessed
                if j != i:
                    
                    #initialize temporal variable of node coincidences
                    num_coincidences = 0
                    nodes = np.array([])
                    # check how many nodes are the same
                    for node in self.cells[j]:
            
                        if node in self.cells[i]:
                            num_coincidences += 1
                            nodes = np.append(nodes, node)

                    if num_coincidences >= 2:
                        idx1 = np.where(self.cells[i] == nodes[0])[0]
                        idx2 = np.where(self.cells[i] == nodes[1])[0]
                        if (idx1 == 0 or idx1 == 1)  and (idx2 == 0 or idx2 == 1):
                            neighbours[0] = j
                        elif (idx1 == 2 or idx1 == 1)  and (idx2 == 2 or idx2 == 1):
                            neighbours[1] = j
                        elif (idx1 == 0 or idx1 == 2)  and (idx2 == 0 or idx2 == 2):
                            neighbours[2] = j   

            # loop to check if there are boundary conditions and where are located
            for side in range(len(self.bc)):
                num_coincidences = 0
                bc_nodes = np.array([])

                for node in self.bc[side]:
                    if node in self.cells[i]:
                        num_coincidences += 1
                        bc_nodes = np.append(bc_nodes, node)

                if num_coincidences >= 2:
                    idx1 = np.where(self.cells[i] == bc_nodes[0])[0]
                    idx2 = np.where(self.cells[i] == bc_nodes[1])[0]

                    if (idx1 == 0 or idx1 == 1)  and (idx2 == 0 or idx2 == 1):
                        neighbours[0] = -(side+1)
                    elif (idx1 == 2 or idx1 == 1)  and (idx2 == 2 or idx2 == 1):
                        neighbours[1] = -(side+1)
                    elif (idx1 == 0 or idx1 == 2)  and (idx2 == 0 or idx2 == 2):
                        neighbours[2] = -(side+1)

            # assign completed temporal array of neighbours to the class property
            self.neighbours[i] = neighbours