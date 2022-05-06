import numpy as np
from boundary_conditions import neumann_convection, dirichlet_convection

def conv_upwind_order1(mesh = None, fluid_prop = None, bc = None, u = None, w = None, t = None):

    """
    function conv_upwind_order1 : This function calculates the matrix convection of the problem using upwind order 1 method. 
    If fluid flow is entering a cell, the temperature value of that cell is the cell as the one
    where flow is leaving.
    In case that there's a boundary condition, the function distinguises which bc is (upper, lower, right, left), and
    which condition is in it (Neumann or Dirichlet).

    PARAMETERS
    ----------
    bc_type : vector that constains the type of boundary condition for each frontier.
    bc_handler : vector of values at each boundary 'n'.
    
    BC : NxN matrix with convection terms associated to boundary conditions
    
    bc : Nx1 vector with convection terms associated to boundary conditions
    conv : convection term for cell 'i' and face/point j
    i : index of cell where calculation is performed
    j : index of the point where calculation is performed
    x : position in x axis [m]
    y : position in y axis [m]
    t : time [s]
    mesh : class with all the mesh properties:
        - Rn : coordinates of nodes forming the mesh
        - Cells : index of nodes forming each cell
        - Volumes : area (2D mesh) of each cell
        - Neighbours : neighboring cells for each cell and/or boundaries
        - Normals : external normal vector for each face of each cell
        - Areas : length of each face of each cell
        - Faces : coordinates of the center of each face for each cell
        
    fluid_prop : class containing the properties of the fluid.
        - k : thermal conductivity [W/m*K]
        - rho : density [kg/m^3]
        - cv : specific heat [J/kg*K]
        
    OUTPUT
    ----------
    C_matrix: N x N matrix, convection matrix where each row corresponds to a cell.
    BC_i : N x N matrix, matrix associated to the convection matrix of the boundary conditions
    bc_i : N x 1 matrix, Nx1 matrix with independent terms
    """

    N = len(mesh.cells)
    C_matrix = np.zeros((N, N))
    bc_vec = np.zeros((N,1))
    BC = np.zeros((N,N))


    for i in range(N):
        V_i = mesh.V[i]
        Centroids = np.array([mesh.Rc[i, 0], mesh.Rc[i, 1]])
        v = u([Centroids[0]], [Centroids[1]], t)
        for j in range(3):
            k = mesh.neighbours[i, j]
            A_ij = mesh.areas[i, j]
            n_ij = np.array([mesh.normals[i, j, 0], mesh.normals[i, j, 1]])
            vn = np.dot(v, n_ij)
            
            conv = np.dot(-A_ij, vn)
            
            if k >= 0:
                if vn < 0:
                    V_k = mesh.V[k]
                    C_matrix[i, k] = conv/V_i
                    C_matrix[k, k] = C_matrix[k, k]-conv/V_k   
            else:
                conv = -A_ij*vn/V_i
                centroid_face = [mesh.faces[i, j, 0], mesh.faces[i, j, 1]]
                if bc.bc_type[np.absolute(k)-1] == 1:
                    [BC_i, bc_i, _, _] = neumann_convection(
                        bc_handler = bc.bc_handler[np.absolute(k)-1] , BC = BC,
                    bc = bc_vec , conv = conv, idx = i, x = centroid_face[0], y = centroid_face[1], t = t,
                    mesh = mesh, fluid_prop = fluid_prop
                    )
                else:
                    [BC_i, bc_i, _, _] = dirichlet_convection(
                        bc_handler = bc.bc_handler[np.absolute(k)-1] , BC = BC,
                    bc = bc_vec , conv = conv, idx = i, x = centroid_face[0], y = centroid_face[1], t = t,
                    mesh = mesh, fluid_prop = fluid_prop
                    )
   
                if vn < 0:
                    bc_vec[i]=bc_i
                else:
                    BC[i,i]=BC_i
   
    return C_matrix, BC, bc_vec

def conv_cds(mesh = None, fluid_prop = None, bc = None, u = None, w = None, t = None):

    
    """
    function conv_cds : this function calculates the matrix convection of the problem using Central Differencing Scheme (CDS). 
    To calculate the temperature of a cell's face, the mean of the temperatures between the adjoining 
    cells is done, so C_matrix is created.
    In case that there's a boundary condition, the function distinguises which bc is (upper, lower, right, left), and
    which condition is in it (Neumann or Dirichlet).

    PARAMETERS
    ----------
    bc_type : vector that constains the type of boundary condition for each frontier.
    bc_handler : vector of values at each boundary 'n'.
    
    BC : NxN matrix with convection terms associated to boundary conditions
    
    bc : Nx1 vector with convection terms associated to boundary conditions
    conv : convection term for cell 'i' and face/point j
    i : index of cell where calculation is performed
    j : index of the point where calculation is performed
    x : position in x axis [m]
    y : position in y axis [m]
    t : time [s]
    mesh : class with all the mesh properties:
        - Rn : coordinates of nodes forming the mesh
        - Cells : index of nodes forming each cell
        - Volumes : area (2D mesh) of each cell
        - Neighbours : neighboring cells for each cell and/or boundaries
        - Normals : external normal vector for each face of each cell
        - Areas : length of each face of each cell
        - Faces : coordinates of the center of each face for each cell
        
    fluid_prop : class containing the properties of the fluid.
        - k : thermal conductivity [W/m*K]
        - rho : density [kg/m^3]
        - cv : specific heat [J/kg*K]
        
    OUTPUT
    ----------
    C_matrix: N x N matrix, convection matrix where each row corresponds to a cell.
    BC_i : N x N matrix, matrix associated to the convection matrix of the boundary conditions
    bc_i : N x 1 matrix, Nx1 matrix with independent terms
    """

    N = len(mesh.cells)
    C_matrix = np.zeros((N, N))
    bc_vec = np.zeros((N,1))
    BC = np.zeros((N,N))


    for i in range(N):
        V_i = mesh.V[i]
        Centroids = np.array([mesh.Rc[i, 0], mesh.Rc[i, 1]])
        v = u([Centroids[0]], [Centroids[1]], t)
        for j in range(3):
            k = mesh.neighbours[i, j]
            A_ij = mesh.areas[i, j]
            n_ij = np.array([mesh.normals[i, j, 0], mesh.normals[i, j, 1]])
            vn = np.dot(v, n_ij)
                        
            if k >= 0:
                conv = np.dot(-A_ij, vn)/2/V_i
                C_matrix[i, i] = C_matrix[i, i] + conv
                C_matrix[i, k] = conv 
            else:
                conv = -A_ij*vn/V_i
                centroid_face = [mesh.faces[i, j, 0], mesh.faces[i, j, 1]]
                if bc.bc_type[np.absolute(k)-1] == 1:
                    [BC_i, bc_i, _, _] = neumann_convection(
                        bc_handler = bc.bc_handler[np.absolute(k)-1] , BC = BC,
                    bc = bc_vec , conv = conv, idx = i, x = centroid_face[0], y = centroid_face[1], t = t,
                    mesh = mesh, fluid_prop = fluid_prop
                    )
                else:
                    [BC_i, bc_i, _, _] = dirichlet_convection(
                        bc_handler = bc.bc_handler[np.absolute(k)-1] , BC = BC,
                    bc = bc_vec , conv = conv, idx = i, x = centroid_face[0], y = centroid_face[1], t = t,
                    mesh = mesh, fluid_prop = fluid_prop
                    )
   
                if vn < 0:
                    bc_vec[i]=bc_i
                else:
                    BC[i,i]=BC_i
   
    return C_matrix, BC, bc_vec
