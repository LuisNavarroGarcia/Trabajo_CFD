import numpy as np
from boundary_conditions import neumann_diffusion, dirichlet_diffusion

def difusion_cds(mesh = None, fluid_prop = None, bc = None, u = None, w = None, t = None):

     """
    function difusion_cds : This function calculates the difusion matrix of the problem using Central Differencing Scheme (CDS). 
    In this method, the temperature of a face is calculated as the mean between the tempeartures of the adjoining cells.
    
    PARAMETERS
    ----------
    bc_type : vector that constains the type of boundary condition for each frontier.
    bc_handler : vector of values at each boundary 'n'.
    
    BC : NxN matrix with convection terms associated to boundary conditions
    
    bc : Nx1 vector with convection terms associated to boundary conditions
    cond : difusion term for cell 'i' and face/point j
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
    C_matrix: N x N matrix, difusion matrix where each row corresponds to a cell.
    BC_i : N x N matrix, matrix associated to the difusion matrix of the boundary conditions
    bc_i : N x 1 matrix, Nx1 matrix with independent terms
    """

     num_cells = len(mesh.cells)
     K = np.zeros((num_cells, num_cells))
     BC = np.zeros((num_cells, num_cells))
     bc_vec = np.zeros((num_cells,1))

     for i in range(num_cells):

          for j in range(3):
               n_ij = np.array([mesh.normals[i, j, 0], mesh.normals[i, j, 1]])
               A_ij = mesh.areas[i, j]
               V_i = mesh.V[i]

               if mesh.neighbours[i, j] >= 0:
               
                    neighbour_index = mesh.neighbours[i, j]
                    r_i = np.array([mesh.Rc[i, 0], mesh.Rc[i, 1]])
                    r_j = np.array([mesh.Rc[neighbour_index, 0], mesh.Rc[neighbour_index, 1]])
                    d_ij = np.subtract(r_i, r_j)
                    d_norm = np.linalg.norm(d_ij)
                    dn = np.dot(d_ij, n_ij)
                    cond = -A_ij / V_i * fluid_prop.k * (dn / (d_norm * d_norm)) \
                                   / (fluid_prop.cv * fluid_prop.rho)
                    K[i, neighbour_index] = cond
                    K[i, i] = K[i, i] - cond

               else:

                    r_i = np.array([mesh.Rc[i, 0], mesh.Rc[i, 1]])
                    r_face = np.array([mesh.faces[i, j, 0], mesh.faces[i, j, 1]])
                    n_ij = np.array([mesh.normals[i, j, 0], mesh.normals[i, j, 1]])
                    d_ic = np.subtract(r_i, r_face)
                    d_norm = np.linalg.norm(d_ic)
                    dn = np.dot(d_ic, n_ij)
                    cond = -A_ij / V_i * fluid_prop.k * (dn / (d_norm * d_norm)) \
                                   / (fluid_prop.cv * fluid_prop.rho)                  

                    if bc.bc_type[np.absolute(mesh.neighbours[i, j])-1] == 1:
                         [BC_i, bc_i, _, _] = neumann_diffusion(
                         bc_handler = bc.bc_handler[np.absolute(mesh.neighbours[i, j])-1] , BC = BC,
                         bc = bc_vec , conv = cond, idx = i, x = r_face[0], y = r_face[1], t = t,
                         mesh = mesh, fluid_prop = fluid_prop
                         )
                    else:
                         [BC_i, bc_i, _, _] = dirichlet_diffusion(
                         bc_handler = bc.bc_handler[np.absolute(mesh.neighbours[i, j])-1] , BC = BC,
                         bc = bc_vec , conv = cond, idx = i, x = r_face[0], y = r_face[1], t = t,
                         mesh = mesh, fluid_prop = fluid_prop
                         )

                    BC[i, i] = BC_i
                    bc_vec[i] = bc_i

     return K, BC, bc_vec

def difusion_cds_weighted(mesh = None, fluid_prop = None, bc = None, u = None, w = None, t = None):

     """
    function difusion_cds_weighted : This function calculates the difusion matrix of the problem.
    The temperature of the cell varies linearly between the centroid of the cell and another, 
    but also depends on the relative distance to the face.
    
    PARAMETERS
    ----------
    bc_type : vector that constains the type of boundary condition for each frontier.
    bc_handler : vector of values at each boundary 'n'.
    
    BC : NxN matrix with convection terms associated to boundary conditions
    
    bc : Nx1 vector with convection terms associated to boundary conditions
    cond : difusion term for cell 'i' and face/point j
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
    C_matrix: N x N matrix, difusion matrix where each row corresponds to a cell.
    BC_i : N x N matrix, matrix associated to the difusion matrix of the boundary conditions
    bc_i : N x 1 matrix, Nx1 matrix with independent terms
    """

     num_cells = len(mesh.cells)
     K = np.zeros((num_cells, num_cells))
     BC = np.zeros((num_cells, num_cells))
     bc_vec = np.zeros((num_cells, 1))

     for i in range(num_cells):

          for j in range(3):
               n_ij = np.array([mesh.normals[i, j, 0], mesh.normals[i, j, 1]])
               A_ij = mesh.areas[i, j]
               V_i = mesh.V[i]

               if mesh.neighbours[i, j] >= 0:
               
                    neighbour_index = mesh.neighbours[i, j]
                    r_i = np.array([mesh.Rc[i, 0], mesh.Rc[i, 1]])
                    r_j = np.array([mesh.Rc[neighbour_index, 0], mesh.Rc[neighbour_index, 1]])
                    r_face = np.array([mesh.faces[i, j, 0], mesh.faces[i, j, 1]])
                    d_ic = np.linalg.norm(r_face - r_i)
                    d_jc = np.linalg.norm(r_face - r_j)
                    d_tot = d_ic + d_jc
                    f_i = (d_tot - d_ic) / d_tot
                    f_j = (d_tot - d_jc) / d_tot
                    d_ij = np.subtract(r_i, r_j) 
                    d_norm = np.linalg.norm(d_ij)
                    dn = np.dot(d_ij, n_ij)
                    cond = -A_ij / V_i * fluid_prop.k * (dn / (d_norm * d_norm)) \
                              / (fluid_prop.cv * fluid_prop.rho)
                    K[i, neighbour_index] = f_j * cond
                    K[i, i] = K[i, i] - f_i * cond

               else:

                    r_i = np.array([mesh.Rc[i, 0], mesh.Rc[i, 1]])
                    r_face = np.array([mesh.faces[i, j, 0], mesh.faces[i, j, 1]])
                    n_ij = np.array([mesh.normals[i, j, 0], mesh.normals[i, j, 1]])
                    d_ia = r_i - r_face
                    d_norm = np.linalg.norm(d_ia)
                    dn = np.dot(d_ia, n_ij)
                    cond = -A_ij / V_i * fluid_prop.k * (dn / (d_norm * d_norm)) \
                              / (fluid_prop.cv * fluid_prop.rho)  

                    if bc.bc_type[np.absolute(mesh.neighbours[i, j])-1] == 1:
                              [BC_i, bc_i, _, _] = neumann_diffusion(
                              bc_handler = bc.bc_handler[np.absolute(mesh.neighbours[i, j])-1] , BC = BC,
                              bc = bc_vec , conv = cond, idx = i, x = r_face[0], y = r_face[1], t = t,
                              mesh = mesh, fluid_prop = fluid_prop
                              )
                    else:
                         [BC_i, bc_i, _, _] = dirichlet_diffusion(
                         bc_handler = bc.bc_handler[np.absolute(mesh.neighbours[i, j])-1] , BC = BC,
                         bc = bc_vec , conv = cond, idx = i, x = r_face[0], y = r_face[1], t = t,
                         mesh = mesh, fluid_prop = fluid_prop
                         )

                    BC[i, i] = BC_i
                    bc_vec[i] = bc_i
                    
     return K, BC, bc_vec   