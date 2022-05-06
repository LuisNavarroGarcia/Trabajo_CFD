import numpy as np
from boundary_conditions import neumann_convection, dirichlet_convection

def conv_upwind_order1(mesh = None, fluid_prop = None, bc = None, u = None, w = None, t = None):

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
