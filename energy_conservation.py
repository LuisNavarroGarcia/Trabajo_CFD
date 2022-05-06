import numpy as np

def energy_conservation(w, t, u, mesh, fluid_prop, convection_integrator, diffusion_integrator, spatial_discretization, bc):
    A, b = spatial_discretization(mesh = mesh, fluid_prop = fluid_prop, bc = bc, u = u, w = w, t = t,
     diffusion_integrator = diffusion_integrator, convection_integrator = convection_integrator)
    
    """
    function energy conservation : 
    
    PARAMETERS
    ----------
    
    w : state vector at a given time
    
    t : time [s]
    
    u : velocity field of the fluid as a function of time [s] and position [m]

    convection_integrator : integration scheme for convection

    diffusion_integrator : integration scheme for diffusion

    bc : class containing boundary condition information

    spatial_discretization : spatial discretization scheme used

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

    dw_dt : time derivative of state vector

    A : NxN Jacobian matrix

    b : Nx1 array with independant terms

    dA_dt : time derivative of the Jacobian matrix

    db_dt : time derivative of the independent term vector

    
    """

    dw_dt = np.add(np.dot(A,w), b)

    dA_dt = np.zeros(np.shape(A))
    db_dt = np.zeros((len(A),1))

    return dw_dt, A, b, dA_dt, db_dt