def neumann_convection(
        bc_handler = None , BC = None, 
    bc = None, conv = None , idx = None , x = None ,
    y = None, t =None, mesh = None, fluid_prop = None
    ):

    """
    function neumann_convection : calculates the convective matrix and
    vector associated to a boundary where flux of heat is imposed 
    (Neumann condition), in finite volumes being a gradient of temperatures.
    
    PARAMETERS
    ----------
    
    bc_handler : vector of values at each boundary 'n'.
    
    BC : NxN matrix with convection terms associated to boundary conditions
    
    bc : Nx1 vector with convection terms associated to boundary conditions

    conv : convection term for cell 'i'

    idx : index of cell where calculation is performed

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
    
    BC_i : term associated to the cell i of the convective matrix of the boundary conditions

    bc_i : term associated to the cell i of the convective vector of the boundary conditions 
    """

    V_i = mesh.V[0]
    q_bc = bc_handler(x, y, t)
    BC_i = BC[idx, idx]
    bc_i = bc[idx] + q_bc/V_i
    w = None
    q = q_bc

    return BC_i, bc_i, w, q

def neumann_diffusion(
        bc_handler = None , BC = None, 
    bc = None, conv = None , idx = None , x = None ,
    y = None, t =None, mesh = None, fluid_prop = None
    ):

    """
    function neumann_diffusion : calculates the convective matrix and
    vector associated to a boundary where flux of heat is imposed 
    (Neumann condition), in finite volumes being a gradient of temperatures.
    
    PARAMETERS
    ----------
    
    bc_handler : vector of values at each boundary 'n'.
    
    BC : NxN matrix with diffusion terms associated to boundary conditions
    
    bc : Nx1 vector with diffusion terms associated to boundary conditions

    conv : diffusion term for cell 'i'

    idx : index of cell where calculation is performed

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
    
    BC_i : term associated to the cell i of the diffusion matrix of the boundary conditions

    bc_i : term associated to the cell i of the diffusion vector of the boundary conditions 
    """

    V_i = mesh.V[0]
    q_bc = bc_handler(x, y, t)
    BC_i = BC[idx, idx]
    bc_i = bc[idx] + q_bc/(V_i * fluid_prop.cv * fluid_prop.rho)
    w = None
    q = q_bc

    return BC_i, bc_i, w, q

def dirichlet_convection(
        bc_handler = None , BC = None, 
    bc = None, conv = None , idx = None , x = None ,
    y = None, t =None, mesh = None, fluid_prop = None
    ):

    """
    function dirichlet_convection : calculates the convective matrix and
    vector associated to a boundary where the temperature of it is imposed 
    (Dirichlet condition).
    
    PARAMETERS
    ----------
    
    bc_handler : vector of values at each boundary 'n'.
    
    BC : NxN matrix with convection terms associated to boundary conditions
    
    bc : Nx1 vector with convection terms associated to boundary conditions

    conv : convection term for cell 'i'

    idx : index of cell where calculation is performed

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
    
    BC_i : term associated to the cell i of the convective matrix of the boundary conditions

    bc_i : term associated to the cell i of the convective vector of the boundary conditions 
    """

    w_bc = bc_handler(x, y, t)
    BC_i = BC[idx, idx] + conv
    bc_i = bc[idx] + conv * w_bc
    w = w_bc
    q = 0

    return BC_i, bc_i, w, q

def dirichlet_diffusion(
        bc_handler = None , BC = None, 
    bc = None, conv = None , idx = None , x = None ,
    y = None, t =None, mesh = None, fluid_prop = None
    ):

    """
    function dirichlet_diffusion : calculates the diffusion matrix and
    vector associated to a boundary where the temperature of it is imposed 
    (Dirichlet condition).
    
    PARAMETERS
    ----------
    
    bc_handler : vector of values at each boundary 'n'.
    
    BC : NxN matrix with diffusion terms associated to boundary conditions
    
    bc : Nx1 vector with diffusion terms associated to boundary conditions

    conv : diffusion term for cell 'i'

    idx : index of cell where calculation is performed

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
    
    BC_i : term associated to the cell i of the diffusion matrix of the boundary conditions

    bc_i : term associated to the cell i of the diffusion vector of the boundary conditions 
    """
    
    w_bc = bc_handler(x, y, t)
    BC_i = BC[idx, idx] - conv
    bc_i = bc[idx] + conv * w_bc
    w = w_bc
    q = 0

    return BC_i, bc_i, w, q

class BC():
    '''
    class BC : stores the properties related to boundaries:
        - bc_handler : array of length n boundaries with the
         value/function that describes it

        - bc_type : array of length n boundaries with the
         type of boundary condition:
         (1): Neumann
         (else): Dirichlet
    '''
    def __init__(self, bc_handler, bc_type):
        self.bc_handler = bc_handler
        self.bc_type = bc_type

