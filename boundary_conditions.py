def neumann_convection(
        bc_handler = None , BC = None, 
    bc = None, conv = None , iteration = None , x = None ,
    y = None, t =None, mesh = None, fluid_prop = None
    ):

    V_i = mesh.V[0]
    q_bc = bc_handler(x, y, t)
    BC_i = BC[iteration, iteration]
    bc_i = bc[iteration] + q_bc/V_i
    w = None
    q = q_bc

    return BC_i, bc_i, w, q

def neumann_diffusion(
        bc_handler = None , BC = None, 
    bc = None, conv = None , iteration = None , x = None ,
    y = None, t =None, mesh = None, fluid_prop = None
    ):

    V_i = mesh.V[0]
    q_bc = bc_handler(x, y, t)
    BC_i = BC[iteration, iteration]
    bc_i = bc[iteration] + q_bc/(V_i * fluid_prop.cv * fluid_prop.rho)
    w = None
    q = q_bc

    return BC_i, bc_i, w, q

def dirichlet_convection(
        bc_handler = None , BC = None, 
    bc = None, conv = None , iteration = None , x = None ,
    y = None, t =None, mesh = None, fluid_prop = None
    ):

    w_bc = bc_handler(x, y, t)
    BC_i = BC[iteration, iteration] + conv
    bc_i = bc[iteration] + conv * w_bc
    w = w_bc
    q = 0

    return BC_i, bc_i, w, q

def dirichlet_diffusion(
        bc_handler = None , BC = None, 
    bc = None, conv = None , iteration = None , x = None ,
    y = None, t =None, mesh = None, fluid_prop = None
    ):
    
    w_bc = bc_handler(x, y, t)
    BC_i = BC[iteration, iteration] - conv
    bc_i = bc[iteration] + conv * w_bc
    w = w_bc
    q = 0

    return BC_i, bc_i, w, q

