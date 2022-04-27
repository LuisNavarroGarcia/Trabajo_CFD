import numpy as np
from fluid_prop import FluidProp
from initial_conditions import InitCond
from propagators import euler_explicit, euler_implicit, euler_pred_corr, runge_kutta, crank_nicolson
from spatial_discretization import spatial_discretization
from conv_upwind_order1 import conv_upwind_order1
from stopping_criterion import stopping_criterion
from mesh import Mesh

fluid_prop = FluidProp() #change to add arguments
    
u = lambda x, y, t: np.full((x, 2), 0.)

num_cells = 4 

bc_type = {
    'upper': 1,
    'lower': 1,
    'left': 1,
    'right': 1
}

bc_handler = {
    'upper': lambda x, y, t : 300,
    'lower': lambda x, y, t : 600,
    'left': lambda x, y, t : 400,
    'right': lambda x, y, t : 1000
}

init_cond = InitCond(t0 = 0, u = u)

propagator = crank_nicolson

type_storage = 1

spatial_discret = lambda mesh, fluid_prop, diffusion_integrator, convection_integrator, bc, u, w, t : spatial_discretization(
    mesh, fluid_prop, bc, u, w, t , type_storage, diffusion_integrator, convection_integrator
)

diffusion_integrator = None
convection_integrator = conv_upwind_order1

dt_calc = None
dt0 = 0.01

v_criteria = None
v_values = None
v_AndOr = None

activation_plots = None

stop_criteria = lambda wsol, t, iteration : stopping_criterion(
    v_criteria, v_values, v_AndOr, activation_plots, wsol, t, iteration
)

sim_config = None

mesh = Mesh(num_cells)
mesh.preprocess()