import numpy as np
from energy_conservation import energy_conservation
from fluid_prop import FluidProp
from initial_conditions import InitCond
from propagators import euler_explicit, euler_implicit, euler_pred_corr, runge_kutta, crank_nicolson
from spatial_discretization import spatial_discretization
from convection import conv_upwind_order1
from diffusion import difusion_cds, difusion_cds_weighted
from stopping_criterion import stopping_criterion
from mesh import Mesh
from timestep import dt_adaptative, dt_constant, courant, DT
import time
from solver import solver
from sim_config import SimConfig
from boundary_conditions import BC
import argparse
from plots import contour_plot

# parser = argparse.ArgumentParser(description = None)
# parser.add_argument('--num_cells', type = int, required = True, help = '')
# parser.add_argument('--num_cells', type = int, required = True, help = '')
# args = parser.parse_args()

k = 0.025 # Thermal Conductivity [W/(m·K)]
cv = 717.5  # Specific Heat [J/(kg·K)]
rho = 1.225 # Density [kg/m^3]

fluid_prop = FluidProp(k = k, cv = cv, rho = rho)
    
u = lambda x, y, t: np.full((len(x), 2), np.array([0., 0.075]))

num_cells = 4

bc_type = np.array([2, 2, 2, 2])

bc_handler = np.array([
    lambda x, y, t : 300,
    lambda x, y, t : 600,
    lambda x, y, t : 400,
    lambda x, y, t : 1000
])

bc = BC(bc_type = bc_type, bc_handler = bc_handler)

init_cond = InitCond(t0 = 0, u = u)

propagator = euler_explicit

type_storage = 0

spatial_discret = lambda mesh, fluid_prop, diffusion_integrator, convection_integrator, bc, u, w, t : spatial_discretization(
    mesh = mesh, fluid_prop = fluid_prop, bc = bc, u = u, w = w, t = t , type_storage = type_storage,
     diffusion_integrator = diffusion_integrator, convection_integrator = convection_integrator
)

diffusion_integrator = difusion_cds
convection_integrator = conv_upwind_order1

dt_calc = dt_constant
dt0 = 0.01

v_criteria = np.array([0, 0, 0, 1])
v_values = np.array([10, 100, 0.007, 2, 0.005, 5])
v_AndOr = np.array([0, 0, 0, 0])

activation_plots = np.array([0, 0, 5])

stop_criteria = lambda wsol, t, iteration : stopping_criterion(
    v_criteria, v_values, v_AndOr, activation_plots, wsol, t, iteration
)

sim_config = SimConfig(courant = 10, t_final = v_values[0])

mesh = Mesh(num_cells)
mesh.preprocess()

if dt_calc == dt_constant:
    dt_calc = dt_constant(dt0)

problem = lambda w, t: energy_conservation(
    w = w, t = t, u = u, mesh = mesh, fluid_prop = fluid_prop,
    convection_integrator = convection_integrator, diffusion_integrator = diffusion_integrator, bc = bc,
    spatial_discretization = spatial_discret
)

dt_courant = courant(mesh, u, sim_config)
dt = DT(maximum= sim_config.tfinal, dt_calc = dt_calc, dt0= dt0, courant = dt_courant)

w0 = np.transpose(init_cond.T(mesh.Rc[:, 0], mesh.Rc[:, 1]))

start = time.time()

w, t, criteria = solver(
    w0 = w0, t0 = init_cond.t0, sim_config = sim_config, problem = problem,
    propagator = propagator, dt = dt, stop_criteria = stop_criteria 
)

end = time.time()

simulation_time = end - start

print(f'The simulation time has been: {simulation_time}')

contour_plot(w = w[:, -1], mesh = mesh, num_map = 2)