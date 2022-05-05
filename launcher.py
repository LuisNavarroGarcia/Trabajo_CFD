import numpy as np
from energy_conservation import energy_conservation
from fluid_prop import FluidProp
from initial_conditions import InitCond
from propagators import euler_explicit, euler_implicit, euler_pred_corr, runge_kutta4, crank_nicolson
from spatial_discretization import spatial_discretization
from convection import conv_upwind_order1, conv_cds
from diffusion import difusion_cds, difusion_cds_weighted
from stopping_criterion import stopping_criterion
from mesh import Mesh
from timestep import dt_adaptative, dt_constant, courant, DT
import time
from solver import solver
from sim_config import SimConfig
from boundary_conditions import BC
from plots import contour_plot
from plots import ErrorPlot
import argparse
import configparser

parser = argparse.ArgumentParser(description = 'Launch calculations given a configuration file for the energy equation solver')
parser.add_argument('--config_file', type = str, default = 'config_default.txt',
    help = 'Specify configuration file to run the simulation (.txt extension)')
args = parser.parse_args()

config = configparser.ConfigParser()
config.read(args.config_file)

k = np.float32(config['fluid properties']['k']) # Thermal Conductivity [W/(m·K)]
cv = np.float32(config['fluid properties']['cv'])  # Specific Heat [J/(kg·K)]
rho = np.float32(config['fluid properties']['rho']) # Density [kg/m^3]

fluid_prop = FluidProp(k = k, cv = cv, rho = rho)

u = eval(config['initial conditions']['u'])

num_cells = np.int32(config['mesh']['num_cells'])
num_bc = np.int32(config['mesh']['num_bc'])

bc_type = np.array([], dtype = np.int8)
for value in config['type boundary conditions'].values():
    if value == 'neumann':
        value = 1
    elif value == 'dirichlet':
        value = 2
    else:
        print('Type of boundary not recognized')
        break
        
    bc_type = np.append(bc_type, np.int8(value))

if len(bc_type) != num_bc:
    raise Exception("Number of boundary conditions and number of types of boundary conditions doesn\'t match")

bc_handler = np.array([])
for value in config['value boundary conditions'].values():
    bc_handler = np.append(bc_handler, eval(value))

if len(bc_handler) != num_bc:
    raise Exception("Number of boundary conditions and number of values of boundary conditions doesn\'t match")

bc = BC(bc_type = bc_type, bc_handler = bc_handler)

init_cond = InitCond(t0 = np.float32(config['initial conditions']['t0']), u = u)

propagator = eval(config['propagators']['propagator'])

spatial_discret = lambda mesh, fluid_prop, diffusion_integrator, convection_integrator, bc, u, w, t : spatial_discretization(
    mesh = mesh, fluid_prop = fluid_prop, bc = bc, u = u, w = w, t = t , diffusion_integrator = diffusion_integrator,
    convection_integrator = convection_integrator
)

diffusion_integrator = eval(config['integrators']['diffusion_integrator'])
convection_integrator = eval(config['integrators']['convection_integrator'])

dt_calc =  eval(config['time configuration']['dt_calc'])
dt0 = np.float32(config['time configuration']['dt0'])

v_criteria = np.array([
    np.int32(config['stopping criteria']['stop_at_time_flag']),
    np.int32(config['stopping criteria']['stop_at_iteration_flag']),
    np.int32(config['stopping criteria']['stop_at_max_error_flag']),
    np.int32(config['stopping criteria']['stop_at_mean_error_flag'])
    ])
v_values = np.array([
    np.float32(config['stopping criteria']['stop_at_time']), 
    np.float32(config['stopping criteria']['stop_at_iteration']), 
    np.float64(config['stopping criteria']['stop_at_max_error']), 
    np.float32(config['stopping criteria']['evaluate_max_error_at_n_last_iterations']), 
    np.float64(config['stopping criteria']['stop_at_mean_error']), 
    np.float32(config['stopping criteria']['evaluate_mean_error_at_n_last_iterations'])
    ])
v_AndOr = np.array([
    np.int8(config['stopping criteria']['OR_AND_condition_for_time_condition']), 
    np.int8(config['stopping criteria']['OR_AND_condition_for_iteration_condition']), 
    np.int8(config['stopping criteria']['OR_AND_condition_for_max_error_condition']), 
    np.int8(config['stopping criteria']['OR_AND_condition_for_mean_error_condition'])
    ])

activation_plots = np.array([
    np.int8(config['plot configuration']['active_convergence_plot']),
    np.int8(config['plot configuration']['plot_max_error']), 
    np.int8(config['plot configuration']['plot_mean_error']), 
    np.int8(config['plot configuration']['plot_frecuency'])
    ])

if activation_plots[0] == 1:
    error_plot = ErrorPlot()
else:
    error_plot = None

stop_criteria = lambda wsol, t, iteration : stopping_criterion(
    v_criteria, v_values, v_AndOr, activation_plots, wsol, t, iteration, error_plot
)

sim_config = SimConfig(courant = np.float32(config['time configuration']['courant']), t_final = v_values[0])

mesh = Mesh(num_cells, num_bc)
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
    propagator = propagator, dt = dt, stop_criteria = stop_criteria,
    activation_plots = activation_plots
)

end = time.time()

simulation_time = end - start

print(f'The simulation time has been: {simulation_time}')

colourmap = config['plot configuration']['contour_colormap'] # Inferno theme for contout plot
interp_method = config['plot configuration']['contour_interpolation_method'] #Cubic interpolation method
contour_plot(w = w[:, -1], mesh = mesh, colourmap = colourmap, interp_method = interp_method)