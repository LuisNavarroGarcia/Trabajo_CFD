import numpy as np

def energy_conservation(w, t, u, mesh, fluid_prop, convection_integrator, diffusion_integrator, spacial_discretization, bc):
    A, b = spacial_discretization(mesh = mesh, fluid_prop = fluid_prop, bc = bc, u = u, w = w, t = t,
     diffusion_integrator = diffusion_integrator, convection_integrator = convection_integrator)
    dw_dt = A*w + b
    dA_dt = np.zeros(np.shape(A))
    db_dt = np.zeros(len(A))

    return dw_dt, A, b, dA_dt, db_dt