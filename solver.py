import numpy as np

def solver(w0, t0, sim_config, problem, propagator, dt, stop_criteria):
    iteration = 0
    w = w0
    wsol = w0
    t = t0
    tsol = np.array([])
    tsol[iteration] = t
    iteration += 1
    stop = False
    updated_dt = dt.dt0

    while not stop:
        w = propagator(w, t, updated_dt, problem)
        t = t + updated_dt

        wsol = np.array([wsol, w])

        stop, criteria = stop_criteria(wsol, t[-1], iteration)

        tsol[iteration] = t
        iteration += 1
        updated_dt = dt.calc(wsol, updated_dt, dt)
    
    return wsol, tsol, criteria