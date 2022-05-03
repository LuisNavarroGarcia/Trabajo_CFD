import numpy as np

def solver(w0, t0, sim_config, problem, propagator, dt, stop_criteria):
    iteration = 0
    w = w0
    wsol = w0
    t = t0
    tsol = np.array([])
    tsol = np.append(tsol, t)
    iteration += 1
    stop = False
    updated_dt = dt.dt0

    while not stop:
        w = propagator(w, t, updated_dt, problem)
        t += updated_dt

        wsol = np.append(wsol, w, axis = 1)

        stop, criteria = stop_criteria(wsol, t, iteration)
  
        tsol = np.append(tsol, t)
        iteration += 1
        try:
            updated_dt = dt.calc(wsol, updated_dt, dt)
        except Exception:
            updated_dt = dt.dt0

    return wsol, tsol, criteria