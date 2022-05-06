import numpy as np
import matplotlib.pyplot as plt

def solver(w0, t0, sim_config, problem, propagator, dt, stop_criteria, activation_plots):
    
    """
    function solver : integrates a probelm in time, updating
    the state vector until it reaches a certain sstopping condition.
    
    PARAMETERS
    ----------
    
    w0 : initial column state vector

    t0 : initial time [s]

    sim_config : class containing these parameters:
        - courant
        - simulation time [s]

    problem : function depending of the state vector and time,
    returning the time derivative of the state vector, its Jacobian matrix
    its independant terms, the time derivative of the Jacobian matrix for a given time

    propagator : time integrator. It's a function of the state vector,
    time, time step and a function that describes the problem. Returns
    the integrated state vector for a given time step.

    dt : class containing:
        - maximum : maximum simulation time
        - calc : static method to calculate the time step given
            the state vector and a given time
        - dt0 : initial time step
        - courant : courant constant

    stop_criteria : the stopping criteria. The calculation runs
    until this function returns stop = True. It depends on the state
    vector and time

    activation_plots : array containing properties that describe the
    generation of the error convergence plots   


    OUTPUT
    ----------
    
    wsol : matrix containing the state vector for each time iteration

    tsol : array containg the time of each iteration

    criteria : value of the final criteria reached
    """    
    
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

        if stop and activation_plots[0] == 1: 
            plt.ioff()
  
        tsol = np.append(tsol, t)
        iteration += 1
        
        try:
            updated_dt = dt.calc(wsol, updated_dt, dt)
        except Exception:
            updated_dt = dt.dt0

    return wsol, tsol, criteria