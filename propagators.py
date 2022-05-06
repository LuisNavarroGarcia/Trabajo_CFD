import numpy as np

def euler_explicit(w, t, dt, problem):

    """
    function euler_explicit : function that integrates a time step with Euler explicit
    method.
    
    PARAMETERS
    ----------
    x : position in x axis [m]
    y : position in y axis [m]
    t : time [s]
    dt: time step [s]
    w: state vector [T]
    problem: handler of the energy conservation function dependent on w and t
        
    OUTPUT
    ----------
    w_new: state vector evaluated at t+dt. Nx1 vector
    """
    dw_dt, _, _, _, _ = problem(w, t) # Returns the first output of problem()
    
    w_new = np.dot(dw_dt, dt) + w
    
    return w_new

def euler_implicit(w, t, dt, problem):

    """
    function euler_implicit : function that integrates a time step with Euler implicit
    method.
    
    PARAMETERS
    ----------
    x : position in x axis [m]
    y : position in y axis [m]
    t : time [s]
    dt: time step [s]
    w: state vector [T]
    problem: handler of the energy conservation function dependent on w and t
        
    OUTPUT
    ----------
    w_new: state vector evaluated at t+dt. Nx1 vector
    """
    dw_dt, A, _, _, _ = problem(w, t) # Returns the first and second outputs of problem()
    
    # Solves the linear system of equations
    w_new = np.linalg.solve(np.identity(len(A)) - np.dot(dt, A), \
      w + np.dot(dt, dw_dt - np.dot(A, w)))

    return w_new

def euler_pred_corr(w, t, dt, problem):

    """
    function euler_pred_corr : function that integrates time step with Second order Corrector Predictor
    method. With it, first is calculated a prediction with euler explicit and corrected with euler implicit
    
    PARAMETERS
    ----------
    x : position in x axis [m]
    y : position in y axis [m]
    t : time [s]
    dt: time step [s]
    w: state vector [T]
    problem: handler of the energy conservation function dependent on w and t
        
    OUTPUT
    ----------
    w_new: state vector evaluated at t+dt. Nx1 vector
    """
    # Actual instant n
    dw_dt, _, _, _, _ = problem(w, t) # Returns the first output of problem()
    
    # Prediction with Euler explicit
    pred_w_1 = np.dot(dw_dt, dt) + w
    pred_dw_dt_1, _, _, _, _ = problem(pred_w_1, t + dt) # Returns the first output of problem()

    # Correction with Euler implicit
    corr_w_1 = np.dot(pred_dw_dt_1, dt) + w

    w_new = corr_w_1

    return w_new

def crank_nicolson(w, t, dt, problem):

    """
    function euler_pred_corr : function that integrates time step with Crank Nicolson method.
        
    PARAMETERS
    ----------
    x : position in x axis [m]
    y : position in y axis [m]
    t : time [s]
    dt: time step [s]
    w: state vector [T]
    problem: handler of the energy conservation function dependent on w and t
        
    OUTPUT
    ----------
    w_new: state vector evaluated at t+dt. Nx1 vector
    """
    dw_dt, A, _, _, _  = problem(w, t) # Returns the first and second outputs of problem()

    # Solves the linear system of equations
    w_new = np.linalg.solve(np.identity(len(A)) - np.dot(dt, A / 2), \
      w + np.dot(dt, (2 * dw_dt - np.dot(A, w)) / 2))

    return w_new

def runge_kutta4(w, t, dt, problem):

    """
    function euler_pred_corr : function that integrates time step with Runge Kutta method. 
    First a prediction is calculated at t+dt/2 with euler explicit. Then, it is corrected with
    Euler implicit. After that, a prediction is done with Euler explicit at t+dt and the result is corrected
    with Simpson's rule.
        
    PARAMETERS
    ----------
    x : position in x axis [m]
    y : position in y axis [m]
    t : time [s]
    dt: time step [s]
    w: state vector [T]
    problem: handler of the energy conservation function dependent on w and t
        
    OUTPUT
    ----------
    w_new: state vector evaluated at t+dt. Nx1 vector
    """
    # At t
    dw_dt, _, _, _, _ = problem(w, t) # Returns the first output of problem()

    # At t + dt / 2 . Prediction with Euler explicit
    pred_w_1_2 = np.dot(dw_dt, dt /2) + w
    pred_dw_dt_1_2_pred, _, _, _, _ = problem(pred_w_1_2, t + dt / 2) # Returns the first output of problem()

    # At t + dt / 2 . Correctoin with Euler implicit
    corr_w_1_2 = np.dot(pred_dw_dt_1_2_pred, dt /2) + w
    corr_dw_dt_1_2, _, _, _, _ = problem(corr_w_1_2, t + dt / 2) # Returns the first output of problem()

    # At t + dt . Prediction with Euler implicit
    pred_w_1 = np.dot(corr_dw_dt_1_2, dt) + w
    pred_dw_dt_1, _, _, _, _ = problem(pred_w_1, t + dt) # Returns the first output of problem()

    # Correction with Simpson's rule
    corr_w_1 = w + np.dot(dt / 6, pred_dw_dt_1 + 2 * pred_dw_dt_1_2_pred + 2 * corr_dw_dt_1_2 + dw_dt)

    w_new = corr_w_1

    return w_new