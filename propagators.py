import numpy as np

def euler_explicit(w, t, dt, problem):

    dw_dt, _, _, _, _ = problem(w, t) # Returns the first output of problem()
    
    w_new = np.dot(dw_dt, dt) + w
    
    return w_new

def euler_implicit(w, t, dt, problem):

    dw_dt, A, _, _, _ = problem(w, t) # Returns the first and second outputs of problem()
    
    # Solves the linear system of equations
    w_new = np.linalg.solve(np.identity(len(A)) - np.dot(dt, A), \
      w + np.dot(dt, dw_dt - np.dot(A, w)))

    return w_new

def euler_pred_corr(w, t, dt, problem):

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

    dw_dt, A, _, _, _  = problem(w, t) # Returns the first and second outputs of problem()

    # Solves the linear system of equations
    w_new = np.linalg.solve(np.identity(len(A)) - np.dot(dt, A / 2), \
      w + np.dot(dt, (2 * dw_dt - np.dot(A, w)) / 2))

    return w_new

def runge_kutta(w, t, dt, problem):

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