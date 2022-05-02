import numpy as np

def courant(mesh, u, sim_config):

  m, n = np.shape(mesh.cells)
  DeltaX_min = 100
    
  for i in range(m):
    xs = np.array([])
    ys = np.array([])

    for j in range(n):
      xs = np.append(xs, mesh.Rn[mesh.cells[i, j]-1][0])
      ys = np.append(ys, mesh.Rn[mesh.cells[i, j]-1][1])

    dxs = np.array([xs[0] - xs[1], xs[0] - xs[2], xs[2] - xs[1]])
    dys = np.array([ys[0] - ys[1], ys[0] - ys[2], ys[2] - ys[1]])
    dist = np.sqrt(dxs**2 + dys**2)
    min_actual = min(dist)
    DeltaX_min = min(min_actual, DeltaX_min)

  # Max. Velocity. Sampling all the function is needed
  lmda = np.max(u(np.linspace(min(mesh.Rn[:, 0]), max(mesh.Rn[:, 0]), 20), \
  np.linspace(min(mesh.Rn[:, 1]), max(mesh.Rn[:, 1]), 20), \
  np. linspace(0, sim_config.tfinal, 20)))
              

  # Applying Courant's stability condition
  C = sim_config.courant
  try:
    dt_courant = C * (DeltaX_min / lmda)
  except ZeroDivisionError:
    dt_courant = ((DeltaX_min*C)/abs(DeltaX_min*C))*np.inf
  return dt_courant

def dt_adaptative(w, d_t, dt):


  TOL = 0.005
  dt_max = 1
  dt_min = 1E-4
  err = abs(np.linalg.norm(w[:, -1] - w[:, -2]) \
        / np.linalg.norm(w[:, -2]))
  dt_opt = (TOL / err) * d_t

  if err > TOL:

    dt_n = min(max(max(dt_opt, d_t * 0.5), dt_min), dt.courant)
  
  else:

    dt_n = min(min(dt_opt, d_t * 1.5), dt_max)

  return dt_n

def dt_constant(dt0):
  ''' dt_constant  Time step calculator
  Updates the required time step from the state vector solution 
  and the previous time step.

  Inputs: dt0: Used time step [s].
  Outputs: dt: Time step [s]. '''
  
  dt = dt0

  return dt

class DT():
  def __init__(self, maximum , dt_calc, dt0, courant):
    self.max = maximum
    self.calc = dt_calc
    self.dt0 = dt0
    self.courant = courant