import numpy as np

def courant(mesh, u, sim_config):

  ''' 
  function courant : calculated the maximun time step to fulfill
  the condition of courant stability.

  PARAMETERS
  ----------

  mesh : class with all the mesh properties:
      - Rn : coordinates of nodes forming the mesh
      - Cells : index of nodes forming each cell
      - Volumes : area (2D mesh) of each cell
      - Neighbours : neighboring cells for each cell and/or boundaries
      - Normals : external normal vector for each face of each cell
      - Areas : length of each face of each cell
      - Faces : coordinates of the center of each face for each cell

  u : velocity of the flow, as a function of time [s] and position [m]

  sim_config : class that stores several properties
      - courant : courant constant
      - tfinal : simulation time [s]

  OUTPUT
  ----------

  dt_courant : maximum time step given by Courant[s]. 
  
  '''

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

  ''' 
  function dt_adaptative : updates the time step for each iteration
  for a given error.

  PARAMETERS
  ----------

  w : state vector at a given time

  t : time [s]

  dt : time step [s]

  OUTPUT
  ----------

  dt_n : adapted time step [s]. 
  
  '''

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
  ''' 
  function dt_constant : returns a constant time step given
  by the initial one specified.

  PARAMETERS
  ----------

  dt0 : Initial time step [s].

  OUTPUT
  ----------

  dt : Time step [s]. 
  
  '''
  
  dt = dt0

  return dt

class DT():

  """
        class DT : stores several properties related to time
        
        PROPERTIES
        ----------
        - maximum : maximum simulation time
        - calc : static method to calculate the time step given
            the state vector and a given time
        - dt0 : initial time step
        - courant : courant constant
        
        """

  def __init__(self, maximum , dt_calc, dt0, courant):
    self.max = maximum
    self.calc = dt_calc
    self.dt0 = dt0
    self.courant = courant