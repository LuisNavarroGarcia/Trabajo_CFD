def dt_constant(dt0):
  ''' dt_constant  Time step calculator
  Updates the required time step from the state vector solution 
  and the previous time step.

  Inputs: dt0: Used time step [s].
  Outputs: dt: Time step [s]. '''
  
  dt = dt0

  return dt