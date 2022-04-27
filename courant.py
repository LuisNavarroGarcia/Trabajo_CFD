def courant(mesh, u, sim_config):

  [m, n] = np.shape(mesh.cells)
  DeltaX_min = 100
    
    '''CUIDADO INDENTADO'''
    for i in range(m):
      xs = []
      ys = []

        for j in range(n):
          xs = np.array(xs, mesh.Rc[mesh.cells[i, j], 0])
          ys = np.array(ys, mesh.Rc[mesh.cells[i, j], 1])

      dxs = np.array(xs[0] - xs[1], xs[0] - xs[2], xs[2] - xs[1])
      dys = np.array(ys[0] - ys[1], ys[0] - ys[2], ys[2] - ys[1])
      dist = np.sqrt(dxs**2 + dys**2)
      min_actual = min(dist)
      DeltaX_min = min(min_actual, DeltaX_min)

    # Max. Velocity. Sampling all the function is needed
    eigenvalue = max(max(u(np.linspace(min(mesh.Rn[:, 0]), max(mesh.Rn[:, 0]), 21), \
    np.linspace(min(mesh.Rn[:, 1]), max(mesh.Rn[:, 1]), 21), \
    np. linspace(0, sim_config.tfinal, 21))))
                

    # Applying Courant's stability condition
    C = sim_config.courant
    dt_courant = C * (DeltaX_min / eigenvalue)
      
    return dt_courant