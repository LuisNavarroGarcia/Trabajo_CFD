def dt_adaptative(w, d_t, t):


  TOL = 0.01
  dt_max = 1E1
  dt_min = 1E-2
  err = np.linalg.norm(w[-1] - w[-2]) \
        / np.linalg.norm(w[-2])
  dt_opt = (TOL / err) * d_t

  if err > TOL:

    dt_n = min(max(max(dt_opt, d_t * 0.5), dt_min), dt.courant)
  
  else:

    dt_n = min(min(dt_opt, d_t * 1.5), dt_max)

  return dt_n