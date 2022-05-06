import numpy as np

def time_step(time, time_max):
    
    """time_step : function to define if the simmulation 
    has converged, but only in case it has 
    reached the maximum time specified.
    
    PARAMETERS
    ----------

    time : last moment in which the calculation has been carried out

    time_max : maximum time over which the simmulation cannot continue. 
    It has to be specified by the user if this converging 
    criteria is active.  

    OUTPUTS 
    -------

    stop : indicates wheter the calculation must stop or continue.
            Returns 0 if the caltulation must continue
            Returns 1 if the calculation must stop

    time : Returns the simmulation time. It is the same value as the input. 
    """
    stop = 0
    if time >= time_max:
        stop = 1
   
    return stop, time

def max_iterations(iteration_number, iteration_max):
    
    """max_iterations : time_step : function to define if the simmulation 
    has converged, but only in case it has 
    reached the maximum time specified.
    
    PARAMETERS
    ----------

    iteration_number : last iteration in which the calculation has been 
        carried out.

    iteration_max : maximum number of iterations over which the 
    simmulation cannot continue. It has to be specified by the user
    if this converging criteria is active.  

    OUTPUTS 
    -------

    stop : indicates wheter the calculation must stop or continue.
            Returns 0 if the caltulation must continue.
            Returns 1 if the calculation must stop.

    iteration_number :Returns the simmulation time. It is the same 
    value as the input. 
    """
    stop = 0
    if iteration_number >= iteration_max:
        stop = 1
    
    return stop, iteration_number    

def mean_error(temperatures, error_mean_limit,
                             num_iterations):
    
    """
    mean_error : function to define if the simmulation 
    has converged, but only using the criterion of the mean error 
    allowed in the last iterations specified.
    
    PARAMETERS
    ----------
    
    temperatures : solution matrix of the problem, containig in each column 
    the state vector correspongding to each iteration. 
    
    error_mean_limit : mean error allowed in the last ierations. 
    It has to be specified by the user if the criteria is active.  
    
    num_iterations : number of iterations where the maximum error is evaluated. 
    It has to be specified by the user if the criteria is active.  
    
    
    OUTPUTS 
    -------
    
    stop : indicates wheter the calculation must stop or continue
            Returns 0 if the caltulation must continue
            Returns 1 if the calculation must stop
            
    error_mean : mean error of the calculations in the last iterations
    specified by the user

    """
    
    b = np.shape(temperatures)[1]

    if b < num_iterations:
        stop = 0
        error_mean = 1
    
    else: 
        mean_matrix = temperatures[:, (b - int(num_iterations)):b]
        maximum_terms = np.max(mean_matrix, axis=1)
        minimum_terms = np.min(mean_matrix, axis=1)
        error_mean = np.mean(np.abs(maximum_terms - minimum_terms)/np.abs(minimum_terms))

        if error_mean <= error_mean_limit:
            stop = 1 #Stop calculation is mean error is low enougth

        else:
            stop = 0 #Continue if mean error is still large

    return stop, error_mean

def max_error(temperatures, error_max_limit,
                             num_iterations):
    
    """
    max_error : function to define if the simmulation 
    has converged, but only using the criterion of the maximum error 
    allowed in the last iterations specified.
    
    PARAMETERS
    ----------
    
    temperatures : solution matrix of the problem, containig in each column 
    the state vector correspongding to each iteration. 
    
    error_max_limit : maximum error allowed in the last ierations. 
    It has to be specified by the user if the criteria is active.  
    
    num_iterations : number of iterations where the mean error is evaluated. 
    It has to be specified by the user if the criteria is active.
    
    
    OUTPUTS 
    -------
    
    stop : indicates wheter the calculation must stop or continue
            Returns 0 if the caltulation must continue
            Returns 1 if the calculation must stop
            
    error_max : maximum error of the calculations in the last iterations
    specified by the user
    """

    b = np.shape(temperatures)[1]

    if b < num_iterations:
        stop = 0
        error_max = 1
    else: 
        max_error_matrix = temperatures[:, (b - int(num_iterations)):b]
        maximum_terms = np.max(max_error_matrix, axis=1)
        minimum_terms = np.min(max_error_matrix, axis=1)
        error_max = np.max(np.abs(maximum_terms - minimum_terms)/np.abs(minimum_terms))

    if error_max <= error_max_limit:
        stop = 1 #Stop calculation is maximum error is low enougth

    else:
        stop = 0 #Continue if maximum error is still large

    return stop, error_max