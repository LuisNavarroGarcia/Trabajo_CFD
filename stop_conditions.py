import numpy as np

def time_step(time, time_max):
    
    """stop_condition_tphysical : it determines if the function
    has converged, only in case it has reached the physical 
    maximum time
    
    
    PARAMETERS
    ----------

    time : last moment in which the calculation has been carried out

    time_max : maximum physical time 

    RETURNS 
    -------

    stop : indicates wheter the calculation must stop or continue
            Returns 0 if the caltulation must continue
            Returns 1 if the calculation must stop

    time : Returns the same variable as the input
        """
    stop = 0
    if time >= time_max:
        stop = 1
   
    return stop, time

def max_iterations(iteration_number, iteration_max):
    
    """stop_condition_iterations : it determines if the function
    has converged, only in case it has reached the number of iterations
    IMPOSED
    
    
    PARAMETERS
    ----------

    iteration_number : last iteration in which the calculation has been 
        carried out

    iteration_max : number of maximum iterations DETERMIENED BY THE USER 

    RETURNS 
    -------

    stop : indicates wheter the calculation must stop or continue
            Returns 0 if the caltulation must continue
            Returns 1 if the calculation must stop

    iteration_number : Returns the same variable as the input
        """
    stop = 0
    if iteration_number >= iteration_max:
        stop = 1
    
    return stop, iteration_number    

def mean_error(temperatures, error_mean_limit,
                             num_iterations):
    
    """stop_condition_error_mean : it determines if the function
    has converged,using the criterion of the mean error allowed in 
    the last iterations specified
    
    PARAMETERS
    ----------
    
    temperatures : solution matrix of the problem
    
    error_mean_limit : mean error allowed in the last ierations 
    
    num_iterations : number of iterations where the maximum error 
                        is evaluated. SPECIFIED BY THE USER
    
    
    RETURNS 
    -------
    
    stop : indicates wheter the calculation must stop or continue
            Returns 0 if the caltulation must continue
            Returns 1 if the calculation must stop
            
    error_mean : mean error of the calculations in the last iterations
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
    
    """stop_condition_error_max : it determines if the function
    has converged,using the criterion of the maximum error allowed in 
    the last iterations specified
    
    PARAMETERS
    ----------
    
    temperatures : solution matrix of the problem
    
    error_max_limit : maximum error allowed in the last ierations 
    
    num_iterations : number of iterations where the maximum error 
                        is evaluated. SPECIFIED BY THE USER
    
    
    RETURNS 
    -------
    
    stop : indicates wheter the calculation must stop or continue
            Returns 0 if the caltulation must continue
            Returns 1 if the calculation must stop
            
    error_max : maximum error of the calculations in the last iterations
    """
    import numpy as np

    b = np.shape(temperatures)[1]

    if b < num_iterations:
        stop = 0
        error_max = 1
    else: 
        max_error_matrix = temperatures[:, (b - int(num_iterations)):b]
        maximum_terms = np.max(max_error_matrix, axis=1)
        minimum_terms = np.min(max_error_matrix, axis=1)
        error_max = np.mean(np.abs(maximum_terms - minimum_terms)/np.abs(minimum_terms))

    if error_max <= error_max_limit:
        stop = 1 #Stop calculation is maximum error is low enougth

    else:
        stop = 0 #Continue if maximum error is still large

    return stop, error_max