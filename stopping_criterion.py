import numpy as np
from stop_conditions import time_step, max_error, max_iterations, mean_error


def stopping_criterion(v_criteria, v_values, v_AndOr, 
                      activation_plots, wsol, t, iteration, error_plot):

    """function stopping_criterion : function to detrmine if the summulation 
            has converged.  It takes into account the convergence ccriterion desired 
            the user, if the criteria are active. 
            

    PARAMETERS
    ----------

    v_criteria : vector array with 4 components
        It indicates the converging criteria that are active
        1 indicates active criterion
        0 indicates inactive criterion
        Component 0: simmulation time
        Component 1: number of maximum iterations
        Component 2: maximum error in the last defined iterations, lower 
        than a certain threshold defined in the vector array v_values 
        Component 3: mean error in the last defined iterations, lower 
        than a certain threshold defined in the vector array v_values 
        

    v_values : vector array with 6 components 
        It gathers the most important values of the different 
        converging criteria. They will have to be defined by the user
        if the corresponding criterion defined in v_criteria is active 
        Component 0: refers to the simmulation time. It indicates the
        maximum physical time of the simulation, in seconds
        Component 1: refers to the number of iterations. It indicates
        the maximum number of iterations
        Component 2: refers to the maximum error. It indicates the 
        value of the maximum error desired (the threshold value)
        Component 3: refers to the maximum error. It indicates the 
        number of iterations desired for the calculation of the maximum error
        Component 4: refers to the mean error. It indicates the 
        value of the mean error desired (the threshold value)
        Component 5: refers to the mean error. It indicates the 
        number of iterations desired for the calculation of the mean error

    v_AndOr : vector with the same number of elements as the numer of 
        sopping criteria
        Each element can be:
        0 if it's an "or" boolean condition
        1 if it's an "and" boolean condition

    activation_plots : configuration of the plot error
        (0): activation of the plot. 
            0 if it is not activated. The converging error plot will not appear
            during the simmulation
            1 if it is activaned. The converging error plot will appear
            during the simmulation
        (1): activation of the maximum error plot
            0 if the maximum error is not active
            1 if the maximum error is active
        (2): activation of the mean error plot
            0 if the mean error is not active
            1 if the mean error is active
        (3): sample frecuency in iterations 

    w_sol: vector with as many elemens as cells there are in the domain 
        The temperature result in [k] corresponding to each cell is
        gathered in each element of the vector

    t: variable indicating the simmulation time

    iteration : variable indicating the itertion of the porblem. It changes
        with each iteration

    RETURNS 
    -------        

    stop : indicates if any of the converging criteria is satisfied
        Returns 0 if the calculation must continue,, according to the 
        criteria defined by the user.
        Returns 1 if the claculation is finished because the condition
        is satisfied, according to the criteria defined by the user.


    stop_condition : returns a value indicating which convergence 
        criterion has been accomplished. 
        Returns 1 if the criterion acomplished is the simmulation time
        Returns 2 if the criterion acomplished is the number of iterations
        Returns 3 if the criterion acomplished is the maximum error 
        Returns 4 if the criterion acomplished is the mean error 
        """ 

    stop_condition = np.zeros(4)
    calculated_value = np.zeros(4)

    for i in range(np.size(v_criteria)):
        if v_criteria[i] == 1:
            
            if i == 0:
                stop_condition[i], calculated_value[i] = time_step(t, v_values[i])

            elif i == 1: 
                stop_condition[i], calculated_value[i] = max_iterations(iteration, v_values[i])
        
            elif i == 2: 
                stop_condition[i], calculated_value[i] = max_error(wsol, v_values[i], v_values[i+1])
                if activation_plots[0] == 1 and v_criteria[2]==1 and activation_plots[1]==1 and (iteration%activation_plots[3]) == 0:
                    error_plot(it = iteration, point = calculated_value[2], error_value = v_values[2], plot_type = 0, act_plot = activation_plots)

            elif i == 3:
                stop_condition[i], calculated_value[i] = mean_error(wsol, v_values[i+1], v_values[i+2])
                if activation_plots[0] == 1 and v_criteria[3]==1 and activation_plots[2]==1 and (iteration%activation_plots[3]) == 0: 
                    error_plot(it = iteration, point = calculated_value[3], error_value = v_values[4], plot_type = 1, act_plot = activation_plots)

    not_v_criteria = 1 - v_criteria
    not_v_AndOr = 1 - v_AndOr

    stop = 0
    for i in range(np.size(v_criteria)):
        if v_AndOr[i] == 0 and stop_condition[i] == 1:
                stop = 1
                break 
                
    if stop == 0: 
        if np.product(not_v_criteria + not_v_AndOr*not_v_criteria + v_AndOr*stop_condition) > 0: 
            stop =1

    return stop, stop_condition

