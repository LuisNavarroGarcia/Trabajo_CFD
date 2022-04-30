from stop_conditions import time_step, max_error, max_iterations, mean_error
from plots import ErrorPlot

def stopping_criterion(v_criteria, v_values, v_AndOr, 
                      activation_plots, wsol, t, iteration):

    """function stopping_criterion : Determines if the problem has 
            coverged. It takes into account the convergence ccriterion desired 
            the user

    PARAMETERS
    ----------

    v_criteria : vector with 4 components
        It indicates the converging criteria that are active
        1 indicates active criterion
        0 indicates inactive criterion
        Component 0: physical time
        Component 1: number of iterations
        Component 2: maximum error in the last A iterations, lower 
        than v_values(i) 
        Component 3: mean error in the last B iterations, lower 
        than v_values(i)
        

    v_values : vector with 6 components 
        It gathers the most important values of the different 
        converging criteria, DEFINED BY THE USER
        Component 0: refers to the physical time. It indicates the
        maximum physical time of the simulation [s]
        Component 1: refers to the number of iterations. It indicates
        the maximum number of iterations
        Component 2: refers to the maximum error. It indicates the 
        value of the maximum error desired
        Component 3: refers to the maximum error. It indicates the 
        number of iterations desired for the maximum error
        Component 4: refers to the mean error. It indicates the 
        value of the mean error desired
        Component 5: refers to the mean error. It indicates the 
        number of iterations desired for the mean error 

    v_AndOr : vector with the same number of elements as the numer of 
        sopping criteria
        Each element can be:
        0 if it's an "or" boolean condition
        1 if it's an "and" boolean condition

    activation_plots : configuration of the plot error
        (0): activation of the plot. 
            0 if it is not activated
            1 if it is activaned
        (1): type of error to plot.
            0 if it is the maximum error
            1 if ir is the mean error
        (2): sample frecuency in iterations 

    w_sol: vector with as many elemens as cells there are in the domain 
        The temperature result in [k] corresponding to each cell is
        gathered in each element of the vector

    t: variable indicating the physical time

    iteration : variable indicating the itertion of the porblem. It changes
        with each iteration

    RETURNS 
    -------        

    stop : indicates if any of the converging criteria is satisfied
        Returns 0 if the calculation must continue 
        Returns 1 if the claculation is finished because the condition
        is satisfied


    stop_condition : returns a value indicating which convergence 
        criterion has been accomplished. 
        Returns 1 if the criterion acomplished is the physical time
        Returns 2 if the criterion acomplished is the number of iterations
        Returns 3 if the criterion acomplished is the maximum time 
        Returns 4 if the criterion acomplished is the mean time 
             """

    import numpy as np

    stop_condition = np.zeros(4)
    calculated_value = np.zeros(4)
    # print(stop_condition)
    # print(calculated_value)

    for i in range(np.size(v_criteria)):
        if v_criteria[i] == 1:
            
            if i == 0:
                stop_condition[i], calculated_value[i] = time_step(t, v_values[i])

            if i == 1: 
                stop_condition[i], calculated_value[i] = max_iterations(iteration, v_values[i])
        
            if i == 2: 
                stop_condition[i], calculated_value[i] = max_error(wsol, v_values[i],
                v_values[i+1])

            if i == 3:
                stop_condition[i], calculated_value[i] = mean_error(wsol, v_values[i+1]
                , v_values[i+2])
    
        
    error_plot = ErrorPlot()

    if activation_plots[0] == 1:
        if i==2 and v_criteria[2]==1 and activation_plots[1]==0 and (iteration%activation_plots[2]) == 0:
            error_plot(it = iteration, point = calculated_value[2], error_value = v_values[2], 
                       labelplot = 'Maximum error [-]' )
        if i==3 and v_criteria[3]==1 and activation_plots[1]==1 and (iteration% activation_plots[2]) == 0: 
            error_plot(it = iteration, point = calculated_value[3], error_value = v_values[4], 
                       labelplot = 'Mean error [-]' )

    not_v_criteria = 1 - v_criteria

    stop = 0
    for i in range(np.size(v_criteria)):
        if v_AndOr[i] == 0 and stop_condition[i] == 1:
                stop = 1
                break 
                
    if stop == 0: 
        stop = int(np.product(not_v_criteria + v_AndOr*stop_condition) > 0)

    return stop, stop_condition