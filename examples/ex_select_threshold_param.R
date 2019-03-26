# First, load and define the data, make the simulation
controls = penda::penda_data_ctrl[1:10, 1:10]
simu_data = penda::penda_data_ctrl[1:10,11:16]
simulation = penda::simplified_simulation(simu_data, fraction = 0.3, threshold = 60)
# Second, lower and higher lists are computed
L_H_list = penda::compute_lower_and_higher_lists(controls, threshold = 0.99, s_max = 50)
# Third, choose the vector of threshold to test
threshold_values = c(0.1, 0.2, 0.3, 0.4, 0.5)
# Fourth, make the test on simulation for different thresholds
which_threshold = penda::choose_threshold(controls,
                                          L_H_list,
                                          iterations = 20,
                                          simulation,
                                          threshold_values,
                                          quant_test = 0,
                                          factor_test = 1)
# Fifth, use the function to automatically choose the better threshold
best_threshold = penda::select_threshold_param(which_threshold, FDR_goal = 0.05)
