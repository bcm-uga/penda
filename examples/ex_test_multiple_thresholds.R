# First, load and define the data, make the simulation
controls = penda::data_ctrl[1:10, 1:10]
simu_data = penda::data_ctrl[1:10,11:16]
simulation = penda::simplified_simulation(simu_data, fraction = 0.3, threshold = 60)
# Second, down and up list are computed
D_U_list = penda::compute_down_and_up_list(controls, threshold = 0.99, s_max = 50)
# Third, choose the vector of threshold to test
threshold_values = c(0.1, 0.2, 0.3, 0.4, 0.5)
# Fourth, make the test on simulation for different thresholds
res_penda = penda::test_multiple_thresholds(controls,
                                D_U_list,
                                simulation,
                                threshold_values,
                                iterations = 20,
                                quant_test = 0,
                                factor_test = 1)
