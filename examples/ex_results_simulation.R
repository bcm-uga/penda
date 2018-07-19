# First, load and define the data, make the simulation
controls = penda::data_ctrl[1:10, 1:10]
simu_data = penda::data_ctrl[1:10,11:16]
simulation = penda::simplified_simulation(simu_data, fraction = 0.3, threshold = 60)
# Second, down and up list are computed
D_U_list = penda::compute_down_and_up_list(controls, threshold = 0.99, s_max = 50)
# Third, test is the expression is deregulated for each gene of a sample
res_penda = penda::penda_test(samples,
                              controls,
                              iterations = 20,
                              D_U_list,
                              threshold = 0.03,
                              quant_test = 0,
                              factor_test = 1)
# Fourth, compare the result of the simulation to the reality
results = penda::results_simulation(res_penda$D,
                                    res_penda$U,
                                    simulation)

