# First, load and define the data, make the simulation
controls = penda::penda_data_ctrl[1:10, 1:13]
median_gene = apply(controls, 1, median, na.rm = TRUE)
median_gene = sort(median_gene)
controls = controls[names(median_gene), ]
simu_data = controls[,11:13]
controls = controls[,-(11:13)]
simulation = penda::simplified_simulation(simu_data, proportion = 0.3, threshold = 60)
samples = simulation$simulated_data
# Second, lower and higher lists are computed
L_H_list = penda::compute_lower_and_higher_lists(controls, threshold = 0.99, s_max = 50)
# Third, test is the expression is deregulated for each gene of a sample
res_penda = penda::penda_test(samples,
                              controls,
                              iterations = 20,
                              L_H_list,
                              threshold = 0.03,
                              quant_test = 0,
                              factor_test = 1)
# Fourth, compare the result of the simulation to the reality
results = penda::results_simulation(res_penda$D,
                                    res_penda$U,
                                    simulation)
# Fifth, draw the results
penda::draw_results(results)

