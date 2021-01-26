# First, load and define the data, make the simulation
controls = penda::penda_data_ctrl[1:10, 1:13]
median_gene = apply(controls, 1, median, na.rm = TRUE)
median_gene = sort(median_gene)
controls = controls[names(median_gene), ]
simu_data = controls[,13]
controls = controls[,-13]

simulation = penda::simplified_simulation(simu_data, proportion = 0.3, threshold = 60)
# Second, lower and higher lists are computed
L_H_list = penda::compute_lower_and_higher_lists(controls, threshold = 0.99, s_max = 50)
# Third, choose the vector of threshold to test
threshold_values = c(0.1, 0.2, 0.3, 0.4, 0.5)
# Fourth, make the test on simulation for different thresholds.
#   When this function is called, quantiles are already computed
quant_test = 0
factor_test = 1
quantile_genes = apply(controls, 1, quantile, c(quant_test,(1-quant_test)), na.rm = TRUE)
quantile_genes[1,] = quantile_genes[1,] / factor_test
# quantile_genes[2,] = quantile_genes[2,] * factor_test
res_penda = penda::test_multiple_thresholds(simulation,
                                threshold_values,
                                iterations = 20)
