# First, load and define the data
dataset = penda::make_dataset(
  penda::penda_data_ctrl[1:10, 1:10],  
  penda::penda_data_case[1:10, 1:2]
)
controls = dataset$data_ctrl
sample = dataset$data_case[,1]
# Second, lower and higher lists are computed
L_H_list = penda::compute_lower_and_higher_lists(controls, threshold = 0.99, s_max = 50)
# Third, test is the expression is deregulated for each gene of a sample.
#When this function is called, quantiles are already computed
quant_test = 0
factor_test = 1
quantile_genes = apply(controls, 1, quantile, c(quant_test,(1-quant_test)), na.rm = TRUE)
quantile_genes[1,] = quantile_genes[1,] / factor_test
quantile_genes[2,] = quantile_genes[2,] * factor_test

res_sample = penda::sample_test(sample,
                                iterations = 20,
                                threshold = 0.03)
