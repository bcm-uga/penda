# First, load and define the data
dataset = penda::make_dataset(
  penda::penda_data_ctrl[1:10, 1:10], 
  penda::penda_data_case[1:10, 1:2]
)
controls = dataset$data_ctrl
sample = dataset$data_case[,1]
gene = names(sample)[1]
# Second, lower and higher list are computed
L_H_list = penda::compute_lower_and_higher_lists(controls, threshold = 0.99, s_max = 50)
# Third, test is the expression is deregulated for a given gene.
#   When this function is called, quantiles are already computed
quant_test = 0
factor_test = 1
quantile_genes = apply(controls, 1, quantile, c(quant_test,(1-quant_test)), na.rm = TRUE)
quantile_genes[1,] = quantile_genes[1,] / factor_test
quantile_genes[2,] = quantile_genes[2,] * factor_test

expression = penda::regulation_test(gene,
                                    L_H_list,
                                    sample,
                                    threshold = 0.03)
