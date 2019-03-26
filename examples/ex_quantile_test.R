# First, load and define the data
controls = penda::penda_data_ctrl[1:10, 1:10]
samples = penda::penda_data_case[1:10, 1:3]
# Second, test is the expression is deregulated for each gene of a sample
res_quantile = penda::quantile_test(controls,
                                samples,
                                quant = 0.03,
                                factor = 1.4)
