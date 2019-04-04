# First, load and define the data
controls = penda::penda_data_ctrl[1:10, 1:10]
samples = penda::penda_data_case[1:10, 1:3]
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
