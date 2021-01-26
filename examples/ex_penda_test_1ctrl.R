# First, load and define the data
control = penda::penda_data_ctrl[1:10, 1]
samples = penda::penda_data_case[1:10, 1:3]
# Second, lower and higher lists are computed
L_H_list = penda::compute_lower_and_higher_lists_1ctrl(control, s_max = 50)
# Third, test is the expression is deregulated for each gene of a sample
res_penda = penda::penda_test_1ctrl(samples,
                              iterations = 20,
                              L_H_list = L_H_list,
                              threshold = 0.03)
