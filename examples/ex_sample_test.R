# First, load and define the data
controls = penda::data_ctrl[1:10, 1:10]
sample = penda::data_case[1:10, 1]
# Second, down and up list are computed
D_U_list = penda::compute_down_and_up_list(controls, threshold = 0.99, s_max = 50)
# Third, test is the expression is deregulated for each gene of a sample
res_sample = penda::sample_test(sample,
                                controls,
                                iterations = 20,
                                D_U_list, 
                                    threshold = 0.03, 
                                    quant_test = 0, 
                                    factor_test = 1)
