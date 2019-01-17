# First, data are filtering to remove low expression values
controls = penda::data_ctrl[1:10, 1:10]
# Second, lower and higher lists are computed
D_U_list = penda::compute_lower_and_higher_lists(controls, threshold = 0.99, s_max = 50)
