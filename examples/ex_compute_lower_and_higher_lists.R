# First, data are filtering to remove low expression values
controls = penda::penda_data_ctrl[1:10, 1:10]
# Second, lower and higher lists are computed
L_H_list = penda::compute_lower_and_higher_lists(controls, threshold = 0.99, s_max = 50)
