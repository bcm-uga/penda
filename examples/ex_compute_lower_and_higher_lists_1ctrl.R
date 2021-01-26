# First, data are ranking by median
control = penda::penda_data_ctrl[1:10, 1]
control = sort(control)
# Second, lower and higher lists are computed
L_H_list = penda::compute_lower_and_higher_lists_1ctrl(control, s_max = 50)
