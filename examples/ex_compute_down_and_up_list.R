# First, data are filtering to remove low expression values
controls = penda::data_ctrl[1:10, 1:10]
# Second, down and up list are computed
D_U_list = penda::compute_down_and_up_list(controls, threshold = 0.99, s_max = 50)
