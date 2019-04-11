# First, data are ranking by median
controls = penda::penda_data_ctrl[1:10, 1:10]
median_gene = apply(controls, 1, median, na.rm = TRUE)
median_gene = sort(median_gene)
controls = controls[names(median_gene), ]
# Second, lower and higher lists are computed
L_H_list = penda::compute_lower_and_higher_lists(controls, threshold = 0.99, s_max = 50)
