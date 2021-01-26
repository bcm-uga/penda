# First, load and define the data
dataset = penda::make_dataset_1ctrl(
  penda::penda_data_ctrl[1:10, 1],
  penda::penda_data_case[1:10, 1:3]
)
controls = dataset$data_ctrl
sample = dataset$data_case[,1]
gene = names(sample)[1]
# Second, lower and higher list are computed
L_H_list = penda::compute_lower_and_higher_lists_1ctrl(controls, s_max = 50)
# Third, test is the expression is deregulated for a given gene.

expression = penda::regulation_test_1ctrl(gene,
                                    L_H_list,
                                    sample,
                                    threshold = 0.03)
