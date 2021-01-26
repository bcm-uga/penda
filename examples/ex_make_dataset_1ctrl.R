# Example for make_dataset_1ctrl function
data_ctrl = penda::penda_data_ctrl[1:10, 1]
data_case = penda::penda_data_case[1:10, 1:3]
dataset = penda::make_dataset_1ctrl(control = data_ctrl,
                                    data_case = data_case)
data_ctrl = dataset$data_ctrl
data_case = dataset$data_case
