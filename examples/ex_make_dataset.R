# Example for detect_zero_value function
data_ctrl = penda::penda_data_ctrl
data_case = penda::penda_data_case
dataset = penda::make_dataset(controls = data_ctrl,
                              cancer_data = data_case,
                              detectlowvalue = TRUE,
                              detectNA = TRUE,
                              threshold = 0.99,
                              val_min = NA)
data_ctrl = dataset$data_ctrl
data_case = dataset$data_case
