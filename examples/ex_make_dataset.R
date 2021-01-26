# Example for make_dataset function
data_ctrl = penda::penda_data_ctrl
data_case = penda::penda_data_case
dataset = penda::make_dataset(controls = data_ctrl,
                              data_case = data_case,
                              detectlowvalue = TRUE,
                              detectNA = TRUE,
                              threshold = 0.99,
                              val_min = NA,
                              bimod = TRUE)
data_ctrl = dataset$data_ctrl
data_case = dataset$data_case
