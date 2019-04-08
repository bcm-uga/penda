# Example for detect_zero_value function
data_ctrl = penda::penda_data_ctrl
data_case = penda::penda_data_case
null_values = penda::detect_zero_value(data_ctrl, data_case, threshold = 0.99, min = 10)
data_ctrl_f = data_case[!null_values,]
data_case_f = data_case[!null_values,]
