# Example for detect_zero_value function
data_ctrl = penda::data_ctrl
data_case = penda::data_case
null_values = detect_zero_value(data_ctrl, data_case, threshold = 0.99, min = 10)
data_ctrl_filtered = data_case[!null_values,]
data_case_filtered = data_case[!null_values,]
