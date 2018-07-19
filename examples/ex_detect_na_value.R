# Example for detect_na_value function
data_ctrl = penda::data_ctrl
data_case = penda::data_case
na_values = detect_na_value(data_ctrl, data_case, threshold = 0.99)
data_ctrl_f = data_case[!na_values,]
data_case_f = data_case[!na_values,]
