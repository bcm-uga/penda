# Example for detect_zero_value function
null_values = detect_zero_value(ctrl_data, simu_data, threshold = 0.8, min = 10)
ctrl_data = ctrl_data[!null_values,]
simu_data = simu_data[!null_values,]