# First, load data
data = penda::penda_data_ctrl[1:10, 1:10]

# Second, simulated the dysregulation
simulation = penda::simplified_simulation(data,
                                          proportion = 0.3,
                                          threshold = 60,
                                          modifier = 30,
                                          factor = 4)
