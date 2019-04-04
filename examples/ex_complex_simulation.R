# First, load and define the data
controls = penda::penda_data_ctrl[1:1000, 1:10]
samples = penda::penda_data_case[1:1000, 1:3]
simu_data = penda::penda_data_ctrl[1:1000,11:16]
# Second, make the simulation
simulation = penda::complex_simulation(controls,
                                       samples,
                                       simu_data,
                                       size_grp = 80,
                                       quant = 0.05)
