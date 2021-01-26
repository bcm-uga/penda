# First, load and define the data, make the simulation
controls = penda::penda_data_ctrl[1:10, 1:10]
simu_data = penda::penda_data_ctrl[1:10,11:16]
simulation = penda::simplified_simulation(simu_data, proportion = 0.3, threshold = 60)
# Second, choose the vectors of quantiles and factors to test
quantile_values = c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.30, 0.35, 0.4, 0.45)
factor_values = c(1, 1.2, 1.4, 1.6, 1.8, 2)
# Third, make the test on simulation for different quantiles
which_quantile = penda::choose_quantile(controls,
                                simulation,
                                factor_values,
                                quantile_values)
