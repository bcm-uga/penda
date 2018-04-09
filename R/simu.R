#'simplified_simulation
#'
#' This function simulates the dysregulation (up and down) of datas.
#'
#'@param data The vector or matrix of original data to edit.
#'@param fraction The fraction of data to modify.
#'@param threshold The expression value under which data modification is +- modifier
#'Above the threshold, dysregulation is expression * or / factor.
#'@param modifier Under the threshold, the up-regulation is to add at the expression this number,
#'and the down-regulation to substract this number.
#'@param factor Below the threshold, the up-regulation is to multiply expression by this number,
#'and the down-regulation to divide.
#'
#'@return This function returns a list of three vectors or matrices :
#'initial_data, datas before the modification
#'simulated_data, datas with up and down modifications
#'changes_idx, index of datas with modifications
#'
#'@example
#'simplified_simulation(simu_data, fraction = 0.3, threshold = 60)
#'
#'@export

simplified_simulation = function(data, fraction, threshold = 60, modifier = 30, factor = 4){
  simu_data = data
  perturb = sample(length(simu_data), fraction*length(simu_data))
  for(i in 1:length(perturb)){
    if (rnorm(1) < 0) {
      if (simu_data[perturb[i]] < threshold) {
        simu_data[perturb[i]] = simu_data[perturb[i]] - modifier
        if (simu_data[perturb[i]] <= 0){
          simu_data[perturb[i]] = 0
        }
      } else {
        simu_data[perturb[i]] = simu_data[perturb[i]] / factor
      }
    } else if (simu_data[perturb[i]] < threshold) {
      simu_data[perturb[i]] = simu_data[perturb[i]] + modifier
    } else {
      simu_data[perturb[i]] = simu_data[perturb[i]] * factor
    }
  }
  return(list(initial_data = data, simulated_data = simu_data, changes_idx = perturb))
}


#'get_random_variable
#'
#' This function gives a random values following a given density function.
#'
#'@param F The density function.
#'@param x A vector with the definition area of the function.
#'@param nb_r The number of random values wanted
#'
#'@return This function returns a vector with nb_r values of F.
#'
#'@example
#'x = seq(-20, 20, 0.05)
#'F = function(x, u, a){
#'  (1 / (a * sqrt(2 * pi))) * exp(-(((x - u)^2) / (2 * a^2)))
#'}
#'get_random_variable(F, x, nb_r = 10000, u = 0, a = 1)

get_random_variable = function (F, x, nb_r,...){
  y = F(x,...)
  cumul = cumsum(y) * diff(x)[1]
  qf = approxfun(cumul, x) #quantile function
  rf = function(n) qf(runif(n)) #random variable
  return(rf(nb_r))
}


#'Heaviside
#'
#' The heaviside function value is zero for argument under the step, and one for argument above.
#'
#'@param x The number to evaluate
#'@param step The limit between heaviside = 0 and heaviside = 1
#'
#'@return This function returns 0 if the argument is under the step, 1 if not.
#'
#'@example
#'Heaviside(1+98-103, step = 0)

Heaviside = function(x, step = 0){
  return((sign(x - step) + 1) / 2)
}


#'P1
#'
#' P1 models the density function of the expression multiplicative factor for the
#' up-regulation of a gene with an expression above 100.
#'
#'@param x vector of input values
#'
#'@return This function returns the value of p1 for x (between 0 and 1).
#'
#'@example

p1 = function(x){
  p1 = Heaviside(x - 0.63) * (1 + tanh((x - 1.4) / 0.4)) * exp(-x / 0.9)
  p1 = p1 / (sum(p1) * diff(x)[1])
  return(p1)
}


#'P2
#'
#' P2 models the density function of the expression divisor for the
#' down-regulation of a gene with an expression above 100.
#'
#'@param x vector of input values
#'
#'@return This function returns the value of p2 for x (between 0 and 1).

p2 = function(x){
  p2 = Heaviside(x - 0.63) * (1 + tanh((x - 1.4) / 0.5)) * exp(-x / 1.3)
  p2 = p2 / (sum(p2) * diff(x)[1])
  return(p2)
}


#'P3
#'
#' P2 models the density function of the expression modifier (to add or substract) for the
#' dysregulation of a gene with an expression under 100.
#'
#'@param x vector of input values
#'
#'@return This function returns the value of p3 for x (between 0 and 1).

p3 = function(x){
  p3 = Heaviside(x - 20) * (exp(-x * 0.007)) * exp((x - 20) * exp(-x * 0.07))
  p3 = p3 / (sum(p3) * diff(x)[1])
  return(p3)
}


#'complex_simulation
#'
#' This function simulates the dysregulation (up and down) of datas, using the distribution of
#' dysregulation in real cancer datas.
#'
#'@param data vector or matrix of original data to edit.
#'
#'@return This function returns a list of three vectors or matrices :
#'initial_data, datas before the modification
#'simulated_data, datas with up and down modifications
#'changes_idx, index of datas with modifications
#'
#'@example
#'complex_simulation(simu_data)
#'
#'@export

complex_simulation = function(data){
  simu_data = data
  simu_data = sapply(simu_data, function(g){
    #if the expression is > 100
    if (g > 100){
      #22% chance of dysregulation
      if (runif(1) <= 0.22) {
        #30% chance of up-regulation
        if (runif(1) <= 0.30) {
          g = g * 2^(get_random_variable(p1, seq(0,100,0.1), 1))
        }
        #else, down-regulation
        else {
          g = g / 2^(get_random_variable(p2, seq(0,100,0.1), 1))
        }
      }
    #if expression < 100, 30% chance of dysregulation
    } else if (runif(1) <= 0.30) {
      #50% chance of down-regulation, 50% chance of up-regulation. While expression < 0, retry.
      while(g < 0){
        if (runif(1) <= 0.50){
          g = g - (get_random_variable(p3, seq(0, 500, 0.5), 1))
        } else {
          g = g + (get_random_variable(p3, seq(0, 500, 0.5), 1))
        }
      }
    }
    return(g)
  })
  return(list(initial_data = data, simulated_data = simu_data, changes_idx = which(data != simu_data)))
}


#'results_simulation
#'
#' This function computes the number of false positive, true positive, false negative
#' and false positive results after the test on a simulation.
#'
#'@param D_matrix The matrix of down-regulated genes.
#'@param U_matrix The matrix of up-regulated genes.
#'@param simulation The list of initial data $initial_data and modified data in $simulated_data
#'
#'@return This function returns a list of four integers :
#'TP, the number of true postive results
#'FP, the number of false postive results
#'FN, the numbler of false negative results
#'TN, the numbler of true negative results
#'
#'@example
#'D_U = find_D_U_ctrl(ctrl_data, quant = 0.001, factor = 4, threshold = 0.99)
#'simulation = simplified_simulation(simu_data, fraction = 0.3, threshold = 60)
#'results_simulation(D_matrix = D_U$D, U_matrix = D_U$U, initial_data = simulation$initial_data, simulated_data = simulation$simulated_data)
#'
#'@export

results_simulation = function(D_matrix, U_matrix, simulation){
  initial_data = simulation$initial_data
  simu_data = simulation$simulated_data
  down = initial_data > simu_data
  up = initial_data < simu_data

  #True positive = real down results and real up results
  TP = (length (which(down == TRUE & D_matrix == TRUE))
        + length (which(up == TRUE & U_matrix == TRUE)))

  #False positive = down or up results, but not changed in reality
  FP = (length (which(down == 0 & D_matrix == 1))
        + length (which(up == 0 & U_matrix == 1)))

  #False negative = not changed results, but down or up in reality
  FN = (length (which((up == 1 | down == 1) & D_matrix == 0 & U_matrix == 0)))

  #TN = real not changed results
  TN = (length (which(down == 0 & up == 0 & D_matrix == 0 & U_matrix == 0)))

  return(list(TP = TP, FP = FP, FN = FN, TN = TN ))
}


#'draw_results
#'
#' This function computes the the FDR and the TPR and makes a ggplot of the results.
#'
#'@param results The list TP, FP, FN and TN return by "results_simulation"
#'
#'@return This function returns a barplot of results, FDR and TPR.
#'
#'@example
#'D_U = find_D_U_ctrl(ctrl_data, quant = 0.001, factor = 4, threshold = 0.99)
#'simulation = simplified_simulation(simu_data, fraction = 0.3, threshold = 60)
#'res = results_simulation(D_matrix = D_U$D, U_matrix = D_U$U, initial_data = simulation$initial_data, simulated_data = simulation$simulated_data)
#'draw_results(res)
#'
#'@export

draw_results = function(results){
  FP = results$FP
  df= data.frame(type = c("TP","FP","FN","TN", "FDR", "TPR"), results=c(unlist(results), (results$FP / (results$TP + results$FP)), results$TP / (results$TP + results$FN)))
  library(ggplot2)
  ggplot(df, aes(x=type, y=results, fill=type)) + geom_bar(stat="identity") + geom_text(aes(label=results), vjust=1)
}


#'DU_rearrangement
#'
#' This function reorganizes matrices resulting from the multiple test. It takes a
#' "DUDUDUDU" matrix and make two independant D and U matrices.
#'
#'@param multiple_tests The results matrix, which has a form of "DUDUDUDU".
#'@param simu_data The matrix of initial_data or simulated_data in the simulation, used for dimnames.
#'@param multiple_values The vector with different values of the test.
#'
#'@return This function returns a list of two logical matrices. The D matrix with TRUE
#'for genes down-regulated, and the U matrix with TRUE for genes up-regulated.
#'
#'@example
#'D_U = find_D_U_ctrl(ctrl_data, quant = 0.001, factor = 4, threshold = 0.99)
#'simulation = simplified_simulation(simu_data, fraction = 0.3, threshold = 60)
#'multiple_values = c(0.01, 0.02, 0.03)
#'multiple_test = sapply(multiple_values, function(value){
#'  patient_test(ctrl_data, simulated_data, quant_0 = 0.03, iterations = 10, D_U, threshold = value)
#'})
#'sorted_matrix = DU_rearrangement(multiple_test, simulation$initial_data, multiple_values)

DU_rearrangement = function(multiple_tests, simu_data, multiple_values){
  #Definition of D and U matrices for all the conditions
  D_simu = matrix(data = NA, nrow = length(simu_data)
                  , ncol = length(multiple_values)
                  , dimnames = list(names(simu_data), multiple_values))
  U_simu = matrix(data = NA, nrow = length(simu_data)
                  , ncol = length(multiple_values)
                  , dimnames = list(names(simu_data), multiple_values))

  #D and U are put in the matrices for all the conditions
  for (r in 1:length(multiple_values)){
    D_simu[,r] = unlist(multiple_tests[(r*2-1)])
    U_simu[,r] = unlist(multiple_tests[r*2])
  }
  return(list(D = D_simu, U = U_simu))
}


#'test_multiple_thresholds
#'
#' This function makes patient_test for different values of the threshold.
#'
#'@param ctrl_data A matrix with genes expressions in controls for all the patients.
#'@param D_U_ctrl The list of Down and Up-expressed genes matrices in the control.
#'@param simulation The list of initial data $initial_data and modified data in $simulated_data
#'@param threshold_values The vector of values of the threshold for regulation_test.
#'@param iterations The maximal number of iterations for the test. If the dysregulation list
#'no longer changes, iterations are stop before.
#'@param quantile The quantile of the expression distribution for the step_0.
#'
#'@return This function returns a list of two logical matrices. The D matrix with for each threshold
#'TRUE for genes down-regulated, and the U matrix with TRUE for genes up-regulated.
#'
#'@example
#'D_U = find_D_U_ctrl(ctrl_data, quant = 0.001, factor = 4, threshold = 0.99)
#'simulation = simplified_simulation(simu_data, fraction = 0.3, threshold = 60)
#'multiple_values = c(0.01, 0.02, 0.03)
#'test_threshold = test_multiple_thresholds(ctrl_data, D_U, simulation, threshold_values = multiple_values, iterations = 10, quantile = 0.03)
#'
#'@export

test_multiple_thresholds = function(ctrl_data, D_U_ctrl, simulation, threshold_values, iterations, quantile){
  #Test for all threshold values
  multiple_test = sapply(threshold_values, function(t){
    print(c("Threshold",t))
    patient_test(ctrl_data, simulation$simulated_data, quantile, iterations, D_U_ctrl, t)
  })
  sorted_test = DU_rearrangement(multiple_test, simulation$initial_data, threshold_values )
  return (sorted_test)
}

test_multiple_thresholds_sl0 = function(ctrl_data, D_U_ctrl, simulation, threshold_values, iterations, quant_test, factor_test){
  #Test for all threshold values
  multiple_test = sapply(threshold_values, function(t){
    print(c("Threshold",t))
    patient_test_sansl0(ctrl_data, simulation$simulated_data, iterations, D_U_ctrl, t, quant_test, factor_test)
  })
  sorted_test = DU_rearrangement(multiple_test, simulation$initial_data, threshold_values )
  return (sorted_test)
}


#'test_multiple_quantiles
#'
#' This function makes patient_test for different values of the quantile.
#'
#'@param ctrl_data A matrix with genes expressions in controls for all the patients.
#'@param D_U_ctrl The list of Down and Up-expressed genes matrices in the control.
#'@param simulation The list of initial data $initial_data and modified data in $simulated_data
#'@param threshold The value of the threshold for regulation_test.
#'@param iterations The maximal number of iterations for the test. If the dysregulation list
#'no longer changes, iterations are stop before.
#'@param quantile_values A vector of values for the quantile of the expression distribution for the step_0.
#'
#'@return This function returns a list of two logical matrices. The D matrix with for each quantile value
#'TRUE for genes down-regulated, and the U matrix with TRUE for genes up-regulated.
#'
#'@example
#'D_U = find_D_U_ctrl(ctrl_data, quant = 0.001, factor = 4, threshold = 0.99)
#'simulation = simplified_simulation(simu_data, fraction = 0.3, threshold = 60)
#'multiple_values = c(0.01, 0.02, 0.03)
#'test_quantile = test_multiple_quantiles(ctrl_data, D_U, simulation, threshold = 0.2, iterations = 10, quantile_values = multiple_values)
#'
#'@export

test_multiple_quantiles = function(ctrl_data, D_U_ctrl, simulation, threshold, iterations, quantile_values){
  multiple_test = sapply(quantile_values, function(q){
    print(c("Quantile",q))
    patient_test(ctrl_data, simulation$simulated_data, q, iterations, D_U_ctrl, threshold)
  })
  sorted_test = DU_rearrangement(multiple_test, simulation$initial_data, quantile_values )
  return (sorted_test)
}


#'multiple_tests
#'
#' This function makes for each patient the test with multiple values for the treshold or the quantile.
#' After that, it computes FDR and TPR.
#'
#'@param ctrl_data A matrix with genes expressions in controls for all the patients.
#'@param D_U_ctrl The list of Down and Up-expressed genes matrices in the control.
#'@param iterations The maximal number of iterations for the test. If the dysregulation list
#'no longer changes, iterations are stop before.
#'@param simulation The list of initial data $initial_data and modified data in $simulated_data
#'@param multiple_values A vector with different values to test for the chosen parameter (threshold or quantile).
#'@param unique_value The value of the parameter that does not vary.
#'@param threshold_change TRUE if you want to vary the threshold, FALSE if you want to vary the quantile.
#'
#'@return This function returns a matrix with 4 columns : the patient number, the value of the changing parameter,
#'the FDR and the TPR of the test.
#'
#'@example
#'D_U = find_D_U_ctrl(ctrl_data, quant = 0.001, factor = 4, threshold = 0.99)
#'simulation = simplified_simulation(simu_data, fraction = 0.3, threshold = 60)
#'multiple_values = c(0.01, 0.02, 0.03)
#'multiple_tests = multiple_tests(ctrl_data, D_U, iterations = 10, simulation, multiple_values, unique_value = 0.3, threshold_change = FALSE)
#'
#'@export

multiple_tests = function(ctrl_data, D_U_ctrl, iterations, simulation, multiple_values, unique_value, threshold_change = TRUE, quant_test = 0, factor_test = 1){
  results = c()
  #For each patient
  for(p in 1:ncol(simulation$initial_data)){
    print(c("Patient number",p))
    simulation_p = list(initial_data = simulation$initial_data[,p], simulated_data = simulation$simulated_data[,p])

    #The test is made for all the values, of threshold if threshold_change = TRUE or of quantile
    if (threshold_change == TRUE){
      test = test_multiple_thresholds(ctrl_data, D_U_ctrl, simulation_p, multiple_values, iterations, unique_value)
    } else {
      test = test_multiple_quantiles(ctrl_data, D_U_ctrl, simulation_p, unique_value, iterations, multiple_values)
    }

    #For each value,
    for(value in 1:ncol(test$D)){
      #Computing of FP, TP, FN, TN
      results_simu = results_simulation(test$D[,value], test$U[,value], simulation_p)
      #Computing of FDR and TPR
      FDR = results_simu$FP / (results_simu$TP + results_simu$FP)
      TPR  = results_simu$TP / (results_simu$TP + results_simu$FN)
      results = rbind(results, c(p, colnames(test$D)[value], FDR, TPR))
    }
  }
  if (threshold_change == TRUE){
    colnames(results) = c("patient", "threshold", "FDR", "TPR")
  } else {
    colnames(results) = c("patient", "quantile", "FDR", "TPR")
  }
  return(results)
}


multiple_tests_sl0 = function(ctrl_data, D_U_ctrl, iterations, simulation, multiple_values, unique_value, threshold_change = TRUE,  quant_test = 0, factor_test = 1){
  results = c()
  #For each patient
  for(p in 1:ncol(simulation$initial_data)){
    print(c("Patient number",p))
    simulation_p = list(initial_data = simulation$initial_data[,p], simulated_data = simulation$simulated_data[,p])

    #The test is made for all the values, of threshold if threshold_change = TRUE or of quantile
    if (threshold_change == TRUE){
      test = test_multiple_thresholds_sl0(ctrl_data, D_U_ctrl, simulation_p, multiple_values, iterations, quant_test, factor_test)
    } else {
      test = test_multiple_quantiles(ctrl_data, D_U_ctrl, simulation_p, unique_value, iterations, multiple_values)
    }

    #For each value,
    for(value in 1:ncol(test$D)){
      #Computing of FP, TP, FN, TN
      results_simu = results_simulation(test$D[,value], test$U[,value], simulation_p)
      #Computing of FDR and TPR
      FDR = results_simu$FP / (results_simu$TP + results_simu$FP)
      TPR  = results_simu$TP / (results_simu$TP + results_simu$FN)
      results = rbind(results, c(p, colnames(test$D)[value], FDR, TPR))
    }
  }
  if (threshold_change == TRUE){
    colnames(results) = c("patient", "threshold", "FDR", "TPR")
  } else {
    colnames(results) = c("patient", "quantile", "FDR", "TPR")
  }
  return(results)
}

