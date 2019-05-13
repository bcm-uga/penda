# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Reorganization of results matrix after the multiple test.
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

DU_rearrangement = function(multiple_tests, simu_data, multiple_values){

  nvalues = length(multiple_values)
  #Definition of D and U matrices for all the conditions.
  D_simu = matrix(data = NA, nrow = length(simu_data)
                  , ncol = nvalues
                  , dimnames = list(names(simu_data), multiple_values))
  U_simu = matrix(data = NA, nrow = length(simu_data)
                  , ncol = nvalues
                  , dimnames = list(names(simu_data), multiple_values))

  #D and U are put in the matrices for all the conditions.
  for (r in seq_len(nvalues)){
    D_loc = unlist(multiple_tests[(r * 2 - 1)])
    U_loc = unlist(multiple_tests[r * 2])
    if(length(D_loc) > 1 & length(U_loc) > 1){
      D_simu[,r] = D_loc
      U_simu[,r] = U_loc
    }
  }
  return(list(D = D_simu, U = U_simu))
}


# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#'Test gene expression variations for multiple thresholds
#'
#' This function makes the Penda dysregulation test for different values of the threshold.
#'
#'@param simulation The simulation with the list of initial data $initial_data and modified data in $simulated_data.
#'@param threshold_values The vector of values of the threshold for the regulation test.
#'@param iterations The maximal number of iterations for the test. If the dysregulation list
#'no longer changes, iterations are stop before.
#'
#'@return This function returns a list of two logical matrices. The D matrix with for each threshold
#'TRUE for genes down-regulated, and the U matrix with TRUE for genes up-regulated.
#'
#'@example examples/ex_test_multiple_thresholds.R
#'
#'@export

test_multiple_thresholds = function(simulation, threshold_values, iterations){
  #Make the deregulation test for all threshold values.
  multiple_test = sapply(threshold_values, function(t){
    print(paste0("Threshold ", t))
    penda::sample_test(simulation$simulated_data, iterations, t)
  })
  sorted_test = DU_rearrangement(multiple_test, simulation$initial_data, threshold_values)
  return(sorted_test)
}


# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Test gene expression variations for multiple thresholds in a serie of samples
#'
#' This function makes for each patient the Penda test with different values for the treshold.
#' It computes FDR, TPR and FPR for each patient and each threshold value.
#'
#'@param controls A matrix with genes expressions in controls for all the patients.
#'@param L_H_list The list of lower and higher expressed genes matrices in the control.
#'@param iterations The maximal number of iterations for the test. If the dysregulation list no longer changes, iterations are stop before.
#'@param simulation The simulation with the list of initial data $initial_data and modified data in $simulated_data.
#'@param threshold_values A vector with different values to test for the threshold.
#'@param quant_test The quantile for the test. When L and H lists are void, we use the naive method.
#'@param factor_test The factor for the test. The limit L will be quant(gene)/factor, and the limit H quant(gene)*factor.
#'
#'@return This function returns a matrix with 5 columns : the patient number, the value of the threshold tested,
#'the FDR, the TPR and the FPR of the test.
#'
#'@example examples/ex_choose_threshold.R
#'
#'@export

choose_threshold = function(controls, L_H_list, iterations, simulation, threshold_values, quant_test = 0, factor_test = 1){

  controls <<- controls
  L_H_list <<- L_H_list
  quantile_genes <<- apply(controls, 1, quantile, c(quant_test,(1-quant_test)), na.rm = TRUE)
  quantile_genes[1,] = quantile_genes[1,] / factor_test
  quantile_genes[2,] = quantile_genes[2,] * factor_test

  results = c()

  #If simulation is a matrix, we do the test and compute errors for each patient
  if (is.matrix(simulation$initial_data)) {
    for(p in 1:ncol(simulation$initial_data)){
      print(paste0("Patient ",p))
      simulation_p = list(initial_data = simulation$initial_data[,p], simulated_data = simulation$simulated_data[,p])

      #The test is made for all the values of threshold.
      test = penda::test_multiple_thresholds(simulation_p, threshold_values, iterations)

      #For each value,
      for(value in 1:ncol(test$D)){
        #Computing of FP, TP, FN, TN
        results_simu = penda::results_simulation(test$D[,value], test$U[,value], simulation_p)
        #Computing of FDR, TPR, FPR
        FDR = results_simu$FP / (results_simu$TP + results_simu$FP)
        TPR  = results_simu$TP / (results_simu$TP + results_simu$FN)
        FPR = results_simu$FP / (results_simu$TN + results_simu$FP)
        results = rbind(results, c(p, colnames(test$D)[value], FDR, TPR, FPR, results_simu$TP, results_simu$FP, results_simu$TN, results_simu$FN))
      }
    }

  } else {
    test = penda::test_multiple_thresholds(simulation, threshold_values, iterations)
    for(value in 1:ncol(test$D)){
      #Computing of FP, TP, FN, TN
      results_simu = penda::results_simulation(test$D[,value], test$U[,value], simulation)
      #Computing of FDR, TPR, FPR
      FDR = results_simu$FP / (results_simu$TP + results_simu$FP)
      TPR  = results_simu$TP / (results_simu$TP + results_simu$FN)
      FPR = results_simu$FP / (results_simu$TN + results_simu$FP)
      results = rbind(results, c(1, colnames(test$D)[value], FDR, TPR, FPR, results_simu$TP, results_simu$FP, results_simu$TN, results_simu$FN))
    }
  }

  colnames(results) = c("patient", "threshold", "FDR", "TPR", "FPR", "TP", "FP", "TN", "FN")
  return(results)
}


# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Select the better threshold to reach a given FDR.
#'
#' This function use the result of choose_threshold to give the parameter which allows to reach
#' the given FDR. If there are multiple patients, it uses the median
#'
#'@param which_threshold The matrix result of choose_threshold, with patient, FDR and threshold.
#'@param FDR_max The maximum FDR wanted.
#'
#'@return This function returns a list of 2 parameters : the best threshold and the associated FDR.
#'
#'@example examples/ex_select_threshold_param.R
#'
#'@export

select_threshold_param = function(which_threshold, FDR_max = 0.05) {

  which_threshold = apply(which_threshold, 2, as.numeric)

  if(length(unique(which_threshold[,1])) > 1){
    results_threshold = c()
    for(value in unique(which_threshold[,2])){
      sum_value = colSums(which_threshold[which_threshold[,"threshold"] == value, ])
      results_threshold = rbind(results_threshold, c(value, sum_value[c(6, 7, 8, 9)]))
    }
  } else {
    results_threshold = which_threshold[, c(2, 6, 7, 8, 9)]
  }
  results_threshold = as.data.frame(results_threshold)

  results_threshold$FDR = results_threshold$FP / (results_threshold$TP + results_threshold$FP)
  results_threshold$TPR = results_threshold$TP / (results_threshold$TP + results_threshold$FN)
  results_threshold$FPR = results_threshold$FP / (results_threshold$TN + results_threshold$FP)

  small_FDR = which(results_threshold$FDR <= FDR_max)

  if(length(small_FDR) != 0 & length(small_FDR) != sum(is.na(small_FDR))){
    # If only one FDR under the FDR max,
    if (length(small_FDR) == 1){
      print(paste0("Best threshold is ",
                    results_threshold[small_FDR, 1],
                    " which has a FDR of ",
                    results_threshold[small_FDR, "FDR"]))
      return(list(threshold = results_threshold[small_FDR, 1],
                  FDR = results_threshold[small_FDR, "FDR"],
                  TPR = results_threshold[small_FDR, "TPR"],
                  FPR = results_threshold[small_FDR, "FPR"]
                  ))
    # If more than one FDR under the FDR max,
    } else {
      small_results = results_threshold[small_FDR, ]
      idx = which.max(small_results[, "TPR"])
      print(paste0("Best threshold is ", small_results[idx, 1]
               , " which has a FDR of ", small_results[idx, "FDR"]))
      return(list(threshold = small_results[idx, 1],
                  FDR = small_results[idx, "FDR"],
                  TPR = small_results[idx, "TPR"],
                  FPR = small_results[idx, "FPR"]
                  ))
    }
  } else {
    warning("Your FDR is not reachable, we return the threshold of the smallest FDR.")
    idx = which.min(results_threshold[,"FDR"])
    return(list(threshold = results_threshold[idx, 1],
                FDR = results_threshold[idx, "FDR"],
                TPR = results_threshold[idx, "TPR"],
                FPR = results_threshold[idx, "FPR"]
                ))
  }
}


# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#'Test gene expression variations with the quantile test for multiple quantile values in a serie of samples
#'
#' This function makes the quantile test for each patient with different values of quantile and factor.
#' for each set of parameters, it computes FDR, TPR and FPR.
#'
#'@param controls A matrix with genes expressions in controls for all the patients.
#'@param simulation The list of initial data $initial_data and modified data in $simulated_data
#'@param factor_values A vector with different values of factor. The down-regulated limit will be quantmin/factor, and the up-regulated limit quantmax*factor.
#'@param quantile_values A vector with different values to test for the quantile.
#'
#'@return This function returns a matrix with 5 columns : the value of the factor test,
#'the value of the quantile tested, the FDR, the TPR and the FPR of the test.
#'
#'@example examples/ex_choose_quantile.R
#'
#'@export

choose_quantile = function(controls, simulation, factor_values = c(1, 1.2, 1.4, 1.6, 1.8, 2), quantile_values = c(0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2)){

  results = c()
  for (f in factor_values){
    print(paste0("Facteur ", f))
    for (q in quantile_values){
      print(paste0("Quantile ", q))
      test = quantile_test(controls, simulation$simulated_data, q, f)
      results_simu = results_simulation(test$D, test$U, simulation)
      FDR = results_simu$FP / (results_simu$TP + results_simu$FP)
      TPR  = results_simu$TP / (results_simu$TP + results_simu$FN)
      FPR = results_simu$FP / (results_simu$TN + results_simu$FP)
      results = rbind(results, c(f, q, FDR, TPR, FPR))
    }
  }
  colnames(results) = c("factor", "quantile", "FDR", "TPR", "FPR")
  return(results)
}



# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Select the better factor and quantile to reach a given FDR.
#'
#' This function use the result of choose_quantile to give parameters which allow to reach
#' the given FDR with the best TPR.
#'
#'@param which_quantile The matrix result of choose_quantile, with factor, quantile, FDR and TPR.
#'@param FDR_max The maximum FDR wanted.
#'
#'@return This function returns a list of 4 parameters : the best factor, the best quantile,
#' and the associated FDR and TPR.
#'
#'@example examples/ex_select_quantile_param.R
#'
#'@export

select_quantile_param = function(which_quantile, FDR_max = 0.15){

#Compute the small FDR.
  small_FDR = which_quantile[which_quantile[,"FDR"] <= FDR_max & !is.na(which_quantile[,"FDR"]) ,]
  if(length(small_FDR) != 0){
    if (is.vector(small_FDR)){
      print(paste0("Best parameters are a factor of ", small_FDR[1],
                   " and a quantile of ", small_FDR[2],
                   " to obtain a FDR of ", small_FDR[3]))
      return(list(factor = small_FDR[1], quantile = small_FDR[2],
                  FDR = small_FDR[3], TPR = small_FDR[4]))

    } else {
      idx = which.max(small_FDR[,"TPR"])
      print(paste0("Best parameters are a factor of ", small_FDR[idx, 1],
                   " and a quantile of ", small_FDR[idx, 2],
                   " to obtain a FDR of ", small_FDR[idx, 3]))

      return(list(factor = small_FDR[idx, 1], quantile = small_FDR[idx, 2],
                  FDR = small_FDR[idx, 3], TPR = small_FDR[idx, 4]))
    }
  } else {
    warning("Your FDR is not reachable, we return the quantile and the factor of the smallest FDR.")
    idx = which.min(which_quantile[,"FDR"])
    return(list(factor = which_quantile[idx, 1], quantile = which_quantile[idx, 2],
                FDR = which_quantile[idx, 3], TPR = which_quantile[idx, 4]))
  }
}

