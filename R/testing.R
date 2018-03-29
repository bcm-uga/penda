#'regulation_test
#'
#' This function makes the test of expression dysregulation for a patient compared to control.
#'
#'@param gene The name of the gene to analyze.
#'@param D_U_ctrl The list of Down and Up-expressed genes matrices in the control.
#'@param patient_genes The vector of genes expressions for one patient.
#'@param threshold If Dd/Du and Uu/ud are under this threshold, the expression not change
#'@param ctrl_data The matrix with genes expressions in controls for all the patients.
#'
#'@return This function return 0 if the gene expression has not changed,
#' 1 if the gene is up-regulated and -1 if the gene is down-regulated.
#'
#'@example
#'simulated_data = ctrl_data[,ncol(ctrl_data)]
#'simulated_data = simplified_simulation(simulated_data, fraction = 0.3, threshold = 60)
#'ctrl_data = ctrl_data[,-ncol(ctrl_data)]
#'D_U_ctrl = find_D_U_ctrl(ctrl_data, quant = 0.001, factor = 4, threshold = 0.99)
#'regulation_test(gene = rownames(ctrl_data)[1], D_U_ctrl, patient_genes = simulated_data$simulated_data, threshold = 0.1, ctrl_data)
#'
#'@export

regulation_test = function(gene, D_U_ctrl, patient_genes, threshold, ctrl_data, quant = 0, factor = 1){

  down_ctrl = D_U_ctrl$D[,gene]
  up_ctrl = D_U_ctrl$U[,gene]
  changement = numeric()
  quantile_gene = quantile(ctrl_data[gene,], c(quant,(1-quant)))
  quantile_gene[1] = quantile_gene[1] / factor
  quantile_gene[2] = quantile_gene[2] * factor

  Du = (patient_genes[down_ctrl == TRUE] > patient_genes[gene])
  Du = sum (Du)
  Ud = (patient_genes[up_ctrl == TRUE] < patient_genes[gene])
  Ud = sum (Ud)
  Dd = sum(down_ctrl) - Du
  Uu = sum(up_ctrl) - Ud

  #If exists both up and down-regulated lists
  if (sum(down_ctrl) != 0 & sum(up_ctrl) != 0){
    if ((Ud / sum(up_ctrl) < threshold) & ((Du / sum(down_ctrl) < threshold))){
      changement = 0
    } else if ((Dd + Ud) < sum(down_ctrl)){
      changement = -1
    } else {
      changement = 1
    }
    return(changement)
  }
  #If exists only up-regulated list
  else if (sum(down_ctrl) == 0 & sum(up_ctrl) != 0){
    if ((Ud / sum(up_ctrl) < threshold) & (patient_genes[gene] >= quantile_gene[1])){
      changement = 0
    } else if (Ud / sum(up_ctrl) > threshold){
      changement = 1
    } else {
      changement = -1
    }
    return(changement)
  }
  #if exists only down-regulated list
  else if (sum(down_ctrl) != 0 & sum(up_ctrl) == 0){
    if ((Du / sum(down_ctrl) < threshold) & (patient_genes[gene] < quantile_gene[2])){
      changement = 0
    } else if (Du / sum(down_ctrl) > threshold){
      changement = -1
    } else {
      changement = 1
    }
    return(changement)
  }
  #If there is no regulation list
  else if (patient_genes[gene] < quantile_gene[1]){
    return(-1)
  } else if (patient_genes[gene] > quantile_gene[2]){
    return(1)
  }  else {
    return(0)
  }
}

original_test = function(gene, D_U_ctrl, patient_genes, threshold, ctrl_data, quant = 0){

  down_ctrl = D_U_ctrl$D[,gene]
  up_ctrl = D_U_ctrl$U[,gene]
  changement = numeric()

  Du = (patient_genes[down_ctrl == TRUE] > patient_genes[gene])
  Du = sum (Du)
  Ud = (patient_genes[up_ctrl == TRUE] < patient_genes[gene])
  Ud = sum (Ud)
  Dd = sum(down_ctrl) - Du
  Uu = sum(up_ctrl) - Ud

  #If exists both up and down-regulated lists
  if (sum(down_ctrl) != 0 & sum(up_ctrl) != 0){
    if ((Ud / sum(up_ctrl) < threshold) & ((Du / sum(down_ctrl) < threshold))){
      changement = 0
    } else if ((Dd + Ud) < sum(down_ctrl)){
      changement = -1
    } else {
      changement = 1
    }
    return(changement)
  }
  #If exists only up-regulated list
  else if (sum(down_ctrl) == 0 & sum(up_ctrl) != 0){
    if ((Ud / sum(up_ctrl) < threshold)){
      changement = 0
    } else {
      changement = 1
    }
    return(changement)
  }
  #if exists only down-regulated list
  else if (sum(down_ctrl) != 0 & sum(up_ctrl) == 0){
    if ((Du / sum(down_ctrl) < threshold)){
      changement = 0
    } else {
      changement = -1
    }
    return(changement)
  }
  #If there is no regulation list
  else {
    return(0)
  }
}

#'step0
#'
#' This function makes the list of genes at the extremities of expression distribution.
#'
#'@param ctrl_data The matrix with genes expressions in controls for all the patients.
#'@param patient_genes The vector of genes expressions for one patient.
#'@param quant The quantile for the expression distribution limit.
#'
#'@return This function return a logical vector with TRUE for genes outside of the distribution limit.
#'
#'@example
#'simulated_data = ctrl_data[,ncol(ctrl_data)]
#'simulated_data = simplified_simulation(simulated_data, fraction = 0.4, threshold = 60)
#'ctrl_data = ctrl_data[,-ncol(ctrl_data)]
#'step0(ctrl_data, simulated_data, quant = 0.03)
#'
#'@export

step0 = function (ctrl_data, patient_genes, quant){
  quantile_gene = apply(ctrl_data, 1, quantile, c(quant,(1-quant)))
  outsiders = ((patient_genes < quantile_gene[1,]) | (patient_genes > quantile_gene[2,]))
  return(outsiders)
}


#'patient_test
#'
#' This function makes the test of the dysregulation of all genes for one patient.
#'@param ctrl_data A matrix with genes expressions in controls for all the patients.
#'@param patient_genes A vector with all the genes expressions for one patient.
#'@param quant_0 The quantile of the expression distribution for the step_0.
#'@param iterations The maximal number of iterations for the test. If the dysregulation list
#'no longer changes, iterations are stop before.
#'@param D_U_ctrl The list of Down and Up-expressed genes matrices in the control.
#'@param threshold The threshold for regulation_test
#'
#'@return This function return a list with two vector :
#'D, with TRUE for genes down-regulated
#'U, with TRUE for genes up-regulated.
#'
#'@example
#'simulated_data = ctrl_data[,ncol(ctrl_data)]
#'simulated_data = simplified_simulation(simulated_data, fraction = 0.3, threshold = 60)
#'ctrl_data = ctrl_data[,-ncol(ctrl_data)]
#'D_U_ctrl = find_D_U_ctrl(ctrl_data, quant = 0.001, factor = 4, threshold = 0.99)
#'patient_test(ctrl_data, simulated_data, quant_0 = 0.03, iterations = 10, D_U_ctrl, threshold = 0.2)
#'
#'@export

patient_test = function (ctrl_data, patient_genes, quant_0, iterations, D_U_ctrl, threshold, factor_test = 1){

  #Suspicious genes with an expression out of bounds are removed of D_U.
  l0 = step0(ctrl_data, patient_genes, quant_0)
  D_U_ctrl$D = D_U_ctrl$D * !l0
  D_U_ctrl$U = D_U_ctrl$U * !l0
  l1 = rep(FALSE, nrow(ctrl_data))
  l1n1 = rep(FALSE, nrow(ctrl_data))
  l1n2 = rep(FALSE, nrow(ctrl_data))
  print("Begining of iterations")
  #For each iteration
  for (i in 1:iterations){
    print(i)
    #Genes dysregulated at the previous iteration are removed of D_U.
    D_U_ctrl_tmp = list()
    D_U_ctrl_tmp$D = D_U_ctrl$D * !abs(l1)
    D_U_ctrl_tmp$U = D_U_ctrl$U * !abs(l1)
    #The dysregulation of each genes is compute.
    l1 = sapply(names(patient_genes), function(gene){
      expression = regulation_test(gene, D_U_ctrl_tmp, patient_genes, threshold, ctrl_data, factor = factor_test)
      return(expression)
    })
    if((sum(l1 == l1n1) == length(l1)) | (sum(l1 == l1n2) == length(l1))){
      print ("Stabilisation of deregulation")
      break
    } else {
      l1n2 = l1n1
      l1n1 = l1
    }
  }
  U_genes = (l1 == 1 | l1n1 == 1)
  D_genes = (l1 == -1 | l1n1 == -1)
  return(list(D = D_genes, U = U_genes))
  gc()
}


original_patient_test = function (ctrl_data, patient_genes, quant_0, iterations, D_U_ctrl, threshold){

  #Suspicious genes with an expression out of bounds are removed of D_U.
  l0 = step0(ctrl_data, patient_genes, quant_0)
  D_U_ctrl$D = D_U_ctrl$D * !l0
  D_U_ctrl$U = D_U_ctrl$U * !l0
  l1 = rep(FALSE, nrow(ctrl_data))
  l1n1 = rep(FALSE, nrow(ctrl_data))
  l1n2 = rep(FALSE, nrow(ctrl_data))
  print("Begining of iterations")
  #For each iteration
  for (i in 1:iterations){
    print(i)
    #Genes dysregulated at the previous iteration are removed of D_U.
    D_U_ctrl_tmp = list()
    D_U_ctrl_tmp$D = D_U_ctrl$D * !abs(l1)
    D_U_ctrl_tmp$U = D_U_ctrl$U * !abs(l1)
    #The dysregulation of each genes is compute.
    l1 = sapply(names(patient_genes), function(gene){
      expression = original_test(gene, D_U_ctrl_tmp, patient_genes, threshold, ctrl_data)
      return(expression)
    })
    if((sum(l1 == l1n1) == length(l1)) | (sum(l1 == l1n2) == length(l1))){
      print ("Stabilisation of deregulation")
      break
    } else {
      l1n2 = l1n1
      l1n1 = l1
    }
  }
  U_genes = (l1 == 1 | l1n1 == 1)
  D_genes = (l1 == -1 | l1n1 == -1)
  return(list(D = D_genes, U = U_genes))
  gc()
}

