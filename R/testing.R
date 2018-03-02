#'regulation_test
#'
#' This function makes the test of expression dysregulation for a patient compared to control
#'
#'@param gene The name of the gene to analyze
#'@param D_U_control The list of Down and Up-expressed genes matrix in the control
#'@param cancer_data A matrix of gene expression in the cancer
#'@param threshold If Dd/Du and Uu/ud are under this threshold, the expression not change
#'
#'@return This function return 0 if the gene expression has not changed,
#' 1 if the gene is up-regulated and -1 if the gene is down-regulated
#'
#'@examples
#'data_regulate = ctrl_data[,ncol(ctrl_data)]
#'data_regulate = dysreg_simulation(data_regulate, fraction = 0.3, threshold = 60)
#'ctrl_data = ctrl_data[,-ncol(ctrl_data)]
#'D_U_ctrl = find_D_U_ctrl(ctrl_data, quant = 0.01, factor = 30, threshold = 0.8)
#'regulation_test(rownames(ctrl_data)[1], D_U_ctrl, data_regulate, 0.1)
#'
#'@export
#'

regulation_test = function(gene, D_U_control, cancer_data, threshold){

  down_ctrl = D_U_control$D[,gene]
  up_ctrl = D_U_control$U[,gene]
  changement = numeric()

  Du = (cancer_data[down_ctrl==TRUE] > cancer_data[gene])
  Du = sum (Du)
  Ud = (cancer_data[up_ctrl==TRUE] < cancer_data[gene])
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
    if (Ud / sum(up_ctrl) < threshold){
      changement = 0
    } else {
      changement = 1
    }
    return(changement)
  }
  #if exists only down-regulated list
  else if (sum(down_ctrl) != 0 & sum(up_ctrl) == 0){
    if (Du / sum(down_ctrl) < threshold){
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
#'@param genes_data a matrix with genes expressions in controls for all the patients
#'@param patient_genes all the genes expression for one patient
#'@param quant the quantile for the expression distribution
#'
#'@return This function return a vector with all the genes outside of the distribution limit
#'
#'@examples
#'data_regulate = ctrl_data[,ncol(ctrl_data)]
#'data_regulate = dysreg_simulation(data_regulate, fraction = 0.3, threshold = 60)
#'ctrl_data = ctrl_data[,-ncol(ctrl_data)]
#'step0(ctrl_data, data_regulate, quant = 0.05)
#'
#'@export

step0 = function (genes_data, patient_genes, quant){
  quantile_gene = apply(genes_data, 1, quantile, c(quant,(1-quant)))
  outsiders = ((patient_genes < quantile_gene[1,]) | (patient_genes > quantile_gene[2,]))
  return(which(outsiders==1))
}

#'test_individu
#'
#' This function makes the test of the dysregulation of all genes for one patient.
#'@param genes_data A matrix with genes expressions in controls for all the patients
#'@param patient_genes All the genes expression for one patient
#'@param quant_0 The quantile for the expression distribution for the step_0
#'@param iterations The number of iterations for the test
#'@param D_U_ctrl The list of Down and Up-expressed genes matrix in the control
#'@param threshold The threshold for regulation_test
#'
#'@return This function return a list with two vector
#'
#'@examples
#'data_regulate = ctrl_data[,ncol(ctrl_data)]
#'data_regulate = dysreg_simulation(data_regulate, fraction = 0.3, threshold = 60)
#'ctrl_data = ctrl_data[,-ncol(ctrl_data)]
#'D_U_ctrl = find_D_U_ctrl(ctrl_data, quant = 0.01, factor = 30, threshold = 0.8)
#'test_individu(ctrl_data, data_regulate, quant_0 = 0.05, iterations = 3, D_U_ctrl, threshold = 0.1)
#'
#'@export


test_individu = function (genes_data, patient_genes, quant_0, iterations, D_U_ctrl, threshold){
  #Suspicious genes are removed
  l0 = rep(TRUE, nrow(genes_data))
  idx_0 = step0(genes_data, patient_genes, quant_0)
  l0[idx_0] = FALSE

  D_U_ctrl$D = D_U_ctrl$D * l0
  D_U_ctrl$U = D_U_ctrl$U * l0

  l1 = rep(FALSE, nrow(genes_data))
  #For each iteration
  for (i in 1:iterations){
    D_U_ctrl_tmp = list()
    D_U_ctrl_tmp$D = D_U_ctrl$D * !abs(l1)
    D_U_ctrl_tmp$U = D_U_ctrl$U * !abs(l1)
    #For each gene
    l1 = sapply(names(patient_genes), function(gene){
      expression = regulation_test(gene, D_U_ctrl_tmp, patient_genes, threshold)
      return(expression)
    })
  }
  U_genes = (l1 == 1)
  D_genes = (l1 == -1)
  return(list(D = D_genes, U = U_genes))
  gc()
}
