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
#'@example simu = simplified_simulation(simu_data, fraction = 0.3, threshold = 60)
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


#'group_genes
#'
#' This function makes group of genes with similar expression. For each group, the function computes the
#' proportion of dysregulated genes in cancer and the difference delta.
#'
#'@param ctrl_data A matrix with genes expressions in controls for all the patients.
#'@param cancer_data A matrix with dysregulated genes expressions for all the patients.
#'@param size_grp The size of each group of genes.
#'@param quant The quantile of gene expression in control. Cancer genes outside this limit are considerd dysregulated.
#'
#'@return This function returns a matrix of four columns :
#'Limit min, the first value of the group,
#'Limit max, the last value of the group,
#'Proportion of dysregulation, the proportion of cancer genes outside the quantile limit.
#'Deltas dysregul, a list with all the value of cancer genes - control genes.
#'
#'@example group = group_genes(ctrl_data, cancer_data)
#'

group_genes = function(ctrl_data, cancer_data, size_grp = 100, quant = 0.05){

  results = c()
  #Vector with all the control genes sorted
  all_ctrl = as.vector(ctrl_data)
  names(all_ctrl) = rep(rownames(ctrl_data), times = ncol(ctrl_data))
  all_ctrl = sort(all_ctrl)
  nb_grp = floor(length(all_ctrl) / size_grp)

  #For each group
  for (grp in 1:nb_grp){
    #We define limits and genes of the group
    if (grp < nb_grp){
      limits = c(1 + (grp - 1) * size_grp, grp * size_grp)
    } else {
      limits = c(1 + (grp - 1) * size_grp, length(all_ctrl))
    }
    genes_grp = all_ctrl[limits[1] : limits[2]]
    delta_gen = c()
    delta_cancer = c()

    #For all genes of the group
    for (g in 1:length(genes_grp)){
      gene_name = names(genes_grp)[g]
      #Difference between the gene expression in control and this gene
      delta = ctrl_data[gene_name, ] - rep(genes_grp[g], times = ncol(ctrl_data))
      delta = delta[-which(delta==0)[1]]
      delta_gen = rbind(delta_gen, delta)
      #Difference between the gene expression in cancer and this gene
      delta = cancer_data[gene_name, ] - rep(genes_grp[g], times = ncol(cancer_data))
      delta_cancer = rbind(delta_cancer, delta)
    }

    #We compute the proportion of genes outside of limits
    q = quantile(delta_gen, c(quant,(1-quant)))
    prop = length(which(delta_cancer < q[1] | delta_cancer > q[2])) / length(delta_cancer)
    delta_deregul = delta_cancer[which(delta_cancer < q[1] | delta_cancer > q[2])]
    results = rbind(results, c(all_ctrl[limits[1]], all_ctrl[limits[2]], prop, list(delta_deregul)))
  }
  colnames(results) = c("Limit min", "Limit max", "Proportion of dysregul", "Deltas dysregul")
  return(results)
}


#'complex_simulation
#'
#' This function simulated the dysregulation of datas, using the real distribution of genes, proportion of
#' dysregulation and difference between control and cancer.
#'
#'@param ctrl_data A matrix with genes expressions in controls for all the patients.
#'@param cancer_data A matrix with dysregulated genes expressions for all the patients.
#'@param data vector or matrix of original data to edit.
#'@param size_grp The size of each group of genes for the grouping.
#'@param quant The quantile of gene expression for the grouping.
#'
#'@return This function returns a list of three vectors or matrices :
#'initial_data, datas before the modification
#'simulated_data, datas with up and down modifications
#'changes_idx, index of datas with modifications
#'
#'@example simu = complex_simulation(ctrl_data, cancer_data, simu_data)
#'
#'@export

complex_simulation = function(ctrl_data, cancer_data, data, size_grp = 100, quant = 0.05){
  simu_data = data
  print("Computing genes groups")
  group = group_genes(ctrl_data, cancer_data, size_grp, quant)
  limits = rbind(as.numeric(unlist(group[,1])), as.numeric(unlist(group[,2])))
  prop = unlist(group[,3])

  print ("Simulating dysregulation")
  for (p in 1:ncol(simu_data)){
    for (g in 1:nrow(simu_data)){
      gene = simu_data[g, p]
      group_gene = max(which(gene >= limits[1,]))
      prop_gene = prop[group_gene]
      if (runif(1) <= prop_gene) {
        all_delta = unlist(group[,4][group_gene])
        all_delta = all_delta > ()
        delta = sample(all_delta, 1)
        while((simu_data[g,p] + delta) < 0){
          delta = sample(all_delta, 1)
        }
        simu_data[g, p] = simu_data[g,p] + delta

      }
    }
  }
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
#'@examples
#' D_U = find_D_U_ctrl(ctrl_data, quant = 0.001, factor = 4, threshold = 0.99)
#' simulation = simplified_simulation(simu_data, fraction = 0.3, threshold = 60)
#' results_simulation(D_matrix = D_U$D, U_matrix = D_U$U, initial_data = simulation$initial_data, simulated_data = simulation$simulated_data)
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
#'@examples
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

