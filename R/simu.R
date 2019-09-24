# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#'Simplified simulation of data deregulation.
#'
#'This function simulates the deregulation of expression (up and down) in the data.
#'
#'@param data The vector or matrix of original data to edit.
#'@param proportion The proportion of data to modify.
#'@param threshold The expression value under which data modification is +- modifier. Above the threshold, dysregulation is expression * or / factor.
#'@param modifier Under the threshold, the deregulation is to add or remove this number.
#'@param factor Below the threshold, the deregulation is to multiply or divide by this number.
#'
#'@return This function returns a list of three vectors or matrices :
#'initial_data, data before the modification
#'simulated_data, data with simulated deregulation
#'changes_idx, index of modifications
#'
#'@example examples/ex_simplified_simulation.R
#'
#'@export

simplified_simulation = function(data, proportion, threshold = 60, modifier = 30, factor = 4){

  simu_data = data
  #Sampling of deregulated genes
  perturb = sample(length(simu_data), proportion*length(simu_data))
  #For each deregulated gene,
  for(g in perturb){
    #Random deregulation (up or down)
    if (rnorm(1) < 0) {
      if (simu_data[g] < threshold) {
        simu_data[g] = simu_data[g] - modifier
        if (simu_data[g] <= 0){
          simu_data[g] = 0
        }
      } else {
        simu_data[g] = simu_data[g] / factor
      }
    } else if (simu_data[g] < threshold) {
      simu_data[g] = simu_data[g] + modifier
    } else {
      simu_data[g] = simu_data[g] * factor
    }
  }
  return(list(initial_data = data, simulated_data = simu_data, changes_idx = sort(perturb)))
}


# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Group genes with similar expression
#'
#' This function makes group of genes with similar expression. For each group,
#' the proportion of dysregulated genes in cancer and the extent of variation are computed.
#'
#'@param controls The matrix of gene expression in controls.
#'@param cancer_data The matrix of gene expression in cases (ex: cancer).
#'@param size_grp The size of the groups of genes.
#'@param quant The quantile of gene expression in control. Cancer genes outside this limit are considered dysregulated.
#'
#'@return This function returns a matrix of four columns :
#'Limit min, the first value of the group,
#'Limit max, the last value of the group,
#'Proportion of dysregulation, the proportion of cancer genes outside the quantile limit.
#'Deltas dysregul, a list with all the value of cancer genes - control genes.

group_genes = function(controls, cancer_data, size_grp = 100, quant = 0.05){

  results = c()
  #Sorting all controls in a vector.
  all_ctrl = as.vector(controls)
  names(all_ctrl) = rep(rownames(controls), times = ncol(controls))
  all_ctrl = sort(all_ctrl)
  nb_grp = floor(length(all_ctrl) / size_grp)

  #For each group,
  for(grp in 1:nb_grp){
    #We define limits and genes of the group.
    if (grp < nb_grp){
      limits = c(1 + (grp - 1) * size_grp, grp * size_grp)
    } else {
      limits = c(1 + (grp - 1) * size_grp, length(all_ctrl))
    }

    genes_grp = all_ctrl[limits[1]:limits[2]]
    delta_gen = c()
    delta_cancer = c()
    #For all genes of the group,
    for (n in 1:length(genes_grp)){
      g = genes_grp[n]
      gene_name = names(g)
      #Computing differences between the gene expression in control and this gene.
      delta = controls[gene_name, ] - g
      delta = delta[-which(delta == 0)[1]]
      delta_gen = rbind(delta_gen, delta)
      #Computing differences between the gene expression in cancer and this gene.
      delta = cancer_data[gene_name, ] - g
      delta_cancer = rbind(delta_cancer, delta)
    }

    #Computing the proportion of genes outside of limits.
    q = quantile(delta_gen, c(quant,(1 - quant)), na.rm = TRUE)
    prop = sum(delta_cancer < q[1] | delta_cancer > q[2]) / length(delta_cancer)
    delta_deregul = delta_cancer[delta_cancer < q[1] | delta_cancer > q[2]]
    results = rbind(results, c(all_ctrl[limits[1]], all_ctrl[limits[2]], prop, list(delta_deregul)))
  }
  colnames(results) = c("Limit min", "Limit max", "Proportion of dysregul", "Deltas dysregul")
  return(results)
}


# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#'Simulation of dysregulation that takes in account real datas
#'
#' This function simulated the dysregulation of datas, using the real distribution of genes, proportion of
#' dysregulation and difference between control and cancer.
#'
#'@param controls A matrix with genes expressions in controls for all the patients.
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
#'@example examples/ex_complex_simulation.R
#'
#'@export

complex_simulation = function(controls, cancer_data, data, size_grp = 100, quant = 0.05){
  simu_data = data
  print("Computing genes groups")
  group = group_genes(controls, cancer_data, size_grp, quant)
  limits = rbind(as.numeric(unlist(group[,1])), as.numeric(unlist(group[,2])))
  prop = unlist(group[,3])

  if(is.matrix(simu_data)){
    print(paste0("Simulating dysregulation of ", ncol(simu_data), " patients."))
    for (p in 1:ncol(simu_data)){
      for (g in 1:nrow(simu_data)){
        gene = simu_data[g, p]
        if (!is.na(gene)){
          sup_group = which(gene >= limits[1,])
          if(length(sup_group) > 0){
            group_gene = max(sup_group)
          } else {
            group_gene = 1
          }
          prop_gene = prop[group_gene]
          if (runif(1) <= prop_gene) {
            all_delta = unlist(group[,4][group_gene])
            delta = sample(all_delta, 1)
            count = 0
            while((simu_data[g,p] + delta) < 0){
              delta = sample(all_delta, 1)
              count = count + 1
              if(count > 500){
                break;
              }
            }
            simu_data[g, p] = simu_data[g,p] + delta
          }
        }
      }
    }
  } else if (is.vector(simu_data)){
    print(paste0("Simulating dysregulation of one patient."))
    for (g in 1:length(simu_data)){
      gene = simu_data[g]
      if (!is.na(gene)){
        sup_group = which(gene >= limits[1,])
        if(length(sup_group) > 0){
          group_gene = max(sup_group)
        } else {
          group_gene = 1
        }
        prop_gene = prop[group_gene]
        if (runif(1) <= prop_gene) {
          all_delta = unlist(group[,4][group_gene])
          delta = sample(all_delta, 1)
          count = 0
          while((simu_data[g] + delta) < 0){
            delta = sample(all_delta, 1)
            count = count + 1
            if(count > 500){
              break;
            }
          }
          simu_data[g] = simu_data[g] + delta
        }
      }
    }
  } else {
    stop("Simu data must be a vector or a matrix.")
  }
  return(list(initial_data = data, simulated_data = simu_data, changes_idx = which(data != simu_data)))
}

# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#'Compare the results of the test on the simulation to the reality.
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
#'@example examples/ex_results_simulation.R
#'
#'@export

results_simulation = function(D_matrix, U_matrix, simulation){

  initial_data = simulation$initial_data
  simu_data = simulation$simulated_data
  down = initial_data > simu_data
  up = initial_data < simu_data

  #True positive = real down results and real up results.
  TP = (length (which(down == TRUE & D_matrix == TRUE))
        + length (which(up == TRUE & U_matrix == TRUE)))

  #False positive = down or up results, but not changed in reality.
  FP = (length (which(down == 0 & D_matrix == 1))
        + length (which(up == 0 & U_matrix == 1)))

  #False negative = not changed results, but down or up in reality.
  FN = (length (which((up == 1 | down == 1) & D_matrix == 0 & U_matrix == 0)))

  #TN = real not changed results.
  TN = (length (which(down == 0 & up == 0 & D_matrix == 0 & U_matrix == 0)))

  return(list(TP = TP, FP = FP, FN = FN, TN = TN ))
}


# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Graphic representation of the simulation results
#'
#' This function computes the the FDR and the TPR and makes a ggplot of the results.
#'
#'@param results The list TP, FP, FN and TN return by "results_simulation"
#'
#'@return This function returns a barplot of results, FDR and TPR.
#'
#'@example examples/ex_draw_results.R
#'
#'@export

draw_results = function(results){
  FP = results$FP
  df= data.frame(type = c("TP","FP","FN","TN", "FDR", "TPR"), results=c(unlist(results), (results$FP / (results$TP + results$FP)), results$TP / (results$TP + results$FN)))
  ggplot2::ggplot(df, ggplot2::aes(x=type, y=results, fill=type)) + ggplot2::geom_bar(stat="identity") + ggplot2::geom_text(ggplot2::aes(label=results), vjust=1)
}
