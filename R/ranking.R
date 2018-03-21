
#'detect_zero_value
#'
#' This function detects genes with more than threshold percent of expression values under the min value.
#'
#'@param ctrl_data A matrix with data to analyze
#'@param cancer_data A matrix with other data to analyze
#'@param threshold The maximum proportion of expression under min tolerated for each gene
#'@param min The minimum value accepted
#'
#'@return This function return a true false vector with true for the values to exclude.
#'
#'@example
#'null_values = detect_zero_value(ctrl_data, simu_data, threshold = 0.8, min = 10)
#'ctrl_data = ctrl_data[!null_values,]
#'simu_data = simu_data[!null_values,]

#'@export

detect_zero_value = function(ctrl_data, cancer_data, threshold, min = 0) {
  binded_data = cbind(ctrl_data, cancer_data)
  idx_ctrl = 1:ncol(ctrl_data)
  idx_lusc = 1:ncol(cancer_data) + ncol(ctrl_data)

  values0 = apply(binded_data, 1, function(l) {
    #Computing the proportion of genes values  < min
    percentLUSC = sum(l[idx_lusc] <= min) / length(idx_lusc)
    percentCtrl = sum(l[idx_ctrl] <= min) / length(idx_ctrl)
    #If the proportion is above the threshold, we return true
    if (percentLUSC >= threshold & percentCtrl >= threshold){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
}


#'find_D_U_ctrl
#'
#' This function ranks each gene and finds the genes which are more or less exprimed.
#'
#'@param ctrl_data A matrix with the gene expressions for each patient
#'@param quant A quantile to delimit the gene expression
#'@param factor A factor to delimit the study limit : between quantile min / factor and quantile max * factor
#'@param threshold The proportion of expression that must be in the conditions
#'
#'@return This function returns a list of two TRUE-FALSE matrix :
#'the D matrix, with TRUE if the row gene has a lower expression than the column gene,
#'and the U Matrix with TRUE if the row gene has a higher expression than the column gene.
#'
#'@example
#'find_D_U_ctrl(ctrl_data, quant = 0.001, factor = 4, threshold = 0.99)
#'
#'@export

find_D_U_ctrl = function (ctrl_data, quant, factor, threshold){

  transposee = t(ctrl_data)
  matrice_u = matrix( nrow = nrow(ctrl_data), ncol=nrow(ctrl_data)
                      , dimnames=list(rownames(ctrl_data),rownames(ctrl_data)))
  matrice_d = matrix( nrow = nrow(ctrl_data), ncol=nrow(ctrl_data)
                      , dimnames=list(rownames(ctrl_data),rownames(ctrl_data)))

  #For each gene,
  matrice_d = apply(ctrl_data, 1, function (g){

    quantile_gene = quantile(g, c(quant,(1-quant)))
    #d_genes are between (quantile min / factor) and gene expression
    d_genes = ((transposee < g) & (transposee >= round(quantile_gene[1] / factor)))
    sum_d_genes = apply(d_genes,2,sum)
    selected_d_genes = (sum_d_genes > (length(g) * threshold))
    return(selected_d_genes)

  })

  matrice_u = apply(ctrl_data, 1, function (g){

    quantile_gene = quantile(g, c(quant,(1-quant)))
    #u_genes are between gene expression and (quantile max * factor)
    u_genes = ((transposee > g) & (transposee <= round(factor * quantile_gene[2])))
    sum_u_genes = apply(u_genes,2,sum)
    selected_u_genes = (sum_u_genes > (length(g) * threshold))
    return(selected_u_genes)

  })
  return(list(D = matrice_d, U = matrice_u))
}


#'check_D_U
#'
#' This function checks if not exists genes both up and down expressed, and for each gene :
#' if the maximal down gene is less expressed in threshold% of cases,
#' if the minimal up gene is more expressed in threshold % of cases.
#'
#'@param D_U_ctrl The list of the two matrix D and U for each genes
#'@param ctrl_data The matrix of gene expression used to make D and U matrix
#'@param threshold The threshold used to make D and U matrix. It's the proportion of genes that
#'must meet the conditions
#'
#'@return This function returns an error message if D or U list don't mette the conditions. In normal case, the function return "OK" message.
#'
#'@example
#' D_U = find_D_U_ctrl(ctrl_data, quant = 0.001, factor = 4, threshold = 0.99)
#' check_D_U(D_U_ctrl = D_U, ctrl_data, threshold = 0.9)
#'
#'@export

check_D_U = function (D_U_ctrl, ctrl_data, threshold){
  #test if exists genes both down and up
  if (length(which(D_U_ctrl$D == 1 & D_U_ctrl$U == 1)) != 0){
    return( list(error = "Warning, genes should not be both down and up regulated "
                 , nb_errors = length(which(D_U_ctrl$D == 1 & D_U_ctrl$U == 1)))
            )
  }

  else {
    #test if maximal down gene is less expressed, and minimal up gene is more expressed
    pb = sapply(colnames(D_U_ctrl$D), function(gene){
      pb_d = FALSE
      pb_u = FALSE
      genes_names_D = names(which(D_U_ctrl$D[,gene]==1))
      genes_names_U = names(which(D_U_ctrl$U[,gene]==1))

      #Research of the maximal down-expressed gene
      if (length(genes_names_D)>1){
        mean_genes_D = apply(ctrl_data[genes_names_D,], 1, mean)
        max_D = names(which(mean_genes_D == max(mean_genes_D), useNames = TRUE))
      } else {
        max_D = genes_names_D
      }
      gene_max_D = ctrl_data[max_D,]

      #Research of the minimal up-expressed gene
      if (length(genes_names_U)>1){
        mean_genes_U = apply(ctrl_data[genes_names_U,], 1, mean)
        max_U = names(which(mean_genes_U == max(mean_genes_U), useNames = TRUE))
      } else {
        max_U = genes_names_U
      }
      gene_max_U = ctrl_data[max_U,]

      #Test if exist problem for down or up
      if (sum(gene_max_D < ctrl_data[gene,]) < (threshold*length(gene_max_D))){
        pb_d = TRUE
      }
      if (sum(gene_max_U > ctrl_data[gene,]) < (threshold*length(gene_max_U))){
        pb_u = TRUE
      }
      return(list(d=pb_d, u=pb_u))
    })
    if (sum(pb$d) > 0){
      return(list(error = "Warning, down genes sometimes overexpressed", nb_errors = sum(pb$d)))
    } else if (sum(pb$u) > 0){
       return(list(error = "Warning, up genes sometimes underexpressed", nb_errors = sum(pb$u)))
    } else {
        return ("OK")
    }
  }
}
