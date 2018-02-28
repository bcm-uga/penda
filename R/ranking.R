
#'detect_zero_value
#'
#' This function detects genes with more than threshold percent of expression values equals to zero.
#'
#'@param ctrl_data A matrix with data to analyze
#'@param cancer_data A matrix with other data to analyze
#'@param threshold The maximum proportion of zero tolerated for each gene
#'
#'@return This function return a true false vector with true for the values to exclude.
#'
#'@examples
#'null_values = detect_zero_value(ctrl_data, simu_data, threshold = 0.8)
#'ctrl_data = ctrl_data[!null_values,]
#'simu_data = simu_data[!null_values,]

#'@export

detect_zero_value = function(ctrl_data, cancer_data, threshold) {
  binded_data = cbind(ctrl_data, cancer_data)
  idx_ctrl = 1:ncol(ctrl_data)
  idx_lusc = 1:ncol(cancer_data) + ncol(ctrl_data)

  values0 = apply(binded_data, 1, function(l) {
    #Computing the proportion of genes values  > 0
    percentLUSC = sum(l[idx_lusc] > 0) / length(idx_lusc)
    percentCtrl = sum(l[idx_ctrl] > 0) / length(idx_ctrl)
    #If these values are under the threshold, we return true
    if (percentLUSC < threshold & percentCtrl < threshold){
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
#'@param genes_data A matrix with the gene expressions for each patient
#'@param quant A quantile to delimit the gene expression
#'@param factor A factor to delimit the study limit : between quantile min / factor and quantile max * factor
#'@param threshold The proportion of expression that must be in the conditions
#'
#'@return This function returns a list of two TRUE-FALSE matrix :
#'the D matrix, with TRUE if the row gene has a lower expression than the column gene,
#'and the U Matrix with TRUE if the row gene has a higher expression than the column gene.
#'
#'@examples
#'find_D_U_ctrl(ctrl_data, quant = 0.01, factor = 30, threshold = 0.8)

#'
#'@export

find_D_U_ctrl = function (genes_data, quant, factor, threshold){

  transposee = t(genes_data)
  matrice_u = matrix( nrow = nrow(genes_data), ncol=nrow(genes_data), dimnames=list(rownames(genes_data),rownames(genes_data)))
  matrice_d = matrix( nrow = nrow(genes_data), ncol=nrow(genes_data), dimnames=list(rownames(genes_data),rownames(genes_data)))

  #For each gene,
  matrice_d = apply(genes_data, 1, function (g){

    quantile_gene = quantile(g, c(quant,(1-quant)))
    #d_genes are between (quantile min / factor) and gene expression
    d_genes = ((transposee < g) & (transposee >= round(quantile_gene[1] / factor)))
    sum_d_genes = apply(d_genes,2,sum)
    selected_d_genes = (sum_d_genes >= (length(g) * threshold))
    return(selected_d_genes)

  })

  matrice_u = apply(genes_data, 1, function (g){

    quantile_gene = quantile(g, c(quant,(1-quant)))
    #u_genes are between gene expression and (quantile max * factor)
    u_genes = ((transposee > g) & (transposee <= round(factor * quantile_gene[2])))
    sum_u_genes = apply(u_genes,2,sum)
    selected_u_genes = (sum_u_genes >= (length(g) * threshold))
    return(selected_u_genes)

  })
  return(list(D = matrice_d, U = matrice_u))
}
