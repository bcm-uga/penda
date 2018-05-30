# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#'detect_zero_value
#'
#' This function detects genes with more than threshold percent of expression values under the min value.
#'
#'@param controls A matrix with datas to analyze.
#'@param cancer_data A matrix with other conditions datas to analyze.
#'@param threshold The maximum proportion of expression under min tolerated for each gene.
#'@param min The minimum value accepted.
#'
#'@return This function returns a true false vector with true for the values to exclude.
#'
#'@example examples/ex_detect_zero_value.R
#'
#'@export

detect_zero_value = function(controls, cancer_data, threshold, min = 0) {

  binded_data = cbind(controls, cancer_data)
  idx_ctrl = 1:ncol(controls)
  idx_lusc = 1:ncol(cancer_data) + ncol(controls)

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
  print(paste0(sum(values0), " genes have less than ", min, " counts in ", threshold*100, " % of the samples."))
  return(values0)
}


# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' compute_down_and_up_list
#'
#' This function ranks each gene and finds the genes which are more or less exprimed.
#' It's the same than "find_D_U_ctrl_size" but faster because uses rcpp to rank.
#'
#'@param controls A matrix with the gene expressions for each patient.
#'@param threshold The proportion of expression that must be in the conditions.
#'@param s_max The maximum number of down and up-expressed gene for each genes.
#'
#'@return This function returns a list of two logical matrices :
#'the D matrix, with TRUE if the row gene has a lower expression than the column gene,
#'and the U Matrix with TRUE if the row gene has a higher expression than the column gene.
#'
#'@example examples/ex_compute_down_and_up_list.R
#'
#'@export

compute_down_and_up_list = function (controls, threshold, s_max = 50){

  print("Computing down and up-expressed genes")
  #Using DU_rcpp to compute down and up-expressed genes.
  DU = penda::compute_DU_cpp(controls, threshold)
  genes_U = unlist(DU$U)
  dimnames(genes_U) = list(DU$n, DU$n)
  genes_D = unlist(DU$D)
  dimnames(genes_D) = list(DU$n, DU$n)

  median_gene = apply(controls, 1, median)

  print("Size restriction")

  #For each gene, if D or U list are too big, we select the closer to g.
  for (i in 1:ncol(genes_D)){
    d_genes = median_gene[genes_D[,i]==1]
    u_genes = median_gene[genes_U[,i]==1]
    if (length(d_genes) > s_max){
      sort_median = sort(d_genes)
      sort_median = sort_median[(length(d_genes) - (s_max-1)) : length(d_genes)]
      genes_D[,i] = FALSE
      genes_D[names(sort_median),i] = TRUE
    }
    if (length(u_genes) > s_max){
      sort_median = sort(u_genes)
      sort_median = sort_median[1 : s_max]
      genes_U[,i] = FALSE
      genes_U[names(sort_median),i] = TRUE
    }
  }
  gc()
  return(list(D = genes_D, U = genes_U))
}


