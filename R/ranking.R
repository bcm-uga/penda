# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Compute L and H list in control samples.
#'
#' For each gene, this function computes two lists.
#' The L list, with genes with a Lower expression in `threshold*100` \% of controls,
#' and the H list with genes with an Higher expression in `threshold*100` \% of controls.
#' These lists can be used in Penda test to compare the gene rank in samples.
#'
#'@param controls A data matrix with the expressions of control samples.
#'@param threshold The proportion of expression that must be in the conditions.
#'@param s_max The maximum number of L and H genes for each gene.
#'
#'@return This function returns a list of two numeric matrices :
#'- the L matrix, with for each row gene the id of the L genes in columns,
#'- the H Matrix, with for each row gene the id of the H genes in columns.
#'
#'@example examples/ex_compute_lower_and_higher_lists.R
#'
#'@export

compute_lower_and_higher_lists = function (controls, threshold, s_max = 50){

  print("Computing genes with lower and higher expression")
  median_gene = sort(apply(controls, 1, median, na.rm = TRUE))
  if(sum(names(median_gene) != rownames(controls)) != 0){
    warning("Genes has to be ordered by their median in controls.")
  }

  #Using DU_rcpp to compute genes with lower and higher expression.
  if(anyNA(controls)){
    LH = penda::compute_LH_cppNA(controls, threshold, rowSums(is.na(controls)), s_max)
  } else {
    LH = penda::compute_LH_cpp(controls, threshold, s_max)
  }
  genes_L = unlist(LH$L)
  rownames(genes_L) = LH$n
  genes_H = unlist(LH$H)
  rownames(genes_H) = LH$n

  gc()
  return(list(L = genes_L, H = genes_H))
}

