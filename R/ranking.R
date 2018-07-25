# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#'Find low values
#'
#' This function detects genes with more than threshold percent of expression values under the min value.
#' NA expression values are not considered.
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
  idx_cancer = 1:ncol(cancer_data) + ncol(controls)

  #If any NA, we don't consider these values.
  if(anyNA(controls)){
    values0 = apply(binded_data, 1, function(l) {
      idx_cancer_sans_na = idx_cancer[!is.na(l[idx_cancer])]
      idx_ctrl_sans_na = idx_ctrl[!is.na(l[idx_ctrl])]
      #Computing the proportion of genes values  < min.
      percentCancer = sum(l[idx_cancer_sans_na] <= min) / length(idx_cancer_sans_na)
      percentCtrl = sum(l[idx_ctrl_sans_na] <= min) / length(idx_ctrl_sans_na)
      #If the proportion is above the threshold, we return true.
      if (percentCancer >= threshold & percentCtrl >= threshold){
        return(TRUE)
      } else {
        return(FALSE)
      }
    })
  } else {

    values0 = apply(binded_data, 1, function(l) {
      #Computing the proportion of genes values  < min.
      percentCancer = sum(l[idx_cancer] <= min)/length(idx_cancer)
      percentCtrl = sum(l[idx_ctrl] <= min)/length(idx_ctrl)
      #If the proportion is above the threshold, we return true.
      if (percentCancer >= threshold & percentCtrl >= threshold) {
        return(TRUE)
      }
      else {
        return(FALSE)
      }
    })
  }
  print(paste0(sum(values0), " genes have less than ", min, " counts in ", threshold*100, " % of the samples."))
  return(values0)
}


# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#'Find NA values
#'
#' This function detects probes with more than threshold percent of value undefined (NA).
#'
#'@param controls A matrix with datas to analyze.
#'@param cancer_data A matrix with other conditions datas to analyze.
#'@param threshold The maximum proportion of NA tolerated for each probe.
#'@param probes TRUE if you want to sort for probes, FALSE to sort for patients.
#'
#'@return This function returns a true false vector with true for the values to exclude.
#'
#'@example examples/ex_detect_zero_value.R
#'
#'@export

detect_na_value = function(controls, cancer_data, threshold, probes = TRUE) {

  if (probes == TRUE){

    nactrl = rowSums(is.na(controls))
    nacancer = rowSums(is.na(cancer_data))
    low_na = nactrl >= threshold*ncol(controls) & nacancer >= threshold*ncol(cancer_data)
    print(paste0(sum(low_na), " probes are NA in at least ", threshold*100, " % of the samples."))
    return(low_na)

  } else {

    nactrl = colSums(is.na(controls))
    nacancer = colSums(is.na(cancer_data))
    low_na = c(nactrl >= threshold*nrow(controls), nacancer >= threshold*nrow(cancer_data))
    print(paste0(sum(low_na), " patients have NA for at least ", threshold*100, " % of the probes."))
    return(low_na)

  }
}


# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Compute D and U list in control samples.
#'
#' For each gene, this function computes two lists. The D list, with genes less expressed in threshold% of controls,
#' and the U list with genes more expressed in threshold% of controls. These lists can be used in Penda test to compare
#' the gene rank in samples.
#'
#'@param controls A matrix with the gene expressions for each patient.
#'@param threshold The proportion of expression that must be in the conditions.
#'@param s_max The maximum number of down and up-expressed gene for each genes.
#'
#'@return This function returns a list of two logical matrices :
#'- the D matrix, with TRUE if the row gene has a lower expression than the column gene,
#'- the U Matrix with TRUE if the row gene has a higher expression than the column gene.
#'
#'@example examples/ex_compute_down_and_up_list.R
#'
#'@export

compute_down_and_up_list = function (controls, threshold, s_max = 50){

  print("Computing down and up-expressed genes")
  #Using DU_rcpp to compute down and up-expressed genes.
  if(anyNA(controls)){
    DU = penda::compute_DU_cppNA(controls, threshold, rowSums(is.na(controls)))
  } else {
    DU = penda::compute_DU_cpp(controls, threshold)
  }
  genes_U = unlist(DU$U)
  dimnames(genes_U) = list(DU$n, DU$n)
  genes_D = unlist(DU$D)
  dimnames(genes_D) = list(DU$n, DU$n)

  median_gene = apply(controls, 1, function(sonde){
    median(sonde, na.rm = TRUE)
  })
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
