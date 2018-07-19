# Authors: Clémentine Decamps, UGA and Magali Richard
# clementine.decamps@univ-grenoble-alpes.fr
# magali.richard@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Test gene expression variation for a given gene
#'
#' This function makes the test of expression dysregulation of one gene for a patient compared to control.
#' It's an hybrid method, which use the rank method and the quantile method in addition when there is no up or down list.
#'
#'@param gene The name of the gene to analyze.
#'@param D_U_list The list of Down and Up-expressed genes matrices in the control.
#'@param sample The vector of genes expressions for one patient.
#'@param threshold If Dd/Du and Uu/ud are under this threshold, the expression not change.
#'@param controls The matrix with genes expressions in controls for all the patients.
#'@param quant_test Quant is used to make the quantile of the normal distribution of a gene when the gene has not D or U lists. A down-regulated gene is under the quantile, and an up-regulated is above.
#'@param factor_test The factor modifiy the quantile, the min quantile is divied by factor, and the max is multiplied.
#'
#'@return This function returns 0 if the gene expression has not changed,
#' 1 if the gene is up-regulated and -1 if the gene is down-regulated.
#'
#'@examples examples/ex_regulation_test.R
#'
#'@export

regulation_test = function(gene, D_U_list, sample, threshold, controls, quant_test = 0, factor_test = 1){

  #If the gene is NA in the sample, we can't compute the dysregulation.
  if (is.na(sample[gene])){
    return(0)
  } else {
    #Compute the quantile for the quantile method, if D or U does not exist.
    changement = numeric()
    quantile_gene = quantile(controls[gene,], c(quant_test,(1-quant_test)), na.rm = TRUE)
    quantile_gene[1] = quantile_gene[1] / factor_test
    quantile_gene[2] = quantile_gene[2] * factor_test

    D = sample[D_U_list$D[,gene] == TRUE]
    D = D[!is.na(D)]
    U = sample[D_U_list$U[,gene] == TRUE]
    U = U[!is.na(U)]

    Du = (D > sample[gene])
    Du = sum (Du)
    Ud = (U < sample[gene])
    Ud = sum (Ud)
    Dd = length(D)-Du
    Uu = length(U) - Ud

    #If exists both up and down-regulated lists,
    if (length(D) != 0 & length(U) != 0){
      if ((Ud / length(U) < threshold) & ((Du / length(D) < threshold))){
        changement = 0
      } else if ((Dd + Ud) < length(D)){
        changement = -1
      } else {
        changement = 1
      }
      return(changement)
    }
    #If exists only up-regulated list,
    else if (length(D) == 0 & length(U) != 0){
      if ((Ud / length(U) < threshold) & (sample[gene] >= quantile_gene[1])){
        changement = 0
      } else if (Ud / length(U) >= threshold){
        changement = 1
      } else {
        changement = -1
      }
      return(changement)
    }
    #If exists only down-regulated list,
    else if (length(D) != 0 & length(U) == 0){
      if ((Du / length(D) < threshold) & (sample[gene] < quantile_gene[2])){
        changement = 0
      } else if (Du / length(D)>= threshold){
        changement = -1
      } else {
        changement = 1
      }
      return(changement)
    }
    #If there is no regulation list,
    else if (sample[gene] < quantile_gene[1]){
      return(-1)
    } else if (sample[gene] > quantile_gene[2]){
      return(1)
    }  else {
      return(0)
    }
  }
}


# Authors: Clémentine Decamps, UGA and Magali Richard
# clementine.decamps@univ-grenoble-alpes.fr
# magali.richard@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Test gene expression variations for a given sample
#'
#' This function makes the test of the dysregulation of all genes for one patient.
#'
#'@param controls A matrix with genes expressions in controls for all the patients.
#'@param sample A vector with all the genes expressions for one patient.
#'@param iterations The maximal number of iterations for the test. If the dysregulation list
#'no longer changes, iterations are stop before.
#'@param D_U_list The list of Down and Up-expressed genes matrices in the control.
#'@param threshold The threshold for regulation_test.
#'@param quant_test The quantile for the test. When D and U lists are void, we use the naive method.
#'@param factor_test The factor for the test. The limit D will be quant(gene)/factor, and the limit U quant(gene)*factor.
#'
#'@return This function return a list with two vectors :
#'D, with TRUE for genes down-regulated
#'U, with TRUE for genes up-regulated.
#'
#'@example examples/ex_sample_test.R
#'
#'@export

sample_test = function (sample, controls, iterations, D_U_list, threshold, quant_test = 0, factor_test = 1){

  l1 = rep(FALSE, nrow(controls))
  l1n1 = rep(FALSE, nrow(controls))
  l1n2 = rep(FALSE, nrow(controls))
  print("Begining of iterations")

  #For each iteration
  for (i in 1:iterations){
    print(i)
    #Genes dysregulated at the previous iteration are removed of D_U.
    D_U_list_tmp = list()
    D_U_list_tmp$D = D_U_list$D * !abs(l1)
    D_U_list_tmp$U = D_U_list$U * !abs(l1)
    #The dysregulation of each gene is compute.
    l1 = sapply(names(sample), function(gene){
      expression = regulation_test(gene, D_U_list_tmp, sample, threshold, controls, quant_test = quant_test, factor_test = factor_test)
      return(expression)
    })
    #Stabilization is checked.
    if((sum(l1 == l1n1) == length(l1)) | (sum(l1 == l1n2) == length(l1))){
      print ("Stabilization of the dysregulated genes list.")
      break
    } else if (i == iterations) {
      print ("Maximum iterations without stabilization.")
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


# Author: Magali Richard, UGA
# magali.richrad@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Test gene expression variations in a serie of samples
#'
#' This function applies the function sample_test to a serie of samples.
#'
#'@param controls A matrix with genes expressions in controls for all the patients.
#'@param samples A vector or a matrix with all the genes expressions for each sample.
#'@param iterations The maximal number of iterations for the test. If the dysregulation list
#'no longer changes, iterations are stop before.
#'@param D_U_list The list of Down and Up-expressed genes matrices in the control.
#'@param threshold The threshold for regulation_test.
#'@param quant_test The quantile for the test. When D and U lists are void, we use the naive method.
#'@param factor_test The factor for the test. The limit D will be quant(gene)/factor, and the limit U quant(gene)*factor.
#'
#'@return This function return a list of two matrices (down_genes and up_genes)
#'
#'@example examples/ex_penda_test.R
#'
#'@export

penda_test = function(samples, controls, iterations, D_U_list, threshold, quant_test = 0, factor_test = 1){
if (is.null(dim(samples)[2])){
  print("compute DE list for one sample")
  res = penda::sample_test(sample = samples,
                            controls = controls,
                            threshold = threshold,
                            iterations =  iterations,
                            D_U_list =  D_U_list,
                            quant_test =  quant_test,
                           factor_test = factor_test)
  down_genes = res$D
  up_genes = res$U
} else if (dim(samples)[2] > 1){
  print(paste0("compute DE list for ", dim(samples)[2], " samples"))
  res = apply(samples, 2, penda::sample_test,
              controls = controls,
              threshold = threshold,
              iterations =  iterations,
              D_U_list =  D_U_list,
              quant_test =  quant_test,
              factor_test = factor_test)
  down_genes = sapply(res, "[[", "D")
  up_genes = sapply(res, "[[", "U")
} else {
  print("Error, samples should correspond to a matrix or a vector")
  down_genes = NULL
  up_genes = NULL
}
return(list(down_genes = down_genes, up_genes = up_genes))
}


# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Simple test of gene expression variations in a serie of samples
#'
#' This function makes the quantile test of the dysregulation of all genes. The patient's gene is compared at the gene expression
#' in controls sample.
#'
#'@param controls A matrix with genes expressions in controls for all the patients.
#'@param samples A vector or a matrix with all the genes expressions for each sample.
#'@param quant The quantile of the control genes expression, quantile(c(quant, 1-quant)).
#'@param factor The factor for the quantile. The D limit will be quantmin/factor, and the U limit quantmax*factor.
#'
#'@return This function return a list with two vectors or matrix :
#'D, with TRUE for genes down-regulated
#'U, with TRUE for genes up-regulated.
#'
#'@example examples/ex_quantile_test.R
#'
#'@export

quantile_test = function(controls, samples, quant, factor){

  quantile_gene = apply(controls, 1, quantile, c(quant,(1-quant)))
  quantile_gene[1,] = quantile_gene[1,] / factor
  quantile_gene[2,] = quantile_gene[2,] * factor

  D_matrix = samples < quantile_gene[1,]
  U_matrix = samples > quantile_gene[2,]
  return(list(D = D_matrix, U = U_matrix))
}
