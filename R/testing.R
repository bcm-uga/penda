# Authors: Clémentine Decamps, UGA and Magali Richard
# clementine.decamps@univ-grenoble-alpes.fr
# magali.richard@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Test gene expression variation for a given gene
#'
#' This function makes the test of expression dysregulation of one gene for a
#' patient compared to control. It's an hybrid method, which use the rank method
#' and the quantile method in addition when there is no significant lower or higher list.
#'
#'@param gene The name of the gene to analyze.
#'@param L_H_list The list of lower and higher expressed genes matrices in the control.
#'@param sample The vector of genes expressions for one patient.
#'@param threshold If Ll/Lh and Hh/Hl are under this threshold, the expression not change.
#'@param n_min The number minimal of genes in the L H list to switch to the quantile method.
#'
#'@return This function returns 0 if the gene expression has not changed,
#' 1 if the gene is up-regulated and -1 if the gene is down-regulated.
#'
#'@example examples/ex_regulation_test.R
#'
#'@export

regulation_test = function(gene, L_H_list, sample, threshold,  n_min = 10){

  #If the gene is NA in the sample, we can't compute the dysregulation.
  if (is.na(sample[gene])){
    return(0)
  } else {
    #Compute the quantile for the quantile method, if L or H does not exist.
    changement = numeric()
    quantile_gene = quantile_genes[ ,gene]

    L = sample[L_H_list$L[gene, ]]
    L = L[!is.na(L)]
    H = sample[L_H_list$H[gene, ]]
    H = H[!is.na(H)]

    Lh = sum(L > sample[gene])
    Hl = sum(H < sample[gene])
    Ll = length(L) - Lh
    Hh = length(H) - Hl

    #If exists both lower and higher expressed lists,
    if (length(L) >= n_min & length(H) >= n_min){
      if (xor((Hl / length(H) >= threshold), ((Lh / length(L) >= threshold)))){
        if (Lh / length(L) >= threshold){
          changement = -1
        } else {
          changement = 1
        }
      } else {
        changement = 0
      }
      return(changement)
    }
    #If exists only the higher expressed list,
    else if (length(L) < n_min & length(H) >= n_min){
      if (xor((Hl / length(H) >= threshold), (sample[gene] < quantile_gene[1]))){
        if (Hl / length(H) >= threshold){
          changement = 1
        } else {
          changement = -1
        }
      } else {
        changement = 0
      }
      return(changement)
    }
    #If exists only the lower expressed list,
    else if (length(L) >= n_min & length(H) < n_min){
      if (xor((Lh / length(L) >= threshold), (sample[gene] > quantile_gene[2]))){
        if (Lh / length(L) >= threshold){
          changement = -1
        } else {
          changement = 1
        }
      } else {
        changement = 0
      }
      return(changement)
    }
    #If there is no expression list,
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
#'@param sample A vector with all the genes expressions for one patient.
#'@param iterations The maximal number of iterations for the test. If the dysregulation list no longer changes, iterations are stop before.
#'@param threshold The threshold for regulation_test.
#'
#'@return This function return a list with two vectors :
#'D, with TRUE for genes down-regulated
#'U, with TRUE for genes up-regulated.
#'
#'@example examples/ex_sample_test.R
#'
#'@export

sample_test = function (sample, iterations, threshold){

  l1 = rep(FALSE, length(sample))
  l1n1 = rep(FALSE, length(sample))
  l1n2 = rep(FALSE, length(sample))
  print("Begining of iterations")

  #For each iteration
  for (i in seq_len(iterations)){
    print(paste0("Iteration ", i))
    #Genes dysregulated at the previous iteration are removed of L_H.
    L_H_list_tmp = L_H_list
    id_l1 = which(l1 != 0)
    L_H_list_tmp$L = apply(L_H_list_tmp$L, 2, function(g){
      g[g %in% id_l1] = 0
      return(g)
    })
    L_H_list_tmp$H = apply(L_H_list_tmp$H, 2, function(g){
      g[g %in% id_l1] = 0
      return(g)
    })

    #The dysregulation of each gene is compute.
    l1 = sapply(names(sample), function(gene){
      expression = penda::regulation_test(gene, L_H_list_tmp, sample, threshold)
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
  lf = l1 + l1n1
  U_genes = lf > 0
  D_genes = lf < 0
  return(list(D = D_genes, U = U_genes))
  gc()
}


# Authors: Magali Richard and Clémentine Decamps, UGA
# magali.richard@univ-grenoble-alpes.fr
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Test gene expression variations in a serie of samples
#'
#' This function applies the function sample_test to a serie of samples.
#'
#'@param samples A vector or a matrix with all the genes expressions for each sample.
#'@param controls A matrix with genes expressions in controls for all the patients.
#'@param iterations The maximal number of iterations for the test. If the dysregulation list no longer changes, iterations are stop before.
#'@param L_H_list The list of lower and higher expressed genes matrices in the control.
#'@param threshold If Ll/Lh and Hh/Hl are under this threshold, the expression not change.
#'@param quant_test The quantile for the test. When L and H lists are void, we use the naive method.
#'@param factor_test The factor for the test. The limit L will be quant(gene)/factor, and the limit H quant(gene)*factor.
#'
#'@return This function return a list of two matrices (down_genes and up_genes)
#'
#'@example examples/ex_penda_test.R
#'
#'@export

penda_test = function(samples, controls, iterations, L_H_list, threshold,
                      quant_test = 0, factor_test = 1){

  controls <<- controls
  L_H_list <<- L_H_list
  quantile_genes <<- apply(controls, 1, quantile,
                           c(quant_test, (1 - quant_test)), na.rm = TRUE)
  quantile_genes[1,] = quantile_genes[1, ] / factor_test
  quantile_genes[2,] = quantile_genes[2, ] * factor_test

  if (is.null(dim(samples)[2])){
    print("Compute DE list for one sample")
    res = penda::sample_test(sample = samples,
                         threshold = threshold,
                         iterations =  iterations)
    down_genes = res$D
    up_genes = res$U
  } else if (dim(samples)[2] > 1){
    print(paste0("Compute DE list for ", dim(samples)[2], " samples"))
    nbpatient = 0
    res = apply(samples, 2, function(s){
      nbpatient <<- nbpatient + 1
      print(paste0("Patient ", nbpatient, "/", dim(samples)[2]))
      penda::sample_test(s, threshold = threshold, iterations =  iterations)
    })
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
#' This function makes the quantile test of the dysregulation of all genes.
#' The patient's gene is compared at the gene expression in controls sample.
#'
#'@param controls A matrix with genes expressions in controls for all the patients.
#'@param samples A vector or a matrix with all the genes expressions for each sample.
#'@param quant The quantile of the control genes expression, quantile(c(quant, 1-quant)).
#'@param factor The factor for the quantile. The down-regulated limit will be quantmin/factor, and the up-regulated limit quantmax*factor.
#'
#'@return This function return a list with two vectors or matrix :
#'D, with TRUE for genes down-regulated
#'U, with TRUE for genes up-regulated.
#'
#'@example examples/ex_quantile_test.R
#'
#'@export

quantile_test = function(controls, samples, quant, factor){

  quantile_gene = apply(controls, 1, quantile, c(quant,(1-quant)), na.rm = TRUE)
  quantile_gene[1,] = quantile_gene[1,] / factor
  quantile_gene[2,] = quantile_gene[2,] * factor

  D_matrix = samples < quantile_gene[1,]
  U_matrix = samples > quantile_gene[2,]
  return(list(D = D_matrix, U = U_matrix))
}
