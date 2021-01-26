# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Making the Penda Dataset for one control samples
#'
#' This function makes the Penda dataset with control and cases
#' sorted by value in control.
#'
#'@param controls The vector with reference data to analyze (ex: control sample).
#'@param data_case The matrix with datas to analyze (ex: tumors samples).
#'
#'@return This function return a list with preprocessed data_ctrl and data_case, and
#'the vector `info` with the different parameters.
#'
#'@example examples/ex_make_dataset_1ctrl.R
#'
#'@export

make_dataset_1ctrl = function(control, data_case){
  resume = c(length(control), 1, ncol(data_case))
  names(resume) = c("init_nb_genes", "init_nb_ctrls", "init_nb_cases")
  median_gene = sort(control)
  control = control[names(median_gene)]
  data_case = data_case[names(median_gene), ]
  return(list(data_ctrl = control, data_case = data_case,
              info = round(resume, 2)))
}

# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Compute L and H list in control sample.
#'
#' For each gene, this function computes two lists.
#' The L list, with s_max genes with the closer Lower expression in control,
#' and the H list with s_max genes with the closer Higher expression in control.
#' These lists can be used in Penda test to compare the gene rank in samples.
#'
#'@param control A vector with the expressions in the control sample.
#'@param s_max The maximum number of L and H genes for each gene.
#'
#'@return This function returns a list of two numeric matrices :
#'- the L matrix, with for each row gene the id of the L genes in columns,
#'- the H Matrix, with for each row gene the id of the H genes in columns.
#'
#'@example examples/ex_compute_lower_and_higher_lists_1ctrl.R
#'
#'@export

compute_lower_and_higher_lists_1ctrl = function (control, s_max = 500){

  print("Computing genes with lower and higher expression")
  median_gene = sort(control)
  if(sum(names(median_gene) != names(control)) != 0){
    warning("Genes has to be ordered by their value in control.")
  }

  controls = cbind(control, control)
  threshold = 1

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

# Authors: Clémentine Decamps, UGA and Magali Richard
# clementine.decamps@univ-grenoble-alpes.fr
# magali.richard@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Test gene expression variation for a given gene
#'
#' This function makes the test of expression dysregulation of one gene for a
#' patient compared to control. It's an adapted method for the datasets with only
#' 1 control which use only the rank method. If you want the normal Penda method used penda::regulation_test
#'
#'@param gene The name of the gene to analyze.
#'@param L_H_list The list of lower and higher expressed genes matrices in the control.
#'@param sample The vector of genes expressions for one patient.
#'@param threshold If Ll/Lh and Hh/Hl are under this threshold, the expression not change.
#'
#'@return This function returns 0 if the gene expression has not changed,
#' 1 if the gene is up-regulated and -1 if the gene is down-regulated.
#'
#'@example examples/ex_regulation_test_1ctrl.R
#'
#'@export

regulation_test_1ctrl = function(gene, L_H_list, sample, threshold){
  #If the gene is NA in the sample, we can't compute the dysregulation.
  if (is.na(sample[gene])){
    return(0)
  } else {

    changement = numeric()
    L = sample[L_H_list$L[gene, ]]
    L = L[!is.na(L)]
    H = sample[L_H_list$H[gene, ]]
    H = H[!is.na(H)]

    Lh = sum(L > sample[gene])
    Hl = sum(H < sample[gene])
    Ll = length(L) - Lh
    Hh = length(H) - Hl

    #If exists both lower and higher expressed lists,
    if (length(L) != 0 & length(H) != 0){
      if ((Hl / length(H) < threshold) & ((Lh / length(L) < threshold))){
        changement = 0
      } else if ((Ll + Hl) < length(L)){
        changement = -1
      } else {
        changement = 1
      }
      return(changement)
    }
    #If exists only the higher expressed list,
    else if (length(L) == 0 & length(H) != 0){
      if (Hl / length(H) >= threshold){
        changement = 1
      } else {
        changement = 0
      }
      return(changement)
    }
    #If exists only the lower expressed list,
    else if (length(L) != 0 & length(H) == 0){
      if (Lh / length(L) >= threshold){
        changement = -1
      } else {
        changement = 0
      }
      return(changement)
    }
    #If there is no expression list,
    else{
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
#' It's the adapted method for the datasets with only 1 control,
#' if you want the normal Penda method used penda::sample_test
#'
#'
#'@param sample A vector with all the genes expressions for one patient.
#'@param iterations The maximal number of iterations for the test. If the dysregulation list no longer changes, iterations are stop before.
#'@param threshold The threshold for regulation_test.
#'
#'@return This function return a list with two vectors :
#'D, with TRUE for genes down-regulated
#'U, with TRUE for genes up-regulated.
#'
#'@example examples/ex_sample_test_1ctrl.R
#'
#'@export
sample_test_1ctrl = function (sample, iterations, threshold){

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
      expression = regulation_test_1ctrl(gene, L_H_list_tmp, sample, threshold)
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

# Authors: Magali Richard and Clémentine Decamps, UGA
# magali.richard@univ-grenoble-alpes.fr
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Test gene expression variations in a serie of samples
#'
#' This function applies the function sample_test_1ctrl to a serie of samples.
#' It's a special version of the Penda method for the datasets with only 1 control,
#' if you want the normal Penda method used penda::penda_test
#'
#'@param samples A vector or a matrix with all the genes expressions for each sample.
#'@param iterations The maximal number of iterations for the test. If the dysregulation list no longer changes, iterations are stop before.
#'@param L_H_list The list of lower and higher expressed genes matrices in the control.
#'@param threshold If Ll/Lh and Hh/Hl are under this threshold, the expression not change.
#'
#'@return This function return a list of two matrices (down_genes and up_genes)
#'
#'@example examples/ex_penda_test_1ctrl.R
#'
#'@export

penda_test_1ctrl = function(samples, iterations, L_H_list, threshold){
  L_H_list <<- L_H_list
  if (is.null(dim(samples)[2])){
    print("Compute DE list for one sample")
    res = sample_test_1ctrl(sample = samples,
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
      sample_test_1ctrl(s, threshold = threshold, iterations =  iterations)
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
