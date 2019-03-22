# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#'Find low values
#'
#' This function detects genes with more than threshold percent of expression
#' values under the min value. NA expression values are not considered.
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
  print(paste0(sum(values0), " genes have less than ", min, " counts in ",
               threshold*100, " % of the samples."))
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
#'@example examples/ex_detect_na_value.R
#'
#'@export

detect_na_value = function(controls, cancer_data, threshold, probes = TRUE) {

  if (probes == TRUE){

    nactrl = rowSums(is.na(controls))
    nacancer = rowSums(is.na(cancer_data))
    low_na = nactrl >= threshold*ncol(controls) & nacancer >= threshold * ncol(cancer_data)
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
#' Making the Penda Dataset
#'
#' This function makes a dataset with controls and cancer_data pre-filtred and sorted by median in controls.
#'
#'@param controls A matrix with datas to analyze.
#'@param cancer_data A matrix with other conditions datas to analyze.
#'@param detectlowvalue If detectlowvalue is true, we apply the function "detect_zero_value".
#'@param detectNA If detectnavalue is true, we apply the function "detect_na_value" on probes and on patients.
#'@param threshold The maximum proportion of expression under val_min or NA tolerated for each gene.
#'@param val_min The minimum value accepted. If val_min is NA, we compute this value with mixtools.
#'
#'@return This function return a list with data_ctrl and data_case.
#'
#'@example examples/ex_make_dataset.R
#'
#'@export

make_dataset = function(controls, cancer_data, detectlowvalue = TRUE, detectNA = TRUE, threshold = 0.99, val_min = NA) {

  if(detectNA == TRUE){

    naprobes = penda::detect_na_value(controls, cancer_data, threshold, probes = TRUE)
    napatient = penda::detect_na_value(controls,  cancer_data, threshold, probes = FALSE)
    controls = controls[!naprobes, !napatient[1:ncol(controls)]]
    cancer_data = cancer_data[!naprobes, !napatient[(ncol(controls)+1):length(napatient)]]

  }
  if(detectlowvalue == TRUE){

    if(is.na(val_min)){
      print("Computing of the low threshold")
      all_data = c((log2(controls + 1)), (log2(cancer_data + 1)))
      mod = mixtools::normalmixEM(na.omit(all_data))
      mu = min(mod$mu)
      sig = mod$sigma[1]
      val_min = 2^qnorm(0.80, mu, sig)
      paste0("The low threshold is ", val_min)
    }
    low_values = penda::detect_zero_value(controls, cancer_data, threshold = threshold, min = val_min)
    controls = controls[!low_values, ]
    cancer_data = cancer_data[!low_values, ]
  }

  median_gene = apply(controls, 1, median, na.rm = TRUE)
  median_gene = sort(median_gene)
  controls = controls[names(median_gene), ]
  cancer_data = cancer_data[names(median_gene), ]

  return(list(data_ctrl = controls, data_case = cancer_data))
}


# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Compute L and H list in control samples.
#'
#' For each gene, this function computes two lists.
#' The L list, with genes with a Lower expression in threshold% of controls,
#' and the H list with genes with an Higher expression in threshold% of controls.
#' These lists can be used in Penda test to compare the gene rank in samples.
#'
#'@param controls A matrix with the gene expressions for each patient.
#'@param threshold The proportion of expression that must be in the conditions.
#'@param s_max The maximum number of L and H genes for each gene.
#'
#'@return This function returns a list of two logical matrices :
#'- the L matrix, with TRUE if the row gene has a lower expression than the column gene,
#'- the H Matrix with TRUE if the row gene has a higher expression than the column gene.
#'
#'@example examples/ex_compute_lower_and_higher_lists.R
#'
#'@export

compute_lower_and_higher_lists = function (controls, threshold, s_max = 50){

  print("Computing genes with lower and higher expression")
  median_gene = sort(apply(controls, 1, median, na.rm = TRUE))
  if( sum(names(median_gene) != rownames(controls)) !=0){
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

