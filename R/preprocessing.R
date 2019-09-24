# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Find low values
#'
#' This function detects genes with more than `threshold` percent of expression
#' values under the min value. NA expression values are not considered.
#'
#'@param controls The first matrix with datas to analyze (ex: controls samples).
#'@param data_case The second matrix with datas to analyze (ex: tumors samples).
#'@param threshold The maximum proportion of expression under `min` tolerated for each gene.
#'@param min The minimum value accepted.
#'
#'@return This function returns a True/False vector with True for the values to exclude.
#'
#'@example examples/ex_detect_zero_value.R
#'
#'@export

detect_zero_value = function(controls, data_case, threshold, min = 0) {

  #Computing for each gene the number of samples under the minimum value
  n_ctrl = rowSums(controls <= min, na.rm = T)
  n_cancer = rowSums(data_case <= min, na.rm = T)

  #Computing genes with too many low values
  if(anyNA(controls)){
    #If any NA, the maximum number of samples depends of the number of NA for each genes
    t_ctrl = rowSums(!is.na(controls)) * threshold
    t_cancer = rowSums(!is.na(data_case)) * threshold
    results = (n_ctrl >= t_ctrl & n_cancer >= t_cancer)
  } else {
    #If no NA, the maximum number of samples is the same for all the genes
    results = (n_ctrl >= ncol(controls)*threshold & n_cancer >= ncol(data_case)*threshold)
  }

  print(paste0(sum(results), " genes have less than ", min, " counts in ",
               threshold*100, " % of the samples."))
  return(results)
}


# Authors: Clémentine Decamps, UGA
# clementine.decamps@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Find NA values
#'
#' This function detects probes or sample with more than threshold percent of
#' value undefined (NA).
#'
#'@param controls The first matrix with datas to analyze (ex: controls samples).
#'@param data_case The second matrix with datas to analyze (ex: tumors samples).
#'@param threshold The maximum proportion of NA tolerated for each probe (or sample).
#'@param probes TRUE if you want to analyze probes, FALSE to analyze samples.
#'
#'@return This function returns a True/False vector with True for the values to exclude.
#'
#'@example examples/ex_detect_na_value.R
#'
#'@export

detect_na_value = function(controls, data_case, threshold, probes = TRUE) {

  if (probes == TRUE){

    nactrl = rowSums(is.na(controls))
    nacancer = rowSums(is.na(data_case))
    low_na = nactrl >= threshold*ncol(controls) & nacancer >= threshold * ncol(data_case)
    print(paste0(sum(low_na), " probes are NA in at least ", threshold*100, " % of the samples."))
    return(low_na)

  } else {

    nactrl = colSums(is.na(controls))
    nacancer = colSums(is.na(data_case))
    low_na = c(nactrl >= threshold*nrow(controls), nacancer >= threshold*nrow(data_case))
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
#' This function makes the Penda dataset with controls and cases pre-filtred and
#' sorted by median in controls.
#'
#'@param controls The first matrix with datas to analyze (ex: controls samples).
#'@param data_case The second matrix with datas to analyze (ex: tumors samples).
#'@param detectlowvalue If detectlowvalue is True, genes with a low values in more than `threshold*100` \% of samples are removed.
#'@param detectNA If detectNA is True, genes and samples with more than `threshold*100` \% of NA values are removed.
#'@param threshold The maximum proportion of expression under val_min or NA tolerated for each gene or sample.
#'@param val_min The minimum value accepted. If val_min is NA, we compute this value as specified by the parameter `bimod`.
#'@param bimod If bimod is True and val_min NA, val_min is computed by mixtools::normalmixEM for a bimodal distribution, to search the value of the first peak. If bimod is False and val_min NA, val_min is computed for an unimodal distribution with a quantile of 0.1.
#'
#'@return This function return a list with preprocessed data_ctrl and data_case, and
#'the vector `info` with the different parameters.
#'
#'@example examples/ex_make_dataset.R
#'
#'@export

make_dataset = function(controls, data_case, detectlowvalue = TRUE, detectNA = TRUE,
                        threshold = 0.99, val_min = NA, bimod = TRUE) {

  resume = c(nrow(controls), ncol(controls), ncol(data_case), threshold,
             0, 0, val_min, 0)
  names(resume) = c("init_nb_genes", "init_nb_ctrls", "init_nb_cases", "threshold",
                    "nb_patients_NA", "nb_genes_NA", "val_min", "nb_genes_0")

  if(detectNA == TRUE){

    print("Removing rows and columns with too much NA")
    naprobes = penda::detect_na_value(controls, data_case, threshold, probes = TRUE)
    napatient = penda::detect_na_value(controls,  data_case, threshold, probes = FALSE)
    controls = controls[!naprobes, !napatient[1:ncol(controls)]]
    data_case = data_case[!naprobes, !napatient[(resume[2] + 1):length(napatient)]]
    resume["nb_patients_NA"] = sum(napatient)
    resume["nb_genes_NA"] = sum(naprobes)

  }

  if(detectlowvalue == TRUE){

    if(is.na(val_min)){

      if (bimod == TRUE){
        print("Computing of the low threshold for a bimodal distribution")
        all_data = c((log2(controls + 1)), (log2(data_case + 1)))
        mod = mixtools::normalmixEM(na.omit(all_data))
        mu = min(mod$mu)
        sig = mod$sigma[1]
        val_min = 2^qnorm(0.80, mu, sig)
        paste0("The low threshold is ", val_min)
        resume["val_min"] = val_min
      } else {
        print("Computing of the low threshold for an unimodal distribution")
        all_data = c(controls, data_case)
        val_min = quantile(all_data, 0.1)
        resume["val_min"] = val_min
      }

    }

    low_values = penda::detect_zero_value(controls, data_case,
                                          threshold = threshold, min = val_min)
    controls = controls[!low_values, ]
    data_case = data_case[!low_values, ]
    resume["nb_genes_0"] = sum(low_values)

  }

  median_gene = apply(controls, 1, median, na.rm = TRUE)
  median_gene = sort(median_gene)
  controls = controls[names(median_gene), ]
  data_case = data_case[names(median_gene), ]
  return(list(data_ctrl = controls, data_case = data_case, info = round(resume, 2)))

}
