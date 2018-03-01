#'dysreg_simulation
#'
#' This function simulates the dysregulation (up and down) of datas.
#'
#'@param simu_data vector of original data to be modified
#'@param fraction fraction of data to modify
#'@param threshold under this threshold, data modification is +-15. Above, it is *2 or /2.
#'
#'@return This function returns a vector of modified data which simulates up and down regulation.
#'
#'@examples
#'dysreg_simulation(simu_data, fraction = 0.3, threshold = 60)
#'
#'@export
#'

dysreg_simulation = function(simu_data, fraction, threshold){

  perturb = round(runif((fraction*length(simu_data)), min=1, max=length(simu_data)))
  perturb = perturb[-which(duplicated(perturb))]
  for(i in 1:length(perturb)){
    if (rnorm(1) < 0) {
      if (simu_data[perturb[i]] < threshold) {
        simu_data[perturb[i]] = simu_data[perturb[i]] -30
        if (simu_data[perturb[i]] <=0){
          simu_data[perturb[i]]=0
        }
      } else {
        simu_data[perturb[i]] = simu_data[perturb[i]] / 4
      }
    } else if (simu_data[perturb[i]] < threshold) {
      simu_data[perturb[i]] = simu_data[perturb[i]] +30
    } else {
      simu_data[perturb[i]] = simu_data[perturb[i]] * 4
    }
  }
  return(simu_data)
}
