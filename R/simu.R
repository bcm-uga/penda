#'simulation
#'
#' This function simulates up and down-regulation for datas
#'
#'@param data_simu vector of original data to be modified
#'@param fraction fraction of data to modify
#'@param threshold under this threshold, data modification is +-15. Above, it is *2 or /2.
#'
#'@return This function return a vector of modified data which simulate up and down regulation.
#'
#'
#'@export
#'

simulation = function(data_simu, fraction, threshold){

  perturb = round(runif((fraction*length(data_simu)), min=1, max=length(data_simu)))
  perturb = perturb[-which(duplicated(perturb))]
  for(i in 1:length(perturb)){
    if (rnorm(1) < 0) {
      if (data_simu[perturb[i]] < threshold) {
        data_simu[perturb[i]] = data_simu[perturb[i]] -15
        if (data_simu[perturb[i]] <=0){
          data_simu[perturb[i]]=0
        }
      } else {
        data_simu[perturb[i]] = data_simu[perturb[i]] / 2
      }
    } else if (data_simu[perturb[i]] < threshold) {
      data_simu[perturb[i]] = data_simu[perturb[i]] +15
    } else {
      data_simu[perturb[i]] = data_simu[perturb[i]] * 2
    }
  }
  return(data_simu)
}
