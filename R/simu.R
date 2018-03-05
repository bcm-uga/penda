#'dysreg_simulation
#'
#' This function simulates the dysregulation (up and down) of datas.
#'
#'@param data vector or matrix of original data to edit
#'@param fraction fraction of data to modify
#'@param threshold under this threshold, data modification is +-15. Above, it is *2 or /2.
#'
#'@return This function returns a list of three vectors or matrix :
#'initial_data, datas before the modification
#'simulated_data, datas with up and down modifications
#'changes_index, index of datas with modifications
#'
#'@examples
#'dysreg_simulation(simu_data, fraction = 0.3, threshold = 60)
#'
#'@export
#'

dysreg_simulation = function(data, fraction, threshold){
  simu_data = data
  perturb = sample(length(simu_data), fraction*length(simu_data))
  for(i in 1:length(perturb)){
    if (rnorm(1) < 0) {
      if (simu_data[perturb[i]] < threshold) {
        simu_data[perturb[i]] = simu_data[perturb[i]] - 30
        if (simu_data[perturb[i]] <= 0){
          simu_data[perturb[i]] = 0
        }
      } else {
        simu_data[perturb[i]] = simu_data[perturb[i]] / 4
      }
    } else if (simu_data[perturb[i]] < threshold) {
      simu_data[perturb[i]] = simu_data[perturb[i]] + 30
    } else {
      simu_data[perturb[i]] = simu_data[perturb[i]] * 4
    }
  }
  return(list(initial_data = data, simulated_data = simu_data, changes_index = perturb))
}
