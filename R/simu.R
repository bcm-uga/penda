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


#'results_simulation
#'
#' This function computes the number of false positive, true positive, false negative
#' and false positive results after the test on a simulation.
#'
#'@param D_matrix The matrix of down-regulated genes
#'@param fraction The matrix of up-regulated genes
#'@param simulation The simulation data : the list with the iniatial_data,
#'simulated_data and changes_index.
#'
#'@return This function returns a list of four integer :
#'TP, the number of true postive results
#'FP, the number of false postive results
#'FN, the numbler of false negative results
#'TN, the numbler of true negative results
#'
#'@examples
#'D_U = find_D_U_ctrl(ctrl_data, quant = 0.001, factor = 4, threshold = 0.99)
#'simulation = dysreg_simulation(simu_data, fraction = 0.3, threshold = 60)
#'results_simulation(D_matrix = D_U$D, U_matrix = D_U$U, simulation  = simulation
#'
#'@export
#'

results_simulation = function(D_matrix, U_matrix, simulation){
  initial_data = simulation$initial_data
  simu_data = simulation$simulated_data
  down = initial_data > simu_data
  up = initial_data < simu_data

  #True positive = real down results and real up results
  TP = (length (which(down == 1 & D_matrix == 1))
        + length (which(up == 1 & U_matrix == 1)))

  #False positive = down or up results, but not changed in reality
  FP = (length (which(down == 0 & D_matrix == 1))
        + length (which(up == 0 & U_matrix == 1)))

  #False negative = not changed results, but down or up in reality
  FN = (length (which((up == 1 | down == 1) & D_matrix == 0 & U_matrix == 0)))

  #TN = real not changed results
  TN = (length (which(down == 0 & up == 0 & D_matrix == 0 & U_matrix == 0)))

  return(list(TP=TP, FP = FP, FN = FN, TN = TN ))
}

#'get_random_variable
#'
#' This function allows to get random values on a density function.
#'
#'@param F The density function.
#'@param x A vector with values of the function.
#'@param nb_r The number of random values wanted
#'
#'@return This function returns a vector with nb_r values of F.
#'
#'@examples
#'x = seq(-20, 20, 0.05)
#'F = function(x, u, a){
#'  (1/(a*sqrt(2*pi))) * exp(-(((x-u)^2)/(2*a^2)))
#'}
#'get_random_variable(F, x, nb_r = 10000, u = 0, a = 1)
#'
#'@export
#'

get_random_variable = function (F, x, nb_r,...){
  y = F(x,...)
  p1 <- cumsum(y)*diff(x)[1] #cumul
  qf <- approxfun(p1,x) #quantile function
  rf <- function(n) qf(runif(n)) #random variable
  return(rf(nb_r))
}

#'Heaviside
#'
#' The heaviside function value is zero for argument under the step, and one for argument above.
#'
#'@param x The number to evaluate
#'@param step The limit between heaviside = 0 and heaviside = 1
#'
#'@return This function returns 0 if the argument is under the step, 1 if not.
#'
#'@examples
#'Heaviside(1+98-103, step = 0)
#'
#'@export
#'

Heaviside = function(x, step = 0){
  return ((sign(x-step)+1)/2)
}

#P1 - distribution of up regulation factor for genes with an expression above 100
p1 = function(x){
  p1 = Heaviside(x-0.63)*(1+tanh((x-1.4)/0.4))*exp(-x/0.9)
  p1 = p1/(sum(p1)*diff(x)[1])
  return(p1)
}

#P2 - distribution of down regulation factor for genes with an expression above 100
p2 = function(x){
  p2 = Heaviside(x-0.63)*(1+tanh((x-1.4)/0.5))*exp(-x/1.3)
  p2 = p2/(sum(p2)*diff(x)[1])
  return(p2)
}

#P3 - distribtuon of expression difference for dysregulated genes with an expression under 100
p3 = function(x){
  p3 = Heaviside(x-20)*(exp(-x*0.007))*exp((x-20)*exp(-x*0.07))
  p3 = p3/(sum(p3)*diff(x)[1])
  return(p3)
}

#'complex_simulation
#'
#' This function simulates the dysregulation (up and down) of datas, using the distribution of
#' dysregulation in real datas.
#'
#'@param data vector or matrix of original data to edit
#'
#'@return This function returns a list of three vectors or matrix :
#'initial_data, datas before the modification
#'simulated_data, datas with up and down modifications
#'changes_index, index of datas with modifications
#'
#'@examples
#'complex_simulation(simu_data)
#'
#'@export
#'
complex_simulation = function(data){
  simu_data = data
  simu_data = sapply(simu_data, function(g){
    #if the expression is > 100
    if (g > 100){
      #22% chance of dysregulation
      if (sample(1:100, 1) <= 22) {
        #30% chance of up-regulation
        if (sample(1:100, 1) <= 30) {
          g = g * 2^(get_random_variable(p1, seq(0,100,0.1), 1))
        }
        #else, down-regulation
        else {
          g = g / 2^(get_random_variable(p2, seq(0,100,0.1), 1))
        }
      }
    #if expression < 100, 30% chance of dysregulation
    } else if (sample(1:100, 1) <= 30) {
      #50% chance of down-regulation, 50% chance of up-regulation. While expression < 0, retry.
      while(g < 0){
        if (sample (1:100,1) <= 50){
          g = g - (get_random_variable(p3, seq(0, 500, 0.5), 1))
        } else {
          g = g + (get_random_variable(p3, seq(0, 500, 0.5), 1))
        }
      }
    }
    return(g)
  })
  return(list(initial_data = data, simulated_data = simu_data, changes_idx = which(data != simu_data)))
}



