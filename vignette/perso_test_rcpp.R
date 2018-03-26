

Rcpp::sourceCpp('R/test_rcpp.cpp')

print("DU 1")
tddu1 = Sys.time()
D_U_ctrl = find_D_U_ctrl(ctrl_data, 0.001, 4, 0.99)
tfdu1 = Sys.time()
print("DU 2")
tddu2 = Sys.time()
D_U_ctrl = find_D_U_ctrl(ctrl_data, 0.001, 4, 0.99)
tfdu2 = Sys.time()
print("DU 3")
tddu3 = Sys.time()
D_U_ctrl = find_D_U_ctrl(ctrl_data, 0.001, 4, 0.99)
tfdu3 = Sys.time()

print("DUrcpp 1")
tddurcpp1 = Sys.time()
D_U_ctrlrcpp = find_D_U_ctrl_rcpp(ctrl_data, 0.001, 4, 0.99)
tfdurcpp1 = Sys.time()
print("DUrcpp 2")
tddurcpp2 = Sys.time()
D_U_ctrlrcpp = find_D_U_ctrl_rcpp(ctrl_data, 0.001, 4, 0.99)
tfdurcpp2 = Sys.time()
print("DUrcpp 3")
tddurcpp3 = Sys.time()
D_U_ctrlrcpp = find_D_U_ctrl_rcpp(ctrl_data, 0.001, 4, 0.99)
tfdurcpp3 = Sys.time()

dimnames(D_U_ctrlrcpp$D) = (list(rownames(ctrl_data), rownames(ctrl_data)))
dimnames(D_U_ctrlrcpp$U) = (list(rownames(ctrl_data), rownames(ctrl_data)))

table(D_U_ctrl$D !=D_U_ctrlrcpp$D)
table(D_U_ctrl$U !=D_U_ctrlrcpp$U)


tdtest = Sys.time()
D_matrix = simuv1$simulated_data
U_matrix = simuv1$simulated_data
for(p in 1:5){
  print(c("Patient number",p))
  D_U = patient_test(ctrl_data, LUSC_data[,3], 0.03, 20, D_U_ctrlrcpp, 0.2)
  D_matrix[,p] = unlist(D_U[1])
  U_matrix[,p] = unlist(D_U[2])
}
tftest = Sys.time()

microbenchmark::microbenchmark(  D_U = patient_test(ctrl_data, LUSC_data[,3], 0.03, 20, D_U_ctrlrcpp, 0.2)
, times = 10)
