#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//rcpp function to use the R quantile
NumericVector computeQuant(NumericVector data, NumericVector qvalue) {
  Environment stats("package:stats");
  Function quantile = stats["quantile"];
  int qsize = qvalue.size();
  NumericVector quantiles(qsize);
  for(int i=0; i<qsize; i++){
    quantiles[i] = as<double>(quantile(data, qvalue[i]));
  }
  return quantiles;
}

// [[Rcpp::export]]
List find_D_U_ctrl_rcpp(NumericMatrix ctrl_data, double quant, double factor, double threshold){

  int size = ctrl_data.nrow();
  LogicalMatrix matrice_u(size, size);
  LogicalMatrix matrice_d(size, size);
  NumericVector quantiles = NumericVector::create(quant, 1-quant);
  //For each gene
  for (int i = 0 ; i < size ; i++){
    NumericVector gene = ctrl_data(i,_);
    NumericVector quantile_gene = computeQuant(gene, quantiles);
    LogicalVector u_genes(size, false);
    LogicalVector d_genes(size, false);
    //For each gene (two-to-two comparison)
    for (int j = 0; j < size; j++){
      NumericVector gene_ctrl = ctrl_data(j,_);
      int nbu = 0;
      int nbd = 0;
      //For each patient
      for (int p = 0 ; p < ctrl_data.ncol() ; p++){
        if ((gene_ctrl[p] > gene[p]) & (gene_ctrl[p] <= (round(quantile_gene[1]*factor)))){
          nbu++;
        }
        if ((gene_ctrl[p] < gene[p]) & (gene_ctrl[p] >= (round(quantile_gene[0]/factor)))){
          nbd++;
        }
      }
      if(nbu > (ctrl_data.ncol() * threshold)){
        u_genes[j] = true;
      }
      if (nbd > (ctrl_data.ncol() * threshold)){
        d_genes[j] = true;
      }
    }
    matrice_u(_,i) = u_genes;
    matrice_d(_,i) = d_genes;
  }
  List total;
  total["U"] = matrice_u;
  total["D"] = matrice_d;
  return(total);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
# ctrl_data = readRDS('~/projects/perso_DA/data_ctrl_process.rds')[1:1000,]
# D_U_ctrl = find_D_U_ctrl(ctrl_data, 0.001, 4, 0.99)
# D_U_ctrlrcpp = find_D_U_ctrl_rcpp(ctrl_data, 0.001, 4, 0.99)
# dimnames(D_U_ctrlrcpp$D) = (list(rownames(ctrl_data), rownames(ctrl_data)))
# dimnames(D_U_ctrlrcpp$U) = (list(rownames(ctrl_data), rownames(ctrl_data)))
#
# table(D_U_ctrl$D != D_U_ctrlrcpp$D)
# table(D_U_ctrl$U != D_U_ctrlrcpp$U)

*/
