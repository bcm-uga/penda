#include <Rcpp.h>
using namespace Rcpp;

//' compute_quant
//'
//' This function uses the R function "quantile" to compute quantiles in a rcpp function.
//'
//' @param data A numeric vector with the gene expressions for each patient.
//' @param qvalue A numeric vector with the probabilities for the quantile.
//'
//' @return This function returns a numeric vector with the quantiles.

NumericVector compute_quant(NumericVector data, NumericVector qvalue) {

  Environment stats("package:stats");
  Function quantile = stats["quantile"];
  int qsize = qvalue.size();
  NumericVector quantiles(qsize);
  for(int i=0; i<qsize; i++){
    quantiles[i] = as<double>(quantile(data, qvalue[i]));
  }
  return quantiles;
}


//' find_D_U_ctrl_rcpp
//'
//' This function ranks each gene and finds the genes which are more or less exprimed.
//'
//' @param ctrl_data A numeric matrix with the genes expressions for each patient.
//' @param quant A quantile to delimit the gene expression.
//' @param factor A factor to delimit the study limit : between quantile min / factor and quantile max * factor.
//' @param threshold The proportion of expression that must be in the conditions.
//'
//' @return This function returns a list of two logical matrices :
//' the D matrix, with TRUE if the row gene has a lower expression than the column gene,
//' and the U Matrix with TRUE if the row gene has a higher expression than the column gene.
//'
//' @export

// [[Rcpp::export]]
List find_D_U_ctrl_rcpp(NumericMatrix ctrl_data, double quant, double factor, double threshold){

  int size = ctrl_data.nrow();
  LogicalMatrix matrice_u(size, size);
  LogicalMatrix matrice_d(size, size);
  NumericVector quantiles = NumericVector::create(quant, 1-quant);

  //For each gene
  for (int i = 0 ; i < size ; i++){
    NumericVector gene = ctrl_data(i,_);
    NumericVector quantile_gene = compute_quant(gene, quantiles);
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


//' DU_rcpp
//'
//' This function makes the matrix D and U for all the genes.
//'
//' @param ctrl_data A numeric matrix with the genes expressions for each patient.
//' @param threshold The proportion of genes that must be under or above the gene.
//'
//' @return This function returns a list of two logical matrices :
//' the D matrix, with TRUE if the row gene has a lower expression than the column gene,
//' and the U Matrix with TRUE if the row gene has a higher expression than the column gene.
//' And a string vector with the genes names.
//'
//' @export

// [[Rcpp::export]]
List DU_rcpp(NumericMatrix ctrl_data, double threshold){

  int size = ctrl_data.nrow();
  LogicalMatrix matrice_u(size, size);
  LogicalMatrix matrice_d(size, size);

  List dimnames = ctrl_data.attr("dimnames");
  StringVector gene_names = dimnames[0];


  //For each gene
  for (int i = 0 ; i < size ; i++){
    NumericVector gene = ctrl_data(i,_);
    LogicalVector d_genes(size, false);
    LogicalVector u_genes(size, false);

    //For each gene (two-to-two comparison)
    for (int j = 0; j < size; j++){
      NumericVector gene_ctrl = ctrl_data(j,_);
      int nbu = 0;
      int nbd = 0;

      //For each patient, we compute the number of down and up-expressed genes
      for (int p = 0 ; p < ctrl_data.ncol() ; p++){
        if (gene_ctrl[p] > gene[p]){
          nbu++;
        }
        if (gene_ctrl[p] < gene[p]){
          nbd++;
        }
      }
      //We compare to the threshold
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
  total["n"] = gene_names;
  return(total);
}
