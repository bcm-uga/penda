#include <Rcpp.h>
using namespace Rcpp;

//' compute_DU_cpp
//'
//' This function makes the matrix D and U for all the genes.
//'
//' @param data_ctrl A numeric matrix with the genes expressions for each patient.
//' @param threshold The proportion of genes that must be under or above the gene.
//'
//' @return This function returns a list of two logical matrices :
//' the D matrix, with TRUE if the row gene has a lower expression than the column gene,
//' and the U Matrix with TRUE if the row gene has a higher expression than the column gene.
//' And a string vector with the genes names.
//'
//' @export
// [[Rcpp::export]]
List compute_DU_cpp(NumericMatrix data_ctrl, double threshold){

  int size = data_ctrl.nrow();
  LogicalMatrix matrice_u(size, size);
  LogicalMatrix matrice_d(size, size);

  List dimnames = data_ctrl.attr("dimnames");
  StringVector gene_names = dimnames[0];

  int i, j, p, nbu, nbd, n = data_ctrl.ncol();
  bool possible_u, possible_d;
  double thr = n * threshold;

  //For each gene
  for (i = 0 ; i < size ; i++){

    //For each gene (two-to-two comparison)
    for (j = 0; j < size; j++){

      nbu = nbd = n;
      possible_u = matrice_u(j, i) = possible_d = matrice_d(j, i) = true;

      //For each patient, we compute the number of down and up-expressed genes
      for (p = 0 ; p < n; p++){

        if (possible_u && (data_ctrl(j, p) <= data_ctrl(i, p))){
          nbu--;
          if (nbu <= thr) matrice_u(j, i) = false;
          possible_u = false;
        }
        if (possible_d && (data_ctrl(j, p) >= data_ctrl(i, p))){
          nbd--;
          if (nbd <= thr) matrice_d(j, i) = false;
          possible_d = false;
        }

        if (!possible_u && !possible_d) break;
      }
    }
  }

  List total;
  total["U"] = matrice_u;
  total["D"] = matrice_d;
  total["n"] = gene_names;
  return(total);
}

//' compute_DU_cppNA
//'
//' This function makes the matrix D and U for all the genes without count NA values.
//'
//' @param data_ctrl A numeric matrix with the genes expressions for each patient.
//' @param threshold The proportion of genes that must be under or above the gene.
//' @param nbNA The number of NA for each gene.
//'
//' @return This function returns a list of two logical matrices :
//' the D matrix, with TRUE if the row gene has a lower expression than the column gene,
//' and the U Matrix with TRUE if the row gene has a higher expression than the column gene.
//' And a string vector with the genes names.
//'
//' @export
// [[Rcpp::export]]
List compute_DU_cppNA(NumericMatrix data_ctrl, double threshold, NumericVector nbNA){

  int size = data_ctrl.nrow();
  LogicalMatrix matrice_u(size, size);
  LogicalMatrix matrice_d(size, size);

  List dimnames = data_ctrl.attr("dimnames");
  StringVector gene_names = dimnames[0];


  //For each gene
  for (int i = 0 ; i < size ; i++){
    NumericVector gene = data_ctrl(i,_);
    LogicalVector d_genes(size, false);
    LogicalVector u_genes(size, false);

    //For each gene (two-to-two comparison)
    for (int j = 0; j < size; j++){
      NumericVector gene_ctrl = data_ctrl(j,_);
      int nbu = 0;
      int nbd = 0;

      //For each patient, we compute the number of down and up-expressed genes
      for (int p = 0 ; p < data_ctrl.ncol() ; p++){
        if (gene_ctrl[p] > gene[p]){
          nbu++;
        }
        if (gene_ctrl[p] < gene[p]){
          nbd++;
        }
      }
      //We compare to the threshold
      if(nbu > ((data_ctrl.ncol()-nbNA[j]) * threshold)){
        u_genes[j] = true;
      }
      if (nbd > ((data_ctrl.ncol()-nbNA[j]) * threshold)){
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

//' compute_DU_cpp_size
//'
//' This function makes the matrix D and U for all the genes.
//'
//' @param data_ctrl A numeric matrix with the genes expressions for each patient.
//' @param threshold The proportion of genes that must be under or above the gene.
//'
//' @return This function returns a list of two logical matrices :
//' the D matrix, with TRUE if the row gene has a lower expression than the column gene,
//' and the U Matrix with TRUE if the row gene has a higher expression than the column gene.
//' And a string vector with the genes names.
//'
//' @export
// [[Rcpp::export]]
List compute_DU_cpp_size(NumericVector data_gene, NumericMatrix data_ctrl, double threshold){

  int size = data_ctrl.nrow();
  List dimnames = data_ctrl.attr("dimnames");
  StringVector gene_names = dimnames[0];

  LogicalVector d_genes(size, false);
  LogicalVector u_genes(size, false);
  //For each gene (two-to-two comparison)
  for (int j = 0; j < size; j++){
    NumericVector gene_ctrl = data_ctrl(j,_);
    int nbu = 0;
    int nbd = 0;
    //For each patient, we compute the number of down and up-expressed genes
    for (int p = 0 ; p < data_ctrl.ncol() ; p++){
      if (gene_ctrl[p] > data_gene[p]){
        nbu++;
      }
      if (gene_ctrl[p] < data_gene[p]){
        nbd++;
      }
    }
    //We compare to the threshold
    if(nbu > (data_ctrl.ncol() * threshold)){
      u_genes[j] = true;
    }
    if (nbd > (data_ctrl.ncol() * threshold)){
      d_genes[j] = true;
    }
  }

  List total;
  total["U"] = u_genes;
  total["D"] = d_genes;
  total["n"] = gene_names;
  return(total);
}
