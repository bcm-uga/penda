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
    for (j = 0; j < size ; j++){

      nbu = nbd = n;
      possible_u = matrice_u(j, i) = possible_d = matrice_d(j, i) = true;

      //For each patient, we compute the number of down and up-expressed genes
      for (p = 0 ; p < n ; p++){

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

  int i, j, p, nbu, nbd, n = data_ctrl.ncol();
  bool possible_u, possible_d;

  //For each gene
  for (i = 0 ; i < size ; i++){

    //For each gene (two-to-two comparison)
    for (j = 0; j < size ; j++){

      nbu = nbd = n;
      possible_u = matrice_u(j, i) = possible_d = matrice_d(j, i) = true;

      //For each patient, we compute the number of down and up-expressed genes
      for (p = 0 ; p < n; p++){

        if (possible_u && (data_ctrl(j, p) <= data_ctrl(i, p))){
          nbu--;
          if (nbu <= ((data_ctrl.ncol()-nbNA[j]) * threshold)){
            matrice_u(j, i) = false;
            possible_u = false;
          }
        }
        if (possible_d && (data_ctrl(j, p) >= data_ctrl(i, p))){
          nbd--;
          if (nbd <= ((data_ctrl.ncol()-nbNA[j]) * threshold)){
            matrice_d(j, i) = false;
            possible_d = false;

          }
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

//' compute_DU_cpp_V2
//'
//' This function makes the matrix D and U for all the genes with a new method.
//'
//' @param data_ctrl A numeric matrix with the genes expressions for each patient. Must be sort by the median.
//' @param threshold The proportion of genes that must be under or above the gene.
//' @param size_max The maximum number of down and up-expressed gene for each genes.
//'
//' @return This function returns a list of two logical matrices :
//' the D matrix, with the id of closest genes with a lower expression,
//' and the U Matrix with the id of closest genes with a higher expression.
//' And a string vector with the genes names.
//'
//' @export
// [[Rcpp::export]]
List compute_DU_cpp_V2(NumericMatrix data_ctrl, double threshold, double size_max){

  int size = data_ctrl.nrow();

  NumericMatrix matrice_u(size, size_max);
  NumericMatrix matrice_d(size, size_max);

  List dimnames = data_ctrl.attr("dimnames");
  StringVector gene_names = dimnames[0];

  int i, j, p, nbu, nbd, n = data_ctrl.ncol();
  bool possible_u, possible_d;

  //For each gene
  for (i = 0 ; i < size ; i++){

    int nb_u = 0;
    int nb_d = 0;
    NumericVector u_genes(size_max);
    NumericVector d_genes(size_max);

    //Compute U genes
    if(i < (size-1)){
      for (j = (i+1) ; j < size ; j++){
        nbu = n;
        possible_u = true;
        //For each patient, we compute the number of down and up-expressed genes
        for (p = 0 ; p < n ; p++){
          if (possible_u && (data_ctrl(j, p) <= data_ctrl(i, p))){
            nbu--;
            if (nbu <= (n * threshold)){
              possible_u = false;
            }
          }
          if (!possible_u) break;
          if(p == (n-1)){
            u_genes[nb_u] = (j+1);
            nb_u ++;
          }
        }
        if(nb_u == 30) break;
      }
    }
    matrice_u(i,_) = u_genes;

    //Compute D genes
    if(i > 1){
      for (j = (i-1) ; j > 0 ; j--){
        nbd = n;
        possible_d = true;
        //For each patient, we compute the number of down and up-expressed genes
        for (p = 0 ; p < n ; p++){
          if (possible_d && (data_ctrl(j, p) >= data_ctrl(i, p))){
            nbd--;
            if (nbd <= (n * threshold)){
              possible_d = false;
            }
          }
          if (!possible_d) break;
          if(p == (n-1)){
            d_genes[nb_d] = (j+1);
            nb_d ++;
          }
        }
        if(nb_d == 30) break;
      }
    }
    matrice_d(i,_) = d_genes;
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
List compute_DU_cppNA_V2(NumericMatrix data_ctrl, double threshold, NumericVector nbNA, double size_max){

  int size = data_ctrl.nrow();
  NumericMatrix matrice_u(size, size_max);
  NumericMatrix matrice_d(size, size_max);

  List dimnames = data_ctrl.attr("dimnames");
  StringVector gene_names = dimnames[0];

  int i, j, p, nbu, nbd, n = data_ctrl.ncol();
  bool possible_u, possible_d;

  //For each gene
  for (i = 0 ; i < size ; i++){

    int nb_u = 0;
    int nb_d = 0;
    NumericVector u_genes(size_max);
    NumericVector d_genes(size_max);

    //Compute U
    if(i < (size-1)){
      for (j = (i+1) ; j < size ; j++){
        nbu = n;
        possible_u = true;

        //For each patient, we compute 30 up-expressed genes
        for (p = 0 ; p < n; p++){
          if (possible_u && (data_ctrl(j, p) <= data_ctrl(i, p))){
            nbu--;
            if (nbu <= ((data_ctrl.ncol()-nbNA[j]) * threshold)){
              possible_u = false;
            }
          }
          if (!possible_u) break;
          if(p == (n-1)){
            u_genes[nb_u] = (j+1);
            nb_u ++;
          }
        }
        if(nb_u == 30) break;
      }
    }
    matrice_u(i,_) = u_genes;

    //Compute D
    if(i > 1){
      for (j = (i-1); j > 0 ; j--){
      nbd = n;
      possible_d = true;

      //For each patient, we compute 30 down-expressed genes
      for (p = 0 ; p < n; p++){

        if (possible_d && (data_ctrl(j, p) >= data_ctrl(i, p))){
          nbd--;
          if (nbd <= ((data_ctrl.ncol()-nbNA[j]) * threshold)){
            possible_d = false;
          }
        }
        if (!possible_d) break;
        if(p == (n-1)){
          d_genes[nb_d] = (j+1);
          nb_d ++;
        }
      }
      if(nb_d == 30) break;
      }
    }
    matrice_d(i,_) = d_genes;
  }

  List total;
  total["U"] = matrice_u;
  total["D"] = matrice_d;
  total["n"] = gene_names;
  return(total);
}
