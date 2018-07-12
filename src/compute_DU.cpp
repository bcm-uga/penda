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

  double gene_p;
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
