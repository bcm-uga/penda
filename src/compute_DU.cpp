#include <Rcpp.h>
using namespace Rcpp;


//' compute_LH_cpp
//'
//' This function makes the matrices L and H for all the genes with a new method.
//'
//' @param data_ctrl A numeric matrix with the genes expressions for each patient. Must be sort by the median.
//' @param threshold The proportion of genes that must be under or above the gene.
//' @param size_max The maximum number of L and H genes for each gene..
//'
//' @return This function returns a list of two logical matrices :
//' the L matrix, with the id of closest genes with a lower expression,
//' the H Matrix with the id of closest genes with a higher expression,
//' And the n string vector with the genes names.
//'
//' @export
// [[Rcpp::export]]
List compute_LH_cpp(NumericMatrix data_ctrl, double threshold, double size_max){

  int size = data_ctrl.nrow();

  NumericMatrix matrice_l(size, size_max);
  NumericMatrix matrice_h(size, size_max);

  List dimnames = data_ctrl.attr("dimnames");
  StringVector gene_names = dimnames[0];

  int i, j, p, nbl, nbh, n = data_ctrl.ncol();
  bool possible_l, possible_h;

  //For each gene
  for (i = 0 ; i < size ; i++){

    int nb_l = 0;
    int nb_h = 0;
    NumericVector l_genes(size_max);
    NumericVector h_genes(size_max);

    //Compute lower genes
    if(i > 1){
      for (j = (i-1) ; j > 0 ; j--){
        nbl = n;
        possible_l = true;
        //For each patient, we compute the number of genes with a lower expression
        for (p = 0 ; p < n ; p++){
          if (possible_l && (data_ctrl(j, p) >= data_ctrl(i, p))){
            nbl--;
            if (nbl <= (n * threshold)){
              possible_l = false;
            }
          }
          if (!possible_l) break;
          if(p == (n-1)){
            l_genes[nb_l] = (j+1);
            nb_l ++;
          }
        }
        if(nb_l == size_max) break;
      }
    }
    matrice_l(i,_) = l_genes;

    //Compute H genes
    if(i < (size-1)){
      for (j = (i+1) ; j < size ; j++){
        nbh = n;
        possible_h = true;
        //For each patient, we compute the number of genes with an higher expression
        for (p = 0 ; p < n ; p++){
          if (possible_h && (data_ctrl(j, p) <= data_ctrl(i, p))){
            nbh--;
            if (nbh <= (n * threshold)){
              possible_h = false;
            }
          }
          if (!possible_h) break;
          if(p == (n-1)){
            h_genes[nb_h] = (j+1);
            nb_h ++;
          }
        }
        if(nb_h == size_max) break;
      }
    }
    matrice_h(i,_) = h_genes;

  }

  List total;
  total["L"] = matrice_l;
  total["H"] = matrice_h;
  total["n"] = gene_names;
  return(total);
}

//' compute_LH_cppNA
//'
//' This function makes the matrices L and H for all the genes without count NA values.
//'
//' @param data_ctrl A numeric matrix with the genes expressions for each patient.
//' @param threshold The proportion of genes that must be under or above the gene.
//' @param nbNA The number of NA for each gene.
//' @param size_max The maximum number of L and H genes for each gene..
//'
//' @return This function returns a list of two logical matrices :
//' the L matrix, with the id of closest genes with a lower expression,
//' the H Matrix with the id of closest genes with a higher expression,
//' And the n string vector with the genes names.
//'
//' @export
// [[Rcpp::export]]
List compute_LH_cppNA(NumericMatrix data_ctrl, double threshold, NumericVector nbNA, double size_max){

  int size = data_ctrl.nrow();
  NumericMatrix matrice_l(size, size_max);
  NumericMatrix matrice_h(size, size_max);

  List dimnames = data_ctrl.attr("dimnames");
  StringVector gene_names = dimnames[0];

  int i, j, p, nbl, nbh, n = data_ctrl.ncol();
  bool possible_l, possible_h;

  //For each gene
  for (i = 0 ; i < size ; i++){

    int nb_l = 0;
    int nb_h = 0;
    NumericVector l_genes(size_max);
    NumericVector h_genes(size_max);

    //Compute L
    if(i > 1){
      for (j = (i-1); j > 0 ; j--){
        nbl = n;
        possible_l = true;

        //For each patient, we compute genes with a lower expression
        for (p = 0 ; p < n; p++){

          if (possible_l && (data_ctrl(j, p) >= data_ctrl(i, p))){
            nbl--;
            if (nbl <= ((data_ctrl.ncol()-nbNA[j]) * threshold)){
              possible_l = false;
            }
          }
          if (!possible_l) break;
          if(p == (n-1)){
            l_genes[nb_l] = (j+1);
            nb_l ++;
          }
        }
        if(nb_l == size_max) break;
      }
    }
    matrice_l(i,_) = l_genes;

    //Compute H
    if(i < (size-1)){
      for (j = (i+1) ; j < size ; j++){
        nbh = n;
        possible_h = true;

        //For each patient, we compute genes with an higher expression
        for (p = 0 ; p < n; p++){
          if (possible_h && (data_ctrl(j, p) <= data_ctrl(i, p))){
            nbh--;
            if (nbh <= ((data_ctrl.ncol()-nbNA[j]) * threshold)){
              possible_h = false;
            }
          }
          if (!possible_h) break;
          if(p == (n-1)){
            h_genes[nb_h] = (j+1);
            nb_h ++;
          }
        }
        if(nb_h == size_max) break;
      }
    }
    matrice_h(i,_) = h_genes;

  }

  List total;
  total["L"] = matrice_l;
  total["H"] = matrice_h;
  total["n"] = gene_names;
  return(total);
}
