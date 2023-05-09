#ifndef SCPLOTSIMPLEX_H
#define SCPLOTSIMPLEX_H

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

mat euclidean_dense(mat query, mat target);

NumericMatrix euclidean_dense2(NumericMatrix query, NumericMatrix target);

mat euclidean_sparse(sp_mat query, sp_mat target);

NumericVector euclideanColNorms(arma::mat x);

mat cosine_dense(mat query, mat target);

mat cosine_sparse(sp_mat query, sp_mat target);

#endif
