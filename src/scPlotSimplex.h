#ifndef SCPLOTSIMPLEX_H
#define SCPLOTSIMPLEX_H

#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_BLAS
#define ARMA_USE_LAPACK

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

mat euclidean_dense(mat query, mat target);

NumericMatrix euclidean_dense2(NumericMatrix query, NumericMatrix target);

sp_mat euclidean_sparse(sp_mat query, sp_mat target);

NumericVector euclideanColNorms(arma::mat x);

mat cosine_dense(mat query, mat target);

sp_mat cosine_sparse(sp_mat query, sp_mat target);

#endif
