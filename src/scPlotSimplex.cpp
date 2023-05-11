#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;


// Distance Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//[[Rcpp::export]]
arma::mat euclidean_dense(arma::mat& query, arma::mat& target) {
    const unsigned int nq=query.n_cols, nt=target.n_cols, nr=query.n_rows;
    arma::mat output = zeros<mat>(nq, nt);
    arma::colvec xt(nr), xq(nr);
    double accum = 0.0;
    arma::colvec outCol = zeros<vec>(nq);
    for (unsigned int j=0; j<nt; j++) {
        xt = target.col(j);
        for (unsigned int i=0; i<nq; i++) {
            xq = query.col(i);
            xq -= xt;
            accum = accu(xq % xq);
            outCol(i) = sqrt(accum);
        }
        output.col(j) = outCol;
    }
    return output;
}


//[[Rcpp::export]]
arma::mat euclidean_sparse(arma::sp_mat query, arma::sp_mat target) {
    const unsigned int nq=query.n_cols, nt=target.n_cols, nr=query.n_rows;
    arma::mat output = zeros<mat>(nq, nt);
    arma::colvec xt(nr), xq(nr);
    double accum = 0.0;
    arma::colvec outCol = zeros<vec>(nq);
    for (unsigned int j=0; j<nt; j++) {
        xt = target.col(j);
        for (unsigned int i=0; i<nq; i++) {
            xq = query.col(i);
            xq -= xt;
            accum = accu(xq % xq);
            outCol(i) = sqrt(accum);
        }
        output.col(j) = outCol;
    }
    return output;
}

arma::vec colFrobNorms_dense(arma::mat x) {
    int nc = x.n_cols;
    arma::vec norms = zeros<vec>(nc);
    for (int j = 0; j < nc; j++) {
        norms(j) = arma::norm(x.col(j), "fro");
    }
    return norms;
}

//[[Rcpp::export]]
arma::mat cosine_dense(arma::mat query, arma::mat target) {
    const uword nq=query.n_cols, nt=target.n_cols;
    arma::mat output(nq,nt);
    output = query.t()*target;
    arma::vec qnorm = colFrobNorms_dense(query);
    arma::vec tnorm = colFrobNorms_dense(target);

    for (uword j = 0; j<nt; j++) {
        output.col(j) /= qnorm;
        output.col(j) /= tnorm(j);
    }
    double tol = 1e-12;
    for (arma::mat::iterator it = output.begin(); it != output.end(); it++) {
        if (*it < -(1-tol)) *it = -1;
        else if (*it > (1-tol)) *it = 1;
        *it = acos(*it) * 180 / M_PI;
    }
    return output;
}

arma::vec colFrobNorms_sparse(arma::sp_mat x) {
    const uword nc = x.n_cols;
    uword j = 0;
    arma::vec norms = zeros<vec>(nc);
    for (sp_mat::iterator it = x.begin(); it != x.end(); ++it) {
        j = it.col();
        norms(j) += (*it)*(*it);
    }
    norms = sqrt(norms);
    return norms;
}

//[[Rcpp::export]]
arma::mat cosine_sparse(arma::sp_mat query, arma::sp_mat target) {
    const uword nq=query.n_cols, nt=target.n_cols;
    arma::mat output(nq,nt);
    output = query.t()*target;
    arma::vec qnorm = colFrobNorms_sparse(query);
    arma::vec tnorm = colFrobNorms_sparse(target);
    for (uword j = 0; j<nt; j++) {
        output.col(j) /= qnorm;
        output.col(j) /= tnorm(j);
    }
    double tol = 1e-12;
    for (arma::mat::iterator it = output.begin(); it != output.end(); it++) {
        if (*it < -(1-tol)) *it = -1;
        else if (*it > (1-tol)) *it = 1;
        *it = acos(*it) * 180 / M_PI;
    }
    return output;
}


// Other Utilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//[[Rcpp::export]]
NumericVector rowMeans_sparse(arma::sp_mat x) {
    NumericVector rm = (x.n_rows);
    for(sp_mat::iterator it = x.begin(); it != x.end(); ++it) {
        rm(it.row()) += *it;
    }
    rm = rm / x.n_cols;
    return rm;
}
