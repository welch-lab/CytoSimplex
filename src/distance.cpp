#include "scPlotSimplex.h"

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
arma::sp_mat euclidean_sparse(arma::sp_mat query, arma::sp_mat target) {
    const int nq=query.n_cols, nt=target.n_cols;
    arma::sp_mat output(nq,nt);
    return output;
}

//[[Rcpp::export]]
NumericVector euclideanColNorms(arma::mat x) {
    int nc = x.n_cols;
    NumericVector norms(nc, 0.0);
    for (int j = 0; j < nc; j++) {
        arma::colvec v = x.col(j);
        norms(j) = sum(square(v));
    }
    norms = sqrt(norms);
    return norms;
}

//[[Rcpp::export]]
arma::mat cosine_dense(arma::mat query, arma::mat target) {
    const int nq=query.n_cols, nt=target.n_cols;
    arma::mat output(nq,nt);
    output = query.t()*target;
    NumericVector qnorm = euclideanColNorms(query);
    NumericVector tnorm = euclideanColNorms(target);
    arma::colvec v(nq);
    for (unsigned j = 0; j<nt; j++) {
        v = output.col(j);
        v /= qnorm;
        v /= tnorm(j);
        output.col(j) = v;
    }
    double tol = 1e-12;
    unsigned int n_clamped = 0;
    for (arma::mat::iterator it = output.begin(); it != output.end(); it++) {
        if (*it < -(1-tol)) {
            if (*it < -(1+tol)) n_clamped++;
            *it = -1;
        }
        else if (*it > (1-tol)) {
            if (*it > (1+tol)) n_clamped++;
            *it = 1;
        }
        *it = acos(*it) * 180 / M_PI;
    }
    return output;
}

//[[Rcpp::export]]
arma::sp_mat cosine_sparse(arma::sp_mat query, arma::sp_mat target) {
    const int nq=query.n_cols, nt=target.n_cols;
    arma::sp_mat output(nq,nt);
    return output;
}
