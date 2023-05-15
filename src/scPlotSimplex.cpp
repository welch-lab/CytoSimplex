#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Rcpp;
using namespace arma;

// Utilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//[[Rcpp::export]]
arma::mat colNormalize_dense(arma::mat x, arma::vec colsums) {
    arma::mat output=zeros<mat>(x.n_rows, x.n_cols);
    for (unsigned int c=0; c<x.n_cols; c++) {
        output.col(c)=x.col(c)/colsums(c);
    }
    return output;
}

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

// WILCOXON TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// [[Rcpp::export]]
Rcpp::List rankMatrix_dense(arma::mat& X) {
    // sizes of tied groups
    arma::inplace_trans(X);
    vector<list<float> > ties(X.n_cols);

    std::vector<pair<float, size_t> > v_sort(X.n_rows);
    for (unsigned c = 0; c < X.n_cols; c++) {
        for (size_t i = 0; i < X.n_rows; i++) {
            v_sort[i] = make_pair(X.col(c)[i], i);
        }
        sort(v_sort.begin(), v_sort.end());

        float rank_sum = 0, n = 1;
        size_t i;
        for (i = 1U; i < v_sort.size(); i++) {
            if (v_sort[i].first != v_sort[i - 1].first) {
                // if current val != prev val
                // set prev val to something
                for (unsigned j = 0; j < n; j++) {
                    X.col(c)[v_sort[i - 1 - j].second] = (rank_sum / n) + 1;
                }
                // restart count ranks
                rank_sum = i;
                if (n > 1) ties[c].push_back(n);
                n = 1;
            } else {
                // if curr val is a tie,
                // don't set anything yet, start computing mean rank
                rank_sum += i;
                n++;
            }
        }
        // set the last element(s)
        for (unsigned j = 0; j < n; j++)
            X.col(c)[v_sort[i - 1 - j].second] = (rank_sum / n) + 1;
    }
    return Rcpp::List::create(Named("X_ranked") = X, Named("ties") = ties);
}

std::list<float> cpp_in_place_rank_mean(arma::vec& v_temp, int idx_begin,
                                        int idx_end) {
    std::list<float> ties;

    if (idx_begin > idx_end) return ties;
    std::vector<pair<float, size_t> > v_sort(idx_end - idx_begin + 1);
    for (size_t i = idx_begin; i <= idx_end; i++) {
        v_sort[i - idx_begin] = make_pair(v_temp[i], i - idx_begin);
    }


    sort(v_sort.begin(), v_sort.end());

    float rank_sum = 0, n = 1;
    size_t i;
    for (i = 1U; i < v_sort.size(); i++) {
        if (v_sort[i].first != v_sort[i - 1].first) {
            // if current val != prev val
            // set prev val to something
            for (unsigned j = 0; j < n; j++) {
                v_temp[v_sort[i - 1 - j].second + idx_begin] =
                    (rank_sum / n) + 1;
            }
            // restart count ranks
            rank_sum = i;
            if (n > 1) ties.push_back(n);
            n = 1;
        } else {
            // if curr val is a tie,
            // don't set anything yet, start computing mean rank
            rank_sum += i;
            n++;
        }
    }
    // set the last element(s)
    for (unsigned j = 0; j < n; j++)
        v_temp[v_sort[i - 1 - j].second + idx_begin] = (rank_sum / n) + 1;

    return ties;
}

// [[Rcpp::export]]
std::vector<std::list<float> > cpp_rank_matrix_dgc(
        arma::vec& x, const arma::vec& p, int nrow, int ncol) {
    vector<list<float> > ties(ncol);
    int n_zero;
    for (int i = 0; i < ncol; i++) {
        if (p[i+1] == p[i]) continue;
        n_zero = nrow - (p[i+1] - p[i]);
        ties[i] = cpp_in_place_rank_mean(x, p[i], p[i + 1] - 1);
        ties[i].push_back(n_zero);
        x.rows(p[i], p[i + 1] - 1) += n_zero;
    }
    return ties;
}



// SUM AGGREGATE %%%%


// [[Rcpp::export]]
arma::mat rowAggregateSum_dense(const arma::mat& X,
                                const arma::uvec& groups,
                                unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_cols);
    for (unsigned r = 0; r < X.n_rows; r++) {
        res.row(groups[r]) += sum(X.row(r), 0);
    }
    return res;
}

// [[Rcpp::export]]
arma::mat rowAggregateSum_sparse(arma::sp_mat& X,
                                 const arma::uvec& groups,
                                 unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_cols);
    for (arma::sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
        res(groups[it.row()], it.col()) += *it;
    }
    return res;
}

// [[Rcpp::export]]
arma::mat colAggregateSum_dense(const arma::mat& X,
                                const arma::uvec& groups,
                                unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_rows);
    for (unsigned c = 0; c < X.n_cols; c++) {
        res.row(groups[c]) += sum(X.col(c), 1).t();
    }
    return res;
}

// [[Rcpp::export]]
arma::mat colAggregateSum_sparse(arma::sp_mat& X,
                                const arma::uvec& groups,
                                unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_rows);
    for (arma::sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
        res(groups[it.col()], it.row()) += *it;
    }
    return res;
}

// Non-zero counting aggregate %%%%
// NNZ - Number of Non-Zero

// [[Rcpp::export]]
arma::mat colNNZAggr_dense(const arma::mat& X,
                           const arma::uvec& groups,
                           unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_rows);
    for (unsigned c = 0; c < X.n_cols; c++) {
        for (unsigned r = 0; r < X.n_rows; r++) {
            if (X(r, c) != 0) res(groups[c], r)++;
        }
    }
    return res;
}

// [[Rcpp::export]]
arma::mat colNNZAggr_sparse(arma::sp_mat& X,
                           const arma::uvec& groups,
                           unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_rows);
    for (arma::sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
        if (*it > 0) res(groups[it.col()], it.row())++;
    }
    return res;
}

// [[Rcpp::export]]
arma::mat rowNNZAggr_dense(const arma::mat& X,
                           const arma::uvec& groups,
                           unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_cols);
    for (unsigned c = 0; c < X.n_cols; c++) {
        for (unsigned r = 0; r < X.n_rows; r++) {
            if (X(r, c) != 0) res(groups[r], c)++;
        }
    }
    return res;
}

// [[Rcpp::export]]
arma::mat rowNNZAggr_sparse(arma::sp_mat& X,
                            const arma::uvec& groups,
                            unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_cols);
    for (arma::sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
        if (*it > 0) res(groups[it.row()], it.col())++;
    }
    return res;
}

