#include <RcppArmadillo.h>

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

// [[Rcpp::export]]
bool is_rawCounts_sparse(const arma::sp_mat& x) {
    for(arma::sp_mat::const_iterator it = x.begin(); it != x.end(); ++it)
    {
        double value = *it;
        if (value != 0 && value != static_cast<int>(value)) {
            return false;
        }
    }
    return true;
}

// [[Rcpp::export]]
bool is_rawCounts_dense(const arma::mat& x) {
  for (arma::mat::const_iterator it = x.begin(); it != x.end(); ++it) {
    double value = *it;
    if (value != 0 && value != static_cast<int>(value)) {
      return false;
    }
  }
  return true;
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

// X - feature x cell
// output - cell x featureRank

// [[Rcpp::export]]
Rcpp::List rankMatrix_dense(arma::mat& X) {
    arma::inplace_trans(X);
    vector<list<double> > ties(X.n_cols);
    std::vector<pair<float, size_t> > v_sort(X.n_rows);

    for (unsigned c = 0; c < X.n_cols; c++) {
        // Get the elements sorted, while using the pair to record their
        // original position
        for (size_t i = 0; i < X.n_rows; i++) {
            v_sort[i] = make_pair(X.col(c)[i], i);
        }
        sort(v_sort.begin(), v_sort.end());

        // Use `tieStartRank` to record where a tie starts
        double tieStartRank = 0;
        // Starting from the second element in the sort result, if the value
        // changes, it marks the end of a tie from `tieStartRank` to `i - 1`
        // Instead of literally summing up all ranking from `tieStartRank` to
        // `i - 1` and divide it by n, `(tieStartRank + i - 1) / 2 ` results in
        // a precise equivalent.
        for (int i = 1; i < v_sort.size(); i++) {
            if (v_sort[i].first != v_sort[i - 1].first) {
                double tieRank = (tieStartRank + i - 1) / 2 + 1;
                for (unsigned j = tieStartRank; j < i; j++) {
                    X.col(c)[v_sort[j].second] = tieRank;
                }
                // Record the tie length
                if ((i - tieStartRank) > 1) ties[c].push_back(i - tieStartRank);
                // Reset the start of a tie
                tieStartRank = i;
            }
        }
        // set the last element(s)
        double tieRank = (tieStartRank + v_sort.size() - 1) / 2 + 1;
        for (unsigned j = tieStartRank; j < v_sort.size(); j++)
            X.col(c)[v_sort[j].second] = tieRank;
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
    if (n > 1) ties.push_back(n);
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
        n_zero = nrow - (p[i+1] - p[i]);
        if (p[i+1] == p[i]) {
            ties[i].push_back(n_zero);
            continue;
        }
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
//
// // [[Rcpp::export]]
// arma::mat rowNNZAggr_dense(const arma::mat& X,
//                            const arma::uvec& groups,
//                            unsigned ngroups) {
//     arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_cols);
//     for (unsigned c = 0; c < X.n_cols; c++) {
//         for (unsigned r = 0; r < X.n_rows; r++) {
//             if (X(r, c) != 0) res(groups[r], c)++;
//         }
//     }
//     return res;
// }

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

