/*
 * Copyright (c) 2020-2022, Dany Cajas
 * All rights reserved.
 * This work is licensed under BSD 3-Clause "New" or "Revised" License.
 * License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
 */

#include "cppfunctions.h"
#include <cstdio>

using namespace std;
using namespace arma;
using namespace std::chrono;

/**
 * Calculate duplication matrix of size "n" as shown
 * in Magnus, J. R., & Neudecker, H. (1980). The Elimination Matrix:
 * Some Lemmas and Applications. In SIAM Journal on Algebraic Discrete
 * Methods (Vol. 1, Issue 4, pp. 422–449). Society for Industrial &
 * Applied Mathematics (SIAM). https://doi.org/10.1137/0601049
 *
 * @n size of symmetric matrix.
 */
arma::Mat<double> duplication_matrix(const int &n) {
    int cols = (n*(n+1))/2;
    int col;
    arma::sp_mat out = sp_mat(n*n, cols);
    for (int j = 0; j < n; ++j) {
        for (int i = j; i < n; ++i) {
            arma::sp_mat u = sp_mat(1, cols);
            col = j*n+i-((j+1)*j)/2;
            u(0, col) = 1.0;
            arma::sp_mat T = sp_mat(n,n);
            T(i,j) = 1.0;
            T(j,i) = 1.0;
            out += arma::vectorise(T) * u;
        }
    }
    arma::Mat<double> D = arma::mat(out);
    return D;
}

/**
 * Calculate duplication elimination matrix of size "n" as shown
 * in Magnus, J. R., & Neudecker, H. (1980). The Elimination Matrix:
 * Some Lemmas and Applications. In SIAM Journal on Algebraic Discrete
 * Methods (Vol. 1, Issue 4, pp. 422–449). Society for Industrial &
 * Applied Mathematics (SIAM). https://doi.org/10.1137/0601049
 *
 * @n size of symmetric matrix.
 */
arma::Mat<double> duplication_elimination_matrix(const int &n) {
    int rows = (n*(n+1))/2;
    int row;
    arma::sp_mat out = sp_mat(rows, n*n);
    for (int j = 0; j < n; ++j) {
        arma::sp_mat e_j = sp_mat(1, n);
        e_j(0, j) = 1.0;
        for (int i = j; i < n; ++i) {
            arma::sp_mat u = sp_mat(rows, 1);
            row = j*n+i-((j+1)*j)/2;
            u(row, 0) = 1.0;
            arma::sp_mat e_i = sp_mat(1, n);
            e_i(0, i) = 1.0;
            out += arma::kron(u, arma::kron(e_j, e_i));
        }
    }
    arma::Mat<double> L = arma::mat(out);
    return L;
}

/**
 * Calculates summation matrix of size "n" as shown in Cajas, D. (2022).
 * Convex Optimization of Portfolio Kurtosis. In SSRN Electronic Journal.
 * Elsevier BV. https://doi.org/10.2139/ssrn.4202967
 *
 * @n size of symmetric matrix.
 */
arma::Mat<double> summation_matrix(const int &n) {
    arma::mat D;
    arma::mat L;
    D = duplication_matrix(n);
    L = duplication_elimination_matrix(n);
    arma::sp_mat D0;
    arma::sp_mat L0;
    arma::sp_mat S0;
    D0 = arma::sp_mat(D);
    L0= arma::sp_mat(L);
    S0 = D0.t() * D0 * L0;
    arma::Mat<double> S;
    S = arma::mat(S0);
    return S;
}

/**
 * Calculates kurtosis square matrix as shown in Cajas, D. (2022).
 * Convex Optimization of Portfolio Kurtosis. In SSRN Electronic Journal.
 * Elsevier BV. https://doi.org/10.2139/ssrn.4202967
 *
 * @Y matrix of assets returns.
 */
arma::Mat<double> kurtosis_matrix(const arma::Mat<double> &Y) {
    double T = Y.n_rows;
    arma::mat mu = arma::mean(Y);
    arma::mat mu_mat(arma::size(Y), fill::zeros);
    mu_mat.each_row() = mu; 
    arma::mat X = Y - mu_mat;
    arma::mat ones(1, X.n_cols, fill::ones);
    arma::Mat<double> M4 = arma::kron(ones, X) % arma::kron(X, ones);
    M4 = 1/T * M4.t() * M4;
    return M4;
}

/**
 * Calculates semi-kurtosis square matrix as shown in Cajas, D. (2022).
 * Convex Optimization of Portfolio Kurtosis. In SSRN Electronic Journal.
 * Elsevier BV. https://doi.org/10.2139/ssrn.4202967
 *
 * @Y matrix of assets returns.
 */
arma::Mat<double> semi_kurtosis_matrix(const arma::Mat<double> &Y) {
    double T = Y.n_rows;
    arma::mat mu = arma::mean(Y);
    arma::mat mu_mat(arma::size(Y), fill::zeros);
    arma::mat zeros_mat(arma::size(Y), fill::zeros);
    mu_mat.each_row() = mu;
    arma::mat X = Y - mu_mat;
    X = arma::min(X, zeros_mat);
    arma::mat ones(1, X.n_cols, fill::ones);
    arma::Mat<double> SM4 = arma::kron(ones, X) % arma::kron(X, ones);
    SM4 = 1/T * SM4.t() * SM4;
    return SM4;
}