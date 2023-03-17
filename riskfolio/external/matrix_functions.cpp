/*
 * Copyright (c) 2020-2023, Dany Cajas
 * All rights reserved.
 * This work is licensed under BSD 3-Clause "New" or "Revised" License.
 * License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
 */

#include <iostream>
#include <carma>
#include <armadillo>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

namespace py = pybind11;
using namespace std;
using namespace carma;
using namespace arma;

/**
 * Calculate duplication matrix of size "n" as shown
 * in Magnus, J. R., & Neudecker, H. (1980). The Elimination Matrix:
 * Some Lemmas and Applications. In SIAM Journal on Algebraic Discrete
 * Methods (Vol. 1, Issue 4, pp. 422–449). Society for Industrial &
 * Applied Mathematics (SIAM). https://doi.org/10.1137/0601049
 *
 * @n size of symmetric matrix.
 */
py::array_t<double> cpp_duplication_matrix(const int & n) {
    int cols = (n*(n+1))/2;
    int rows = n*n;
    int col;
    arma::sp_mat out = sp_mat(rows, cols);
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
    return carma::mat_to_arr<double>(D, true);
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
py::array_t<double> cpp_duplication_elimination_matrix(const int &n) {
    int rows = (n*(n+1))/2;
    int cols = n*n;
    int row;
    arma::sp_mat out = sp_mat(rows, cols);
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
    return carma::mat_to_arr<double> (L, true);
}

/**
 * Calculates summation matrix of size "n" as shown in Cajas, D. (2022).
 * Convex Optimization of Portfolio Kurtosis. In SSRN Electronic Journal.
 * Elsevier BV. https://doi.org/10.2139/ssrn.4202967
 *
 * @n size of symmetric matrix.
 */
py::array_t<double> cpp_duplication_summation_matrix(const int &n) {
    arma::mat D;
    arma::mat L;
    D = carma::arr_to_mat(cpp_duplication_matrix(n));
    L = carma::arr_to_mat(cpp_duplication_elimination_matrix(n));
    arma::sp_mat D0;
    arma::sp_mat L0;
    arma::sp_mat S0;
    D0 = arma::sp_mat(D);
    L0= arma::sp_mat(L);
    S0 = D0.t() * D0 * L0;
    arma::mat S;
    S = arma::mat(S0);
    return carma::mat_to_arr<double> (S, true);
}


void bind_duplication_matrix(py::module &m) {
    m.def(
        "cpp_duplication_matrix",
        &cpp_duplication_matrix,
        R"pbdoc(
            Calculate duplication matrix of size "n" as shown in :cite:`d-Magnus1980`.

            Parameters
            ----------
            n : int
                Number of assets.

            Returns
            -------
            D: np.ndarray
                Duplication matrix
        )pbdoc",
        py::arg("n")
        );
}

void bind_duplication_elimination_matrix(py::module &m) {
    m.def(
        "cpp_duplication_elimination_matrix",
        &cpp_duplication_elimination_matrix,
        R"pbdoc(
            Calculate duplication elimination matrix of size "n" as shown in :cite:`d-Magnus1980`.

            Parameters
            ----------
            n : int
                Number of assets.

            Returns
            -------
            L: np.ndarray
                Duplication matrix
        )pbdoc",
        py::arg("n")
    );
}

void bind_duplication_summation_matrix(py::module &m) {
    m.def(
        "cpp_duplication_summation_matrix",
        &cpp_duplication_summation_matrix,
        R"pbdoc(
            Calculate duplication summation matrix of size "n" as shown in :cite:`d-Cajas4`.

            Parameters
            ----------
            n : int
                Number of assets.

            Returns
            -------
            S: np.ndarray
                Duplication summation matrix.
        )pbdoc",
        py::arg("n")
    );
}
