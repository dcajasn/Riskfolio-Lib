/*
 * Copyright (c) 2020-2025, Dany Cajas
 * All rights reserved.
 * This work is licensed under BSD 3-Clause "New" or "Revised" License.
 * License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <map>
#include <Eigen>
#include <Eigen/Core>
#include <Eigen/KroneckerProduct>
#include <Spectra/SymEigsSolver.h>

namespace py = pybind11;
using namespace std;
using namespace Eigen;
using namespace Spectra;

std::vector<int> cumsum(const std::vector<int>& vec) {
    std::vector<int> result(vec.size());
    std::partial_sum(vec.begin(), vec.end(), result.begin());
    return result;
}

/**
 * Calculate duplication matrix of size "n" as shown
 * in Magnus, J. R., & Neudecker, H. (1980). The Elimination Matrix:
 * Some Lemmas and Applications. In SIAM Journal on Algebraic Discrete
 * Methods (Vol. 1, Issue 4, pp. 422–449). Society for Industrial &
 * Applied Mathematics (SIAM). https://doi.org/10.1137/0601049
 *
 * @n size of symmetric matrix.
 */

Eigen::SparseMatrix<double> cpp_duplication_matrix(const int &n, const bool &diag = true) {
    /* Old Version
    Eigen::SparseMatrix<double> out(n*n, (n*(n+1))/2);
    for (int j = 0; j < n; ++j) {
        for (int i = j; i < n; ++i) {
            Eigen::SparseMatrix<double> u(1, (n*(n+1))/2);
            u.insert(0, j*n+i-((j+1)*j)/2) = 1.0;
            Eigen::MatrixXd T0 = Eigen::MatrixXd::Zero(n, n);
            T0(i,j) = 1.0;
            T0(j,i) = 1.0;
            Eigen::SparseMatrix<double> T(n * n, 1);
            T = (T0.reshaped(n * n, 1)).sparseView();
            out += T * u;
        }
    }
    */
    // New version based on https://www.mathworks.com/matlabcentral/answers/473737-efficient-algorithm-for-a-duplication-matrix
    int m = 1;
    int nsq = n * n;
    int r = 0;
    int a = 1;

    std::vector<int> v(nsq, 0);
    std::vector<int> sequence(n - 1); // Create cn which is the cumulative sum of n to 2
    std::iota(sequence.begin(), sequence.end(), 2); // sequence: [2, 3, 4, ... , n]
    std::reverse(sequence.begin(), sequence.end()); // reverse: [n, n-1, n-2, ..., 2]
    std::vector<int> cn = cumsum(sequence); // cumulative sum
    std::vector<int> filtered_rows, filtered_cols;

    for (int i = 1; i <= n; ++i) {
        for (int j = 0; j < i - 1; ++j) {
            v[r + j] = i - n + cn[j];
        }
        r += i - 1;
        for (int j = 0; j < n - i + 1; ++j) {
            v[r + j] = a + j;
        }
        r += n - i + 1;
        a += n - i + 1;
    }

    if (diag==true){
        m = n * (n + 1) / 2;
        for (int i = 0; i < nsq; ++i) {
            filtered_rows.push_back(i);
            filtered_cols.push_back(v[i]-1);
        }
    }else{
        m = n * (n - 1) / 2;
        std::vector<int> rows(nsq), cols(m);
        std::iota(rows.begin(), rows.end(), 0);

        for (size_t i; i < nsq; ++i){
            std::cout << rows[i] <<' ';
        }

        std::unordered_map<int, int> counts;
        for (int i : v) {
            counts[i]++;
        }

        std::unordered_set<int> non_unique_elements;
        for (const auto& pair : counts) {  // Use pair explicitly
            if (pair.second > 1) {
                    non_unique_elements.insert(pair.first);
            }
        }

        std::list<int> new_non_unique_elements;
        for (int i : non_unique_elements) {
            new_non_unique_elements.push_back(i);
        }
        new_non_unique_elements.sort();

        std::unordered_map<int, int> new_cols;
        int currentNumber = 0;
        for (int col : new_non_unique_elements) {
            new_cols[col] = currentNumber++;
        }

        for (int i = 0; i < nsq; ++i) {
            if (non_unique_elements.count(v[i])) {
                filtered_rows.push_back(rows[i]);
                filtered_cols.push_back(new_cols[v[i]]);
            }
        }
    }

    Eigen::SparseMatrix<double> out(nsq, m);
    for (int i = 0; i < filtered_rows.size(); ++i) {
        out.insert(filtered_rows[i], filtered_cols[i]) = 1;
    }
    return out;
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
Eigen::SparseMatrix<double> cpp_duplication_elimination_matrix(const int &n, const bool &diag = true) {
    /* Old Version
    Eigen::SparseMatrix<double> out((n*(n+1))/2, n*n);
    for (int j = 0; j < n; ++j) {
        Eigen::SparseMatrix<double> e_j(1, n);
        e_j.insert(0, j) = 1.0;
        for (int i = j; i < n; ++i) {
            Eigen::SparseMatrix<double> u((n*(n+1))/2, 1);
            u.insert(j*n+i-((j+1)*j)/2, 0) = 1.0;
            Eigen::SparseMatrix<double> e_i(1, n);
            e_i.insert(0, i) = 1.0;
            out += Eigen::kroneckerProduct(u, Eigen::kroneckerProduct(e_j, e_i).eval()).eval();
        }
    }
    */
    // New version based on https://www.mathworks.com/matlabcentral/answers/473737-efficient-algorithm-for-a-duplication-matrix
    
    int m = 1;
    int j = 1;
    int nsq = n * n;
    int r = 0;
    int a = 1;

    std::vector<int> v(nsq, 0);
    std::vector<int> sequence(n); // Create cn which is the cumulative sum of 0 to n-1
    std::iota(sequence.begin(), sequence.end(), 0); // sequence: [0, 1, 2, ... , n-1]
    std::vector<int> cn = cumsum(sequence); // cumulative sum

    if (diag == true){
        m = n * (n + 1) / 2;
    }else{
        m = n * (n - 1) / 2;
        j = 2;
    }

    for (int i = j; i <= n; ++i) {
        r += i - 1;
        // Filling v[r:r + n - i] with a:a + n - i
        for (int j = 0; j < n - i + 1; ++j) {
            v[r + j] = a + j + cn[i-1];
        }
        r += n - i + 1;
        a += n - i + 1;
    }
    v.erase(std::remove(v.begin(), v.end(), 0),v.end());

    Eigen::SparseMatrix<double> out(m, nsq);
    for (int i = 0; i < m; ++i) {
        out.insert(i, v[i] - 1) = 1;
    }
    return out;
}

/**
 * Calculates summation matrix of size "n" as shown in Cajas, D. (2022).
 * Convex Optimization of Portfolio Kurtosis. In SSRN Electronic Journal.
 * Elsevier BV. https://doi.org/10.2139/ssrn.4202967
 *
 * @n size of symmetric matrix.
 */
Eigen::SparseMatrix<double> cpp_duplication_summation_matrix(const int &n, const bool &diag = true) {
    /* Old version
    std::vector<int> a = {0};
    for (int k = 0; k < n-1; ++k) {
        a.push_back(a.back() + n - k);
    }
    int rows = (n*(n+1))/2;
    Eigen::SparseMatrix<double> out(rows, n*n);
    for (int j = 0; j < n; ++j) {
        Eigen::SparseMatrix<double> e_j(1, n);
        e_j.insert(0, j) = 1.0;
        for (int i = j; i < n; ++i) {
            Eigen::SparseMatrix<double> u(rows, 1);
            int row = j*n+i-((j+1)*j)/2;
            auto it = std::find(a.begin(), a.end(), row);
            if (it != a.end()) {
                u.insert(row, 0) = 1.0;
            } else {
                u.insert(row, 0) = 2.0;
            }
            Eigen::SparseMatrix<double> e_i(1, n);
            e_i.insert(0, i) = 1.0;
            out += Eigen::kroneckerProduct(u, Eigen::kroneckerProduct(e_j, e_i).eval()).eval();
        }
    }
    */
    // New version based on https://www.mathworks.com/matlabcentral/answers/473737-efficient-algorithm-for-a-duplication-matrix
    int m = 1;
    int j = 1;
    int nsq = n * n;
    int r = 0;
    int a = 1;

    std::vector<int> v(nsq, 0);
    std::vector<int> v2(nsq, 0);
    std::vector<int> rows2(nsq, 0);
    std::vector<int> sequence(n); // Create cn which is the cumulative sum of 0 to n-1
    std::iota(sequence.begin(), sequence.end(), 0); // sequence: [0, 1, 2, ... , n-1]
    std::vector<int> cn = cumsum(sequence); // cumulative sum

    if (diag == true){
        m = n * (n + 1) / 2;
    }else{
        m = n * (n - 1) / 2;
        j = 2;
    }

    for (int i = j; i <= n; ++i) {
        r += i - 1;
        // Filling v[r:r + n - i] with a:a + n - i
        for (int j = 0; j < n - i + 1; ++j) {
            v[r + j] = a + j + cn[i-1];
        }
        for (int j = 1; j < n - i + 1; ++j) {
            v2[r + j] = a + j + cn[i-1];
            rows2[r + j] = a + j;
        }
        r += n - i + 1;
        a += n - i + 1;
    }

    v.erase(std::remove(v.begin(), v.end(), 0), v.end());
    v2.erase(std::remove(v2.begin(), v2.end(), 0), v2.end());
    rows2.erase(std::remove(rows2.begin(), rows2.end(), 0), rows2.end());

    Eigen::SparseMatrix<double> out1(m, nsq);
    Eigen::SparseMatrix<double> out2(m, nsq);
    Eigen::SparseMatrix<double> out;

    if (diag == true){
        for (int i = 0; i < m; ++i) {
            out1.insert(i, v[i] - 1) = 1;
        }
        for (int i = 0; i < m-n; ++i) {
            out2.insert(rows2[i]-1, v2[i] - 1) = 1;
        }
        out = out1 + out2;
    } else {
        for (int i = 0; i < m; ++i) {
            out1.insert(i, v[i] - 1) = 2;
        }
        out = out1;
    }
    
    return out;
}

/**
 * Calculates the commutation matrix for a matrix of size T x n.
 *
 * @T number of rows.
 * @n number of columns.
 */
Eigen::SparseMatrix<double> cpp_commutation_matrix(const int &T, const int &n) {
    std::vector<int> i(T * n);
    std::iota (std::begin(i), std::end(i), 0);
    std::vector<int> j;
    for (int k = 0; k < T; ++k) {
        for (int l = 0; l < n; ++l) {
            j.push_back(T * l + k);
        }
    }
    Eigen::SparseMatrix<double> out(T * n, T * n);
    for (int idx = 0; idx < T * n; ++idx) {
        out.insert(i[idx], j[idx]) = 1.0;
    }
    return out;
}

/**
 * Calculates coskewness tensor of size "n" as shown in Cajas, D. (2022).
 * Convex Optimization of Portfolio Kurtosis. In SSRN Electronic Journal.
 * Elsevier BV. https://doi.org/10.2139/ssrn.4202967
 *
 * @Y array of size n_samples x n_features.
 * @semi flag that indicates if calculates lower semi skewness rectangular matrix.
 */
Eigen::MatrixXd cpp_coskewness_matrix(Eigen::MatrixXd Y, const bool &semi=false) {
    double T = Y.rows();
    Eigen::VectorXd mu = Y.colwise().mean();
    Eigen::MatrixXd X = Y - mu.transpose().replicate(T, 1);
    if (semi == true){
        X = X.cwiseMin(0);
    }
    Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(1, X.cols());
    Eigen::MatrixXd V1 = Eigen::kroneckerProduct(ones, X).eval();
    Eigen::MatrixXd V2 = Eigen::kroneckerProduct(X, ones).eval();
    Eigen::MatrixXd V = (V1.array() * V2.array()).matrix();
    Eigen::MatrixXd M3 = X.transpose() * V / T;
    return M3;
}

/**
 * Calculates cokurtosis square tensor of size "n" as shown in Cajas, D. (2022).
 * Convex Optimization of Portfolio Kurtosis. In SSRN Electronic Journal.
 * Elsevier BV. https://doi.org/10.2139/ssrn.4202967
 *
 * @Y array of size n_samples x n_features.
 * @semi flag that indicates if calculates lower semi cokurtosis square matrix.
 */
Eigen::MatrixXd cpp_cokurtosis_matrix(Eigen::MatrixXd Y, const bool &semi=false) {
    double T = Y.rows();
    Eigen::VectorXd mu = Y.colwise().mean();
    Eigen::MatrixXd X = Y - mu.transpose().replicate(T, 1);
    if (semi == true){
            X = X.cwiseMin(0);
    }
    Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(1, X.cols());
    Eigen::MatrixXd V1 = Eigen::kroneckerProduct(ones, X).eval();
    Eigen::MatrixXd V2 = Eigen::kroneckerProduct(X, ones).eval();
    Eigen::MatrixXd V = (V1.array() * V2.array()).matrix();
    Eigen::MatrixXd S4 = V.transpose() * V / T;
    return S4;
}

/**
 * Calculates first k eigenvalues and eigenvectors of a symmetric matrix.
 * 
 * @M symmetric matrix.
 * @k number of largest eigenvalues.
 */
std::tuple<Eigen::VectorXd, Eigen::MatrixXd>  cpp_k_eigh(Eigen::MatrixXd M, const int &k) {
    DenseSymMatProd<double> op(M);
    SymEigsSolver<DenseSymMatProd<double>> eigs(op, k, 2 * k);
    eigs.init();
    eigs.compute(SortRule::LargestAlge);
    Eigen::VectorXd eigen_values;
    Eigen::MatrixXd eigen_vectors;
    if(eigs.info() == CompInfo::Successful){
        eigen_values = eigs.eigenvalues();
        eigen_vectors = eigs.eigenvectors();
    }
    return std::make_tuple(eigen_values, eigen_vectors);
}

/**
 * Calculates the euclidean distance matrix of matrix X.
 * 
 * @X A matrix.
 */
Eigen::MatrixXd cpp_pdist(Eigen::MatrixXd& X) {
    int n = X.rows();
    Eigen::MatrixXd distances(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            distances(i, j) = (X.row(i) - X.row(j)).norm();
            distances(j, i) = distances(i, j);
        }
    }
    return distances;
}

/**
 * Calculates the distance correlation of matrices X and Y.
 * 
 * @X A matrix.
 * @Y A matrix.
 */
double cpp_dcorr(Eigen::MatrixXd X, Eigen::MatrixXd Y){
    Eigen::MatrixXd a = cpp_pdist(X);
    Eigen::MatrixXd b = cpp_pdist(Y);

    int n = a.rows();

    // Calculate mean and center the matrices A and B
    Eigen::MatrixXd A = a - (a.colwise().mean()).replicate(n, 1) - (a.rowwise().mean()).replicate(1, n) + a.mean() * Eigen::MatrixXd::Ones(n,n);
    Eigen::MatrixXd B = b - (b.colwise().mean()).replicate(n, 1) - (b.rowwise().mean()).replicate(1, n) + b.mean() * Eigen::MatrixXd::Ones(n,n);

    // Calculate dcov2_xy, dcov2_xx, dcov2_yy
    double dcov2_xy = (A.array() * B.array()).sum() / (n * n);
    double dcov2_xx = (A.array() * A.array()).sum() / (n * n);
    double dcov2_yy = (B.array() * B.array()).sum() / (n * n);

    // Calculate the distance correlation value
    double value = std::sqrt(dcov2_xy) / std::sqrt(std::sqrt(dcov2_xx) * std::sqrt(dcov2_yy));

    return value;
}

/**
 * Calculates the distance correlation matrix of a matrix of variables Y.
 * 
 * @Y A matrix which each column represents a variable.
 */
Eigen::MatrixXd cpp_dcorr_matrix(Eigen::MatrixXd Y){
    int n = Y.cols();
    Eigen::MatrixXd corr = Eigen::MatrixXd::Ones(n,n);

    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            corr(i, j) = cpp_dcorr(Y.col(i), Y.col(j));
            corr(j, i) = corr(i, j);
        }
    }

    return corr;
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
            diag : bool
                Include diagonal elements or not.

            Returns
            -------
            D: np.ndarray
                Duplication matrix
        )pbdoc",
        py::arg("n"),
        py::arg("diag")
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
            diag : bool
                Include diagonal elements or not.

            Returns
            -------
            L: np.ndarray
                Duplication matrix
        )pbdoc",
        py::arg("n"),
        py::arg("diag")
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
            diag : bool
                Include diagonal elements or not.

            Returns
            -------
            S: np.ndarray
                Duplication summation matrix.
        )pbdoc",
        py::arg("n"),
        py::arg("diag")
    );
}

void bind_commutation_matrix(py::module &m) {
    m.def(
        "cpp_commutation_matrix",
        &cpp_commutation_matrix,
        R"pbdoc(
            Calculate commutation matrix of size T x n.

            Parameters
            ----------
            T : int
                Number of rows.
            n : int
                Number of columns.

            Returns
            -------
            K: np.ndarray
                Commutation matrix.
        )pbdoc",
        py::arg("T"),
        py::arg("n")
    );
}

void bind_coskewness_matrix(py::module &m) {
    m.def(
        "cpp_coskewness_matrix",
        &cpp_coskewness_matrix,
        R"pbdoc(
            Calculate coskewness and lower semi coskewness rectangular matrix as shown in :cite:`d-Cajas4`.

            Parameters
            ----------
            Y : ndarray or dataframe
                Returns series of shape n_samples x n_features.

            semi: bool
                Flag that indicates if we calculate coskewness or lower semi coskewness rectangular matrix.

            Returns
            -------
            M3: np.ndarray
                Coskewness or lower semi coskewness rectangular matrix.
        )pbdoc",
        py::arg("n"),
        py::arg("semi")=false
    );
}

void bind_cokurtosis_matrix(py::module &m) {
    m.def(
        "cpp_cokurtosis_matrix",
        &cpp_cokurtosis_matrix,
        R"pbdoc(
            Calculate cokurtosis and lower semi cokurtosis square matrix as shown in :cite:`d-Cajas4`.

            Parameters
            ----------
            Y : ndarray or dataframe
                Returns series of shape n_samples x n_features.

            semi: bool
                Flag that indicates if we calculate cokurtosis or lower semi cokurtosis square matrix.

            Returns
            -------
            S4: np.ndarray
                Cokurtosis or lower semi cokurtosis square matrix.
        )pbdoc",
        py::arg("n"),
        py::arg("semi")=false
    );
}

void bind_k_eigh(py::module &m) {
    m.def(
        "cpp_k_eigh",
        &cpp_k_eigh,
        R"pbdoc(
            Calculate the first largest "k" eigenvalues and eigenvectors.

            Parameters
            ----------
            M : ndarray or dataframe
                A symmetric square matrix.
            k : int
                number of eigenvalues and eigenvectors calculated.

            Returns
            -------
            value: tuple
                Tuple which first element is the eigenvalues vector of M and second element is the eigenvectors matrix of M.
        )pbdoc",
        py::arg("M"),
        py::arg("k")
    );
}

void bind_dcorr(py::module &m) {
    m.def(
        "cpp_dcorr",
        &cpp_dcorr,
        R"pbdoc(
            Calculate the distance correlation.

            Parameters
            ----------
            X : ndarray
                A matrix of variables.
            Y : ndarray
                A matrix of variables.

            Returns
            -------
            value: float
                Distance correlation.
        )pbdoc",
        py::arg("X"),
        py::arg("Y")
    );
}

void bind_dcorr_matrix(py::module &m) {
    m.def(
        "cpp_dcorr_matrix",
        &cpp_dcorr_matrix,
        R"pbdoc(
            Calculate the distance correlation matrix.

            Parameters
            ----------
            Y : ndarray
                A matrix of variables.

            Returns
            -------
            corr: ndarray
                Distance correlation matrix.
        )pbdoc",
        py::arg("Y")
    );
}


