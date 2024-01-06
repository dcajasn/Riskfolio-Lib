/*
 * Copyright (c) 2020-2024, Dany Cajas
 * All rights reserved.
 * This work is licensed under BSD 3-Clause "New" or "Revised" License.
 * License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
 */

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include "matrix_functions.cpp"

namespace py = pybind11;

PYBIND11_MODULE(functions, m) {
    bind_duplication_matrix(m);
    bind_duplication_elimination_matrix(m);
    bind_duplication_summation_matrix(m);
    bind_commutation_matrix(m);
    bind_coskewness_matrix(m);
    bind_cokurtosis_matrix(m);
    bind_k_eigh(m);
    bind_dcorr(m);
    bind_dcorr_matrix(m);
}
