/*  carma/typechecks.h: Template type checks
 *  Copyright (c) 2020 Ralph Urlus <rurlus.dev@gmail.com>
 *  All rights reserved. Use of this source code is governed by a
 *  Apache-2.0 license that can be found in the LICENSE file.
 */
#ifndef INCLUDE_CARMA_BITS_TYPECHECK_H_
#define INCLUDE_CARMA_BITS_TYPECHECK_H_
#include <armadillo>

#include <memory>
#include <utility>
#include <type_traits>

namespace carma {

using aconf =  arma::arma_config;

// Mat catches Row and Col as well
template <typename T>
struct is_convertible {
    static const bool value = (arma::is_Mat<T>::value || arma::is_Cube<T>::value);
};

template <typename T>
struct p_is_Vec {
    static const bool value = (arma::is_Row<T>::value || arma::is_Col<T>::value);
};

// for reference see: https://www.fluentcpp.com/2019/08/23/how-to-make-sfinae-pretty-and-robust/
template <typename armaT>
using is_Cube = std::enable_if_t<arma::is_Cube<armaT>::value, int>;
template <typename armaT>
using is_Vec = std::enable_if_t<p_is_Vec<armaT>::value, int>;
template <typename armaT>
using is_Mat = std::enable_if_t<arma::is_Mat<armaT>::value, int>;
template <typename armaT>
using is_Mat_only = std::enable_if_t<arma::is_Mat_only<armaT>::value, int>;

template <typename T>
struct is_mat : std::false_type {};
template <typename T>
struct is_mat<arma::Mat<T>> : std::true_type {};

template <typename T>
struct is_col : std::false_type {};
template <typename T>
struct is_col<arma::Col<T>> : std::true_type {};

template <typename T>
struct is_row : std::false_type {};
template <typename T>
struct is_row<arma::Row<T>> : std::true_type {};

template <typename T>
struct is_cube : std::false_type {};
template <typename T>
struct is_cube<arma::Cube<T>> : std::true_type {};

}  // namespace carma
#endif  // INCLUDE_CARMA_BITS_TYPECHECK_H_
