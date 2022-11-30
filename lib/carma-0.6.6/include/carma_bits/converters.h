/*  carma/converters.h: Coverter of Numpy arrays and Armadillo matrices
 *  Copyright (c) 2020 Ralph Urlus <rurlus.dev@gmail.com>
 *  All rights reserved. Use of this source code is governed by a
 *  Apache-2.0 license that can be found in the LICENSE file.
 *
 *  Adapated from:
 *
 *      pybind11/eigen.h: Transparent conversion for dense and sparse Eigen matrices
 *      Copyright (c) 2016 Wolf Vollprecht <w.vollprecht@gmail.com>
 *                         Wenzel Jakob <wenzel.jakob@epfl.ch>
 *      All rights reserved. Use of this source code is governed by a
 *      BSD-style license that can be found in the pybind11/LICENSE file.
 *
 *      arma_wrapper/arma_wrapper.h:
 *      Copyright (C) 2019 Paul Sangrey governed by Apache 2.0 License
 */
#ifndef INCLUDE_CARMA_BITS_CONVERTERS_H_
#define INCLUDE_CARMA_BITS_CONVERTERS_H_
#include <pybind11/buffer_info.h>  // NOLINT
#include <pybind11/detail/common.h>  // NOLINT
#include <pybind11/numpy.h>  // NOLINT
#include <pybind11/pybind11.h>  // NOLINT

/* carma headers */
#include <carma_bits/debug.h>  // NOLINT
#include <carma_bits/typecheck.h>  // NOLINT
#include <carma_bits/cnumpy.h>  // NOLINT
#include <carma_bits/nparray.h>  // NOLINT
#include <carma_bits/numpytoarma.h>  // NOLINT
#include <carma_bits/armatonumpy.h>  // NOLINT
#include <carma_bits/config.h> // NOLINT
#include <carma_bits/exceptions.h> // NOLINT

#include <armadillo>  // NOLINT

#include <memory>
#include <type_traits>
#include <utility>

namespace py = pybind11;


namespace carma {

namespace details {
using aconf = arma::arma_config;
}  // namespace details

/*****************************************************************************************
 *                                   Numpy to Armadillo                                   *
 *****************************************************************************************/

/* Convert numpy array to Armadillo Matrix with copy
 * If the array is 1D we create a column oriented matrix (N, 1) */
template <typename T>
inline arma::Mat<T> arr_to_mat(const py::array_t<T>& src) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_mat<T>(info);
    // copy and ensure fortran order
    data = details::steal_copy_array<T>(src.ptr());
    return details::arr_to_mat(info, data, true, false);
} /* arr_to_mat */

/* Convert numpy array to Armadillo Matrix as view
 * If the array is 1D we create a column oriented matrix (N, 1) */
template <typename T>
inline const arma::Mat<T> arr_to_mat_view(const py::array_t<T>& src) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_mat<T>(info);
    PyObject* obj = src.ptr();
    if (!well_conditioned(obj)) {
        data = details::steal_copy_array<T>(obj);
        return details::arr_to_mat(info, data, true, true);
    }
    return details::arr_to_mat(info, data, false, true);
} /* arr_to_mat_view */

/* Convert numpy array to Armadillo Matrix by stealing the data
 *
 * We copy the array if:
 * - ndim == 2 && (optional) not F contiguous memory
 * - memory is not aligned
 * - writeable is false
 * - owndata is false (optional)
 *
 * If the array is 1D we create a column oriented matrix (N, 1) */
template <typename T>
inline arma::Mat<T> arr_to_mat(py::array_t<T>&& src) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_mat<T>(info);
    // steal memory and copy if needed
    data = details::steal_andor_copy(src.ptr(), data);
    return details::arr_to_mat(info, data, true, false);
} /* arr_to_mat */

/* Convert numpy array to Armadillo Matrix
 *
 * The default behaviour is to borrow the array, we copy if:
 * - copy is true
 * - writeable is false
 * - memory is not aligned
 * - ndim == 2 && (optional) not F contiguous memory
 * - owndata is false (optional)
 *
 * Note that the user set behaviour is overridden is one of the above conditions
 * is true
 *
 * If the array is 1D we create a column oriented matrix (N, 1) */
template <typename T>
inline arma::Mat<T> arr_to_mat(py::array_t<T>& src, bool copy = false) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_mat<T>(info);
    PyObject* obj = src.ptr();
#ifdef CARMA_EXTRA_DEBUG
    if (!well_behaved(obj)) {
        debug::print_copy_of_data<T>(data, obj);
    }
#endif
    if (copy) {
        // copy and ensure fortran order
        data = details::steal_copy_array<T>(obj);
        return details::arr_to_mat(info, data, true, false);
    }
    if (!well_conditioned(obj)) {
        // copy and ensure fortran order and swap with src array
        data = details::swap_copy_array<T>(obj);
    }
    return details::arr_to_mat(info, data, false, true);
} /* arr_to_mat */


// #########################################################################
//                                   COL
// #########################################################################

/* Convert numpy array to Armadillo Column with copy */
template <typename T>
inline arma::Col<T> arr_to_col(const py::array_t<T>& src) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_col<T>(info);
    // copy and ensure fortran order
    data = details::steal_copy_array<T>(src.ptr());
    return details::arr_to_col(info, data, true, false);
} /* arr_to_col */

/* Convert numpy array to Armadillo Col as view
 * If the array is 1D we create a column oriented matrix (N, 1) */
template <typename T>
inline const arma::Col<T> arr_to_col_view(const py::array_t<T>& src) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_col<T>(info);
    PyObject* obj = src.ptr();
    if (!well_conditioned(obj)) {
        data = details::steal_copy_array<T>(obj);
        return details::arr_to_col(info, data, true, true);
    }
    return details::arr_to_col(info, data, false, true);
} /* arr_to_col_view */

/* Convert numpy array to Armadillo Column by stealing the data
 *
 * We copy the array if:
 * - writeable is false
 * - owndata is false
 * - memory is not aligned
 * - NOTE if platform is windows
 */
template <typename T>
inline arma::Col<T> arr_to_col(py::array_t<T>&& src) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_col<T>(info);
    // steal memory and copy if needed
    data = details::steal_andor_copy(src.ptr(), data);
    return details::arr_to_col(info, data, true, false);
} /* arr_to_col */

/* Convert numpy array to Armadillo Col
 *
 * The default behaviour is to borrow the array, we copy if:
 * - copy is true
 * - writeable is false
 * - owndata is false
 * - memory is not aligned
 * Note that the user set behaviour is overridden is one of the above conditions
 * is true
 */
template <typename T>
inline arma::Col<T> arr_to_col(py::array_t<T>& src, bool copy = false) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_col<T>(info);
    PyObject* obj = src.ptr();
#ifdef CARMA_EXTRA_DEBUG
    if (!well_behaved(obj)) {
        debug::print_copy_of_data<T>(data, obj);
    }
#endif
    if (copy) {
        // copy and ensure fortran order
        data = details::steal_copy_array<T>(obj);
        return details::arr_to_col(info, data, true, false);
    }
    if (!well_conditioned(obj)) {
        // copy and ensure fortran order and swap with src array
        data = details::swap_copy_array<T>(obj);
    }
    return details::arr_to_col(info, data, false, true);
} /* arr_to_col */

// #########################################################################
//                                   ROW
// #########################################################################

/* Convert numpy array to Armadillo Row with copy */
template <typename T>
inline arma::Row<T> arr_to_row(const py::array_t<T>& src) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_row<T>(info);
    // copy and ensure fortran order
    data = details::steal_copy_array<T>(src.ptr());
    return details::arr_to_row(info, data, true, false);
} /* arr_to_row */

/* Convert numpy array to Armadillo Row as view
 * If the array is 1D we create a column oriented matrix (N, 1) */
template <typename T>
inline const arma::Row<T> arr_to_row_view(const py::array_t<T>& src) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_row<T>(info);
    PyObject* obj = src.ptr();
    if (!well_conditioned(obj)) {
        data = details::steal_copy_array<T>(obj);
        return details::arr_to_row(info, data, true, true);
    }
    return details::arr_to_row(info, data, false, true);
} /* arr_to_row_view */

/* Convert numpy array to Armadillo Row by stealing the data
 *
 * We copy the array if:
 * - writeable is false
 * - owndata is false
 * - memory is not aligned
 * - NOTE if platform is windows
 */
template <typename T>
inline arma::Row<T> arr_to_row(py::array_t<T>&& src) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_row<T>(info);
    // steal memory and copy if needed
    data = details::steal_andor_copy(src.ptr(), data);
    return details::arr_to_row(info, data, true, false);
} /* arr_to_row */

/* Convert numpy array to Armadillo Row
 *
 * The default behaviour is to borrow the array, we copy if:
 * - copy is true
 * - writeable is false
 * - owndata is false
 * - memory is not aligned
 * Note that the user set behaviour is overridden is one of the above conditions
 * is true
 */
template <typename T>
inline arma::Row<T> arr_to_row(py::array_t<T>& src, bool copy = false) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_row<T>(info);
    PyObject* obj = src.ptr();
#ifdef CARMA_EXTRA_DEBUG
    if (!well_behaved(obj)) {
        debug::print_copy_of_data<T>(data, obj);
    }
#endif
    if (copy) {
        // copy and ensure fortran order
        data = details::steal_copy_array<T>(obj);
        return details::arr_to_row(info, data, true, false);
    }
    if (!well_conditioned(obj)) {
        // copy and ensure fortran order and swap with src array
        data = details::swap_copy_array<T>(obj);
    }
    return details::arr_to_row(info, data, false, true);
} /* arr_to_row */

// #########################################################################
//                                   Cube
// #########################################################################

/* Convert numpy array to Armadillo Cube with copy */
template <typename T>
inline arma::Cube<T> arr_to_cube(const py::array_t<T>& src) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_cube<T>(info);
    // copy and ensure fortran order
    data = details::steal_copy_array<T>(src.ptr());
    return details::arr_to_cube(info, data, true, false);
} /* arr_to_cube */

/* Convert numpy array to Armadillo Cube as view
 * If the array is 1D we create a column oriented matrix (N, 1) */
template <typename T>
inline const arma::Cube<T> arr_to_cube_view(const py::array_t<T>& src) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_cube<T>(info);
    PyObject* obj = src.ptr();
    if (!well_conditioned(obj)) {
        data = details::steal_copy_array<T>(obj);
        return details::arr_to_cube(info, data, true, true);
    }
    return details::arr_to_cube(info, data, false, true);
} /* arr_to_cube_view */

/* Convert numpy array to Armadillo Cube by stealing
 *
 * We copy if:
 * - (optional) not F contiguous memory
 * - writeable is false
 * - owndata is false
 * - memory is not aligned
 * - NOTE if platform is windows
 * Note that the user set behaviour is overridden is one of the above conditions
 * is true
 */
template <typename T>
inline arma::Cube<T> arr_to_cube(py::array_t<T>&& src) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_cube<T>(info);
    // steal memory and copy if needed
    data = details::steal_andor_copy(src.ptr(), data);
    return details::arr_to_cube(info, data, true, false);
} /* arr_to_cube */

/* Convert numpy array to Armadillo Cube
 *
 * The default behaviour is to borrow the array, we copy if:
 * - copy is true
 * - (optional) not F contiguous memory
 * - writeable is false
 * - owndata is false
 * - memory is not aligned
 * Note that the user set behaviour is overridden is one of the above conditions
 * is true
 */
template <typename T>
inline arma::Cube<T> arr_to_cube(py::array_t<T>& src, bool copy = false) {
    py::buffer_info info = src.request();
    T* data = details::validate_from_array_cube<T>(info);
    PyObject* obj = src.ptr();
#ifdef CARMA_EXTRA_DEBUG
    if (!well_behaved(obj)) {
        debug::print_copy_of_data<T>(data, obj);
    }
#endif
    if (copy) {
        // copy and ensure fortran order
        data = details::steal_copy_array<T>(obj);
        return details::arr_to_cube(info, data, true, false);
    }
    if (!well_conditioned(obj)) {
        // copy and ensure fortran order and swap with src array
        data = details::swap_copy_array<T>(obj);
    }
    return details::arr_to_cube(info, data, false, true);
} /* arr_to_cube */

// #########################################################################
//                                   TO_ARMA
// #########################################################################
/* The below functor approach is ported from:
 *     Arma_Wrapper - Paul Sangrey 2019
 *     Apache 2.0 License
 * This is a templated functor that has overloads that convert the various
 * types that I want to pass from Python to C++.
 */
template <typename returnT, typename SFINAE = std::true_type>
struct to_arma {
    static_assert(!SFINAE::value, "The general case is not defined.");
    template <typename innerT>
    static returnT from(innerT&&);
}; /* to_arma */

template <typename returnT>
struct to_arma<returnT, typename is_row<returnT>::type> {
    /* Overload concept on return type; convert to row */
    static returnT from(const py::array_t<typename returnT::elem_type>& arr) {
        return arr_to_row<typename returnT::elem_type>(arr, true, false);
    }
    static returnT from(py::array_t<typename returnT::elem_type>& arr, bool copy) {
        return arr_to_row<typename returnT::elem_type>(arr, copy);
    }
    static returnT from(py::array_t<typename returnT::elem_type>&& arr) {
        return arr_to_row<typename returnT::elem_type>(std::move(arr));
    }
}; /* to_arma */

template <typename returnT>
struct to_arma<returnT, typename is_col<returnT>::type> {
    /* Overload concept on return type; convert to col */
    static returnT from(const py::array_t<typename returnT::elem_type>& arr) {
        return arr_to_col<typename returnT::elem_type>(arr, true, false);
    }
    static returnT from(py::array_t<typename returnT::elem_type>& arr, bool copy) {
        return arr_to_col<typename returnT::elem_type>(arr, copy);
    }
    static returnT from(py::array_t<typename returnT::elem_type>&& arr) {
        return arr_to_col<typename returnT::elem_type>(std::move(arr));
    }
}; /* to_arma */

template <typename returnT>
struct to_arma<returnT, typename is_mat<returnT>::type> {
    /* Overload concept on return type; convert to matrix */
    static returnT from(const py::array_t<typename returnT::elem_type>& arr) {
        return arr_to_mat<typename returnT::elem_type>(arr, true, false);
    }
    static returnT from(py::array_t<typename returnT::elem_type>& arr, bool copy) {
        return arr_to_mat<typename returnT::elem_type>(arr, copy);
    }
    static returnT from(py::array_t<typename returnT::elem_type>&& arr) {
        return arr_to_mat<typename returnT::elem_type>(std::move(arr));
    }
}; /* to_arma */

template <typename returnT>
struct to_arma<returnT, typename is_cube<returnT>::type> {
    /* Overload concept on return type; convert to cube */
    static returnT from(const py::array_t<typename returnT::elem_type>& arr) {
        return arr_to_cube<typename returnT::elem_type>(arr, true, false);
    }
    static returnT from(py::array_t<typename returnT::elem_type>& arr, bool copy) {
        return arr_to_cube<typename returnT::elem_type>(arr, copy);
    }
    static returnT from(py::array_t<typename returnT::elem_type>&& arr) {
        return arr_to_cube<typename returnT::elem_type>(std::move(arr));
    }
}; /* to_arma */

/*****************************************************************************************
 *                                   Armadillo to Numpy                                   *
 *****************************************************************************************/

/* ######################################## Row ######################################## */

template <typename T>
inline py::array_t<T> row_to_arr(const arma::Row<T>& src) {
    /* Convert armadillo row to numpy array */
    auto data = new arma::Row<T>(src);
    return details::construct_array<T>(data);
} /* row_to_arr */

template <typename T>
inline py::array_t<T> row_to_arr(arma::Row<T>&& src) {
    /* Convert armadillo row to numpy array */
    auto data = new arma::Row<T>(std::move(src));
    return details::construct_array<T>(data);
} /* row_to_arr */

template <typename T>
inline py::array_t<T> row_to_arr(arma::Row<T>& src, int copy = false) {
    /* Convert armadillo row to numpy array */
    arma::Row<T>* data;
    if (!copy) {
        data = new arma::Row<T>(std::move(src));
    } else {
        data = new arma::Row<T>(src.memptr(), src.n_elem, true);
    }
    return details::construct_array<T>(data);
} /* row_to_arr */

template <typename T>
inline py::array_t<T> row_to_arr(arma::Row<T>* src, int copy = false) {
    /* Convert armadillo row to numpy array */
    arma::Row<T>* data;
    if (!copy) {
        data = new arma::Row<T>(std::move(*src));
    } else {
        data = new arma::Row<T>(src->memptr(), src->n_elem, true);
    }
    return details::construct_array<T>(data);
} /* row_to_arr */

/* ######################################## Col ######################################## */

template <typename T>
inline py::array_t<T> col_to_arr(const arma::Col<T>& src) {
    /* Convert armadillo col to numpy array */
    auto data = new arma::Col<T>(src);
    return details::construct_array<T>(data);
} /* col_to_arr */

template <typename T>
inline py::array_t<T> col_to_arr(arma::Col<T>&& src) {
    /* Convert armadillo col to numpy array */
    auto data = new arma::Col<T>(std::move(src));
    return details::construct_array<T>(data);
} /* col_to_arr */

template <typename T>
inline py::array_t<T> col_to_arr(arma::Col<T>& src, int copy = false) {
    /* Convert armadillo col to numpy array */
    arma::Col<T>* data;
    if (!copy) {
        data = new arma::Col<T>(std::move(src));
    } else {
        data = new arma::Col<T>(src.memptr(), src.n_elem, true);
    }
    return details::construct_array<T>(data);
} /* col_to_arr */

template <typename T>
inline py::array_t<T> col_to_arr(arma::Col<T>* src, int copy = 0) {
    /* Convert armadillo col to numpy array */
    arma::Col<T>* data;
    if (!copy) {
        data = new arma::Col<T>(std::move(*src));
    } else {
        data = new arma::Col<T>(src->memptr(), src->n_elem, true);
    }
    return details::construct_array<T>(data);
} /* col_to_arr */

/* ######################################## Mat ######################################## */

template <typename T>
inline py::array_t<T> mat_to_arr(const arma::Mat<T>& src) {
    return details::construct_array<T>(new arma::Mat<T>(src));
} /* mat_to_arr */

template <typename T>
inline py::array_t<T> mat_to_arr(arma::Mat<T>&& src) {
    return details::construct_array<T>(new arma::Mat<T>(std::move(src)));
} /* mat_to_arr */

template <typename T>
inline py::array_t<T> mat_to_arr(arma::Mat<T>& src, int copy = 0) {
    if (!copy) {
        return details::construct_array<T>(new arma::Mat<T>(std::move(src)));
    }
    return details::construct_array<T>(
        new arma::Mat<T>(src.memptr(), src.n_rows, src.n_cols, true)
    );
} /* mat_to_arr */

template <typename T>
inline py::array_t<T> mat_to_arr(arma::Mat<T>* src, int copy = 0) {
    if (!copy) {
        return details::construct_array<T>(new arma::Mat<T>(std::move(*src)));
    }
    return details::construct_array(
        new arma::Mat<T>(src->memptr(), src->n_rows, src->n_cols, true)
    );
} /* mat_to_arr */

/* ######################################## Cube ######################################## */

template <typename T>
inline py::array_t<T> cube_to_arr(const arma::Cube<T>& src) {
    auto data = new arma::Cube<T>(src);
    return details::construct_array<T>(data);
} /* cube_to_arr */

template <typename T>
inline py::array_t<T> cube_to_arr(arma::Cube<T>&& src) {
    auto data = new arma::Cube<T>(std::move(src));
    return details::construct_array<T>(data);
} /* cube_to_arr */

template <typename T>
inline py::array_t<T> cube_to_arr(arma::Cube<T>& src, int copy = 0) {
    arma::Cube<T>* data;
    if (!copy) {
        data = new arma::Cube<T>(std::move(src));
    } else {
        data = new arma::Cube<T>(src.memptr(), src.n_rows, src.n_cols, src.n_slices, true);
    }
    return details::construct_array<T>(data);
} /* cube_to_arr */

template <typename T>
inline py::array_t<T> cube_to_arr(arma::Cube<T>* src, int copy = 0) {
    arma::Cube<T>* data;
    if (!copy) {
        data = new arma::Cube<T>(std::move(*src));
    } else {
        data = new arma::Cube<T>(src->memptr(), src->n_rows, src->n_cols, src->n_slices, true);
    }
    return details::construct_array<T>(data);
} /* cube_to_arr */

/* ---------------------------------- to_numpy ---------------------------------- */
template <typename armaT, typename T = typename armaT::elem_type, is_Cube<armaT> = 0>
inline py::array_t<T> to_numpy(const armaT& src) {
    auto data = new arma::Cube<T>(src);
    return details::construct_array<T>(data);
} /* cube_to_arr */

template <typename armaT, typename T = typename armaT::elem_type, is_Cube<armaT> = 0>
inline py::array_t<T> to_numpy(armaT&& src) {
    auto data = new arma::Cube<T>(std::forward<arma::Cube<T>>(src));
    return details::construct_array<T>(data);
} /* cube_to_arr */

template <typename armaT, typename T = typename armaT::elem_type, is_Cube<armaT> = 0>
inline py::array_t<T> to_numpy(armaT& src, int copy = 0) {
    arma::Cube<T>* data;
    if (!copy) {
        data = new arma::Cube<T>(std::move(src));
    } else {
        data = new arma::Cube<T>(src.memptr(), src.n_rows, src.n_cols, src.n_slices, true);
    }
    return details::construct_array<T>(data);
} /* cube_to_arr */

template <typename armaT, typename T = typename armaT::elem_type, is_Cube<armaT> = 0>
inline py::array_t<T> to_numpy(armaT* src, int copy = 0) {
    arma::Cube<T>* data;
    if (!copy) {
        data = new arma::Cube<T>(std::move(*src));
    } else {
        data = new arma::Cube<T>(src->memptr(), src->n_rows, src->n_cols, src->n_slices, true);
    }
    return details::construct_array<T>(data);
} /* cube_to_arr */

template <typename armaT, typename T = typename armaT::elem_type, is_Mat<armaT> = 1>
inline py::array_t<T> to_numpy(const armaT& src) {
    // use armadillo copy constructor
    auto data = new armaT(src);
    return details::construct_array<T>(data);
} /* to_numpy */

template <typename armaT, typename T = typename armaT::elem_type, is_Mat<armaT> = 1>
inline py::array_t<T> to_numpy(armaT&& src) {
    // steal mem
    auto data = new armaT(std::forward<armaT>(src));
    return details::construct_array<T>(data);
} /* to_numpy */

template <typename armaT, typename T = typename armaT::elem_type, is_Mat_only<armaT> = 2>
inline py::array_t<T> to_numpy(armaT& src, int copy = 0) {
    // if not copy we steal
    armaT* data;
    if (!copy) {
        data = new armaT(std::move(src));
    } else {
        data = new armaT(src.memptr(), src.n_rows, src.n_cols, true);
    }
    return details::construct_array<T>(data);
} /* to_numpy */

template <typename armaT, typename T = typename armaT::elem_type, is_Mat_only<armaT> = 2>
inline py::array_t<T> to_numpy(armaT* src, int copy = 0) {
    // if not copy we steal
    armaT* data;
    if (!copy) {
        data = new armaT(std::move(*src));
    } else {
        data = new armaT(src->memptr(), src->n_rows, src->n_cols, true);
    }
    return details::construct_array<T>(data);
} /* to_numpy */

template <typename armaT, typename T = typename armaT::elem_type, is_Vec<armaT> = 3>
inline py::array_t<T> to_numpy(armaT& src, int copy = 0) {
    // if not copy we steal
    armaT* data;
    if (!copy) {
        data = new armaT(std::move(src));
    } else {
        data = new armaT(src.memptr(), src.n_elem, true);
    }
    return details::construct_array<T>(data);
} /* to_numpy */

template <typename armaT, typename T = typename armaT::elem_type, is_Vec<armaT> = 3>
inline py::array_t<T> to_numpy(armaT* src, int copy = 0) {
    // if not copy we steal
    armaT* data;
    if (!copy) {
        data = new armaT(std::move(*src));
    } else {
        data = new armaT(src->memptr(), src->n_elem, true);
    }
    return details::construct_array<T>(data);
} /* to_numpy */

/* ---------------------------------- to_numpy_view ---------------------------------- */
template <typename armaT, typename T = typename armaT::elem_type>
inline py::array_t<T> to_numpy_view(const armaT& src) {
    const armaT* data;
    if (src.n_elem > details::aconf::mat_prealloc) {
        data = &src;
    } else {
        data = new armaT(src);
    }
    return details::construct_array<T>(data);
} /* to_numpy_view */

}  // namespace carma

namespace pybind11 {
namespace detail {

template <typename armaT>
struct type_caster<armaT, enable_if_t<carma::is_convertible<armaT>::value>> {
    using T = typename armaT::elem_type;

    /* Convert numpy array to Armadillo Matrix
     *
     * The default behaviour is to avoid copying, we copy if:
     * - ndim == 2 && (optional) not F contiguous memory
     * - writeable is false
     * - owndata is false
     * - memory is not aligned
     * Note that the user set behaviour is overridden is one of the above conditions
     * is true
     *
     * If the array is 1D we create a column oriented matrix (N, 1) */
    bool load(handle src, bool) {
        // set as array buffer
        py::array_t<T> buffer = py::array_t<T>::ensure(src);
        if (!buffer) {
            throw carma::ConversionError("CARMA: Input cannot be interpreted as array.");
        }
        // borrow the array
        armaT tmp = carma::to_arma<armaT>::from(buffer, false);
        // meet conditions that allow armadillo to steal the matrix
        arma::access::rw(tmp.n_alloc) = tmp.n_elem;
        arma::access::rw(tmp.mem_state) = 0;
        // move the created matrix into the empty but instantiated default
        value = std::move(tmp);
        // reset settings such that the array is in a borrow state
        arma::access::rw(value.n_alloc) = 0;
        arma::access::rw(value.mem_state) = 2;
        return true;
    }

 private:
    // Cast implementation
    static handle cast_impl(const armaT& src, return_value_policy policy, handle) {
        switch (policy) {
            case return_value_policy::automatic:
                return carma::to_numpy<armaT>(src).release();
            case return_value_policy::copy:
                return carma::to_numpy<armaT>(src).release();
            default:
                throw cast_error("unhandled return_value_policy");
        }
    }

    static handle cast_impl(armaT&& src, return_value_policy policy, handle) {
        switch (policy) {
            case return_value_policy::move:
                return carma::to_numpy<armaT>(std::move(src)).release();
            case return_value_policy::automatic:
                return carma::to_numpy<armaT>(std::move(src)).release();
            case return_value_policy::take_ownership:
                return carma::to_numpy<armaT>(std::move(src)).release();
            case return_value_policy::copy:
                return carma::to_numpy<armaT>(src, true).release();
            default:
                throw cast_error("unhandled return_value_policy");
        }
    }

    static handle cast_impl(armaT* src, return_value_policy policy, handle) {
        switch (policy) {
            case return_value_policy::move:
                return carma::to_numpy<armaT>(src).release();
            case return_value_policy::automatic:
                return carma::to_numpy<armaT>(src).release();
            case return_value_policy::take_ownership:
                return carma::to_numpy<armaT>(src).release();
            case return_value_policy::copy:
                return carma::to_numpy<armaT>(src, true).release();
            default:
                throw cast_error("unhandled return_value_policy");
        }
    }

 public:
    // Normal returned non-reference, non-const value: we steal
    static handle cast(armaT&& src, return_value_policy policy, handle parent) {
        return cast_impl(std::move(src), policy, parent);
    }
    // If you return a non-reference const; we copy
    static handle cast(const armaT&& src, return_value_policy policy, handle parent) {
        policy = return_value_policy::copy;
        return cast_impl(src, policy, parent);
    }
    // lvalue reference return; default (automatic) becomes steal
    static handle cast(armaT& src, return_value_policy policy, handle parent) {
        return cast_impl(&src, policy, parent);
    }
    // const lvalue reference return; default (automatic) becomes copy
    static handle cast(const armaT& src, return_value_policy policy, handle parent) {
        policy = return_value_policy::copy;
        return cast_impl(src, policy, parent);
    }
    // non-const pointer return; we steal
    static handle cast(armaT* src, return_value_policy policy, handle parent) {
        return cast_impl(src, policy, parent);
    }
    // const pointer return; we copy
    static handle cast(const armaT* src, return_value_policy policy, handle parent) {
        policy = return_value_policy::copy;
        return cast_impl(*src, policy, parent);
    }

    PYBIND11_TYPE_CASTER(armaT, _("Numpy.ndarray[") + npy_format_descriptor<T>::name + _("]"));
};
} /* namespace detail */
} /* namespace pybind11 */
#endif  // INCLUDE_CARMA_BITS_CONVERTERS_H_
