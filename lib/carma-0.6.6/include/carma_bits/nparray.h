/*  carma/nparray.h: Condition checks numpy arrays
 *  Copyright (c) 2020 Ralph Urlus <rurlus.dev@gmail.com>
 *  All rights reserved. Use of this source code is governed by a
 *  Apache-2.0 license that can be found in the LICENSE file.
 *
 *  Adapated from:
 *      pybind11/numpy.h: Basic NumPy support, vectorize() wrapper
 *
 *      Copyright (c) 2016 Wenzel Jakob <wenzel.jakob@epfl.ch>
 *      All rights reserved. Use of this source code is governed by a
 *      BSD-style license that can be found in the LICENSE file.
 */
#ifndef INCLUDE_CARMA_BITS_NPARRAY_H_
#define INCLUDE_CARMA_BITS_NPARRAY_H_

#include <pybind11/numpy.h>  // NOLINT
#include <pybind11/pybind11.h>  // NOLINT

#include <carma_bits/cnumpy.h> // NOLINT
#include <carma_bits/config.h> // NOLINT

#include <memory>
#include <type_traits>
#include <utility>

namespace py = pybind11;

namespace carma {

template <typename T>
inline bool is_f_contiguous(const py::array_t<T>& arr) {
    return details::is_f_contiguous(arr.ptr());
}

template <typename T>
inline bool is_c_contiguous(const py::array_t<T>& arr) {
    return details::is_c_contiguous(arr.ptr());
}

template <typename T>
inline bool is_contiguous(const py::array_t<T>& arr) {
    return is_f_contiguous(arr) || is_c_contiguous(arr);
}

template <typename T>
inline bool is_writeable(const py::array_t<T>& arr) {
    return details::is_writeable(arr.ptr());
}

template <typename T>
inline bool is_owndata(const py::array_t<T>& arr) {
    return details::is_owndata(arr.ptr());
}

template <typename T>
inline bool is_aligned(const py::array_t<T>& arr) {
    return details::is_aligned(arr.ptr());
}

template <typename T>
inline bool is_well_behaved(const py::array_t<T>& arr) {
    return well_behaved(arr.ptr());
}

template <typename T>
inline void set_owndata(py::array_t<T>& arr) {
    details::set_owndata(arr.ptr());
}

template <typename T>
inline py::array_t<T> set_owndata(py::array_t<T>&& arr) {
    details::set_owndata(arr.ptr());
    return std::move(arr);
}

template <typename T>
inline void set_not_owndata(py::array_t<T>& arr) {
    details::set_not_owndata(arr.ptr());
}

template <typename T>
inline py::array_t<T> set_not_owndata(py::array_t<T>&& arr) {
    details::set_not_owndata(arr.ptr());
    return std::move(arr);
}

template <typename T>
inline void set_writeable(py::array_t<T>& arr) {
    details::set_writeable(arr.ptr());
}

template <typename T>
inline py::array_t<T> set_writeable(py::array_t<T>&& arr) {
    details::set_writeable(arr.ptr());
    return std::move(arr);
}

template <typename T>
inline void set_not_writeable(py::array_t<T>& arr) {
    details::set_not_writeable(arr.ptr());
}

template <typename T>
inline py::array_t<T> set_not_writeable(py::array_t<T>&& arr) {
    details::set_not_writeable(arr.ptr());
    return std::move(arr);
}

}  // namespace carma

#endif  // INCLUDE_CARMA_BITS_NPARRAY_H_
