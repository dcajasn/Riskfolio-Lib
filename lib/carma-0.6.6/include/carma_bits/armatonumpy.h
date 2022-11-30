/*  carma/armatonumpy.h: Coverter of Armadillo matrices to numpy arrays
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
#ifndef INCLUDE_CARMA_BITS_ARMATONUMPY_H_
#define INCLUDE_CARMA_BITS_ARMATONUMPY_H_
#include <pybind11/pybind11.h>  // NOLINT
#include <pybind11/numpy.h>  // NOLINT
#include <carma_bits/config.h> // NOLINT
#include <carma_bits/nparray.h> // NOLINT

#include <armadillo>  // NOLINT
#include <utility>

namespace py = pybind11;

namespace carma {
namespace details {

template <typename armaT>
inline py::capsule create_capsule(armaT* data) {
    return py::capsule(data, [](void* f) {
        auto mat = reinterpret_cast<armaT*>(f);
#ifdef CARMA_EXTRA_DEBUG
        std::cout << "\n-----------\nCARMA DEBUG\n-----------" << "\n";
        // if in debug mode let us know what pointer is being freed
        std::cout << "Freeing memory @" << mat->memptr() << std::endl;
        std::cout << "-----------" << "\n";
#endif
        delete mat;
    });
} /* create_capsule */

template <typename armaT>
inline py::capsule create_dummy_capsule(const armaT* data) {
#ifdef CARMA_EXTRA_DEBUG
    return py::capsule(data, [](void* f) {
        auto mat = reinterpret_cast<armaT*>(f);
        std::cout << "\n-----------\nCARMA DEBUG\n-----------" << "\n";
        // if in debug mode let us know what pointer is being freed
        std::cout << "Destructing view on memory @" << mat->memptr() << std::endl;
        std::cout << "-----------" << "\n";
    });
#else
    return py::capsule(data, [](void*) {});
#endif
} /* create_capsule */


template <typename T>
inline py::array_t<T> construct_array(arma::Row<T>* data) {
    constexpr auto tsize = static_cast<ssize_t>(sizeof(T));
    auto ncols = static_cast<ssize_t>(data->n_cols);

    py::capsule base = create_capsule<arma::Row<T>>(data);

    return py::array_t<T>(
        {static_cast<ssize_t>(1), ncols},  // shape
        {tsize, tsize},                    // F-style contiguous strides
        data->memptr(),                    // the data pointer
        base                               // numpy array references this parent
    );
} /* construct_array */

template <typename T>
inline py::array_t<T> construct_array(arma::Col<T>* data) {
    constexpr auto tsize = static_cast<ssize_t>(sizeof(T));
    auto nrows = static_cast<ssize_t>(data->n_rows);

    py::capsule base = create_capsule<arma::Col<T>>(data);

    return py::array_t<T>(
        {nrows, static_cast<ssize_t>(1)},  // shape
        {tsize, nrows * tsize},            // F-style contiguous strides
        data->memptr(),                    // the data pointer
        base                               // numpy array references this parent
    );
} /* construct_array */

template <typename T>
inline py::array_t<T> construct_array(arma::Mat<T>* data) {
    constexpr auto tsize = static_cast<ssize_t>(sizeof(T));
    auto nrows = static_cast<ssize_t>(data->n_rows);
    auto ncols = static_cast<ssize_t>(data->n_cols);

    py::capsule base = create_capsule<arma::Mat<T>>(data);

    return py::array_t<T>(
        {nrows, ncols},          // shape
        {tsize, nrows * tsize},  // F-style contiguous strides
        data->memptr(),          // the data pointer
        base                     // numpy array references this parent
    );
} /* construct_array */

template <typename T>
inline py::array_t<T> construct_array(arma::Cube<T>* data) {
    constexpr auto tsize = static_cast<ssize_t>(sizeof(T));
    auto nrows = static_cast<ssize_t>(data->n_rows);
    auto ncols = static_cast<ssize_t>(data->n_cols);
    auto nslices = static_cast<ssize_t>(data->n_slices);

    return py::array_t<T>(
        // shape
        {nrows, ncols, nslices},
        // F-style contiguous strides
        {tsize, nrows * tsize, tsize * nrows * ncols},
        // the data pointer
        data->memptr(),
        // numpy array references this parent
        create_capsule<arma::Cube<T>>(data)
    );
} /* construct_array */

template <typename T>
inline py::array_t<T> construct_array(const arma::Row<T>* data) {
    constexpr auto tsize = static_cast<ssize_t>(sizeof(T));
    auto ncols = static_cast<ssize_t>(data->n_cols);

    py::capsule base = create_dummy_capsule<arma::Row<T>>(data);

    auto arr = py::array_t<T>(
        {static_cast<ssize_t>(1), ncols},  // shape
        {tsize, tsize},                    // F-style contiguous strides
        data->memptr(),                    // the data pointer
        base                               // numpy array references this parent
    );
    carma::set_not_writeable(arr);
    return arr;
} /* construct_array */

template <typename T>
inline py::array_t<T> construct_array(const arma::Col<T>* data) {
    constexpr auto tsize = static_cast<ssize_t>(sizeof(T));
    auto nrows = static_cast<ssize_t>(data->n_rows);

    py::capsule base = create_dummy_capsule<arma::Col<T>>(data);

    auto arr = py::array_t<T>(
        {nrows, static_cast<ssize_t>(1)},  // shape
        {tsize, nrows * tsize},            // F-style contiguous strides
        data->memptr(),                    // the data pointer
        base                               // numpy array references this parent
    );
    carma::set_not_writeable(arr);
    return arr;
} /* construct_array */

template <typename T>
inline py::array_t<T> construct_array(const arma::Mat<T>* data) {
    constexpr auto tsize = static_cast<ssize_t>(sizeof(T));
    auto nrows = static_cast<ssize_t>(data->n_rows);
    auto ncols = static_cast<ssize_t>(data->n_cols);

    auto arr = py::array_t<T>(
        {nrows, ncols},          // shape
        {tsize, nrows * tsize},  // F-style contiguous strides
        data->memptr(),          // the data pointer
        // numpy array references this parent
        create_dummy_capsule<arma::Mat<T>>(data)
    );
    carma::set_not_writeable(arr);
    return arr;
} /* construct_array */

template <typename T>
inline py::array_t<T> construct_array(const arma::Cube<T>* data) {
    constexpr auto tsize = static_cast<ssize_t>(sizeof(T));
    auto nrows = static_cast<ssize_t>(data->n_rows);
    auto ncols = static_cast<ssize_t>(data->n_cols);
    auto nslices = static_cast<ssize_t>(data->n_slices);
    auto arr = py::array_t<T>(
        // shape
        {nrows, ncols, nslices},
        // F-style contiguous strides
        {tsize, nrows * tsize, tsize * nrows * ncols},
        // the data pointer
        data->memptr(),
        // numpy array references this parent
        create_dummy_capsule<arma::Cube<T>>(data)
    );
    carma::set_not_writeable(arr);
    return arr;
} /* construct_array */

}  // namespace details
}  // namespace carma
#endif  // INCLUDE_CARMA_BITS_ARMATONUMPY_H_
