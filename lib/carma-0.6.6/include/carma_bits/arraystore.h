/*  carma/arraystore.h: Store arrays as armadillo matrices
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
#ifndef INCLUDE_CARMA_BITS_ARRAYSTORE_H_
#define INCLUDE_CARMA_BITS_ARRAYSTORE_H_

#include <pybind11/numpy.h>  // NOLINT
#include <pybind11/pybind11.h>  // NOLINT
#include <carma_bits/debug.h>  // NOLINT
#include <carma_bits/config.h> // NOLINT
#include <carma_bits/typecheck.h>  // NOLINT
#include <carma_bits/converters.h>  // NOLINT

#include <utility>

namespace py = pybind11;

namespace carma {

template <typename armaT>
class ArrayStore {
    using T = typename armaT::elem_type;

 protected:
    constexpr static ssize_t tsize = sizeof(T);
    py::capsule base;

 public:
    armaT mat;

 public:
    ArrayStore(py::array_t<T>& arr, bool copy) :
    mat{to_arma<armaT>::from(arr, copy)} {
        base = details::create_dummy_capsule(&mat);
    }

    explicit ArrayStore(const armaT& src) : mat{armaT(src)} {
        base = details::create_dummy_capsule(&mat);
    }

    ArrayStore(arma::Mat<T>& src, bool copy) {
        if (copy) {
            mat = armaT(src.memptr(), src.n_rows, src.n_cols, true);
        } else {
            mat = std::move(src);
        }
        base = details::create_dummy_capsule(&mat);
    }

    ArrayStore(arma::Cube<T>& src, bool copy) {
        if (copy) {
            mat = armaT(src.memptr(), src.n_rows, src.n_cols, src.n_slices, true);
        } else {
            mat = std::move(src);
        }
        base = details::create_dummy_capsule(&mat);
    }

    // SFINAE by adding additional parameter as
    // to avoid shadowing the class template
    template <typename U = armaT>
    ArrayStore(armaT& src, bool copy, is_Vec<U>) {
        if (copy) {
            mat = armaT(src.memptr(), src.n_elem, true);
        } else {
            mat = std::move(src);
        }
        base = details::create_dummy_capsule(&mat);
    }

    explicit ArrayStore(armaT&& src) noexcept : mat{std::move(src)} {
        base = details::create_dummy_capsule(&mat);
    }

    // Function requires different name than set_data
    // as overload could not be resolved without
    void set_array(py::array_t<T>& arr, bool copy) {
        if (copy) {
            mat = to_arma<armaT>::from(arr, true);
        } else {
            mat = to_arma<armaT>::from(std::move(arr));
        }
        base = details::create_dummy_capsule(&mat);
    }

    void set_data(const armaT& src) {
        mat = armaT(src);
        base = details::create_dummy_capsule(&mat);
    }

    void set_data(arma::Mat<T>& src, bool copy) {
        if (copy) {
            mat = armaT(src.memptr(), src.n_rows, src.n_cols, true);
        } else {
            mat = std::move(src);
        }
        base = details::create_dummy_capsule(&mat);
    }

    // SFINAE by adding additional parameter as
    // to avoid shadowing the class template
    template <typename U = armaT>
    void set_data(armaT& src, bool copy, is_Vec<U>) {
        if (copy) {
            mat = armaT(src.memptr(), src.n_elem, true);
        } else {
            mat = std::move(src);
        }
        base = details::create_dummy_capsule(&mat);
    }

    void set_data(arma::Cube<T>& src, bool copy) {
        if (copy) {
            mat = armaT(src.memptr(), src.n_rows, src.n_cols, src.n_slices, true);
        } else {
            mat = std::move(src);
        }
        base = details::create_dummy_capsule(&mat);
    }

    void set_data(armaT&& src) {
        mat = std::move(src);
        base = details::create_dummy_capsule(&mat);
    }

    py::array_t<T> get_view(bool writeable) {
        ssize_t nslices;
        auto nelem = static_cast<ssize_t>(mat.n_elem);
        auto nrows = static_cast<ssize_t>(mat.n_rows);
        auto ncols = static_cast<ssize_t>(mat.n_cols);
        ssize_t rc_elem = nrows * ncols;

        py::array_t<T> arr;

        // detect cubes
        if (rc_elem != nelem) {
            nslices = nelem / rc_elem;
            arr = py::array_t<T>(
                {nrows, ncols, nslices},
                // F-style contiguous strides
                {tsize, nrows * tsize, tsize * nrows * ncols},
                mat.memptr(),                                   // the data pointer
                base                                           // numpy array references this parent
            );
        } else {
            arr = py::array_t<T>(
                {nrows, ncols},          // shape
                {tsize, nrows * tsize},  // F-style contiguous strides
                mat.memptr(),            // the data pointer
                base                    // numpy array references this parent
            );
        }

        // inform numpy it does not own the buffer
        set_not_owndata(arr);

        if (!writeable)
            set_not_writeable(arr);
        return arr;
    }
};

} /* namespace carma */

#endif  // INCLUDE_CARMA_BITS_ARRAYSTORE_H_
