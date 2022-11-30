/*  carma/debug.h: Debug printers
 *  Copyright (c) 2020 Ralph Urlus <rurlus.dev@gmail.com>
 *  All rights reserved. Use of this source code is governed by a
 *  Apache-2.0 license that can be found in the LICENSE file. */
/* External headers */
#ifndef INCLUDE_CARMA_BITS_DEBUG_H_
#define INCLUDE_CARMA_BITS_DEBUG_H_

#include <Python.h>
#include <pymem.h>
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>

#include <pybind11/numpy.h>  // NOLINT
#include <pybind11/pybind11.h>  // NOLINT

#include <carma_bits/config.h> // NOLINT
#include <armadillo>  // NOLINT

#include <string>
#include <memory>
#include <type_traits>
#include <utility>

namespace py = pybind11;

namespace carma {

#ifdef CARMA_EXTRA_DEBUG
namespace debug {
using aconf =  arma::arma_config;

template <typename T>
inline void print_array_info(PyObject* src) {
    PyArrayObject* arr = reinterpret_cast<PyArrayObject*>(src);
    T* data = reinterpret_cast<T*>(PyArray_DATA(arr));
    int ndim = PyArray_NDIM(arr);
    npy_intp * dims = PyArray_DIMS(arr);
    bool first = true;
    std::cout << "\nThe array has shape: ";
    for (int i = 0; i < ndim; i++) {
        std::cout << (first ? "(" : ", ") << dims[i];
        first = false;
    }
    std::cout << ")" << "\n";
    std::cout << "with first element: " << data[0] << "\n";
}  // print_array_info

inline void print_opening() {
    static const std::string m = "\n-----------\nCARMA DEBUG\n-----------\n";
    std::cout << m;
}

inline void print_closing() {
    static const std::string m = "-----------\n";
    std::cout << m;
}

template <typename T>
inline void print_copy_of_data(T* data) {
    print_opening();
    std::cout << "Memory @" << data <<  " will be copied." << "\n";
    std::cout << "It is not well behaved." << "\n";
    print_closing();
}

template <typename T>
inline void print_copy_of_data(T* data, PyObject* obj) {
    print_opening();
    std::cout << "Memory @" << data <<  " will be copied." << "\n";
    std::cout << "It is not well behaved." << "\n";
    print_array_info<T>(obj);
    print_closing();
}

template <typename T>
inline void print_cannot_steal(T* data) {
    print_opening();
    std::cout << "Memory @" << data <<  " cannot be stolen and will be copied." << "\n";
    std::cout << "It is not well behaved." << "\n";
    print_closing();
}

template <typename T>
inline void print_cannot_steal(T* data, PyObject* obj) {
    print_opening();
    std::cout << "Memory @" << data <<  " cannot be stolen and will be copied." << "\n";
    std::cout << "It is not well behaved." << "\n";
    print_array_info<T>(obj);
    print_closing();
}

template <typename T>
inline void print_prealloc(T* data) {
    static constexpr int pl = aconf::mat_prealloc;
    print_opening();
    std::cout << "Memory at @" << data << " will be copied." << "\n";
    std::cout << "It is smaller than armadillo's preallocation limit: " << pl << "\n";
    print_closing();
}

template <typename T>
inline void print_prealloc(T* data, PyObject* obj) {
    static constexpr int pl = aconf::mat_prealloc;
    print_opening();
    std::cout << "Memory at @" << data << " will be copied." << "\n";
    std::cout << "It is smaller than armadillo's preallocation limit: " << pl << "\n";
    print_array_info<T>(obj);
    print_closing();
}

}  // namespace debug
#endif  // CARMA_EXTRA_DEBUG

}  // namespace carma
#endif  // INCLUDE_CARMA_BITS_DEBUG_H_
