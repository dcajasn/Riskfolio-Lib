/*  carma/cnalloc.h: Wrapper around PyDataMem_* for use with Armadillo
 *  Copyright (c) 2021 Ralph Urlus <rurlus.dev@gmail.com>
 *  All rights reserved. Use of this source code is governed by a
 *  Apache-2.0 license that can be found in the LICENSE file.
 */
#ifndef INCLUDE_CARMA_BITS_CNALLOC_H_
#define INCLUDE_CARMA_BITS_CNALLOC_H_

// pybind11 include required even if not explicitly used 
// to prevent link with pythonXX_d.lib on Win32 
// (cf Py_DEBUG defined in numpy headers and https://github.com/pybind/pybind11/issues/1295)
#include <pybind11/pybind11.h>     
#define NPY_NO_DEPRECATED_API NPY_1_14_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>

#include <cstddef>
#ifdef CARMA_DEV_DEBUG
#include <iostream>
#endif

namespace cnalloc {

inline void* npy_malloc(size_t bytes) {
    if (PyArray_API == nullptr) {
        _import_array();
    }
#ifdef CARMA_DEV_DEBUG
    std::cout << "\n-----------\nCARMA DEBUG\n-----------\n";
    std::cout << "Using numpy allocator" << "\n";
    std::cout << "-----------\n";
#endif  // ARMA_EXTRA_DEBUG
    return PyDataMem_NEW(bytes);
} // npy_malloc

inline void npy_free(void* ptr) {
    if (PyArray_API == nullptr) {
        _import_array();
    }
#ifdef CARMA_DEV_DEBUG
    std::cout << "\n-----------\nCARMA DEBUG\n-----------\n";
    std::cout << "Using numpy deallocator\n";
    std::cout << "-----------\n";
#endif  // ARMA_EXTRA_DEBUG
    PyDataMem_FREE(ptr);
} // npy_free

} // namespace cnalloc

#define ARMA_ALIEN_MEM_ALLOC_FUNCTION cnalloc::npy_malloc
#define ARMA_ALIEN_MEM_FREE_FUNCTION cnalloc::npy_free
#ifndef CARMA_ARMA_ALIEN_MEM_FUNCTIONS_SET
  #define CARMA_ARMA_ALIEN_MEM_FUNCTIONS_SET
#endif
#endif  // INCLUDE_CARMA_BITS_CNALLOC_H_
