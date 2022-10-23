// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------



#if !defined(ARMA_WARN_LEVEL)
  #define ARMA_WARN_LEVEL 2
#endif
//// The level of warning messages printed to ARMA_CERR_STREAM.
//// Must be an integer >= 0. The default value is 2.
//// 0 = no warnings; generally not recommended
//// 1 = only critical warnings about arguments and/or data which are likely to lead to incorrect results
//// 2 = as per level 1, and warnings about poorly conditioned systems (low rcond) detected by solve(), spsolve(), etc
//// 3 = as per level 2, and warnings about failed decompositions, failed saving/loading, etc

#if !defined(ARMA_USE_LAPACK)
#cmakedefine ARMA_USE_LAPACK
//// Comment out the above line if you don't have LAPACK or a high-speed replacement for LAPACK,
//// such as OpenBLAS, Intel MKL, or the Accelerate framework.
//// LAPACK is required for matrix decompositions (eg. SVD) and matrix inverse.
#endif

#if !defined(ARMA_USE_BLAS)
#cmakedefine ARMA_USE_BLAS
//// Comment out the above line if you don't have BLAS or a high-speed replacement for BLAS,
//// such as OpenBLAS, Intel MKL, or the Accelerate framework.
//// BLAS is used for matrix multiplication.
//// Without BLAS, matrix multiplication will still work, but might be slower.
#endif

#if !defined(ARMA_USE_NEWARP)
#define ARMA_USE_NEWARP
//// Uncomment the above line to enable the built-in partial emulation of ARPACK.
//// This is used for eigen decompositions of real (non-complex) sparse matrices, eg. eigs_sym(), svds() 
#endif

#if !defined(ARMA_USE_ARPACK)
#cmakedefine ARMA_USE_ARPACK
//// Uncomment the above line if you have ARPACK or a high-speed replacement for ARPACK.
//// ARPACK is required for eigen decompositions of complex sparse matrices
#endif

#if !defined(ARMA_USE_SUPERLU)
#cmakedefine ARMA_USE_SUPERLU
//// Uncomment the above line if you have SuperLU.
//// SuperLU is used for solving sparse linear systems via spsolve()
//// Caveat: only SuperLU version 5.2 can be used!
#endif

#if !defined(ARMA_SUPERLU_INCLUDE_DIR)
#define ARMA_SUPERLU_INCLUDE_DIR ${ARMA_SUPERLU_INCLUDE_DIR}/
//// If you're using SuperLU and want to explicitly include the SuperLU headers,
//// uncomment the above define and specify the appropriate include directory.
//// Make sure the directory has a trailing /
#endif

#if !defined(ARMA_USE_ATLAS)
#cmakedefine ARMA_USE_ATLAS
//// NOTE: support for ATLAS is deprecated and will be removed.
#endif

#cmakedefine ARMA_USE_WRAPPER
//// Comment out the above line if you're getting linking errors when compiling your programs,
//// or if you prefer to directly link with LAPACK, BLAS + etc instead of the Armadillo runtime library.
//// You will then need to link your programs directly with -llapack -lblas instead of -larmadillo

// #define ARMA_BLAS_CAPITALS
//// Uncomment the above line if your BLAS and LAPACK libraries have capitalised function names

#define ARMA_BLAS_UNDERSCORE
//// Uncomment the above line if your BLAS and LAPACK libraries have function names with a trailing underscore.
//// Conversely, comment it out if the function names don't have a trailing underscore.

// #define ARMA_BLAS_LONG
//// Uncomment the above line if your BLAS and LAPACK libraries use "long" instead of "int"

// #define ARMA_BLAS_LONG_LONG
//// Uncomment the above line if your BLAS and LAPACK libraries use "long long" instead of "int"

// #define ARMA_BLAS_NOEXCEPT
//// Uncomment the above line if you require BLAS functions to have the 'noexcept' specification

// #define ARMA_LAPACK_NOEXCEPT
//// Uncomment the above line if you require LAPACK functions to have the 'noexcept' specification

#define ARMA_USE_FORTRAN_HIDDEN_ARGS
//// Comment out the above line to call BLAS and LAPACK functions without using so-called "hidden" arguments.
//// Fortran functions (compiled without a BIND(C) declaration) that have char arguments
//// (like many BLAS and LAPACK functions) also have associated "hidden" arguments.
//// For each char argument, the corresponding "hidden" argument specifies the number of characters.
//// These "hidden" arguments are typically tacked onto the end of function definitions.

// #define ARMA_USE_TBB_ALLOC
//// Uncomment the above line if you want to use Intel TBB scalable_malloc() and scalable_free() instead of standard malloc() and free()

// #define ARMA_USE_MKL_ALLOC
//// Uncomment the above line if you want to use Intel MKL mkl_malloc() and mkl_free() instead of standard malloc() and free()

// #define ARMA_USE_MKL_TYPES
//// Uncomment the above line if you want to use Intel MKL types for complex numbers.
//// You will need to include appropriate MKL headers before the Armadillo header.
//// You may also need to enable or disable the following options:
//// ARMA_BLAS_LONG, ARMA_BLAS_LONG_LONG, ARMA_USE_FORTRAN_HIDDEN_ARGS

#if !defined(ARMA_USE_OPENMP)
// #define ARMA_USE_OPENMP
//// Uncomment the above line to forcefully enable use of OpenMP for parallelisation.
//// Note that ARMA_USE_OPENMP is automatically enabled when a compiler supporting OpenMP 3.1 is detected.
#endif

#if !defined(ARMA_64BIT_WORD)
// #define ARMA_64BIT_WORD
//// Uncomment the above line if you require matrices/vectors capable of holding more than 4 billion elements.
//// Note that ARMA_64BIT_WORD is automatically enabled when std::size_t has 64 bits and ARMA_32BIT_WORD is not defined.
#endif

#if !defined(ARMA_USE_HDF5)
// #define ARMA_USE_HDF5
//// Uncomment the above line to allow the ability to save and load matrices stored in HDF5 format;
//// the hdf5.h header file must be available on your system,
//// and you will need to link with the hdf5 library (eg. -lhdf5)
#endif

#if !defined(ARMA_OPTIMISE_BAND)
  #define ARMA_OPTIMISE_BAND
  //// Comment out the above line if you don't want automatically optimised handling
  //// of band matrices by solve() and chol()
#endif

#if !defined(ARMA_OPTIMISE_SYMPD)
  #define ARMA_OPTIMISE_SYMPD
  //// Comment out the above line if you don't want automatically optimised handling
  //// of symmetric/hermitian positive definite matrices by various functions:
  //// solve(), inv(), pinv(), expmat(), logmat(), sqrtmat(), rcond(), rank()
#endif

#if !defined(ARMA_OPTIMISE_INVEXPR)
  #define ARMA_OPTIMISE_INVEXPR
  //// Comment out the above line if you don't want automatically optimised handling
  //// of inv() and inv_sympd() within compound expressions
#endif

#if !defined(ARMA_CHECK_NONFINITE)
  #define ARMA_CHECK_NONFINITE
  //// Comment out the above line if you don't want automatic checking for nonfinite matrices
#endif

#cmakedefine ARMA_USE_HDF5_CMAKE
#if defined(ARMA_USE_HDF5_CMAKE) && defined(ARMA_USE_WRAPPER)
  #undef  ARMA_USE_HDF5
  #define ARMA_USE_HDF5
  
  #define ARMA_HDF5_INCLUDE_DIR ${ARMA_HDF5_INCLUDE_DIR}/
#endif

#if !defined(ARMA_MAT_PREALLOC)
  #define ARMA_MAT_PREALLOC 16
#endif
//// This is the number of preallocated elements used by matrices and vectors;
//// it must be an integer that is at least 1.
//// If you mainly use lots of very small vectors (eg. <= 4 elements),
//// change the number to the size of your vectors.

#if !defined(ARMA_OPENMP_THRESHOLD)
  #define ARMA_OPENMP_THRESHOLD 320
#endif
//// The minimum number of elements in a matrix to allow OpenMP based parallelisation;
//// it must be an integer that is at least 1.

#if !defined(ARMA_OPENMP_THREADS)
  #define ARMA_OPENMP_THREADS 8
#endif
//// The maximum number of threads to use for OpenMP based parallelisation;
//// it must be an integer that is at least 1.

// #define ARMA_NO_DEBUG
//// Uncomment the above line if you want to disable all run-time checks.
//// This will result in faster code, but you first need to make sure that your code runs correctly!
//// We strongly recommend to have the run-time checks enabled during development,
//// as this greatly aids in finding mistakes in your code, and hence speeds up development.
//// We recommend that run-time checks be disabled _only_ for the shipped version of your program.

// #define ARMA_EXTRA_DEBUG
//// Uncomment the above line if you want to see the function traces of how Armadillo evaluates expressions.
//// This is mainly useful for debugging of the library.


#if defined(ARMA_DEFAULT_OSTREAM)
  #pragma message ("WARNING: support for ARMA_DEFAULT_OSTREAM is deprecated and will be removed;")
  #pragma message ("WARNING: use ARMA_COUT_STREAM and ARMA_CERR_STREAM instead")
#endif


#if !defined(ARMA_COUT_STREAM)
  #if defined(ARMA_DEFAULT_OSTREAM)
    // for compatibility with earlier versions of Armadillo
    #define ARMA_COUT_STREAM ARMA_DEFAULT_OSTREAM
  #else
    #define ARMA_COUT_STREAM std::cout
  #endif
#endif

#if !defined(ARMA_CERR_STREAM)
  #if defined(ARMA_DEFAULT_OSTREAM)
    // for compatibility with earlier versions of Armadillo
    #define ARMA_CERR_STREAM ARMA_DEFAULT_OSTREAM
  #else
    #define ARMA_CERR_STREAM std::cerr
  #endif
#endif


#if !defined(ARMA_PRINT_EXCEPTIONS)
  // #define ARMA_PRINT_EXCEPTIONS
  #if defined(ARMA_PRINT_EXCEPTIONS_INTERNAL)
    #undef  ARMA_PRINT_EXCEPTIONS
    #define ARMA_PRINT_EXCEPTIONS
  #endif
#endif

#if !defined(ARMA_PRINT_HDF5_ERRORS)
// #define ARMA_PRINT_HDF5_ERRORS
#endif

#if defined(ARMA_DONT_USE_LAPACK)
  #undef ARMA_USE_LAPACK
#endif

#if defined(ARMA_DONT_USE_BLAS)
  #undef ARMA_USE_BLAS
#endif

#if defined(ARMA_DONT_USE_NEWARP) || !defined(ARMA_USE_LAPACK)
  #undef ARMA_USE_NEWARP
#endif

#if defined(ARMA_DONT_USE_ARPACK)
  #undef ARMA_USE_ARPACK
#endif

#if defined(ARMA_DONT_USE_SUPERLU)
  #undef ARMA_USE_SUPERLU
  #undef ARMA_SUPERLU_INCLUDE_DIR
#endif

#if defined(ARMA_DONT_USE_ATLAS)
  #undef ARMA_USE_ATLAS
#endif

#if defined(ARMA_DONT_USE_WRAPPER)
  #undef ARMA_USE_WRAPPER
  #undef ARMA_USE_HDF5_CMAKE
#endif

#if defined(ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS)
  #undef ARMA_USE_FORTRAN_HIDDEN_ARGS
#endif

#if !defined(ARMA_DONT_USE_STD_MUTEX)
  // #define ARMA_DONT_USE_STD_MUTEX
  //// Uncomment the above line to disable use of std::mutex
#endif

// for compatibility with earlier versions of Armadillo
#if defined(ARMA_DONT_USE_CXX11_MUTEX)
  #pragma message ("WARNING: support for ARMA_DONT_USE_CXX11_MUTEX is deprecated and will be removed;")
  #pragma message ("WARNING: use ARMA_DONT_USE_STD_MUTEX instead")
  #undef  ARMA_DONT_USE_STD_MUTEX
  #define ARMA_DONT_USE_STD_MUTEX
#endif

#if defined(ARMA_DONT_USE_OPENMP)
  #undef ARMA_USE_OPENMP
#endif

#if defined(ARMA_USE_WRAPPER)
  #if !defined(ARMA_USE_EXTERN_RNG)
    #cmakedefine ARMA_USE_EXTERN_RNG
  #endif
#endif

#if defined(ARMA_DONT_USE_EXTERN_RNG)
  #undef ARMA_USE_EXTERN_RNG
#endif

// for compatibility with earlier versions of Armadillo
#if defined(ARMA_DONT_USE_EXTERN_CXX11_RNG)
  #pragma message ("WARNING: support for ARMA_DONT_USE_EXTERN_CXX11_RNG is deprecated and will be removed;")
  #pragma message ("WARNING: use ARMA_DONT_USE_EXTERN_RNG instead")
  #undef ARMA_USE_EXTERN_RNG
#endif

#if defined(ARMA_32BIT_WORD)
  #undef ARMA_64BIT_WORD
#endif

#if defined(ARMA_DONT_USE_HDF5)
  #undef ARMA_USE_HDF5
  #undef ARMA_USE_HDF5_CMAKE
#endif

#if defined(ARMA_DONT_OPTIMISE_BAND) || defined(ARMA_DONT_OPTIMISE_SOLVE_BAND)
  #undef ARMA_OPTIMISE_BAND
#endif

#if defined(ARMA_DONT_OPTIMISE_SYMPD) || defined(ARMA_DONT_OPTIMISE_SOLVE_SYMPD)
  #undef ARMA_OPTIMISE_SYMPD
#endif

#if defined(ARMA_DONT_OPTIMISE_INVEXPR)
  #undef ARMA_OPTIMISE_INVEXPR
#endif

#if defined(ARMA_DONT_CHECK_NONFINITE)
  #undef ARMA_CHECK_NONFINITE
#endif

#if defined(ARMA_DONT_PRINT_ERRORS)
  #pragma message ("INFO: support for ARMA_DONT_PRINT_ERRORS option has been removed")
  
  #if defined(ARMA_PRINT_EXCEPTIONS)
    #pragma message ("INFO: suggest to use ARMA_WARN_LEVEL and ARMA_DONT_PRINT_EXCEPTIONS options instead")
  #else
    #pragma message ("INFO: suggest to use ARMA_WARN_LEVEL option instead")
  #endif
  
  #pragma message ("INFO: see the documentation for details")
#endif

#if defined(ARMA_DONT_PRINT_EXCEPTIONS)
  #undef ARMA_PRINT_EXCEPTIONS
#endif

#if !defined(ARMA_DONT_ZERO_INIT)
  // #define ARMA_DONT_ZERO_INIT
  //// Uncomment the above line to disable initialising elements to zero during construction of dense matrices and cubes
#endif

#if defined(ARMA_DONT_PRINT_HDF5_ERRORS)
  #undef ARMA_PRINT_HDF5_ERRORS
#endif

#if defined(ARMA_NO_CRIPPLED_LAPACK)
  #undef ARMA_CRIPPLED_LAPACK
#endif


// if Armadillo was installed on this system via CMake and ARMA_USE_WRAPPER is not defined,
// ARMA_AUX_LIBS lists the libraries required by Armadillo on this system, and
// ARMA_AUX_INCDIRS lists the include directories required by Armadillo on this system.
// Do not use these unless you know what you are doing.
#define ARMA_AUX_LIBS ${ARMA_LIBS}
#define ARMA_AUX_INCDIRS ${CMAKE_REQUIRED_INCLUDES}
