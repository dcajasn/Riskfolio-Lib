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


//! \addtogroup typedef_elem
//! @{


#if (defined(ARMA_U8_TYPE) && defined(ARMA_S8_TYPE))
    typedef ARMA_U8_TYPE     u8;
    typedef ARMA_S8_TYPE     s8;
#else
  #if   UCHAR_MAX >= 0xff
    typedef unsigned char    u8;
    typedef          char    s8;
  #elif defined(UINT8_MAX)
    typedef          uint8_t u8;
    typedef           int8_t s8;
  #else
    #error "don't know how to typedef 'u8' on this system"
  #endif
#endif

// NOTE:
// "char" is not guaranteed to be the same as "signed char" 
// https://en.wikipedia.org/wiki/C_data_types


#if   USHRT_MAX >= 0xffff
  typedef unsigned short    u16;
  typedef          short    s16;
#elif defined(UINT16_MAX)
  typedef          uint16_t u16;
  typedef           int16_t s16;
#else
  #error "don't know how to typedef 'u16' on this system"
#endif


#if   UINT_MAX  >= 0xffffffff
  typedef unsigned int      u32;
  typedef          int      s32;
#elif defined(UINT32_MAX)
  typedef          uint32_t u32;
  typedef           int32_t s32;
#else
  #error "don't know how to typedef 'u32' on this system"
#endif


#if   ULLONG_MAX >= 0xffffffffffffffff
  typedef unsigned long long u64;
  typedef          long long s64;
#elif defined(UINT64_MAX)
  typedef          uint64_t  u64;
  typedef           int64_t  s64;
#else
    #error "don't know how to typedef 'u64' on this system"
#endif


// for compatibility with earlier versions of Armadillo
typedef unsigned long ulng_t;
typedef          long slng_t;


#if defined(ARMA_64BIT_WORD)
  typedef u64 uword;
  typedef s64 sword;
  
  typedef u32 uhword;
  typedef s32 shword;

  #define ARMA_MAX_UWORD  0xffffffffffffffff
  #define ARMA_MAX_UHWORD 0xffffffff
#else
  typedef u32 uword;
  typedef s32 sword;

  typedef u16 uhword;
  typedef s16 shword;
  
  #define ARMA_MAX_UWORD  0xffffffff
  #define ARMA_MAX_UHWORD 0xffff
#endif


typedef std::complex<float>  cx_float;
typedef std::complex<double> cx_double;

typedef void* void_ptr;


//


#if   defined(ARMA_BLAS_LONG_LONG)
  typedef long long blas_int;
  #define ARMA_MAX_BLAS_INT 0x7fffffffffffffffULL
#elif defined(ARMA_BLAS_LONG)
  typedef long      blas_int;
  #define ARMA_MAX_BLAS_INT 0x7fffffffffffffffUL
#else
  typedef int       blas_int;
  #define ARMA_MAX_BLAS_INT 0x7fffffffU
#endif


//


#if defined(ARMA_USE_MKL_TYPES)
  // for compatibility with MKL
  typedef MKL_Complex8  blas_cxf;
  typedef MKL_Complex16 blas_cxd;
#else
  // standard BLAS and LAPACK prototypes use "void*" pointers for complex arrays
  typedef void blas_cxf;
  typedef void blas_cxd;
#endif


//


// NOTE: blas_len is the fortran type for "hidden" arguments that specify the length of character arguments;
// NOTE: it varies across compilers, compiler versions and systems (eg. 32 bit vs 64 bit);
// NOTE: the default setting of "size_t" is an educated guess.
// NOTE: ---
// NOTE: for gcc / gfortran:  https://gcc.gnu.org/onlinedocs/gfortran/Argument-passing-conventions.html
// NOTE: gcc 7 and earlier: int
// NOTE: gcc 8 and 9:       size_t
// NOTE: ---
// NOTE: for ifort (intel fortran compiler): 
// NOTE: "Intel Fortran Compiler User and Reference Guides", Document Number: 304970-006US, 2009, p. 301
// NOTE: http://www.complexfluids.ethz.ch/MK/ifort.pdf
// NOTE: the type is unsigned 4-byte integer on 32 bit systems
// NOTE: the type is unsigned 8-byte integer on 64 bit systems
// NOTE: ---
// NOTE: for NAG fortran: https://www.nag.co.uk/nagware/np/r62_doc/manual/compiler_11_1.html#AUTOTOC_11_1
// NOTE: Chrlen = usually int, or long long on 64-bit Windows
// NOTE: ---
// TODO: flang:  https://github.com/flang-compiler/flang/wiki
// TODO: other compilers: http://fortranwiki.org/fortran/show/Compilers

#if !defined(ARMA_FORTRAN_CHARLEN_TYPE)
  #if defined(__GNUC__) && !defined(__clang__)
    #if (__GNUC__ <= 7)
      #define ARMA_FORTRAN_CHARLEN_TYPE int
    #else
      #define ARMA_FORTRAN_CHARLEN_TYPE size_t
    #endif
  #else
    // TODO: determine the type for other compilers
    #define ARMA_FORTRAN_CHARLEN_TYPE size_t
  #endif
#endif

typedef ARMA_FORTRAN_CHARLEN_TYPE blas_len;


//! @}
