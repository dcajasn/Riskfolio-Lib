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



#if defined(ARMA_USE_BLAS)

#if defined(dgemm) || defined(DGEMM)
  #pragma message ("WARNING: detected possible interference with definitions of BLAS functions;")
  #pragma message ("WARNING: include the armadillo header before any other header as a workaround")
#endif


#if defined(ARMA_BLAS_NOEXCEPT)
  #undef  ARMA_NOEXCEPT
  #define ARMA_NOEXCEPT noexcept
#else
  #undef  ARMA_NOEXCEPT
  #define ARMA_NOEXCEPT
#endif


#if !defined(ARMA_BLAS_CAPITALS)
  
  #define arma_sasum sasum
  #define arma_dasum dasum
  
  #define arma_snrm2 snrm2
  #define arma_dnrm2 dnrm2
  
  #define arma_sdot  sdot
  #define arma_ddot  ddot
  
  #define arma_sgemv sgemv
  #define arma_dgemv dgemv
  #define arma_cgemv cgemv
  #define arma_zgemv zgemv
  
  #define arma_sgemm sgemm
  #define arma_dgemm dgemm
  #define arma_cgemm cgemm
  #define arma_zgemm zgemm
  
  #define arma_ssyrk ssyrk
  #define arma_dsyrk dsyrk
  
  #define arma_cherk cherk
  #define arma_zherk zherk
  
#else
  
  #define arma_sasum SASUM
  #define arma_dasum DASUM
  
  #define arma_snrm2 SNRM2
  #define arma_dnrm2 DNRM2
  
  #define arma_sdot  SDOT
  #define arma_ddot  DDOT
  
  #define arma_sgemv SGEMV
  #define arma_dgemv DGEMV
  #define arma_cgemv CGEMV
  #define arma_zgemv ZGEMV
  
  #define arma_sgemm SGEMM
  #define arma_dgemm DGEMM
  #define arma_cgemm CGEMM
  #define arma_zgemm ZGEMM
  
  #define arma_ssyrk SSYRK
  #define arma_dsyrk DSYRK
  
  #define arma_cherk CHERK
  #define arma_zherk ZHERK
  
#endif


// NOTE: "For arguments of CHARACTER type, the character length is passed as a hidden argument at the end of the argument list."
// NOTE: https://gcc.gnu.org/onlinedocs/gfortran/Argument-passing-conventions.html


extern "C"
{
#if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
  
  float  arma_fortran(arma_sasum)(const blas_int* n, const float*  x, const blas_int* incx) ARMA_NOEXCEPT;
  double arma_fortran(arma_dasum)(const blas_int* n, const double* x, const blas_int* incx) ARMA_NOEXCEPT;
  
  float  arma_fortran(arma_snrm2)(const blas_int* n, const float*  x, const blas_int* incx) ARMA_NOEXCEPT;
  double arma_fortran(arma_dnrm2)(const blas_int* n, const double* x, const blas_int* incx) ARMA_NOEXCEPT;
  
  float  arma_fortran(arma_sdot)(const blas_int* n, const float*  x, const blas_int* incx, const float*  y, const blas_int* incy) ARMA_NOEXCEPT;
  double arma_fortran(arma_ddot)(const blas_int* n, const double* x, const blas_int* incx, const double* y, const blas_int* incy) ARMA_NOEXCEPT;
  
  void arma_fortran(arma_sgemv)(const char* transA, const blas_int* m, const blas_int* n, const float*    alpha, const float*    A, const blas_int* ldA, const float*    x, const blas_int* incx, const float*    beta, float*    y, const blas_int* incy, blas_len transA_len) ARMA_NOEXCEPT;
  void arma_fortran(arma_dgemv)(const char* transA, const blas_int* m, const blas_int* n, const double*   alpha, const double*   A, const blas_int* ldA, const double*   x, const blas_int* incx, const double*   beta, double*   y, const blas_int* incy, blas_len transA_len) ARMA_NOEXCEPT;
  void arma_fortran(arma_cgemv)(const char* transA, const blas_int* m, const blas_int* n, const blas_cxf* alpha, const blas_cxf* A, const blas_int* ldA, const blas_cxf* x, const blas_int* incx, const blas_cxf* beta, blas_cxf* y, const blas_int* incy, blas_len transA_len) ARMA_NOEXCEPT;
  void arma_fortran(arma_zgemv)(const char* transA, const blas_int* m, const blas_int* n, const blas_cxd* alpha, const blas_cxd* A, const blas_int* ldA, const blas_cxd* x, const blas_int* incx, const blas_cxd* beta, blas_cxd* y, const blas_int* incy, blas_len transA_len) ARMA_NOEXCEPT;
  
  void arma_fortran(arma_sgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const float*    alpha, const float*    A, const blas_int* ldA, const float*    B, const blas_int* ldB, const float*    beta, float*    C, const blas_int* ldC, blas_len transA_len, blas_len transB_len) ARMA_NOEXCEPT;
  void arma_fortran(arma_dgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const double*   alpha, const double*   A, const blas_int* ldA, const double*   B, const blas_int* ldB, const double*   beta, double*   C, const blas_int* ldC, blas_len transA_len, blas_len transB_len) ARMA_NOEXCEPT;
  void arma_fortran(arma_cgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const blas_cxf* alpha, const blas_cxf* A, const blas_int* ldA, const blas_cxf* B, const blas_int* ldB, const blas_cxf* beta, blas_cxf* C, const blas_int* ldC, blas_len transA_len, blas_len transB_len) ARMA_NOEXCEPT;
  void arma_fortran(arma_zgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const blas_cxd* alpha, const blas_cxd* A, const blas_int* ldA, const blas_cxd* B, const blas_int* ldB, const blas_cxd* beta, blas_cxd* C, const blas_int* ldC, blas_len transA_len, blas_len transB_len) ARMA_NOEXCEPT;
  
  void arma_fortran(arma_ssyrk)(const char* uplo, const char* transA, const blas_int* n, const blas_int* k, const  float* alpha, const  float* A, const blas_int* ldA, const  float* beta,  float* C, const blas_int* ldC, blas_len uplo_len, blas_len transA_len) ARMA_NOEXCEPT;
  void arma_fortran(arma_dsyrk)(const char* uplo, const char* transA, const blas_int* n, const blas_int* k, const double* alpha, const double* A, const blas_int* ldA, const double* beta, double* C, const blas_int* ldC, blas_len uplo_len, blas_len transA_len) ARMA_NOEXCEPT;
  
  void arma_fortran(arma_cherk)(const char* uplo, const char* transA, const blas_int* n, const blas_int* k, const  float* alpha, const blas_cxf* A, const blas_int* ldA, const  float* beta, blas_cxf* C, const blas_int* ldC, blas_len uplo_len, blas_len transA_len) ARMA_NOEXCEPT;
  void arma_fortran(arma_zherk)(const char* uplo, const char* transA, const blas_int* n, const blas_int* k, const double* alpha, const blas_cxd* A, const blas_int* ldA, const double* beta, blas_cxd* C, const blas_int* ldC, blas_len uplo_len, blas_len transA_len) ARMA_NOEXCEPT;
  
#else
  
  // prototypes without hidden arguments
  
  float  arma_fortran(arma_sasum)(const blas_int* n, const float*  x, const blas_int* incx) ARMA_NOEXCEPT;
  double arma_fortran(arma_dasum)(const blas_int* n, const double* x, const blas_int* incx) ARMA_NOEXCEPT;
  
  float  arma_fortran(arma_snrm2)(const blas_int* n, const float*  x, const blas_int* incx) ARMA_NOEXCEPT;
  double arma_fortran(arma_dnrm2)(const blas_int* n, const double* x, const blas_int* incx) ARMA_NOEXCEPT;
  
  float  arma_fortran(arma_sdot)(const blas_int* n, const float*  x, const blas_int* incx, const float*  y, const blas_int* incy) ARMA_NOEXCEPT;
  double arma_fortran(arma_ddot)(const blas_int* n, const double* x, const blas_int* incx, const double* y, const blas_int* incy) ARMA_NOEXCEPT;
  
  void arma_fortran(arma_sgemv)(const char* transA, const blas_int* m, const blas_int* n, const float*    alpha, const float*    A, const blas_int* ldA, const float*    x, const blas_int* incx, const float*    beta, float*    y, const blas_int* incy) ARMA_NOEXCEPT;
  void arma_fortran(arma_dgemv)(const char* transA, const blas_int* m, const blas_int* n, const double*   alpha, const double*   A, const blas_int* ldA, const double*   x, const blas_int* incx, const double*   beta, double*   y, const blas_int* incy) ARMA_NOEXCEPT;
  void arma_fortran(arma_cgemv)(const char* transA, const blas_int* m, const blas_int* n, const blas_cxf* alpha, const blas_cxf* A, const blas_int* ldA, const blas_cxf* x, const blas_int* incx, const blas_cxf* beta, blas_cxf* y, const blas_int* incy) ARMA_NOEXCEPT;
  void arma_fortran(arma_zgemv)(const char* transA, const blas_int* m, const blas_int* n, const blas_cxd* alpha, const blas_cxd* A, const blas_int* ldA, const blas_cxd* x, const blas_int* incx, const blas_cxd* beta, blas_cxd* y, const blas_int* incy) ARMA_NOEXCEPT;
  
  void arma_fortran(arma_sgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const float*    alpha, const float*    A, const blas_int* ldA, const float*    B, const blas_int* ldB, const float*    beta, float*    C, const blas_int* ldC) ARMA_NOEXCEPT;
  void arma_fortran(arma_dgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const double*   alpha, const double*   A, const blas_int* ldA, const double*   B, const blas_int* ldB, const double*   beta, double*   C, const blas_int* ldC) ARMA_NOEXCEPT;
  void arma_fortran(arma_cgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const blas_cxf* alpha, const blas_cxf* A, const blas_int* ldA, const blas_cxf* B, const blas_int* ldB, const blas_cxf* beta, blas_cxf* C, const blas_int* ldC) ARMA_NOEXCEPT;
  void arma_fortran(arma_zgemm)(const char* transA, const char* transB, const blas_int* m, const blas_int* n, const blas_int* k, const blas_cxd* alpha, const blas_cxd* A, const blas_int* ldA, const blas_cxd* B, const blas_int* ldB, const blas_cxd* beta, blas_cxd* C, const blas_int* ldC) ARMA_NOEXCEPT;
  
  void arma_fortran(arma_ssyrk)(const char* uplo, const char* transA, const blas_int* n, const blas_int* k, const  float* alpha, const  float* A, const blas_int* ldA, const  float* beta,  float* C, const blas_int* ldC) ARMA_NOEXCEPT;
  void arma_fortran(arma_dsyrk)(const char* uplo, const char* transA, const blas_int* n, const blas_int* k, const double* alpha, const double* A, const blas_int* ldA, const double* beta, double* C, const blas_int* ldC) ARMA_NOEXCEPT;
  
  void arma_fortran(arma_cherk)(const char* uplo, const char* transA, const blas_int* n, const blas_int* k, const  float* alpha, const  blas_cxf* A, const blas_int* ldA, const  float* beta, blas_cxf* C, const blas_int* ldC) ARMA_NOEXCEPT;
  void arma_fortran(arma_zherk)(const char* uplo, const char* transA, const blas_int* n, const blas_int* k, const double* alpha, const  blas_cxd* A, const blas_int* ldA, const double* beta, blas_cxd* C, const blas_int* ldC) ARMA_NOEXCEPT;
  
#endif
}

#undef ARMA_NOEXCEPT

#endif
