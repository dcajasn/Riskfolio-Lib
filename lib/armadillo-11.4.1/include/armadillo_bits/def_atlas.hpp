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


// TODO: remove support for ATLAS in next major version

#if defined(ARMA_USE_ATLAS)


typedef enum
  {
  atlas_CblasRowMajor = 101,
  atlas_CblasColMajor = 102
  }
  atlas_CBLAS_LAYOUT;

typedef enum
  {
  atlas_CblasNoTrans   = 111,
  atlas_CblasTrans     = 112,
  atlas_CblasConjTrans = 113
  }
  atlas_CBLAS_TRANS;
  
typedef enum
  {
  atlas_CblasUpper = 121,
  atlas_CblasLower = 122
  }
  atlas_CBLAS_UPLO;


extern "C"
  {
  float  arma_wrapper(cblas_sasum)(const int N, const float  *X, const int incX);
  double arma_wrapper(cblas_dasum)(const int N, const double *X, const int incX);
  
  float  arma_wrapper(cblas_snrm2)(const int N, const float  *X, const int incX);
  double arma_wrapper(cblas_dnrm2)(const int N, const double *X, const int incX);
  
  float  arma_wrapper(cblas_sdot)(const int N, const float  *X, const int incX, const float  *Y, const int incY);
  double arma_wrapper(cblas_ddot)(const int N, const double *X, const int incX, const double *Y, const int incY);
  
  void arma_wrapper(cblas_cdotu_sub)(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotu);
  void arma_wrapper(cblas_zdotu_sub)(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotu);
  
  void arma_wrapper(cblas_sgemv)(const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_TRANS TransA, const int M, const int N, const float   alpha, const float  *A, const int lda, const float  *X, const int incX, const float   beta, float  *Y, const int incY);
  void arma_wrapper(cblas_dgemv)(const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_TRANS TransA, const int M, const int N, const double  alpha, const double *A, const int lda, const double *X, const int incX, const double  beta, double *Y, const int incY);
  void arma_wrapper(cblas_cgemv)(const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_TRANS TransA, const int M, const int N, const void   *alpha, const void   *A, const int lda, const void   *X, const int incX, const void   *beta, void   *Y, const int incY);
  void arma_wrapper(cblas_zgemv)(const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_TRANS TransA, const int M, const int N, const void   *alpha, const void   *A, const int lda, const void   *X, const int incX, const void   *beta, void   *Y, const int incY);
  
  void arma_wrapper(cblas_sgemm)(const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_TRANS TransA, const atlas_CBLAS_TRANS TransB, const int M, const int N, const int K, const float   alpha, const float  *A, const int lda, const float  *B, const int ldb, const float   beta, float  *C, const int ldc);
  void arma_wrapper(cblas_dgemm)(const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_TRANS TransA, const atlas_CBLAS_TRANS TransB, const int M, const int N, const int K, const double  alpha, const double *A, const int lda, const double *B, const int ldb, const double  beta, double *C, const int ldc);
  void arma_wrapper(cblas_cgemm)(const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_TRANS TransA, const atlas_CBLAS_TRANS TransB, const int M, const int N, const int K, const void   *alpha, const void   *A, const int lda, const void   *B, const int ldb, const void   *beta, void   *C, const int ldc);
  void arma_wrapper(cblas_zgemm)(const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_TRANS TransA, const atlas_CBLAS_TRANS TransB, const int M, const int N, const int K, const void   *alpha, const void   *A, const int lda, const void   *B, const int ldb, const void   *beta, void   *C, const int ldc);
  
  void arma_wrapper(cblas_ssyrk)(const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_UPLO Uplo, const atlas_CBLAS_TRANS Trans, const int N, const int K, const float  alpha, const float  *A, const int lda, const float  beta, float  *C, const int ldc);
  void arma_wrapper(cblas_dsyrk)(const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_UPLO Uplo, const atlas_CBLAS_TRANS Trans, const int N, const int K, const double alpha, const double *A, const int lda, const double beta, double *C, const int ldc);
  
  void arma_wrapper(cblas_cherk)(const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_UPLO Uplo, const atlas_CBLAS_TRANS Trans, const int N, const int K, const float  alpha, const void *A, const int lda, const float  beta, void *C, const int ldc);
  void arma_wrapper(cblas_zherk)(const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_UPLO Uplo, const atlas_CBLAS_TRANS Trans, const int N, const int K, const double alpha, const void *A, const int lda, const double beta, void *C, const int ldc);
  }


#endif
