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


#if defined(ARMA_USE_ATLAS)


// TODO: remove support for ATLAS in next major version

//! \namespace atlas namespace for ATLAS functions
namespace atlas
  {
  
  template<typename eT>
  inline static const eT& tmp_real(const eT& X)              { return X; }
  
  template<typename T>
  inline static const  T  tmp_real(const std::complex<T>& X) { return X.real(); }
  
  
  
  template<typename eT>
  arma_inline
  eT
  cblas_asum(const int N, const eT* X)
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_float<eT>::value)
      {
      typedef float T;
      return eT( arma_wrapper(cblas_sasum)(N, (const T*)X, 1) );
      }
    else
    if(is_double<eT>::value)
      {
      typedef double T;
      return eT( arma_wrapper(cblas_dasum)(N, (const T*)X, 1) );
      }
    
    return eT(0);
    }
  
  
  
  template<typename eT>
  arma_inline
  eT
  cblas_nrm2(const int N, const eT* X)
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_float<eT>::value)
      {
      typedef float T;
      return eT( arma_wrapper(cblas_snrm2)(N, (const T*)X, 1) );
      }
    else
    if(is_double<eT>::value)
      {
      typedef double T;
      return eT( arma_wrapper(cblas_dnrm2)(N, (const T*)X, 1) );
      }
    
    return eT(0);
    }
  
  
  
  template<typename eT>
  arma_inline
  eT
  cblas_dot(const int N, const eT* X, const eT* Y)
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_float<eT>::value)
      {
      typedef float T;
      return eT( arma_wrapper(cblas_sdot)(N, (const T*)X, 1, (const T*)Y, 1) );
      }
    else
    if(is_double<eT>::value)
      {
      typedef double T;
      return eT( arma_wrapper(cblas_ddot)(N, (const T*)X, 1, (const T*)Y, 1) );
      }
    
    return eT(0);
    }
  
  
  
  template<typename eT>
  arma_inline
  eT
  cblas_cx_dot(const int N, const eT* X, const eT* Y)
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_cx_float<eT>::value)
      {
      typedef typename std::complex<float> T;
      
      T out;    
      arma_wrapper(cblas_cdotu_sub)(N, (const T*)X, 1, (const T*)Y, 1, &out);
      
      return eT(out);
      }
    else
    if(is_cx_double<eT>::value)
      {
      typedef typename std::complex<double> T;
      
      T out;
      arma_wrapper(cblas_zdotu_sub)(N, (const T*)X, 1, (const T*)Y, 1, &out);
      
      return eT(out);
      }
    
    return eT(0);
    }
  
  
  
  template<typename eT>
  inline
  void
  cblas_gemv
    (
    const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_TRANS TransA,
    const int M, const int N,
    const eT alpha,
    const eT *A, const int lda,
    const eT *X, const int incX,
    const eT beta,
    eT *Y, const int incY
    )
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_float<eT>::value)
      {
      typedef float T;
      arma_wrapper(cblas_sgemv)(layout, TransA, M, N, (const T)tmp_real(alpha), (const T*)A, lda, (const T*)X, incX, (const T)tmp_real(beta), (T*)Y, incY);
      }
    else
    if(is_double<eT>::value)
      {
      typedef double T;
      arma_wrapper(cblas_dgemv)(layout, TransA, M, N, (const T)tmp_real(alpha), (const T*)A, lda, (const T*)X, incX, (const T)tmp_real(beta), (T*)Y, incY);
      }
    else
    if(is_cx_float<eT>::value)
      {
      typedef std::complex<float> T;
      arma_wrapper(cblas_cgemv)(layout, TransA, M, N, (const T*)&alpha, (const T*)A, lda, (const T*)X, incX, (const T*)&beta, (T*)Y, incY);
      }
    else
    if(is_cx_double<eT>::value)
      {
      typedef std::complex<double> T;
      arma_wrapper(cblas_zgemv)(layout, TransA, M, N, (const T*)&alpha, (const T*)A, lda, (const T*)X, incX, (const T*)&beta, (T*)Y, incY);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  cblas_gemm
    (
    const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_TRANS TransA,
    const atlas_CBLAS_TRANS TransB, const int M, const int N,
    const int K, const eT alpha, const eT *A,
    const int lda, const eT *B, const int ldb,
    const eT beta, eT *C, const int ldc
    )
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_float<eT>::value)
      {
      typedef float T;
      arma_wrapper(cblas_sgemm)(layout, TransA, TransB, M, N, K, (const T)tmp_real(alpha), (const T*)A, lda, (const T*)B, ldb, (const T)tmp_real(beta), (T*)C, ldc);
      }
    else
    if(is_double<eT>::value)
      {
      typedef double T;
      arma_wrapper(cblas_dgemm)(layout, TransA, TransB, M, N, K, (const T)tmp_real(alpha), (const T*)A, lda, (const T*)B, ldb, (const T)tmp_real(beta), (T*)C, ldc);
      }
    else
    if(is_cx_float<eT>::value)
      {
      typedef std::complex<float> T;
      arma_wrapper(cblas_cgemm)(layout, TransA, TransB, M, N, K, (const T*)&alpha, (const T*)A, lda, (const T*)B, ldb, (const T*)&beta, (T*)C, ldc);
      }
    else
    if(is_cx_double<eT>::value)
      {
      typedef std::complex<double> T;
      arma_wrapper(cblas_zgemm)(layout, TransA, TransB, M, N, K, (const T*)&alpha, (const T*)A, lda, (const T*)B, ldb, (const T*)&beta, (T*)C, ldc);
      }
    }
  
  
  
  template<typename eT>
  inline
  void
  cblas_syrk
    (
    const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_UPLO Uplo, const atlas_CBLAS_TRANS Trans,
    const int N, const int K, const eT alpha,
    const eT* A, const int lda, const eT beta, eT* C, const int ldc
    )
    {
    arma_type_check((is_supported_blas_type<eT>::value == false));
    
    if(is_float<eT>::value)
      {
      typedef float T;
      arma_wrapper(cblas_ssyrk)(layout, Uplo, Trans, N, K, (const T)alpha, (const T*)A, lda, (const T)beta, (T*)C, ldc);
      }
    else
    if(is_double<eT>::value)
      {
      typedef double T;
      arma_wrapper(cblas_dsyrk)(layout, Uplo, Trans, N, K, (const T)alpha, (const T*)A, lda, (const T)beta, (T*)C, ldc);
      }
    }
  
  
  
  template<typename T>
  inline
  void
  cblas_herk
    (
    const atlas_CBLAS_LAYOUT layout, const atlas_CBLAS_UPLO Uplo, const atlas_CBLAS_TRANS Trans,
    const int N, const int K, const T alpha,
    const std::complex<T>* A, const int lda, const T beta, std::complex<T>* C, const int ldc
    )
    {
    arma_type_check((is_supported_blas_type<T>::value == false));
    
    if(is_float<T>::value)
      {
      typedef float                  TT;
      typedef std::complex<float> cx_TT;
      
      arma_wrapper(cblas_cherk)(layout, Uplo, Trans, N, K, (const TT)alpha, (const cx_TT*)A, lda, (const TT)beta, (cx_TT*)C, ldc);
      }
    else
    if(is_double<T>::value)
      {
      typedef double                  TT;
      typedef std::complex<double> cx_TT;
      
      arma_wrapper(cblas_zherk)(layout, Uplo, Trans, N, K, (const TT)alpha, (const cx_TT*)A, lda, (const TT)beta, (cx_TT*)C, ldc);
      }
    }
  
  }

#endif
