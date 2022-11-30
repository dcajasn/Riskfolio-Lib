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



#if defined(ARMA_USE_LAPACK)


//! \namespace lapack namespace for LAPACK functions
namespace lapack
  {
  
  template<typename eT>
  inline
  void
  getrf(blas_int* m, blas_int* n, eT* a, blas_int* lda, blas_int* ipiv, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgetrf)(m, n, (T*)a, lda, ipiv, info); }
    else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgetrf)(m, n, (T*)a, lda, ipiv, info); }
    else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgetrf)(m, n, (T*)a, lda, ipiv, info); }
    else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgetrf)(m, n, (T*)a, lda, ipiv, info); }
    }
    
    
    
  template<typename eT>
  inline
  void
  getrs(char* trans, blas_int* n, blas_int* nrhs, eT* a, blas_int* lda, blas_int* ipiv, eT* b, blas_int* ldb, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgetrs)(trans, n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info, 1); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgetrs)(trans, n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info, 1); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgetrs)(trans, n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info, 1); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgetrs)(trans, n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info, 1); }
    #else
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgetrs)(trans, n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgetrs)(trans, n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgetrs)(trans, n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgetrs)(trans, n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info); }
    #endif
    }
    
    
    
  template<typename eT>
  inline
  void
  getri(blas_int* n,  eT* a, blas_int* lda, blas_int* ipiv, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgetri)(n, (T*)a, lda, ipiv, (T*)work, lwork, info); }
    else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgetri)(n, (T*)a, lda, ipiv, (T*)work, lwork, info); }
    else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgetri)(n, (T*)a, lda, ipiv, (T*)work, lwork, info); }
    else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgetri)(n, (T*)a, lda, ipiv, (T*)work, lwork, info); }
    }
  
  
  
  template<typename eT>
  inline
  void
  trtri(char* uplo, char* diag, blas_int* n, eT* a, blas_int* lda, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_strtri)(uplo, diag, n, (T*)a, lda, info, 1, 1); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dtrtri)(uplo, diag, n, (T*)a, lda, info, 1, 1); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_ctrtri)(uplo, diag, n, (T*)a, lda, info, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_ztrtri)(uplo, diag, n, (T*)a, lda, info, 1, 1); }
    #else
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_strtri)(uplo, diag, n, (T*)a, lda, info); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dtrtri)(uplo, diag, n, (T*)a, lda, info); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_ctrtri)(uplo, diag, n, (T*)a, lda, info); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_ztrtri)(uplo, diag, n, (T*)a, lda, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  geev(char* jobvl, char* jobvr, blas_int* n, eT* a, blas_int* lda, eT* wr, eT* wi, eT* vl, blas_int* ldvl, eT* vr, blas_int* ldvr, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgeev)(jobvl, jobvr, n, (T*)a, lda, (T*)wr, (T*)wi, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info, 1, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgeev)(jobvl, jobvr, n, (T*)a, lda, (T*)wr, (T*)wi, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info, 1, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgeev)(jobvl, jobvr, n, (T*)a, lda, (T*)wr, (T*)wi, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgeev)(jobvl, jobvr, n, (T*)a, lda, (T*)wr, (T*)wi, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  cx_geev(char* jobvl, char* jobvr, blas_int* n, eT* a, blas_int* lda, eT* w, eT* vl, blas_int* ldvl, eT* vr, blas_int* ldvr, eT* work, blas_int* lwork, typename eT::value_type* rwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  T; typedef blas_cxf cx_T; arma_fortran(arma_cgeev)(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)w, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef double T; typedef blas_cxd cx_T; arma_fortran(arma_zgeev)(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)w, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info, 1, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  T; typedef blas_cxf cx_T; arma_fortran(arma_cgeev)(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)w, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info); }
      else if(is_cx_double<eT>::value)  { typedef double T; typedef blas_cxd cx_T; arma_fortran(arma_zgeev)(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)w, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  geevx(char* balanc, char* jobvl, char* jobvr, char* sense, blas_int* n, eT* a, blas_int* lda, eT* wr, eT* wi, eT* vl, blas_int* ldvl, eT* vr, blas_int* ldvr, blas_int* ilo, blas_int* ihi, eT* scale, eT* abnrm, eT* rconde, eT* rcondv, eT* work, blas_int* lwork, blas_int* iwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));

    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgeevx)(balanc, jobvl, jobvr, sense, n, (T*)(a), lda, (T*)(wr), (T*)(wi), (T*)(vl), ldvl, (T*)(vr), ldvr, ilo, ihi, (T*)(scale), (T*)(abnrm), (T*)(rconde), (T*)(rcondv), (T*)(work), lwork, iwork, info, 1, 1, 1, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgeevx)(balanc, jobvl, jobvr, sense, n, (T*)(a), lda, (T*)(wr), (T*)(wi), (T*)(vl), ldvl, (T*)(vr), ldvr, ilo, ihi, (T*)(scale), (T*)(abnrm), (T*)(rconde), (T*)(rcondv), (T*)(work), lwork, iwork, info, 1, 1, 1, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgeevx)(balanc, jobvl, jobvr, sense, n, (T*)(a), lda, (T*)(wr), (T*)(wi), (T*)(vl), ldvl, (T*)(vr), ldvr, ilo, ihi, (T*)(scale), (T*)(abnrm), (T*)(rconde), (T*)(rcondv), (T*)(work), lwork, iwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgeevx)(balanc, jobvl, jobvr, sense, n, (T*)(a), lda, (T*)(wr), (T*)(wi), (T*)(vl), ldvl, (T*)(vr), ldvr, ilo, ihi, (T*)(scale), (T*)(abnrm), (T*)(rconde), (T*)(rcondv), (T*)(work), lwork, iwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  cx_geevx(char* balanc, char* jobvl, char* jobvr, char* sense, blas_int* n, eT* a, blas_int* lda, eT* w, eT* vl, blas_int* ldvl, eT* vr, blas_int* ldvr, blas_int* ilo, blas_int* ihi, typename eT::value_type* scale, typename eT::value_type* abnrm, typename eT::value_type* rconde, typename eT::value_type* rcondv, eT* work, blas_int* lwork, typename eT::value_type* rwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  T; typedef blas_cxf cx_T; arma_fortran(arma_cgeevx)(balanc, jobvl, jobvr, sense, n, (cx_T*)(a), lda, (cx_T*)(w), (cx_T*)(vl), ldvl, (cx_T*)(vr), ldvr, ilo, ihi, (T*)(scale), (T*)(abnrm), (T*)(rconde), (T*)(rcondv), (cx_T*)(work), lwork, (T*)(rwork), info, 1, 1, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef double T; typedef blas_cxd cx_T; arma_fortran(arma_zgeevx)(balanc, jobvl, jobvr, sense, n, (cx_T*)(a), lda, (cx_T*)(w), (cx_T*)(vl), ldvl, (cx_T*)(vr), ldvr, ilo, ihi, (T*)(scale), (T*)(abnrm), (T*)(rconde), (T*)(rcondv), (cx_T*)(work), lwork, (T*)(rwork), info, 1, 1, 1, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  T; typedef blas_cxf cx_T; arma_fortran(arma_cgeevx)(balanc, jobvl, jobvr, sense, n, (cx_T*)(a), lda, (cx_T*)(w), (cx_T*)(vl), ldvl, (cx_T*)(vr), ldvr, ilo, ihi, (T*)(scale), (T*)(abnrm), (T*)(rconde), (T*)(rcondv), (cx_T*)(work), lwork, (T*)(rwork), info); }
      else if(is_cx_double<eT>::value)  { typedef double T; typedef blas_cxd cx_T; arma_fortran(arma_zgeevx)(balanc, jobvl, jobvr, sense, n, (cx_T*)(a), lda, (cx_T*)(w), (cx_T*)(vl), ldvl, (cx_T*)(vr), ldvr, ilo, ihi, (T*)(scale), (T*)(abnrm), (T*)(rconde), (T*)(rcondv), (cx_T*)(work), lwork, (T*)(rwork), info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  syev(char* jobz, char* uplo, blas_int* n, eT* a, blas_int* lda, eT* w,  eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  {  typedef float  T; arma_fortran(arma_ssyev)(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, info, 1, 1); }
      else if(is_double<eT>::value)  {  typedef double T; arma_fortran(arma_dsyev)(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, info, 1, 1); }
    #else
           if( is_float<eT>::value)  {  typedef float  T; arma_fortran(arma_ssyev)(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, info); }
      else if(is_double<eT>::value)  {  typedef double T; arma_fortran(arma_dsyev)(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  heev
    (
    char* jobz, char* uplo, blas_int* n,
    eT* a, blas_int* lda, typename eT::value_type* w,
    eT* work, blas_int* lwork, typename eT::value_type* rwork,
    blas_int* info
    )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  T; typedef blas_cxf cx_T; arma_fortran(arma_cheev)(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, info, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef double T; typedef blas_cxd cx_T; arma_fortran(arma_zheev)(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, info, 1, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  T; typedef blas_cxf cx_T; arma_fortran(arma_cheev)(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, info); }
      else if(is_cx_double<eT>::value)  { typedef double T; typedef blas_cxd cx_T; arma_fortran(arma_zheev)(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  syevd(char* jobz, char* uplo, blas_int* n, eT* a, blas_int* lda, eT* w,  eT* work, blas_int* lwork, blas_int* iwork, blas_int* liwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_ssyevd)(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, iwork, liwork, info, 1, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dsyevd)(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, iwork, liwork, info, 1, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_ssyevd)(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, iwork, liwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dsyevd)(jobz, uplo, n, (T*)a, lda, (T*)w, (T*)work, lwork, iwork, liwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  heevd
    (
    char* jobz, char* uplo, blas_int* n,
    eT* a, blas_int* lda, typename eT::value_type* w,
    eT* work, blas_int* lwork, typename eT::value_type* rwork, 
    blas_int* lrwork, blas_int* iwork, blas_int* liwork,
    blas_int* info
    )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  T; typedef blas_cxf cx_T; arma_fortran(arma_cheevd)(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, lrwork, iwork, liwork, info, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef double T; typedef blas_cxd cx_T; arma_fortran(arma_zheevd)(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, lrwork, iwork, liwork, info, 1, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  T; typedef blas_cxf cx_T; arma_fortran(arma_cheevd)(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, lrwork, iwork, liwork, info); }
      else if(is_cx_double<eT>::value)  { typedef double T; typedef blas_cxd cx_T; arma_fortran(arma_zheevd)(jobz, uplo, n, (cx_T*)a, lda, (T*)w, (cx_T*)work, lwork, (T*)rwork, lrwork, iwork, liwork, info); }
    #endif
    }
  
	
	
  template<typename eT>
  inline
  void
  ggev
    (
    char* jobvl, char* jobvr, blas_int* n,
    eT* a, blas_int* lda, eT* b, blas_int* ldb,
    eT* alphar, eT* alphai, eT* beta,
    eT* vl, blas_int* ldvl, eT* vr, blas_int* ldvr,
    eT* work, blas_int* lwork,
    blas_int* info
    )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sggev)(jobvl, jobvr, n, (T*)a, lda, (T*)b, ldb, (T*)alphar, (T*)alphai, (T*)beta, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info, 1, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dggev)(jobvl, jobvr, n, (T*)a, lda, (T*)b, ldb, (T*)alphar, (T*)alphai, (T*)beta, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info, 1, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sggev)(jobvl, jobvr, n, (T*)a, lda, (T*)b, ldb, (T*)alphar, (T*)alphai, (T*)beta, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dggev)(jobvl, jobvr, n, (T*)a, lda, (T*)b, ldb, (T*)alphar, (T*)alphai, (T*)beta, (T*)vl, ldvl, (T*)vr, ldvr, (T*)work, lwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  cx_ggev
    (
    char* jobvl, char* jobvr, blas_int* n,
    eT* a, blas_int* lda, eT* b, blas_int* ldb,
    eT* alpha, eT* beta,
    eT* vl, blas_int* ldvl, eT* vr, blas_int* ldvr,
    eT* work, blas_int* lwork, typename eT::value_type* rwork,
    blas_int* info
    )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  T; typedef blas_cxf cx_T; arma_fortran(arma_cggev)(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)b, ldb, (cx_T*)alpha, (cx_T*)beta, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef double T; typedef blas_cxd cx_T; arma_fortran(arma_zggev)(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)b, ldb, (cx_T*)alpha, (cx_T*)beta, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info, 1, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  T; typedef blas_cxf cx_T; arma_fortran(arma_cggev)(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)b, ldb, (cx_T*)alpha, (cx_T*)beta, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info); }
      else if(is_cx_double<eT>::value)  { typedef double T; typedef blas_cxd cx_T; arma_fortran(arma_zggev)(jobvl, jobvr, n, (cx_T*)a, lda, (cx_T*)b, ldb, (cx_T*)alpha, (cx_T*)beta, (cx_T*)vl, ldvl, (cx_T*)vr, ldvr, (cx_T*)work, lwork, (T*)rwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  potrf(char* uplo, blas_int* n, eT* a, blas_int* lda, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_spotrf)(uplo, n, (T*)a, lda, info, 1); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dpotrf)(uplo, n, (T*)a, lda, info, 1); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cpotrf)(uplo, n, (T*)a, lda, info, 1); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zpotrf)(uplo, n, (T*)a, lda, info, 1); }
    #else
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_spotrf)(uplo, n, (T*)a, lda, info); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dpotrf)(uplo, n, (T*)a, lda, info); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cpotrf)(uplo, n, (T*)a, lda, info); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zpotrf)(uplo, n, (T*)a, lda, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  potrs(char* uplo, blas_int* n, const blas_int* nrhs, eT* a, blas_int* lda, eT* b, blas_int* ldb, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_spotrs)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info, 1); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dpotrs)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info, 1); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cpotrs)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info, 1); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zpotrs)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info, 1); }
    #else
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_spotrs)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dpotrs)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cpotrs)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zpotrs)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  pbtrf(char* uplo, blas_int* n, blas_int* kd, eT* ab, blas_int* ldab, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_spbtrf)(uplo, n, kd, (T*)ab, ldab, info, 1); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dpbtrf)(uplo, n, kd, (T*)ab, ldab, info, 1); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cpbtrf)(uplo, n, kd, (T*)ab, ldab, info, 1); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zpbtrf)(uplo, n, kd, (T*)ab, ldab, info, 1); }
    #else
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_spbtrf)(uplo, n, kd, (T*)ab, ldab, info); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dpbtrf)(uplo, n, kd, (T*)ab, ldab, info); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cpbtrf)(uplo, n, kd, (T*)ab, ldab, info); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zpbtrf)(uplo, n, kd, (T*)ab, ldab, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  potri(char* uplo, blas_int* n, eT* a, blas_int* lda, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_spotri)(uplo, n, (T*)a, lda, info, 1); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dpotri)(uplo, n, (T*)a, lda, info, 1); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cpotri)(uplo, n, (T*)a, lda, info, 1); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zpotri)(uplo, n, (T*)a, lda, info, 1); }
    #else
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_spotri)(uplo, n, (T*)a, lda, info); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dpotri)(uplo, n, (T*)a, lda, info); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cpotri)(uplo, n, (T*)a, lda, info); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zpotri)(uplo, n, (T*)a, lda, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  geqrf(blas_int* m, blas_int* n, eT* a, blas_int* lda, eT* tau, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgeqrf)(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info); }
    else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgeqrf)(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info); }
    else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgeqrf)(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info); }
    else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgeqrf)(m, n, (T*)a, lda, (T*)tau, (T*)work, lwork, info); }
    }
  
  
  
  template<typename eT>
  inline
  void
  geqp3(blas_int* m, blas_int* n, eT* a, blas_int* lda, blas_int* jpvt, eT* tau, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgeqp3)(m, n, (T*)a, lda, jpvt, (T*)tau, (T*)work, lwork, info); }
    else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgeqp3)(m, n, (T*)a, lda, jpvt, (T*)tau, (T*)work, lwork, info); }
    }
  
  
  
  template<typename eT>
  inline
  void
  cx_geqp3(blas_int* m, blas_int* n, eT* a, blas_int* lda, blas_int* jpvt, eT* tau, eT* work, blas_int* lwork, typename eT::value_type* rwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if( is_cx_float<eT>::value)  { typedef  float T; typedef blas_cxf cx_T; arma_fortran(arma_cgeqp3)(m, n, (cx_T*)a, lda, jpvt, (cx_T*)tau, (cx_T*)work, lwork, (T*)rwork, info); }
    else if(is_cx_double<eT>::value)  { typedef double T; typedef blas_cxd cx_T; arma_fortran(arma_zgeqp3)(m, n, (cx_T*)a, lda, jpvt, (cx_T*)tau, (cx_T*)work, lwork, (T*)rwork, info); }
    }
  
  
  
  template<typename eT>
  inline
  void
  orgqr(blas_int* m, blas_int* n, blas_int* k, eT* a, blas_int* lda, eT* tau, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sorgqr)(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info); }
    else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dorgqr)(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info); }
    }
  
  
  
  template<typename eT>
  inline
  void
  ungqr(blas_int* m, blas_int* n, blas_int* k, eT* a, blas_int* lda, eT* tau, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cungqr)(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info); }
    else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zungqr)(m, n, k, (T*)a, lda, (T*)tau, (T*)work, lwork, info); }
    }
  
  
  
  template<typename eT>
  inline
  void
  gesvd
    (
    char* jobu, char* jobvt, blas_int* m, blas_int* n, eT* a, blas_int* lda,
    eT* s, eT* u, blas_int* ldu, eT* vt, blas_int* ldvt,
    eT* work, blas_int* lwork, blas_int* info
    )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgesvd)(jobu, jobvt, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, info, 1, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgesvd)(jobu, jobvt, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, info, 1, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgesvd)(jobu, jobvt, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgesvd)(jobu, jobvt, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, info); }
    #endif
    }
  
  
  
  template<typename T>
  inline
  void
  cx_gesvd
    (
    char* jobu, char* jobvt, blas_int* m, blas_int* n, std::complex<T>* a, blas_int* lda,
    T* s, std::complex<T>* u, blas_int* ldu, std::complex<T>* vt, blas_int* ldvt, 
    std::complex<T>* work, blas_int* lwork, T* rwork, blas_int* info
    )
    {
    arma_type_check(( is_supported_blas_type<T>::value == false ));
    arma_type_check(( is_supported_blas_type< std::complex<T> >::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<T>::value)  { typedef float  bT; typedef blas_cxf cx_bT; arma_fortran(arma_cgesvd)(jobu, jobvt, m, n, (cx_bT*)a, lda, (bT*)s, (cx_bT*)u, ldu, (cx_bT*)vt, ldvt, (cx_bT*)work, lwork, (bT*)rwork, info, 1, 1); }
      else if(is_double<T>::value)  { typedef double bT; typedef blas_cxd cx_bT; arma_fortran(arma_zgesvd)(jobu, jobvt, m, n, (cx_bT*)a, lda, (bT*)s, (cx_bT*)u, ldu, (cx_bT*)vt, ldvt, (cx_bT*)work, lwork, (bT*)rwork, info, 1, 1); }
    #else
           if( is_float<T>::value)  { typedef float  bT; typedef blas_cxf cx_bT; arma_fortran(arma_cgesvd)(jobu, jobvt, m, n, (cx_bT*)a, lda, (bT*)s, (cx_bT*)u, ldu, (cx_bT*)vt, ldvt, (cx_bT*)work, lwork, (bT*)rwork, info); }
      else if(is_double<T>::value)  { typedef double bT; typedef blas_cxd cx_bT; arma_fortran(arma_zgesvd)(jobu, jobvt, m, n, (cx_bT*)a, lda, (bT*)s, (cx_bT*)u, ldu, (cx_bT*)vt, ldvt, (cx_bT*)work, lwork, (bT*)rwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  gesdd
    (
    char* jobz, blas_int* m, blas_int* n,
    eT* a, blas_int* lda, eT* s, eT* u, blas_int* ldu, eT* vt, blas_int* ldvt,
    eT* work, blas_int* lwork, blas_int* iwork, blas_int* info
    )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgesdd)(jobz, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, iwork, info, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgesdd)(jobz, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, iwork, info, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgesdd)(jobz, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, iwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgesdd)(jobz, m, n, (T*)a, lda, (T*)s, (T*)u, ldu, (T*)vt, ldvt, (T*)work, lwork, iwork, info); }
    #endif
    }
  
  
  
  template<typename T>
  inline
  void
  cx_gesdd
    (
    char* jobz, blas_int* m, blas_int* n,
    std::complex<T>* a, blas_int* lda, T* s, std::complex<T>* u, blas_int* ldu, std::complex<T>* vt, blas_int* ldvt,
    std::complex<T>* work, blas_int* lwork, T* rwork, blas_int* iwork, blas_int* info
    )
    {
    arma_type_check(( is_supported_blas_type<T>::value == false ));
    arma_type_check(( is_supported_blas_type< std::complex<T> >::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<T>::value)  { typedef float  bT; typedef blas_cxf cx_bT; arma_fortran(arma_cgesdd)(jobz, m, n, (cx_bT*)a, lda, (bT*)s, (cx_bT*)u, ldu, (cx_bT*)vt, ldvt, (cx_bT*)work, lwork, (bT*)rwork, iwork, info, 1); }
      else if(is_double<T>::value)  { typedef double bT; typedef blas_cxd cx_bT; arma_fortran(arma_zgesdd)(jobz, m, n, (cx_bT*)a, lda, (bT*)s, (cx_bT*)u, ldu, (cx_bT*)vt, ldvt, (cx_bT*)work, lwork, (bT*)rwork, iwork, info, 1); }
    #else
           if( is_float<T>::value)  { typedef float  bT; typedef blas_cxf cx_bT; arma_fortran(arma_cgesdd)(jobz, m, n, (cx_bT*)a, lda, (bT*)s, (cx_bT*)u, ldu, (cx_bT*)vt, ldvt, (cx_bT*)work, lwork, (bT*)rwork, iwork, info); }
      else if(is_double<T>::value)  { typedef double bT; typedef blas_cxd cx_bT; arma_fortran(arma_zgesdd)(jobz, m, n, (cx_bT*)a, lda, (bT*)s, (cx_bT*)u, ldu, (cx_bT*)vt, ldvt, (cx_bT*)work, lwork, (bT*)rwork, iwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  gesv(blas_int* n, blas_int* nrhs, eT* a, blas_int* lda, blas_int* ipiv, eT* b, blas_int* ldb, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgesv)(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info); }
    else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgesv)(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info); }
    else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgesv)(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info); }
    else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgesv)(n, nrhs, (T*)a, lda, ipiv, (T*)b, ldb, info); }
    }
  
  
  
  template<typename eT>
  inline
  void
  gesvx(char* fact, char* trans, blas_int* n, blas_int* nrhs, eT* a, blas_int* lda, eT* af, blas_int* ldaf, blas_int* ipiv, char* equed, eT* r, eT* c, eT* b, blas_int* ldb, eT* x, blas_int* ldx, eT* rcond, eT* ferr, eT* berr, eT* work, blas_int* iwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgesvx)(fact, trans, n, nrhs, (T*)a, lda, (T*)af, ldaf, ipiv, equed, (T*)r, (T*)c, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info, 1, 1, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgesvx)(fact, trans, n, nrhs, (T*)a, lda, (T*)af, ldaf, ipiv, equed, (T*)r, (T*)c, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info, 1, 1, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgesvx)(fact, trans, n, nrhs, (T*)a, lda, (T*)af, ldaf, ipiv, equed, (T*)r, (T*)c, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgesvx)(fact, trans, n, nrhs, (T*)a, lda, (T*)af, ldaf, ipiv, equed, (T*)r, (T*)c, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info); }
    #endif
    }
  
  
  
  template<typename T, typename eT>
  inline
  void
  cx_gesvx(char* fact, char* trans, blas_int* n, blas_int* nrhs, eT* a, blas_int* lda, eT* af, blas_int* ldaf, blas_int* ipiv, char* equed, T* r, T* c, eT* b, blas_int* ldb, eT* x, blas_int* ldx, T* rcond, T* ferr, T* berr, eT* work, T* rwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cgesvx)(fact, trans, n, nrhs, (cx_T*)a, lda, (cx_T*)af, ldaf, ipiv, equed, (pod_T*)r, (pod_T*)c, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info, 1, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zgesvx)(fact, trans, n, nrhs, (cx_T*)a, lda, (cx_T*)af, ldaf, ipiv, equed, (pod_T*)r, (pod_T*)c, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info, 1, 1, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cgesvx)(fact, trans, n, nrhs, (cx_T*)a, lda, (cx_T*)af, ldaf, ipiv, equed, (pod_T*)r, (pod_T*)c, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zgesvx)(fact, trans, n, nrhs, (cx_T*)a, lda, (cx_T*)af, ldaf, ipiv, equed, (pod_T*)r, (pod_T*)c, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  posv(char* uplo, blas_int* n, blas_int* nrhs, eT* a, blas_int* lda, eT* b, blas_int* ldb, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sposv)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info, 1); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dposv)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info, 1); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cposv)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info, 1); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zposv)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info, 1); }
    #else
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sposv)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dposv)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cposv)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zposv)(uplo, n, nrhs, (T*)a, lda, (T*)b, ldb, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  posvx(char* fact, char* uplo, blas_int* n, blas_int* nrhs, eT* a, blas_int* lda, eT* af, blas_int* ldaf, char* equed, eT* s, eT* b, blas_int* ldb, eT* x, blas_int* ldx, eT* rcond, eT* ferr, eT* berr, eT* work, blas_int* iwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sposvx)(fact, uplo, n, nrhs, (T*)a, lda, (T*)af, ldaf, equed, (T*)s, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info, 1, 1, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dposvx)(fact, uplo, n, nrhs, (T*)a, lda, (T*)af, ldaf, equed, (T*)s, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info, 1, 1, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sposvx)(fact, uplo, n, nrhs, (T*)a, lda, (T*)af, ldaf, equed, (T*)s, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dposvx)(fact, uplo, n, nrhs, (T*)a, lda, (T*)af, ldaf, equed, (T*)s, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info); }
    #endif
    }
  
  
  
  template<typename T, typename eT>
  inline
  void
  cx_posvx(char* fact, char* uplo, blas_int* n, blas_int* nrhs, eT* a, blas_int* lda, eT* af, blas_int* ldaf, char* equed, T* s, eT* b, blas_int* ldb, eT* x, blas_int* ldx, T* rcond, T* ferr, T* berr, eT* work, T* rwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cposvx)(fact, uplo, n, nrhs, (cx_T*)a, lda, (cx_T*)af, ldaf, equed, (pod_T*)s, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info, 1, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zposvx)(fact, uplo, n, nrhs, (cx_T*)a, lda, (cx_T*)af, ldaf, equed, (pod_T*)s, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info, 1, 1, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cposvx)(fact, uplo, n, nrhs, (cx_T*)a, lda, (cx_T*)af, ldaf, equed, (pod_T*)s, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zposvx)(fact, uplo, n, nrhs, (cx_T*)a, lda, (cx_T*)af, ldaf, equed, (pod_T*)s, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  gels(char* trans, blas_int* m, blas_int* n, blas_int* nrhs, eT* a, blas_int* lda, eT* b, blas_int* ldb, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info, 1); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info, 1); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info, 1); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info, 1); }
    #else
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgels)(trans, m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)work, lwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  gelsd(blas_int* m, blas_int* n, blas_int* nrhs, eT* a, blas_int* lda, eT* b, blas_int* ldb, eT* S, eT* rcond, blas_int* rank, eT* work, blas_int* lwork, blas_int* iwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgelsd)(m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)S, (T*)rcond, rank, (T*)work, lwork, iwork, info); }
    else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgelsd)(m, n, nrhs, (T*)a, lda, (T*)b, ldb, (T*)S, (T*)rcond, rank, (T*)work, lwork, iwork, info); }
    }
  
  
  
  template<typename T>
  inline
  void
  cx_gelsd(blas_int* m, blas_int* n, blas_int* nrhs, std::complex<T>* a, blas_int* lda, std::complex<T>* b, blas_int* ldb, T* S, T* rcond, blas_int* rank, std::complex<T>* work, blas_int* lwork, T* rwork, blas_int* iwork, blas_int* info)
    {
    typedef typename std::complex<T> eT;
    
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cgelsd)(m, n, nrhs, (cx_T*)a, lda, (cx_T*)b, ldb, (pod_T*)S, (pod_T*)rcond, rank, (cx_T*)work, lwork, (pod_T*)rwork, iwork, info); }
    else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zgelsd)(m, n, nrhs, (cx_T*)a, lda, (cx_T*)b, ldb, (pod_T*)S, (pod_T*)rcond, rank, (cx_T*)work, lwork, (pod_T*)rwork, iwork, info); }
    }
  
  
  
  template<typename eT>
  inline
  void
  trtrs(char* uplo, char* trans, char* diag, blas_int* n, blas_int* nrhs, const eT* a, blas_int* lda, eT* b, blas_int* ldb, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_strtrs)(uplo, trans, diag, n, nrhs, (T*)a, lda, (T*)b, ldb, info, 1, 1, 1); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dtrtrs)(uplo, trans, diag, n, nrhs, (T*)a, lda, (T*)b, ldb, info, 1, 1, 1); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_ctrtrs)(uplo, trans, diag, n, nrhs, (T*)a, lda, (T*)b, ldb, info, 1, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_ztrtrs)(uplo, trans, diag, n, nrhs, (T*)a, lda, (T*)b, ldb, info, 1, 1, 1); }
    #else
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_strtrs)(uplo, trans, diag, n, nrhs, (T*)a, lda, (T*)b, ldb, info); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dtrtrs)(uplo, trans, diag, n, nrhs, (T*)a, lda, (T*)b, ldb, info); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_ctrtrs)(uplo, trans, diag, n, nrhs, (T*)a, lda, (T*)b, ldb, info); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_ztrtrs)(uplo, trans, diag, n, nrhs, (T*)a, lda, (T*)b, ldb, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  gbtrf(blas_int* m, blas_int* n, blas_int* kl, blas_int* ku, eT* ab, blas_int* ldab, blas_int* ipiv, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgbtrf)(m, n, kl, ku, (T*)ab, ldab, ipiv, info); }
    else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgbtrf)(m, n, kl, ku, (T*)ab, ldab, ipiv, info); }
    else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgbtrf)(m, n, kl, ku, (T*)ab, ldab, ipiv, info); }
    else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgbtrf)(m, n, kl, ku, (T*)ab, ldab, ipiv, info); }
    }
  
  
  
  template<typename eT>
  inline
  void
  gbtrs(char* trans, blas_int* n, blas_int* kl, blas_int* ku, blas_int* nrhs, eT* ab, blas_int* ldab, blas_int* ipiv, eT* b, blas_int* ldb, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgbtrs)(trans, n, kl, ku, nrhs, (T*)ab, ldab, ipiv, (T*)b, ldb, info, 1); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgbtrs)(trans, n, kl, ku, nrhs, (T*)ab, ldab, ipiv, (T*)b, ldb, info, 1); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgbtrs)(trans, n, kl, ku, nrhs, (T*)ab, ldab, ipiv, (T*)b, ldb, info, 1); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgbtrs)(trans, n, kl, ku, nrhs, (T*)ab, ldab, ipiv, (T*)b, ldb, info, 1); }
    #else
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgbtrs)(trans, n, kl, ku, nrhs, (T*)ab, ldab, ipiv, (T*)b, ldb, info); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgbtrs)(trans, n, kl, ku, nrhs, (T*)ab, ldab, ipiv, (T*)b, ldb, info); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgbtrs)(trans, n, kl, ku, nrhs, (T*)ab, ldab, ipiv, (T*)b, ldb, info); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgbtrs)(trans, n, kl, ku, nrhs, (T*)ab, ldab, ipiv, (T*)b, ldb, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  gbsv(blas_int* n, blas_int* kl, blas_int* ku, blas_int* nrhs, eT* ab, blas_int* ldab, blas_int* ipiv, eT* b, blas_int* ldb, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgbsv)(n, kl, ku, nrhs, (T*)ab, ldab, ipiv, (T*)b, ldb, info); }
    else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgbsv)(n, kl, ku, nrhs, (T*)ab, ldab, ipiv, (T*)b, ldb, info); }
    else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgbsv)(n, kl, ku, nrhs, (T*)ab, ldab, ipiv, (T*)b, ldb, info); }
    else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgbsv)(n, kl, ku, nrhs, (T*)ab, ldab, ipiv, (T*)b, ldb, info); }
    }
  
  
  
  template<typename eT>
  inline
  void
  gbsvx(char* fact, char* trans, blas_int* n, blas_int* kl, blas_int* ku, blas_int* nrhs, eT* ab, blas_int* ldab, eT* afb, blas_int* ldafb, blas_int* ipiv, char* equed, eT* r, eT* c, eT* b, blas_int* ldb, eT* x, blas_int* ldx, eT* rcond, eT* ferr, eT* berr, eT* work, blas_int* iwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgbsvx)(fact, trans, n, kl, ku, nrhs, (T*)ab, ldab, (T*)afb, ldafb, ipiv, equed, (T*)r, (T*)c, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info, 1, 1, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgbsvx)(fact, trans, n, kl, ku, nrhs, (T*)ab, ldab, (T*)afb, ldafb, ipiv, equed, (T*)r, (T*)c, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info, 1, 1, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgbsvx)(fact, trans, n, kl, ku, nrhs, (T*)ab, ldab, (T*)afb, ldafb, ipiv, equed, (T*)r, (T*)c, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgbsvx)(fact, trans, n, kl, ku, nrhs, (T*)ab, ldab, (T*)afb, ldafb, ipiv, equed, (T*)r, (T*)c, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info); }
    #endif
    }
  
  
  
  template<typename T, typename eT>
  inline
  void
  cx_gbsvx(char* fact, char* trans, blas_int* n, blas_int* kl, blas_int* ku, blas_int* nrhs, eT* ab, blas_int* ldab, eT* afb, blas_int* ldafb, blas_int* ipiv, char* equed, T* r, T* c, eT* b, blas_int* ldb, eT* x, blas_int* ldx, T* rcond, T* ferr, T* berr, eT* work, T* rwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cgbsvx)(fact, trans, n, kl, ku, nrhs, (cx_T*)ab, ldab, (cx_T*)afb, ldafb, ipiv, equed, (pod_T*)r, (pod_T*)c, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info, 1, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zgbsvx)(fact, trans, n, kl, ku, nrhs, (cx_T*)ab, ldab, (cx_T*)afb, ldafb, ipiv, equed, (pod_T*)r, (pod_T*)c, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info, 1, 1, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cgbsvx)(fact, trans, n, kl, ku, nrhs, (cx_T*)ab, ldab, (cx_T*)afb, ldafb, ipiv, equed, (pod_T*)r, (pod_T*)c, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zgbsvx)(fact, trans, n, kl, ku, nrhs, (cx_T*)ab, ldab, (cx_T*)afb, ldafb, ipiv, equed, (pod_T*)r, (pod_T*)c, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  gtsv(blas_int* n, blas_int* nrhs, eT* dl, eT* d, eT* du, eT* b, blas_int* ldb, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgtsv)(n, nrhs, (T*)dl, (T*)d, (T*)du, (T*)b, ldb, info); }
    else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgtsv)(n, nrhs, (T*)dl, (T*)d, (T*)du, (T*)b, ldb, info); }
    else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgtsv)(n, nrhs, (T*)dl, (T*)d, (T*)du, (T*)b, ldb, info); }
    else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgtsv)(n, nrhs, (T*)dl, (T*)d, (T*)du, (T*)b, ldb, info); }
    }
  
  
  
  template<typename eT>
  inline
  void
  gtsvx(char* fact, char* trans, blas_int* n, blas_int* nrhs, eT* dl, eT* d, eT* du, eT* dlf, eT* df, eT* duf, eT* du2, blas_int* ipiv, eT* b, blas_int* ldb, eT* x, blas_int* ldx, eT* rcond, eT* ferr, eT* berr, eT* work, blas_int* iwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgtsvx)(fact, trans, n, nrhs, (T*)dl, (T*)d, (T*)du, (T*)dlf, (T*)df, (T*)duf, (T*)du2, ipiv, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info, 1, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgtsvx)(fact, trans, n, nrhs, (T*)dl, (T*)d, (T*)du, (T*)dlf, (T*)df, (T*)duf, (T*)du2, ipiv, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info, 1, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgtsvx)(fact, trans, n, nrhs, (T*)dl, (T*)d, (T*)du, (T*)dlf, (T*)df, (T*)duf, (T*)du2, ipiv, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgtsvx)(fact, trans, n, nrhs, (T*)dl, (T*)d, (T*)du, (T*)dlf, (T*)df, (T*)duf, (T*)du2, ipiv, (T*)b, ldb, (T*)x, ldx, (T*)rcond, (T*)ferr, (T*)berr, (T*)work, iwork, info); }
    #endif
    }
  
  
  
  template<typename T, typename eT>
  inline
  void
  cx_gtsvx(char* fact, char* trans, blas_int* n, blas_int* nrhs, eT* dl, eT* d, eT* du, eT* dlf, eT* df, eT* duf, eT* du2, blas_int* ipiv, eT* b, blas_int* ldb, eT* x, blas_int* ldx, T* rcond, T* ferr, T* berr, eT* work, T* rwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cgtsvx)(fact, trans, n, nrhs, (cx_T*)dl, (cx_T*)d, (cx_T*)du, (cx_T*)dlf, (cx_T*)df, (cx_T*)duf, (cx_T*)du2, ipiv, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zgtsvx)(fact, trans, n, nrhs, (cx_T*)dl, (cx_T*)d, (cx_T*)du, (cx_T*)dlf, (cx_T*)df, (cx_T*)duf, (cx_T*)du2, ipiv, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info, 1, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cgtsvx)(fact, trans, n, nrhs, (cx_T*)dl, (cx_T*)d, (cx_T*)du, (cx_T*)dlf, (cx_T*)df, (cx_T*)duf, (cx_T*)du2, ipiv, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zgtsvx)(fact, trans, n, nrhs, (cx_T*)dl, (cx_T*)d, (cx_T*)du, (cx_T*)dlf, (cx_T*)df, (cx_T*)duf, (cx_T*)du2, ipiv, (cx_T*)b, ldb, (cx_T*)x, ldx, (pod_T*)rcond, (pod_T*)ferr, (pod_T*)berr, (cx_T*)work, (pod_T*)rwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  gees(char* jobvs, char* sort, void* select, blas_int* n, eT* a, blas_int* lda, blas_int* sdim, eT* wr, eT* wi, eT* vs, blas_int* ldvs, eT* work, blas_int* lwork, blas_int* bwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgees)(jobvs, sort, (fn_select_s2)select, n, (T*)a, lda, sdim, (T*)wr, (T*)wi, (T*)vs, ldvs, (T*)work, lwork, bwork, info, 1, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgees)(jobvs, sort, (fn_select_d2)select, n, (T*)a, lda, sdim, (T*)wr, (T*)wi, (T*)vs, ldvs, (T*)work, lwork, bwork, info, 1, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgees)(jobvs, sort, (fn_select_s2)select, n, (T*)a, lda, sdim, (T*)wr, (T*)wi, (T*)vs, ldvs, (T*)work, lwork, bwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgees)(jobvs, sort, (fn_select_d2)select, n, (T*)a, lda, sdim, (T*)wr, (T*)wi, (T*)vs, ldvs, (T*)work, lwork, bwork, info); }
    #endif
    }
  
  
  
  template<typename T>
  inline
  void
  cx_gees(char* jobvs, char* sort, void* select, blas_int* n, std::complex<T>* a, blas_int* lda, blas_int* sdim, std::complex<T>* w, std::complex<T>* vs, blas_int* ldvs, std::complex<T>* work, blas_int* lwork, T* rwork, blas_int* bwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<T>::value == false ));
    arma_type_check(( is_supported_blas_type< std::complex<T> >::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<T>::value)  { typedef float  bT; typedef blas_cxf cT; arma_fortran(arma_cgees)(jobvs, sort, (fn_select_c1)select, n, (cT*)a, lda, sdim, (cT*)w, (cT*)vs, ldvs, (cT*)work, lwork, (bT*)rwork, bwork, info, 1, 1); }
      else if(is_double<T>::value)  { typedef double bT; typedef blas_cxd cT; arma_fortran(arma_zgees)(jobvs, sort, (fn_select_z1)select, n, (cT*)a, lda, sdim, (cT*)w, (cT*)vs, ldvs, (cT*)work, lwork, (bT*)rwork, bwork, info, 1, 1); }
    #else
           if( is_float<T>::value)  { typedef float  bT; typedef blas_cxf cT; arma_fortran(arma_cgees)(jobvs, sort, (fn_select_c1)select, n, (cT*)a, lda, sdim, (cT*)w, (cT*)vs, ldvs, (cT*)work, lwork, (bT*)rwork, bwork, info); }
      else if(is_double<T>::value)  { typedef double bT; typedef blas_cxd cT; arma_fortran(arma_zgees)(jobvs, sort, (fn_select_z1)select, n, (cT*)a, lda, sdim, (cT*)w, (cT*)vs, ldvs, (cT*)work, lwork, (bT*)rwork, bwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  trsyl(char* transa, char* transb, blas_int* isgn, blas_int* m, blas_int* n, const eT* a, blas_int* lda, const eT* b, blas_int* ldb, eT* c, blas_int* ldc, eT* scale, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_strsyl)(transa, transb, isgn, m, n, (T*)a, lda, (T*)b, ldb, (T*)c, ldc, (T*)scale,      info, 1, 1); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dtrsyl)(transa, transb, isgn, m, n, (T*)a, lda, (T*)b, ldb, (T*)c, ldc, (T*)scale,      info, 1, 1); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_ctrsyl)(transa, transb, isgn, m, n, (T*)a, lda, (T*)b, ldb, (T*)c, ldc, (float*)scale,  info, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_ztrsyl)(transa, transb, isgn, m, n, (T*)a, lda, (T*)b, ldb, (T*)c, ldc, (double*)scale, info, 1, 1); }
    #else
           if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_strsyl)(transa, transb, isgn, m, n, (T*)a, lda, (T*)b, ldb, (T*)c, ldc, (T*)scale,      info); }
      else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dtrsyl)(transa, transb, isgn, m, n, (T*)a, lda, (T*)b, ldb, (T*)c, ldc, (T*)scale,      info); }
      else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_ctrsyl)(transa, transb, isgn, m, n, (T*)a, lda, (T*)b, ldb, (T*)c, ldc, (float*)scale,  info); }
      else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_ztrsyl)(transa, transb, isgn, m, n, (T*)a, lda, (T*)b, ldb, (T*)c, ldc, (double*)scale, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  gges
    (
    char* jobvsl, char* jobvsr, char* sort, void* selctg, blas_int* n,
    eT* a, blas_int* lda, eT* b, blas_int* ldb, blas_int* sdim,
    eT* alphar, eT* alphai, eT* beta,
    eT* vsl, blas_int* ldvsl, eT* vsr, blas_int* ldvsr,
    eT* work, blas_int* lwork,
    blas_int* bwork,
    blas_int* info
    )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgges)(jobvsl, jobvsr, sort, (fn_select_s3)selctg, n, (T*)a, lda, (T*)b, ldb, sdim, (T*)alphar, (T*)alphai, (T*)beta, (T*)vsl, ldvsl, (T*)vsr, ldvsr, (T*)work, lwork, bwork, info, 1, 1, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgges)(jobvsl, jobvsr, sort, (fn_select_d3)selctg, n, (T*)a, lda, (T*)b, ldb, sdim, (T*)alphar, (T*)alphai, (T*)beta, (T*)vsl, ldvsl, (T*)vsr, ldvsr, (T*)work, lwork, bwork, info, 1, 1, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgges)(jobvsl, jobvsr, sort, (fn_select_s3)selctg, n, (T*)a, lda, (T*)b, ldb, sdim, (T*)alphar, (T*)alphai, (T*)beta, (T*)vsl, ldvsl, (T*)vsr, ldvsr, (T*)work, lwork, bwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgges)(jobvsl, jobvsr, sort, (fn_select_d3)selctg, n, (T*)a, lda, (T*)b, ldb, sdim, (T*)alphar, (T*)alphai, (T*)beta, (T*)vsl, ldvsl, (T*)vsr, ldvsr, (T*)work, lwork, bwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  cx_gges
    (
    char* jobvsl, char* jobvsr, char* sort, void* selctg, blas_int* n,
    eT* a, blas_int* lda, eT* b, blas_int* ldb, blas_int* sdim,
    eT* alpha, eT* beta,
    eT* vsl, blas_int* ldvsl, eT* vsr, blas_int* ldvsr,
    eT* work, blas_int* lwork, typename eT::value_type* rwork,
    blas_int* bwork,
    blas_int* info
    )
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  T; typedef blas_cxf cx_T; arma_fortran(arma_cgges)(jobvsl, jobvsr, sort, (fn_select_c2)selctg, n, (cx_T*)a, lda, (cx_T*)b, ldb, sdim, (cx_T*)alpha, (cx_T*)beta, (cx_T*)vsl, ldvsl, (cx_T*)vsr, ldvsr, (cx_T*)work, lwork, (T*)rwork, bwork, info, 1, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef double T; typedef blas_cxd cx_T; arma_fortran(arma_zgges)(jobvsl, jobvsr, sort, (fn_select_z2)selctg, n, (cx_T*)a, lda, (cx_T*)b, ldb, sdim, (cx_T*)alpha, (cx_T*)beta, (cx_T*)vsl, ldvsl, (cx_T*)vsr, ldvsr, (cx_T*)work, lwork, (T*)rwork, bwork, info, 1, 1, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  T; typedef blas_cxf cx_T; arma_fortran(arma_cgges)(jobvsl, jobvsr, sort, (fn_select_c2)selctg, n, (cx_T*)a, lda, (cx_T*)b, ldb, sdim, (cx_T*)alpha, (cx_T*)beta, (cx_T*)vsl, ldvsl, (cx_T*)vsr, ldvsr, (cx_T*)work, lwork, (T*)rwork, bwork, info); }
      else if(is_cx_double<eT>::value)  { typedef double T; typedef blas_cxd cx_T; arma_fortran(arma_zgges)(jobvsl, jobvsr, sort, (fn_select_z2)selctg, n, (cx_T*)a, lda, (cx_T*)b, ldb, sdim, (cx_T*)alpha, (cx_T*)beta, (cx_T*)vsl, ldvsl, (cx_T*)vsr, ldvsr, (cx_T*)work, lwork, (T*)rwork, bwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  typename get_pod_type<eT>::result
  lange(char* norm, blas_int* m, blas_int* n, eT* a, blas_int* lda, typename get_pod_type<eT>::result* work)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    typedef typename get_pod_type<eT>::result out_T;
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float  pod_T; typedef float    T; return out_T( arma_fortran(arma_slange)(norm, m, n, (T*)a, lda, (pod_T*)work, 1) ); }
      else if(   is_double<eT>::value)  { typedef double pod_T; typedef double   T; return out_T( arma_fortran(arma_dlange)(norm, m, n, (T*)a, lda, (pod_T*)work, 1) ); }
      else if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf T; return out_T( arma_fortran(arma_clange)(norm, m, n, (T*)a, lda, (pod_T*)work, 1) ); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd T; return out_T( arma_fortran(arma_zlange)(norm, m, n, (T*)a, lda, (pod_T*)work, 1) ); }
    #else
           if(    is_float<eT>::value)  { typedef float  pod_T; typedef float    T; return out_T( arma_fortran(arma_slange)(norm, m, n, (T*)a, lda, (pod_T*)work) ); }
      else if(   is_double<eT>::value)  { typedef double pod_T; typedef double   T; return out_T( arma_fortran(arma_dlange)(norm, m, n, (T*)a, lda, (pod_T*)work) ); }
      else if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf T; return out_T( arma_fortran(arma_clange)(norm, m, n, (T*)a, lda, (pod_T*)work) ); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd T; return out_T( arma_fortran(arma_zlange)(norm, m, n, (T*)a, lda, (pod_T*)work) ); }
    #endif
    
    return out_T(0);
    }
  


  template<typename eT>
  inline
  typename get_pod_type<eT>::result
  lansy(char* norm, char* uplo, blas_int* n, eT* a, blas_int* lda, typename get_pod_type<eT>::result* work)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    typedef typename get_pod_type<eT>::result out_T;
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float  pod_T; typedef float    T; return out_T( arma_fortran(arma_slansy)(norm, uplo, n, (T*)a, lda, (pod_T*)work, 1, 1) ); }
      else if(   is_double<eT>::value)  { typedef double pod_T; typedef double   T; return out_T( arma_fortran(arma_dlansy)(norm, uplo, n, (T*)a, lda, (pod_T*)work, 1, 1) ); }
      else if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf T; return out_T( arma_fortran(arma_clansy)(norm, uplo, n, (T*)a, lda, (pod_T*)work, 1, 1) ); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd T; return out_T( arma_fortran(arma_zlansy)(norm, uplo, n, (T*)a, lda, (pod_T*)work, 1, 1) ); }
    #else
           if(    is_float<eT>::value)  { typedef float  pod_T; typedef float    T; return out_T( arma_fortran(arma_slansy)(norm, uplo, n, (T*)a, lda, (pod_T*)work) ); }
      else if(   is_double<eT>::value)  { typedef double pod_T; typedef double   T; return out_T( arma_fortran(arma_dlansy)(norm, uplo, n, (T*)a, lda, (pod_T*)work) ); }
      else if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf T; return out_T( arma_fortran(arma_clansy)(norm, uplo, n, (T*)a, lda, (pod_T*)work) ); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd T; return out_T( arma_fortran(arma_zlansy)(norm, uplo, n, (T*)a, lda, (pod_T*)work) ); }
    #endif
    
    return out_T(0);
    }
  
  
  
  template<typename eT>
  inline
  typename get_pod_type<eT>::result
  lanhe(char* norm, char* uplo, blas_int* n, eT* a, blas_int* lda, typename get_pod_type<eT>::result* work)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    typedef typename get_pod_type<eT>::result out_T;
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf T; return out_T( arma_fortran(arma_clanhe)(norm, uplo, n, (T*)a, lda, (pod_T*)work, 1, 1) ); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd T; return out_T( arma_fortran(arma_zlanhe)(norm, uplo, n, (T*)a, lda, (pod_T*)work, 1, 1) ); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf T; return out_T( arma_fortran(arma_clanhe)(norm, uplo, n, (T*)a, lda, (pod_T*)work) ); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd T; return out_T( arma_fortran(arma_zlanhe)(norm, uplo, n, (T*)a, lda, (pod_T*)work) ); }
    #endif
    
    return out_T(0);
    }
  


  template<typename eT>
  inline
  typename get_pod_type<eT>::result
  langb(char* norm, blas_int* n, blas_int* kl, blas_int* ku, eT* ab, blas_int* ldab, typename get_pod_type<eT>::result* work)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    typedef typename get_pod_type<eT>::result out_T;
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float  pod_T; typedef float    T; return out_T( arma_fortran(arma_slangb)(norm, n, kl, ku, (T*)ab, ldab, (pod_T*)work, 1) ); }
      else if(   is_double<eT>::value)  { typedef double pod_T; typedef double   T; return out_T( arma_fortran(arma_dlangb)(norm, n, kl, ku, (T*)ab, ldab, (pod_T*)work, 1) ); }
      else if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf T; return out_T( arma_fortran(arma_clangb)(norm, n, kl, ku, (T*)ab, ldab, (pod_T*)work, 1) ); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd T; return out_T( arma_fortran(arma_zlangb)(norm, n, kl, ku, (T*)ab, ldab, (pod_T*)work, 1) ); }
    #else
           if(    is_float<eT>::value)  { typedef float  pod_T; typedef float    T; return out_T( arma_fortran(arma_slangb)(norm, n, kl, ku, (T*)ab, ldab, (pod_T*)work) ); }
      else if(   is_double<eT>::value)  { typedef double pod_T; typedef double   T; return out_T( arma_fortran(arma_dlangb)(norm, n, kl, ku, (T*)ab, ldab, (pod_T*)work) ); }
      else if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf T; return out_T( arma_fortran(arma_clangb)(norm, n, kl, ku, (T*)ab, ldab, (pod_T*)work) ); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd T; return out_T( arma_fortran(arma_zlangb)(norm, n, kl, ku, (T*)ab, ldab, (pod_T*)work) ); }
    #endif
    
    return out_T(0);
    }
  
  
  
  template<typename eT>
  inline
  void
  gecon(char* norm, blas_int* n, const eT* a, blas_int* lda, const eT* anorm, eT* rcond, eT* work, blas_int* iwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgecon)(norm, n, (T*)a, lda, (T*)anorm, (T*)rcond, (T*)work, iwork, info, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgecon)(norm, n, (T*)a, lda, (T*)anorm, (T*)rcond, (T*)work, iwork, info, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgecon)(norm, n, (T*)a, lda, (T*)anorm, (T*)rcond, (T*)work, iwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgecon)(norm, n, (T*)a, lda, (T*)anorm, (T*)rcond, (T*)work, iwork, info); }
    #endif
    }
  
  
  
  template<typename T>
  inline
  void
  cx_gecon(char* norm, blas_int* n, const std::complex<T>* a, blas_int* lda, const T* anorm, T* rcond, std::complex<T>* work, T* rwork, blas_int* info)
    {
    typedef typename std::complex<T> eT;
    
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cgecon)(norm, n, (cx_T*)a, lda, (pod_T*)anorm, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info, 1); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zgecon)(norm, n, (cx_T*)a, lda, (pod_T*)anorm, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cgecon)(norm, n, (cx_T*)a, lda, (pod_T*)anorm, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zgecon)(norm, n, (cx_T*)a, lda, (pod_T*)anorm, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  pocon(char* uplo, blas_int* n, const eT* a, blas_int* lda, const eT* anorm, eT* rcond, eT* work, blas_int* iwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_spocon)(uplo, n, (T*)a, lda, (T*)anorm, (T*)rcond, (T*)work, iwork, info, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dpocon)(uplo, n, (T*)a, lda, (T*)anorm, (T*)rcond, (T*)work, iwork, info, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_spocon)(uplo, n, (T*)a, lda, (T*)anorm, (T*)rcond, (T*)work, iwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dpocon)(uplo, n, (T*)a, lda, (T*)anorm, (T*)rcond, (T*)work, iwork, info); }
    #endif
    }
  
  
  
  template<typename T>
  inline
  void
  cx_pocon(char* uplo, blas_int* n, const std::complex<T>* a, blas_int* lda, const T* anorm, T* rcond, std::complex<T>* work, T* rwork, blas_int* info)
    {
    typedef typename std::complex<T> eT;
    
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cpocon)(uplo, n, (cx_T*)a, lda, (pod_T*)anorm, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info, 1); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zpocon)(uplo, n, (cx_T*)a, lda, (pod_T*)anorm, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cpocon)(uplo, n, (cx_T*)a, lda, (pod_T*)anorm, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zpocon)(uplo, n, (cx_T*)a, lda, (pod_T*)anorm, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info); }
    #endif
    }
  


  template<typename eT>
  inline
  void
  trcon(char* norm, char* uplo, char* diag, blas_int* n, const eT* a, blas_int* lda, eT* rcond, eT* work, blas_int* iwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_strcon)(norm, uplo, diag, n, (T*)a, lda, (T*)rcond, (T*)work, iwork, info, 1, 1, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dtrcon)(norm, uplo, diag, n, (T*)a, lda, (T*)rcond, (T*)work, iwork, info, 1, 1, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_strcon)(norm, uplo, diag, n, (T*)a, lda, (T*)rcond, (T*)work, iwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dtrcon)(norm, uplo, diag, n, (T*)a, lda, (T*)rcond, (T*)work, iwork, info); }
    #endif
    }
  
  
  
  template<typename T>
  inline
  void
  cx_trcon(char* norm, char* uplo, char* diag, blas_int* n, const std::complex<T>* a, blas_int* lda, T* rcond, std::complex<T>* work, T* rwork, blas_int* info)
    {
    typedef typename std::complex<T> eT;
    
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_ctrcon)(norm, uplo, diag, n, (cx_T*)a, lda, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info, 1, 1, 1); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_ztrcon)(norm, uplo, diag, n, (cx_T*)a, lda, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info, 1, 1, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_ctrcon)(norm, uplo, diag, n, (cx_T*)a, lda, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_ztrcon)(norm, uplo, diag, n, (cx_T*)a, lda, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  gbcon(char* norm, blas_int* n, blas_int* kl, blas_int* ku, const eT* ab, blas_int* ldab, const blas_int* ipiv, const eT* anorm, eT* rcond, eT* work, blas_int* iwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgbcon)(norm, n, kl, ku, (T*)ab, ldab, ipiv, (T*)anorm, (T*)rcond, (T*)work, iwork, info, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgbcon)(norm, n, kl, ku, (T*)ab, ldab, ipiv, (T*)anorm, (T*)rcond, (T*)work, iwork, info, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sgbcon)(norm, n, kl, ku, (T*)ab, ldab, ipiv, (T*)anorm, (T*)rcond, (T*)work, iwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dgbcon)(norm, n, kl, ku, (T*)ab, ldab, ipiv, (T*)anorm, (T*)rcond, (T*)work, iwork, info); }
    #endif
    }
  
  
  
  template<typename T>
  inline
  void
  cx_gbcon(char* norm, blas_int* n, blas_int* kl, blas_int* ku, const std::complex<T>* ab, blas_int* ldab, const blas_int* ipiv, const T* anorm, T* rcond, std::complex<T>* work, T* rwork, blas_int* info)
    {
    typedef typename std::complex<T> eT;
    
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cgbcon)(norm, n, kl, ku, (cx_T*)ab, ldab, ipiv, (pod_T*)anorm, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info, 1); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zgbcon)(norm, n, kl, ku, (cx_T*)ab, ldab, ipiv, (pod_T*)anorm, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info, 1); }
    #else
           if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf cx_T; arma_fortran(arma_cgbcon)(norm, n, kl, ku, (cx_T*)ab, ldab, ipiv, (pod_T*)anorm, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd cx_T; arma_fortran(arma_zgbcon)(norm, n, kl, ku, (cx_T*)ab, ldab, ipiv, (pod_T*)anorm, (pod_T*)rcond, (cx_T*)work, (pod_T*)rwork, info); }
    #endif
    }
  
  
  
  inline
  blas_int
  laenv(blas_int* ispec, char* name, char* opts, blas_int* n1, blas_int* n2, blas_int* n3, blas_int* n4, blas_len name_len, blas_len opts_len)
    {
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
      return arma_fortran(arma_ilaenv)(ispec, name, opts, n1, n2, n3, n4, name_len, opts_len);
    #else
      arma_ignore(name_len);
      arma_ignore(opts_len);
      return arma_fortran(arma_ilaenv)(ispec, name, opts, n1, n2, n3, n4);  // not advised!
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  lahqr(blas_int* wantt, blas_int* wantz, blas_int* n, blas_int* ilo, blas_int* ihi, eT* h, blas_int* ldh, eT* wr, eT* wi, blas_int* iloz, blas_int* ihiz, eT* z, blas_int* ldz, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_slahqr)(wantt, wantz, n, ilo, ihi, (T*)h, ldh, (T*)wr, (T*)wi, iloz, ihiz, (T*)z, ldz, info); }
    else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dlahqr)(wantt, wantz, n, ilo, ihi, (T*)h, ldh, (T*)wr, (T*)wi, iloz, ihiz, (T*)z, ldz, info); }
    }
  
  
  
  template<typename eT>
  inline
  void
  stedc(char* compz, blas_int* n, eT* d, eT* e, eT* z, blas_int* ldz, eT* work, blas_int* lwork, blas_int* iwork, blas_int* liwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sstedc)(compz, n, (T*)d, (T*)e, (T*)z, ldz, (T*)work, lwork, iwork, liwork, info, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dstedc)(compz, n, (T*)d, (T*)e, (T*)z, ldz, (T*)work, lwork, iwork, liwork, info, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_sstedc)(compz, n, (T*)d, (T*)e, (T*)z, ldz, (T*)work, lwork, iwork, liwork, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dstedc)(compz, n, (T*)d, (T*)e, (T*)z, ldz, (T*)work, lwork, iwork, liwork, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  trevc(char* side, char* howmny, blas_int* select, blas_int* n, eT* t, blas_int* ldt, eT* vl, blas_int* ldvl, eT* vr, blas_int* ldvr, blas_int* mm, blas_int* m, eT* work, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_strevc)(side, howmny, select, n, (T*)t, ldt, (T*)vl, ldvl, (T*)vr, ldvr, mm, m, (T*)work, info, 1, 1); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dtrevc)(side, howmny, select, n, (T*)t, ldt, (T*)vl, ldvl, (T*)vr, ldvr, mm, m, (T*)work, info, 1, 1); }
    #else
           if( is_float<eT>::value)  { typedef float  T; arma_fortran(arma_strevc)(side, howmny, select, n, (T*)t, ldt, (T*)vl, ldvl, (T*)vr, ldvr, mm, m, (T*)work, info); }
      else if(is_double<eT>::value)  { typedef double T; arma_fortran(arma_dtrevc)(side, howmny, select, n, (T*)t, ldt, (T*)vl, ldvl, (T*)vr, ldvr, mm, m, (T*)work, info); }
    #endif
    }
  
  
  
  template<typename eT>
  inline
  void
  gehrd(blas_int* n, blas_int* ilo, blas_int* ihi, eT* a, blas_int* lda, eT* tao, eT* work, blas_int* lwork, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
         if(    is_float<eT>::value)  { typedef float    T; arma_fortran(arma_sgehrd)(n, ilo, ihi, (T*)a, lda, (T*)tao, (T*)work, lwork, info); }
    else if(   is_double<eT>::value)  { typedef double   T; arma_fortran(arma_dgehrd)(n, ilo, ihi, (T*)a, lda, (T*)tao, (T*)work, lwork, info); }
    else if( is_cx_float<eT>::value)  { typedef blas_cxf T; arma_fortran(arma_cgehrd)(n, ilo, ihi, (T*)a, lda, (T*)tao, (T*)work, lwork, info); }
    else if(is_cx_double<eT>::value)  { typedef blas_cxd T; arma_fortran(arma_zgehrd)(n, ilo, ihi, (T*)a, lda, (T*)tao, (T*)work, lwork, info); }
    }
  
  
  
  template<typename eT>
  inline
  void
  pstrf(const char* uplo, const blas_int* n, eT* a, const blas_int* lda, blas_int* piv, blas_int* rank, const typename get_pod_type<eT>::result* tol, const typename get_pod_type<eT>::result* work, blas_int* info)
    {
    arma_type_check(( is_supported_blas_type<eT>::value == false ));
    
    #if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
           if(    is_float<eT>::value)  { typedef float  pod_T; typedef float    T; arma_fortran(arma_spstrf)(uplo, n, (T*)a, lda, piv, rank, (const pod_T*)tol, (pod_T*)work, info, 1); }
      else if(   is_double<eT>::value)  { typedef double pod_T; typedef double   T; arma_fortran(arma_dpstrf)(uplo, n, (T*)a, lda, piv, rank, (const pod_T*)tol, (pod_T*)work, info, 1); }
      else if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf T; arma_fortran(arma_cpstrf)(uplo, n, (T*)a, lda, piv, rank, (const pod_T*)tol, (pod_T*)work, info, 1); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd T; arma_fortran(arma_zpstrf)(uplo, n, (T*)a, lda, piv, rank, (const pod_T*)tol, (pod_T*)work, info, 1); }
    #else
           if(    is_float<eT>::value)  { typedef float  pod_T; typedef float    T; arma_fortran(arma_spstrf)(uplo, n, (T*)a, lda, piv, rank, (const pod_T*)tol, (pod_T*)work, info); }
      else if(   is_double<eT>::value)  { typedef double pod_T; typedef double   T; arma_fortran(arma_dpstrf)(uplo, n, (T*)a, lda, piv, rank, (const pod_T*)tol, (pod_T*)work, info); }
      else if( is_cx_float<eT>::value)  { typedef float  pod_T; typedef blas_cxf T; arma_fortran(arma_cpstrf)(uplo, n, (T*)a, lda, piv, rank, (const pod_T*)tol, (pod_T*)work, info); }
      else if(is_cx_double<eT>::value)  { typedef double pod_T; typedef blas_cxd T; arma_fortran(arma_zpstrf)(uplo, n, (T*)a, lda, piv, rank, (const pod_T*)tol, (pod_T*)work, info); }
    #endif
    }
  
  
  }


#endif
