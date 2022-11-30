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


//! \addtogroup auxlib
//! @{



template<typename eT>
inline
bool
auxlib::inv(Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  if(A.is_empty())  { return true; }
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(A);
    
    blas_int n     = blas_int(A.n_rows);
    blas_int lda   = blas_int(A.n_rows);
    blas_int lwork = (std::max)(blas_int(podarray_prealloc_n_elem::val), n);
    blas_int info  = 0;
    
    podarray<blas_int> ipiv(A.n_rows);
    
    arma_extra_debug_print("lapack::getrf()");
    lapack::getrf(&n, &n, A.memptr(), &lda, ipiv.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    if(n > 16)
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = -1;
      
      arma_extra_debug_print("lapack::getri()");
      lapack::getri(&n, A.memptr(), &lda, ipiv.memptr(), &work_query[0], &lwork_query, &info);
      
      if(info != 0)  { return false; }
      
      blas_int lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
      
      lwork = (std::max)(lwork_proposed, lwork);
      }
    
    podarray<eT> work( static_cast<uword>(lwork) );
    
    arma_extra_debug_print("lapack::getri()");
    lapack::getri(&n, A.memptr(), &lda, ipiv.memptr(), work.memptr(), &lwork, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(A);
    arma_stop_logic_error("inv(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::inv(Mat<eT>& out, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  out = X;
  
  return auxlib::inv(out);
  }



template<typename eT>
inline
bool
auxlib::inv_rcond(Mat<eT>& A, typename get_pod_type<eT>::result& out_rcond)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  out_rcond = T(0);
  
  if(A.is_empty())  { return true; }
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(A);
    
    char     norm_id  = '1';
    blas_int n        = blas_int(A.n_rows);
    blas_int lda      = blas_int(A.n_rows);
    blas_int lwork    = (std::max)(blas_int(podarray_prealloc_n_elem::val), n);
    blas_int info     = 0;
    T        norm_val = T(0);
    
    podarray<T>        junk(1);
    podarray<blas_int> ipiv(A.n_rows);
    
    arma_extra_debug_print("lapack::lange()");
    norm_val = lapack::lange<eT>(&norm_id, &n, &n, A.memptr(), &lda, junk.memptr());
    
    arma_extra_debug_print("lapack::getrf()");
    lapack::getrf(&n, &n, A.memptr(), &lda, ipiv.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    out_rcond = auxlib::lu_rcond<T>(A, norm_val);
    
    if(n > 16)
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = -1;
      
      arma_extra_debug_print("lapack::getri()");
      lapack::getri(&n, A.memptr(), &lda, ipiv.memptr(), &work_query[0], &lwork_query, &info);
      
      if(info != 0)  { return false; }
      
      blas_int lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
      
      lwork = (std::max)(lwork_proposed, lwork);
      }
    
    podarray<eT> work( static_cast<uword>(lwork) );
    
    arma_extra_debug_print("lapack::getri()");
    lapack::getri(&n, A.memptr(), &lda, ipiv.memptr(), work.memptr(), &lwork, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(A);
    arma_stop_logic_error("inv_rcond(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::inv_tr(Mat<eT>& A, const uword layout)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    if(A.is_empty())  { return true; }
    
    arma_debug_assert_blas_size(A);
    
    char     uplo = (layout == 0) ? 'U' : 'L';
    char     diag = 'N';
    blas_int n    = blas_int(A.n_rows);
    blas_int info = 0;
    
    arma_extra_debug_print("lapack::trtri()");
    lapack::trtri(&uplo, &diag, &n, A.memptr(), &n, &info);
    
    if(info != 0)  { return false; }
    
    return true;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(layout);
    arma_stop_logic_error("inv(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::inv_tr_rcond(Mat<eT>& A, typename get_pod_type<eT>::result& out_rcond, const uword layout)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename get_pod_type<eT>::result T;
    
    if(A.is_empty())  { return true; }
    
    out_rcond = auxlib::rcond_trimat(A, layout);
    
    arma_debug_assert_blas_size(A);
    
    char     uplo = (layout == 0) ? 'U' : 'L';
    char     diag = 'N';
    blas_int n    = blas_int(A.n_rows);
    blas_int info = 0;
    
    arma_extra_debug_print("lapack::trtri()");
    lapack::trtri(&uplo, &diag, &n, A.memptr(), &n, &info);
    
    if(info != 0)  { out_rcond = T(0); return false; }
    
    return true;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(out_rcond);
    arma_ignore(layout);
    arma_stop_logic_error("inv(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::inv_sympd(Mat<eT>& A, bool& out_sympd_state)
  {
  arma_extra_debug_sigprint();
  
  out_sympd_state = false;
  
  if(A.is_empty())  { return true; }
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(A);
    
    char     uplo = 'L';
    blas_int n    = blas_int(A.n_rows);
    blas_int info = 0;
    
    // NOTE: for complex matrices, zpotrf() assumes the matrix is hermitian (not simply symmetric)
    
    arma_extra_debug_print("lapack::potrf()");
    lapack::potrf(&uplo, &n, A.memptr(), &n, &info);
    
    if(info != 0)  { return false; }
    
    out_sympd_state = true;
    
    arma_extra_debug_print("lapack::potri()");
    lapack::potri(&uplo, &n, A.memptr(), &n, &info);
    
    if(info != 0)  { return false; }
    
    A = symmatl(A);
    
    return true;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(out_sympd_state);
    arma_stop_logic_error("inv_sympd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::inv_sympd(Mat<eT>& out, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  out = X;
  
  bool sympd_state_junk = false;
  
  return auxlib::inv_sympd(out, sympd_state_junk);
  }



template<typename eT>
inline
bool
auxlib::inv_sympd_rcond(Mat<eT>& A, bool& out_sympd_state, eT& out_rcond, const eT rcond_threshold)
  {
  arma_extra_debug_sigprint();
  
  out_sympd_state = false;
  
  if(A.is_empty())  { return true; }
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename get_pod_type<eT>::result T;
    
    arma_debug_assert_blas_size(A);
    
    char     norm_id  = '1';
    char     uplo     = 'L';
    blas_int n        = blas_int(A.n_rows);
    blas_int info     = 0;
    T        norm_val = T(0);
    
    podarray<T> work(A.n_rows);
    
    arma_extra_debug_print("lapack::lansy()");
    norm_val = lapack::lansy(&norm_id, &uplo, &n, A.memptr(), &n, work.memptr());
    
    arma_extra_debug_print("lapack::potrf()");
    lapack::potrf(&uplo, &n, A.memptr(), &n, &info);
    
    if(info != 0)  { out_rcond = eT(0); return false; }
    
    out_sympd_state = true;
    
    out_rcond = auxlib::lu_rcond_sympd<T>(A, norm_val);
    
    if( arma_isnan(out_rcond) || ((rcond_threshold > eT(0)) && (out_rcond < rcond_threshold)) )  { return false; }
    
    arma_extra_debug_print("lapack::potri()");
    lapack::potri(&uplo, &n, A.memptr(), &n, &info);
    
    if(info != 0)  { return false; }
    
    A = symmatl(A);
    
    return true;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(out_sympd_state);
    arma_ignore(out_rcond);
    arma_ignore(rcond_threshold);
    arma_stop_logic_error("inv_sympd_rcond(): use LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T>
inline
bool
auxlib::inv_sympd_rcond(Mat< std::complex<T> >& A, bool& out_sympd_state, T& out_rcond, const T rcond_threshold)
  {
  arma_extra_debug_sigprint();
  
  out_sympd_state = false;
  
  if(A.is_empty())  { return true; }
  
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_ignore(A);
    arma_ignore(out_sympd_state);
    arma_ignore(out_rcond);
    arma_ignore(rcond_threshold);
    return false;
    }
  #elif defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(A);
    
    char     norm_id  = '1';
    char     uplo     = 'L';
    blas_int n        = blas_int(A.n_rows);
    blas_int info     = 0;
    T        norm_val = T(0);
    
    podarray<T> work(A.n_rows);
    
    arma_extra_debug_print("lapack::lanhe()");
    norm_val = lapack::lanhe(&norm_id, &uplo, &n, A.memptr(), &n, work.memptr());
    
    arma_extra_debug_print("lapack::potrf()");
    lapack::potrf(&uplo, &n, A.memptr(), &n, &info);
    
    if(info != 0)  { out_rcond = T(0); return false; }
    
    out_sympd_state = true;
    
    out_rcond = auxlib::lu_rcond_sympd<T>(A, norm_val);
    
    if( arma_isnan(out_rcond) || ((rcond_threshold > T(0)) && (out_rcond < rcond_threshold)) )  { return false; }
    
    arma_extra_debug_print("lapack::potri()");
    lapack::potri(&uplo, &n, A.memptr(), &n, &info);
    
    if(info != 0)  { return false; }
    
    A = symmatl(A);
    
    return true;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(out_sympd_state);
    arma_ignore(out_rcond);
    arma_ignore(rcond_threshold);
    arma_stop_logic_error("inv_sympd_rcond(): use LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! determinant of a matrix
template<typename eT>
inline
bool
auxlib::det(eT& out_val, Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  if(A.is_empty())  { out_val = eT(1); return true; }
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(A);
    
    podarray<blas_int> ipiv(A.n_rows);
    
    blas_int info   = 0;
    blas_int n_rows = blas_int(A.n_rows);
    blas_int n_cols = blas_int(A.n_cols);
    
    arma_extra_debug_print("lapack::getrf()");
    lapack::getrf(&n_rows, &n_cols, A.memptr(), &n_rows, ipiv.memptr(), &info);
    
    if(info < 0)  { return false; }
    
    // on output A appears to be L+U_alt, where U_alt is U with the main diagonal set to zero
    eT val = A.at(0,0);
    for(uword i=1; i < A.n_rows; ++i)  { val *= A.at(i,i); }
    
    blas_int sign = +1;
    for(uword i=0; i < A.n_rows; ++i)
      {
      // NOTE: adjustment of -1 is required as Fortran counts from 1
      if( blas_int(i) != (ipiv.mem[i] - 1) )  { sign *= -1; }
      }
    
    out_val = (sign < 0) ? eT(-val) : eT(val);
    
    return true;
    }
  #else
    {
    arma_ignore(out_val);
    arma_ignore(A);
    arma_stop_logic_error("det(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! log determinant of a matrix
template<typename eT>
inline
bool
auxlib::log_det(eT& out_val, typename get_pod_type<eT>::result& out_sign, Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  if(A.is_empty())  { out_val  = eT(0); out_sign =  T(1); return true; }
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(A);
    
    podarray<blas_int> ipiv(A.n_rows);
    
    blas_int info   = 0;
    blas_int n_rows = blas_int(A.n_rows);
    blas_int n_cols = blas_int(A.n_cols);
    
    arma_extra_debug_print("lapack::getrf()");
    lapack::getrf(&n_rows, &n_cols, A.memptr(), &n_rows, ipiv.memptr(), &info);
    
    if(info < 0)  { return false; }
    
    // on output A appears to be L+U_alt, where U_alt is U with the main diagonal set to zero
    
    sword sign = (is_cx<eT>::no) ? ( (access::tmp_real( A.at(0,0) ) < T(0)) ? -1 : +1 ) : +1;
    eT    val  = (is_cx<eT>::no) ? std::log( (access::tmp_real( A.at(0,0) ) < T(0)) ? A.at(0,0)*T(-1) : A.at(0,0) ) : std::log( A.at(0,0) );
    
    for(uword i=1; i < A.n_rows; ++i)
      {
      const eT x = A.at(i,i);
      
      sign *= (is_cx<eT>::no) ? ( (access::tmp_real(x) < T(0)) ? -1 : +1 ) : +1;
      val  += (is_cx<eT>::no) ? std::log( (access::tmp_real(x) < T(0)) ? x*T(-1) : x ) : std::log(x);
      }
    
    for(uword i=0; i < A.n_rows; ++i)
      {
      if( blas_int(i) != (ipiv.mem[i] - 1) )  // NOTE: adjustment of -1 is required as Fortran counts from 1
        {
        sign *= -1;
        }
      }
    
    out_val  = val;
    out_sign = T(sign);
    
    return true;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(out_val);
    arma_ignore(out_sign);
    arma_stop_logic_error("log_det(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::log_det_sympd(typename get_pod_type<eT>::result& out_val, Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  if(A.is_empty())  { out_val = T(0); return true; }
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(A);
    
    char     uplo = 'L';
    blas_int n    = blas_int(A.n_rows);
    blas_int info = 0;
    
    arma_extra_debug_print("lapack::potrf()");
    lapack::potrf(&uplo, &n, A.memptr(), &n, &info);
    
    if(info != 0)  { return false; }
    
    T val = T(0);
    
    for(uword i=0; i < A.n_rows; ++i)  { val += std::log( access::tmp_real(A.at(i,i)) ); }
    
    out_val = T(2) * val;
    
    return true;
    }
  #else
    {
    arma_ignore(out_val);
    arma_ignore(A);
    arma_stop_logic_error("log_det_sympd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! LU decomposition of a matrix
template<typename eT, typename T1>
inline
bool
auxlib::lu(Mat<eT>& L, Mat<eT>& U, podarray<blas_int>& ipiv, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  U = X.get_ref();
  
  const uword U_n_rows = U.n_rows;
  const uword U_n_cols = U.n_cols;
  
  if(U.is_empty())  { L.set_size(U_n_rows, 0); U.set_size(0, U_n_cols); ipiv.reset(); return true; }
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(U);
    
    ipiv.set_size( (std::min)(U_n_rows, U_n_cols) );
    
    blas_int info = 0;
    
    blas_int n_rows = blas_int(U_n_rows);
    blas_int n_cols = blas_int(U_n_cols);
    
    arma_extra_debug_print("lapack::getrf()");
    lapack::getrf(&n_rows, &n_cols, U.memptr(), &n_rows, ipiv.memptr(), &info);
    
    if(info < 0)  { return false; }
    
    // take into account that Fortran counts from 1
    arrayops::inplace_minus(ipiv.memptr(), blas_int(1), ipiv.n_elem);
    
    L.copy_size(U);
    
    for(uword col=0; col < U_n_cols; ++col)
      {
      for(uword row=0; (row < col) && (row < U_n_rows); ++row)
        {
        L.at(row,col) = eT(0);
        }
      
      if( L.in_range(col,col) )
        {
        L.at(col,col) = eT(1);
        }
      
      for(uword row = (col+1); row < U_n_rows; ++row)
        {
        L.at(row,col) = U.at(row,col);
        U.at(row,col) = eT(0);
        }
      }
    
    return true;
    }
  #else
    {
    arma_stop_logic_error("lu(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT, typename T1>
inline
bool
auxlib::lu(Mat<eT>& L, Mat<eT>& U, Mat<eT>& P, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  podarray<blas_int> ipiv1;
  const bool status = auxlib::lu(L, U, ipiv1, X);
  
  if(status == false)  { return false; }
  
  if(U.is_empty())
    {
    // L and U have been already set to the correct empty matrices
    P.eye(L.n_rows, L.n_rows);
    return true;
    }
  
  const uword n      = ipiv1.n_elem;
  const uword P_rows = U.n_rows;
  
  podarray<blas_int> ipiv2(P_rows);
  
  const blas_int* ipiv1_mem = ipiv1.memptr();
        blas_int* ipiv2_mem = ipiv2.memptr();
  
  for(uword i=0; i<P_rows; ++i)
    {
    ipiv2_mem[i] = blas_int(i);
    }
  
  for(uword i=0; i<n; ++i)
    {
    const uword k = static_cast<uword>(ipiv1_mem[i]);
    
    if( ipiv2_mem[i] != ipiv2_mem[k] )
      {
      std::swap( ipiv2_mem[i], ipiv2_mem[k] );
      }
    }
  
  P.zeros(P_rows, P_rows);
  
  for(uword row=0; row<P_rows; ++row)
    {
    P.at(row, static_cast<uword>(ipiv2_mem[row])) = eT(1);
    }
  
  if(L.n_cols > U.n_rows)
    {
    L.shed_cols(U.n_rows, L.n_cols-1);
    }
    
  if(U.n_rows > L.n_cols)
    {
    U.shed_rows(L.n_cols, U.n_rows-1);
    }
  
  return true;
  }



template<typename eT, typename T1>
inline
bool
auxlib::lu(Mat<eT>& L, Mat<eT>& U, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  podarray<blas_int> ipiv1;
  const bool status = auxlib::lu(L, U, ipiv1, X);
  
  if(status == false)  { return false; }
  
  if(U.is_empty())
    {
    // L and U have been already set to the correct empty matrices
    return true;
    }
  
  const uword n      = ipiv1.n_elem;
  const uword P_rows = U.n_rows;
  
  podarray<blas_int> ipiv2(P_rows);
  
  const blas_int* ipiv1_mem = ipiv1.memptr();
        blas_int* ipiv2_mem = ipiv2.memptr();
  
  for(uword i=0; i<P_rows; ++i)
    {
    ipiv2_mem[i] = blas_int(i);
    }
  
  for(uword i=0; i<n; ++i)
    {
    const uword k = static_cast<uword>(ipiv1_mem[i]);
    
    if( ipiv2_mem[i] != ipiv2_mem[k] )
      {
      std::swap( ipiv2_mem[i], ipiv2_mem[k] );
      L.swap_rows( static_cast<uword>(ipiv2_mem[i]), static_cast<uword>(ipiv2_mem[k]) );
      }
    }
  
  if(L.n_cols > U.n_rows)
    {
    L.shed_cols(U.n_rows, L.n_cols-1);
    }
    
  if(U.n_rows > L.n_cols)
    {
    U.shed_rows(L.n_cols, U.n_rows-1);
    }
  
  return true;
  }



//! eigen decomposition of general square matrix (real)
template<typename T1>
inline
bool
auxlib::eig_gen
  (
         Mat< std::complex<typename T1::pod_type> >& vals,
         Mat< std::complex<typename T1::pod_type> >& vecs,
  const bool                                         vecs_on,
  const Base<typename T1::pod_type,T1>&              expr
  )
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type T;
    
    Mat<T> X = expr.get_ref();
    
    arma_debug_check( (X.is_square() == false), "eig_gen(): given matrix must be square sized" );
    
    arma_debug_assert_blas_size(X);
    
    if(X.is_empty())  { vals.reset(); vecs.reset(); return true; }
    
    if(arma_config::check_nonfinite && X.has_nonfinite())  { return false; }
    
    vals.set_size(X.n_rows, 1);
    
    Mat<T> tmp(1, 1, arma_nozeros_indicator());
    
    if(vecs_on)
      {
      vecs.set_size(X.n_rows, X.n_rows);
       tmp.set_size(X.n_rows, X.n_rows);
      }
    
    podarray<T> junk(1);
    
    char     jobvl = 'N';
    char     jobvr = (vecs_on) ? 'V' : 'N';
    blas_int N     = blas_int(X.n_rows);
    T*       vl    = junk.memptr();
    T*       vr    = (vecs_on) ? tmp.memptr() : junk.memptr();
    blas_int ldvl  = blas_int(1);
    blas_int ldvr  = (vecs_on) ? blas_int(tmp.n_rows) : blas_int(1);
    blas_int lwork = 64*N;  // lwork_min = (vecs_on) ? (std::max)(blas_int(1), 4*N) : (std::max)(blas_int(1), 3*N)
    blas_int info  = 0;
    
    podarray<T> work( static_cast<uword>(lwork) );
    
    podarray<T> vals_real(X.n_rows);
    podarray<T> vals_imag(X.n_rows);
    
    arma_extra_debug_print("lapack::geev() -- START");
    lapack::geev(&jobvl, &jobvr, &N, X.memptr(), &N, vals_real.memptr(), vals_imag.memptr(), vl, &ldvl, vr, &ldvr, work.memptr(), &lwork, &info);
    arma_extra_debug_print("lapack::geev() -- END");
    
    if(info != 0)  { return false; }
    
    arma_extra_debug_print("reformatting eigenvalues and eigenvectors");
    
    std::complex<T>* vals_mem = vals.memptr();
    
    for(uword i=0; i < X.n_rows; ++i)  { vals_mem[i] = std::complex<T>(vals_real[i], vals_imag[i]); }
    
    if(vecs_on)
      {
      for(uword j=0; j < X.n_rows; ++j)
        {
        if( (j < (X.n_rows-1)) && (vals_mem[j] == std::conj(vals_mem[j+1])) )
          {
          for(uword i=0; i < X.n_rows; ++i)
            {
            vecs.at(i,j)   = std::complex<T>( tmp.at(i,j),  tmp.at(i,j+1) );
            vecs.at(i,j+1) = std::complex<T>( tmp.at(i,j), -tmp.at(i,j+1) );
            }
          
          ++j;
          }
        else
          {
          for(uword i=0; i<X.n_rows; ++i)
            {
            vecs.at(i,j) = std::complex<T>(tmp.at(i,j), T(0));
            }
          }
        }
      }
    
    return true;
    }
  #else
    {
    arma_ignore(vals);
    arma_ignore(vecs);
    arma_ignore(vecs_on);
    arma_ignore(expr);
    arma_stop_logic_error("eig_gen(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! eigen decomposition of general square matrix (complex)
template<typename T1>
inline
bool
auxlib::eig_gen
  (
         Mat< std::complex<typename T1::pod_type> >&     vals,
         Mat< std::complex<typename T1::pod_type> >&     vecs, 
  const bool                                             vecs_on,
  const Base< std::complex<typename T1::pod_type>, T1 >& expr
  )
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type     T;
    typedef typename std::complex<T> eT;
    
    Mat<eT> X = expr.get_ref();
    
    arma_debug_check( (X.is_square() == false), "eig_gen(): given matrix must be square sized" );
    
    arma_debug_assert_blas_size(X);
    
    if(X.is_empty())  { vals.reset(); vecs.reset(); return true; }
    
    if(arma_config::check_nonfinite && X.has_nonfinite())  { return false; }
    
    vals.set_size(X.n_rows, 1);
    
    if(vecs_on)  { vecs.set_size(X.n_rows, X.n_rows); }
    
    podarray<eT> junk(1);
    
    char     jobvl = 'N';
    char     jobvr = (vecs_on) ? 'V' : 'N';
    blas_int N     = blas_int(X.n_rows);
    eT*      vl    = junk.memptr();
    eT*      vr    = (vecs_on) ? vecs.memptr() : junk.memptr();
    blas_int ldvl  = blas_int(1);
    blas_int ldvr  = (vecs_on) ? blas_int(vecs.n_rows) : blas_int(1);
    blas_int lwork = 64*N;  // lwork_min = (std::max)(blas_int(1), 2*N)
    blas_int info  = 0;
    
    podarray<eT>  work( static_cast<uword>(lwork) );
    podarray< T> rwork( static_cast<uword>(2*N)   );
    
    arma_extra_debug_print("lapack::cx_geev() -- START");
    lapack::cx_geev(&jobvl, &jobvr, &N, X.memptr(), &N, vals.memptr(), vl, &ldvl, vr, &ldvr, work.memptr(), &lwork, rwork.memptr(), &info);
    arma_extra_debug_print("lapack::cx_geev() -- END");
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(vals);
    arma_ignore(vecs);
    arma_ignore(vecs_on);
    arma_ignore(expr);
    arma_stop_logic_error("eig_gen(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! eigen decomposition of general square matrix (real, balance given matrix)
template<typename T1>
inline
bool
auxlib::eig_gen_balance
  (
         Mat< std::complex<typename T1::pod_type> >& vals,
         Mat< std::complex<typename T1::pod_type> >& vecs,
  const bool                                         vecs_on,
  const Base<typename T1::pod_type,T1>&              expr
  )
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type T;
    
    Mat<T> X = expr.get_ref();
    
    arma_debug_check( (X.is_square() == false), "eig_gen(): given matrix must be square sized" );
    
    arma_debug_assert_blas_size(X);
    
    if(X.is_empty())  { vals.reset(); vecs.reset(); return true; }
    
    if(arma_config::check_nonfinite && X.has_nonfinite())  { return false; }
    
    vals.set_size(X.n_rows, 1);
    
    Mat<T> tmp(1, 1, arma_nozeros_indicator());
    
    if(vecs_on)
      {
      vecs.set_size(X.n_rows, X.n_rows);
       tmp.set_size(X.n_rows, X.n_rows);
      }
    
    podarray<T> junk(1);
    
    char     bal   = 'B'; 
    char     jobvl = 'N';
    char     jobvr = (vecs_on) ? 'V' : 'N';
    char     sense = 'N';
    blas_int N     = blas_int(X.n_rows);
    T*       vl    = junk.memptr();
    T*       vr    = (vecs_on) ? tmp.memptr() : junk.memptr();
    blas_int ldvl  = blas_int(1);
    blas_int ldvr  = (vecs_on) ? blas_int(tmp.n_rows) : blas_int(1);
    blas_int ilo   = blas_int(0);
    blas_int ihi   = blas_int(0);
    T        abnrm = T(0);
    blas_int lwork = 64*N;  // lwork_min = (vecs_on) ? (std::max)(blas_int(1), 2*N) : (std::max)(blas_int(1), 3*N)
    blas_int info  = blas_int(0);
    
    podarray<T>  scale(X.n_rows);
    podarray<T> rconde(X.n_rows);
    podarray<T> rcondv(X.n_rows);
    
    podarray<T>         work( static_cast<uword>(lwork) );
    podarray<blas_int> iwork( uword(1) );  // iwork not used by lapack::geevx() as sense = 'N'
    
    podarray<T> vals_real(X.n_rows);
    podarray<T> vals_imag(X.n_rows);
    
    arma_extra_debug_print("lapack::geevx() -- START");
    lapack::geevx(&bal, &jobvl, &jobvr, &sense, &N, X.memptr(), &N, vals_real.memptr(), vals_imag.memptr(), vl, &ldvl, vr, &ldvr, &ilo, &ihi, scale.memptr(), &abnrm, rconde.memptr(), rcondv.memptr(), work.memptr(), &lwork, iwork.memptr(), &info);
    arma_extra_debug_print("lapack::geevx() -- END");
    
    if(info != 0)  { return false; }
    
    arma_extra_debug_print("reformatting eigenvalues and eigenvectors");
    
    std::complex<T>* vals_mem = vals.memptr();
    
    for(uword i=0; i < X.n_rows; ++i)  { vals_mem[i] = std::complex<T>(vals_real[i], vals_imag[i]); }
    
    if(vecs_on)
      {
      for(uword j=0; j < X.n_rows; ++j)
        {
        if( (j < (X.n_rows-1)) && (vals_mem[j] == std::conj(vals_mem[j+1])) )
          {
          for(uword i=0; i < X.n_rows; ++i)
            {
            vecs.at(i,j)   = std::complex<T>( tmp.at(i,j),  tmp.at(i,j+1) );
            vecs.at(i,j+1) = std::complex<T>( tmp.at(i,j), -tmp.at(i,j+1) );
            }
          
          ++j;
          }
        else
          {
          for(uword i=0; i<X.n_rows; ++i)
            {
            vecs.at(i,j) = std::complex<T>(tmp.at(i,j), T(0));
            }
          }
        }
      }
    
    return true;
    }
  #else
    {
    arma_ignore(vals);
    arma_ignore(vecs);
    arma_ignore(vecs_on);
    arma_ignore(expr);
    arma_stop_logic_error("eig_gen(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! eigen decomposition of general square matrix (complex, balance given matrix)
template<typename T1>
inline
bool
auxlib::eig_gen_balance
  (
         Mat< std::complex<typename T1::pod_type> >&     vals,
         Mat< std::complex<typename T1::pod_type> >&     vecs, 
  const bool                                             vecs_on,
  const Base< std::complex<typename T1::pod_type>, T1 >& expr
  )
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_extra_debug_print("auxlib::eig_gen_balance(): redirecting to auxlib::eig_gen() due to crippled LAPACK");
    
    return auxlib::eig_gen(vals, vecs, vecs_on, expr);
    }
  #elif defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type     T;
    typedef typename std::complex<T> eT;
    
    Mat<eT> X = expr.get_ref();
    
    arma_debug_check( (X.is_square() == false), "eig_gen(): given matrix must be square sized" );
    
    arma_debug_assert_blas_size(X);
    
    if(X.is_empty())  { vals.reset(); vecs.reset(); return true; }
    
    if(arma_config::check_nonfinite && X.has_nonfinite())  { return false; }
    
    vals.set_size(X.n_rows, 1);
    
    if(vecs_on)  { vecs.set_size(X.n_rows, X.n_rows); }
    
    podarray<eT> junk(1);
    
    char     bal   = 'B';
    char     jobvl = 'N';
    char     jobvr = (vecs_on) ? 'V' : 'N';
    char     sense = 'N';
    blas_int N     = blas_int(X.n_rows);
    eT*      vl    = junk.memptr();
    eT*      vr    = (vecs_on) ? vecs.memptr() : junk.memptr();
    blas_int ldvl  = blas_int(1);
    blas_int ldvr  = (vecs_on) ? blas_int(vecs.n_rows) : blas_int(1);
    blas_int ilo   = blas_int(0);
    blas_int ihi   = blas_int(0);
    T        abnrm = T(0);
    blas_int lwork = 64*N;  // lwork_min = (std::max)(blas_int(1), blas_int(2*N))
    blas_int info  = blas_int(0);
    
    podarray<T>  scale(X.n_rows);
    podarray<T> rconde(X.n_rows);
    podarray<T> rcondv(X.n_rows);
    
    podarray<eT>  work( static_cast<uword>(lwork) );
    podarray< T> rwork( static_cast<uword>(2*N)   );
    
    arma_extra_debug_print("lapack::cx_geevx() -- START");
    lapack::cx_geevx(&bal, &jobvl, &jobvr, &sense, &N, X.memptr(), &N, vals.memptr(), vl, &ldvl, vr, &ldvr, &ilo, &ihi, scale.memptr(), &abnrm, rconde.memptr(), rcondv.memptr(), work.memptr(), &lwork, rwork.memptr(), &info);
    arma_extra_debug_print("lapack::cx_geevx() -- END");
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(vals);
    arma_ignore(vecs);
    arma_ignore(vecs_on);
    arma_ignore(expr);
    arma_stop_logic_error("eig_gen(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! two-sided eigen decomposition of general square matrix (real)
template<typename T1>
inline
bool
auxlib::eig_gen_twosided
  (
         Mat< std::complex<typename T1::pod_type> >&  vals,
         Mat< std::complex<typename T1::pod_type> >& lvecs,
         Mat< std::complex<typename T1::pod_type> >& rvecs,
  const Base<typename T1::pod_type,T1>&              expr
  )
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type T;
    
    Mat<T> X = expr.get_ref();
    
    arma_debug_check( (X.is_square() == false), "eig_gen(): given matrix must be square sized" );
    
    arma_debug_assert_blas_size(X);
    
    if(X.is_empty())  { vals.reset(); lvecs.reset(); rvecs.reset(); return true; }
    
    if(arma_config::check_nonfinite && X.has_nonfinite())  { return false; }
    
    vals.set_size(X.n_rows, 1);
    
    lvecs.set_size(X.n_rows, X.n_rows);
    rvecs.set_size(X.n_rows, X.n_rows);
    
    Mat<T> ltmp(X.n_rows, X.n_rows, arma_nozeros_indicator());
    Mat<T> rtmp(X.n_rows, X.n_rows, arma_nozeros_indicator());
    
    char     jobvl = 'V';
    char     jobvr = 'V';
    blas_int N     = blas_int(X.n_rows);
    blas_int ldvl  = blas_int(ltmp.n_rows);
    blas_int ldvr  = blas_int(rtmp.n_rows);
    blas_int lwork = 64*N;  // lwork_min = (std::max)(blas_int(1), 4*N)
    blas_int info  = 0;
    
    podarray<T> work( static_cast<uword>(lwork) );
    
    podarray<T> vals_real(X.n_rows);
    podarray<T> vals_imag(X.n_rows);
    
    arma_extra_debug_print("lapack::geev() -- START");
    lapack::geev(&jobvl, &jobvr, &N, X.memptr(), &N, vals_real.memptr(), vals_imag.memptr(), ltmp.memptr(), &ldvl, rtmp.memptr(), &ldvr, work.memptr(), &lwork, &info);
    arma_extra_debug_print("lapack::geev() -- END");
    
    if(info != 0)  { return false; }
    
    arma_extra_debug_print("reformatting eigenvalues and eigenvectors");
    
    std::complex<T>* vals_mem = vals.memptr();
    
    for(uword i=0; i < X.n_rows; ++i)  { vals_mem[i] = std::complex<T>(vals_real[i], vals_imag[i]); }
    
    for(uword j=0; j < X.n_rows; ++j)
      {
      if( (j < (X.n_rows-1)) && (vals_mem[j] == std::conj(vals_mem[j+1])) )
        {
        for(uword i=0; i < X.n_rows; ++i)
          {
          lvecs.at(i,j)   = std::complex<T>( ltmp.at(i,j),  ltmp.at(i,j+1) );
          lvecs.at(i,j+1) = std::complex<T>( ltmp.at(i,j), -ltmp.at(i,j+1) );
          rvecs.at(i,j)   = std::complex<T>( rtmp.at(i,j),  rtmp.at(i,j+1) );
          rvecs.at(i,j+1) = std::complex<T>( rtmp.at(i,j), -rtmp.at(i,j+1) );
          }
        ++j;
        }
      else
        {
        for(uword i=0; i<X.n_rows; ++i)
          {
          lvecs.at(i,j) = std::complex<T>(ltmp.at(i,j), T(0));
          rvecs.at(i,j) = std::complex<T>(rtmp.at(i,j), T(0));
          }
        }
      }
    
    return true;
    }
  #else
    {
    arma_ignore(vals);
    arma_ignore(lvecs);
    arma_ignore(rvecs);
    arma_ignore(expr);
    arma_stop_logic_error("eig_gen(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! two-sided eigen decomposition of general square matrix (complex)
template<typename T1>
inline
bool
auxlib::eig_gen_twosided
  (
         Mat< std::complex<typename T1::pod_type> >&      vals,
         Mat< std::complex<typename T1::pod_type> >&     lvecs,
         Mat< std::complex<typename T1::pod_type> >&     rvecs,
  const Base< std::complex<typename T1::pod_type>, T1 >& expr
  )
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type     T;
    typedef typename std::complex<T> eT;
    
    Mat<eT> X = expr.get_ref();
    
    arma_debug_check( (X.is_square() == false), "eig_gen(): given matrix must be square sized" );
    
    arma_debug_assert_blas_size(X);
    
    if(X.is_empty())  { vals.reset(); lvecs.reset(); rvecs.reset(); return true; }
    
    if(arma_config::check_nonfinite && X.has_nonfinite())  { return false; }
    
    vals.set_size(X.n_rows, 1);
    
    lvecs.set_size(X.n_rows, X.n_rows);
    rvecs.set_size(X.n_rows, X.n_rows);
    
    char     jobvl = 'V';
    char     jobvr = 'V';
    blas_int N     = blas_int(X.n_rows);
    blas_int ldvl  = blas_int(lvecs.n_rows);
    blas_int ldvr  = blas_int(rvecs.n_rows);
    blas_int lwork = 64*N;  // lwork_min = (std::max)(blas_int(1), 2*N)
    blas_int info  = 0;
    
    podarray<eT>  work( static_cast<uword>(lwork) );
    podarray< T> rwork( static_cast<uword>(2*N)   );
    
    arma_extra_debug_print("lapack::cx_geev() -- START");
    lapack::cx_geev(&jobvl, &jobvr, &N, X.memptr(), &N, vals.memptr(), lvecs.memptr(), &ldvl, rvecs.memptr(), &ldvr, work.memptr(), &lwork, rwork.memptr(), &info);
    arma_extra_debug_print("lapack::cx_geev() -- END");
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(vals);
    arma_ignore(lvecs);
    arma_ignore(rvecs);
    arma_ignore(expr);
    arma_stop_logic_error("eig_gen(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! two-sided eigen decomposition of general square matrix (real, balance given matrix)
template<typename T1>
inline
bool
auxlib::eig_gen_twosided_balance
  (
         Mat< std::complex<typename T1::pod_type> >&  vals,
         Mat< std::complex<typename T1::pod_type> >& lvecs,
         Mat< std::complex<typename T1::pod_type> >& rvecs,
  const Base<typename T1::pod_type,T1>&              expr
  )
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type T;
    
    Mat<T> X = expr.get_ref();
    
    arma_debug_check( (X.is_square() == false), "eig_gen(): given matrix must be square sized" );
    
    arma_debug_assert_blas_size(X);
    
    if(X.is_empty())  { vals.reset(); lvecs.reset(); rvecs.reset(); return true; }
    
    if(arma_config::check_nonfinite && X.has_nonfinite())  { return false; }
    
    vals.set_size(X.n_rows, 1);
    
    lvecs.set_size(X.n_rows, X.n_rows);
    rvecs.set_size(X.n_rows, X.n_rows);
    
    Mat<T> ltmp(X.n_rows, X.n_rows, arma_nozeros_indicator());
    Mat<T> rtmp(X.n_rows, X.n_rows, arma_nozeros_indicator());
    
    char     bal   = 'B';
    char     jobvl = 'V';
    char     jobvr = 'V';
    char     sense = 'N';
    blas_int N     = blas_int(X.n_rows);
    blas_int ldvl  = blas_int(ltmp.n_rows);
    blas_int ldvr  = blas_int(rtmp.n_rows);
    blas_int ilo   = blas_int(0);
    blas_int ihi   = blas_int(0);
    T        abnrm = T(0);
    blas_int lwork = 64*N;  // lwork_min = (std::max)(blas_int(1), blas_int(3*N))
    blas_int info  = blas_int(0);
    
    podarray<T>  scale(X.n_rows);
    podarray<T> rconde(X.n_rows);
    podarray<T> rcondv(X.n_rows);
    
    podarray<T>         work( static_cast<uword>(lwork) );
    podarray<blas_int> iwork( uword(1) );  // iwork not used by lapack::geevx() as sense = 'N'
    
    podarray<T> vals_real(X.n_rows);
    podarray<T> vals_imag(X.n_rows);
    
    arma_extra_debug_print("lapack::geevx() -- START");
    lapack::geevx(&bal, &jobvl, &jobvr, &sense, &N, X.memptr(), &N, vals_real.memptr(), vals_imag.memptr(), ltmp.memptr(), &ldvl, rtmp.memptr(), &ldvr, &ilo, &ihi, scale.memptr(), &abnrm, rconde.memptr(), rcondv.memptr(), work.memptr(), &lwork, iwork.memptr(), &info);
    arma_extra_debug_print("lapack::geevx() -- END");
    
    if(info != 0)  { return false; }
    
    arma_extra_debug_print("reformatting eigenvalues and eigenvectors");
    
    std::complex<T>* vals_mem = vals.memptr();
    
    for(uword i=0; i < X.n_rows; ++i)  { vals_mem[i] = std::complex<T>(vals_real[i], vals_imag[i]); }
    
    for(uword j=0; j < X.n_rows; ++j)
      {
      if( (j < (X.n_rows-1)) && (vals_mem[j] == std::conj(vals_mem[j+1])) )
        {
        for(uword i=0; i < X.n_rows; ++i)
          {
          lvecs.at(i,j)   = std::complex<T>( ltmp.at(i,j),  ltmp.at(i,j+1) );
          lvecs.at(i,j+1) = std::complex<T>( ltmp.at(i,j), -ltmp.at(i,j+1) );
          rvecs.at(i,j)   = std::complex<T>( rtmp.at(i,j),  rtmp.at(i,j+1) );
          rvecs.at(i,j+1) = std::complex<T>( rtmp.at(i,j), -rtmp.at(i,j+1) );
          }
        ++j;
        }
      else
        {
        for(uword i=0; i<X.n_rows; ++i)
          {
          lvecs.at(i,j) = std::complex<T>(ltmp.at(i,j), T(0));
          rvecs.at(i,j) = std::complex<T>(rtmp.at(i,j), T(0));
          }
        }
      }
    
    return true;
    }
  #else
    {
    arma_ignore(vals);
    arma_ignore(lvecs);
    arma_ignore(rvecs);
    arma_ignore(expr);
    arma_stop_logic_error("eig_gen(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! two-sided eigen decomposition of general square matrix (complex, balance given matrix)
template<typename T1>
inline
bool
auxlib::eig_gen_twosided_balance
  (
         Mat< std::complex<typename T1::pod_type> >&      vals,
         Mat< std::complex<typename T1::pod_type> >&     lvecs,
         Mat< std::complex<typename T1::pod_type> >&     rvecs,
  const Base< std::complex<typename T1::pod_type>, T1 >& expr
  )
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_extra_debug_print("auxlib::eig_gen_twosided_balance(): redirecting to auxlib::eig_gen() due to crippled LAPACK");
    
    return auxlib::eig_gen(vals, lvecs, rvecs, expr);
    }
  #elif defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type     T;
    typedef typename std::complex<T> eT;
    
    Mat<eT> X = expr.get_ref();
    
    arma_debug_check( (X.is_square() == false), "eig_gen(): given matrix must be square sized" );
    
    arma_debug_assert_blas_size(X);
    
    if(X.is_empty())  { vals.reset(); lvecs.reset(); rvecs.reset(); return true; }
    
    if(arma_config::check_nonfinite && X.has_nonfinite())  { return false; }
    
    vals.set_size(X.n_rows, 1);
    
    lvecs.set_size(X.n_rows, X.n_rows);
    rvecs.set_size(X.n_rows, X.n_rows);
    
    char     bal   = 'B';
    char     jobvl = 'V';
    char     jobvr = 'V';
    char     sense = 'N';
    blas_int N     = blas_int(X.n_rows);
    blas_int ldvl  = blas_int(lvecs.n_rows);
    blas_int ldvr  = blas_int(rvecs.n_rows);
    blas_int ilo   = blas_int(0);
    blas_int ihi   = blas_int(0);
    T        abnrm = T(0);
    blas_int lwork = 64*N;  // lwork_min = (std::max)(blas_int(1), blas_int(2*N))
    blas_int info  = blas_int(0);
    
    podarray<T>  scale(X.n_rows);
    podarray<T> rconde(X.n_rows);
    podarray<T> rcondv(X.n_rows);
    
    podarray<eT>  work( static_cast<uword>(lwork) );
    podarray< T> rwork( static_cast<uword>(2*N)   );
    
    arma_extra_debug_print("lapack::cx_geevx() -- START");
    lapack::cx_geevx(&bal, &jobvl, &jobvr, &sense, &N, X.memptr(), &N, vals.memptr(), lvecs.memptr(), &ldvl, rvecs.memptr(), &ldvr, &ilo, &ihi, scale.memptr(), &abnrm, rconde.memptr(), rcondv.memptr(), work.memptr(), &lwork, rwork.memptr(), &info);
    arma_extra_debug_print("lapack::cx_geevx() -- END");
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(vals);
    arma_ignore(lvecs);
    arma_ignore(rvecs);
    arma_ignore(expr);
    arma_stop_logic_error("eig_gen(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! eigendecomposition of general square matrix pair (real)
template<typename T1, typename T2>
inline
bool
auxlib::eig_pair
  (
        Mat< std::complex<typename T1::pod_type> >& vals,
        Mat< std::complex<typename T1::pod_type> >& vecs,
  const bool                                        vecs_on,
  const Base<typename T1::pod_type,T1>&             A_expr,
  const Base<typename T1::pod_type,T2>&             B_expr
  )
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type  T;
    typedef std::complex<T>       eT;
    
    Mat<T> A(A_expr.get_ref());
    Mat<T> B(B_expr.get_ref());
    
    arma_debug_check( ((A.is_square() == false) || (B.is_square() == false)), "eig_pair(): given matrices must be square sized" );
    
    arma_debug_check( (A.n_rows != B.n_rows), "eig_pair(): given matrices must have the same size" );
    
    arma_debug_assert_blas_size(A);
    
    if(A.is_empty())  { vals.reset(); vecs.reset(); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    if(arma_config::check_nonfinite && B.has_nonfinite())  { return false; }
    
    vals.set_size(A.n_rows, 1);
    
    Mat<T> tmp(1, 1, arma_nozeros_indicator());
    
    if(vecs_on)
      {
      vecs.set_size(A.n_rows, A.n_rows);
       tmp.set_size(A.n_rows, A.n_rows);
      }
    
    podarray<T> junk(1);
    
    char     jobvl = 'N';
    char     jobvr = (vecs_on) ? 'V' : 'N';
    blas_int N     = blas_int(A.n_rows);
    T*       vl    = junk.memptr();
    T*       vr    = (vecs_on) ? tmp.memptr() : junk.memptr();
    blas_int ldvl  = blas_int(1);
    blas_int ldvr  = (vecs_on) ? blas_int(tmp.n_rows) : blas_int(1);
    blas_int lwork = 64*N;  // lwork_min = (std::max)(blas_int(1), 8*N)
    blas_int info  = 0;
    
    podarray<T> alphar(A.n_rows);
    podarray<T> alphai(A.n_rows);
    podarray<T>   beta(A.n_rows);
    
    podarray<T> work( static_cast<uword>(lwork) );
    
    arma_extra_debug_print("lapack::ggev()");
    lapack::ggev(&jobvl, &jobvr, &N, A.memptr(), &N,  B.memptr(), &N, alphar.memptr(), alphai.memptr(), beta.memptr(), vl, &ldvl, vr, &ldvr, work.memptr(), &lwork, &info);
    
    if(info != 0)  { return false; }
    
    arma_extra_debug_print("reformatting eigenvalues and eigenvectors");
    
          eT*   vals_mem =   vals.memptr();
    const  T* alphar_mem = alphar.memptr();
    const  T* alphai_mem = alphai.memptr();
    const  T*   beta_mem =   beta.memptr();
    
    bool beta_has_zero = false;
    
    for(uword j=0; j<A.n_rows; ++j)
      {
      const T alphai_val = alphai_mem[j];
      const T   beta_val =   beta_mem[j];
      
      const T re = alphar_mem[j] / beta_val;
      const T im = alphai_val    / beta_val;
      
      beta_has_zero = (beta_has_zero || (beta_val == T(0)));
      
      vals_mem[j] = std::complex<T>(re, im);
      
      if( (alphai_val > T(0)) && (j < (A.n_rows-1)) )
        {
        ++j;
        vals_mem[j] = std::complex<T>(re,-im);  // force exact conjugate
        }
      }
    
    if(beta_has_zero)  { arma_debug_warn_level(1, "eig_pair(): given matrices appear ill-conditioned"); }
    
    if(vecs_on)
      {
      for(uword j=0; j<A.n_rows; ++j)
        {
        if( (j < (A.n_rows-1)) && (vals_mem[j] == std::conj(vals_mem[j+1])) )
          {
          for(uword i=0; i<A.n_rows; ++i)
            {
            vecs.at(i,j)   = std::complex<T>( tmp.at(i,j),  tmp.at(i,j+1) );
            vecs.at(i,j+1) = std::complex<T>( tmp.at(i,j), -tmp.at(i,j+1) );
            }
          
          ++j;
          }
        else
          {
          for(uword i=0; i<A.n_rows; ++i)
            {
            vecs.at(i,j) = std::complex<T>(tmp.at(i,j), T(0));
            }
          }
        }
      }
    
    return true;
    }
  #else
    {
    arma_ignore(vals);
    arma_ignore(vecs);
    arma_ignore(vecs_on);
    arma_ignore(A_expr);
    arma_ignore(B_expr);
    arma_stop_logic_error("eig_pair(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! eigendecomposition of general square matrix pair (complex)
template<typename T1, typename T2>
inline
bool
auxlib::eig_pair
  (
        Mat< std::complex<typename T1::pod_type> >&      vals,
        Mat< std::complex<typename T1::pod_type> >&      vecs,
  const bool                                             vecs_on,
  const Base< std::complex<typename T1::pod_type>, T1 >& A_expr,
  const Base< std::complex<typename T1::pod_type>, T2 >& B_expr
  )
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type     T;
    typedef typename std::complex<T> eT;
    
    Mat<eT> A(A_expr.get_ref());
    Mat<eT> B(B_expr.get_ref());
    
    arma_debug_check( ((A.is_square() == false) || (B.is_square() == false)), "eig_pair(): given matrices must be square sized" );
    
    arma_debug_check( (A.n_rows != B.n_rows), "eig_pair(): given matrices must have the same size" );
    
    arma_debug_assert_blas_size(A);
    
    if(A.is_empty())  { vals.reset(); vecs.reset(); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    if(arma_config::check_nonfinite && B.has_nonfinite())  { return false; }
    
    vals.set_size(A.n_rows, 1);
    
    if(vecs_on)  { vecs.set_size(A.n_rows, A.n_rows); }
    
    podarray<eT> junk(1);
    
    char     jobvl = 'N';
    char     jobvr = (vecs_on) ? 'V' : 'N';
    blas_int N     = blas_int(A.n_rows);
    eT*      vl    = junk.memptr();
    eT*      vr    = (vecs_on) ? vecs.memptr() : junk.memptr();
    blas_int ldvl  = blas_int(1);
    blas_int ldvr  = (vecs_on) ? blas_int(vecs.n_rows) : blas_int(1);
    blas_int lwork = 64*N;  // lwork_min = (std::max)(blas_int(1),2*N)
    blas_int info  = 0;
    
    podarray<eT> alpha(A.n_rows);
    podarray<eT>  beta(A.n_rows);
    
    podarray<eT>  work( static_cast<uword>(lwork) );
    podarray<T>  rwork( static_cast<uword>(8*N)   );
    
    arma_extra_debug_print("lapack::cx_ggev()");
    lapack::cx_ggev(&jobvl, &jobvr, &N, A.memptr(), &N, B.memptr(), &N, alpha.memptr(), beta.memptr(), vl, &ldvl, vr, &ldvr, work.memptr(), &lwork, rwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
          eT*   vals_mem =  vals.memptr();
    const eT*  alpha_mem = alpha.memptr();
    const eT*   beta_mem =  beta.memptr();
    
    const std::complex<T> zero(T(0), T(0));
    
    bool beta_has_zero = false;
    
    for(uword i=0; i<A.n_rows; ++i)
      {
      const eT& beta_val = beta_mem[i];
      
      vals_mem[i] = alpha_mem[i] / beta_val;
      
      beta_has_zero = (beta_has_zero || (beta_val == zero));
      }
    
    if(beta_has_zero)  { arma_debug_warn_level(1, "eig_pair(): given matrices appear ill-conditioned"); }
    
    return true;
    }
  #else
    {
    arma_ignore(vals);
    arma_ignore(vecs);
    arma_ignore(vecs_on);
    arma_ignore(A_expr);
    arma_ignore(B_expr);
    arma_stop_logic_error("eig_pair(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! two-sided eigendecomposition of general square matrix pair (real)
template<typename T1, typename T2>
inline
bool
auxlib::eig_pair_twosided
  (
        Mat< std::complex<typename T1::pod_type> >&  vals,
        Mat< std::complex<typename T1::pod_type> >& lvecs,
        Mat< std::complex<typename T1::pod_type> >& rvecs,
  const Base<typename T1::pod_type,T1>&             A_expr,
  const Base<typename T1::pod_type,T2>&             B_expr
  )
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type  T;
    typedef std::complex<T>       eT;
    
    Mat<T> A(A_expr.get_ref());
    Mat<T> B(B_expr.get_ref());
    
    arma_debug_check( ((A.is_square() == false) || (B.is_square() == false)), "eig_pair(): given matrices must be square sized" );
    
    arma_debug_check( (A.n_rows != B.n_rows), "eig_pair(): given matrices must have the same size" );
    
    arma_debug_assert_blas_size(A);
    
    if(A.is_empty())  { vals.reset(); lvecs.reset(); rvecs.reset(); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    if(arma_config::check_nonfinite && B.has_nonfinite())  { return false; }
    
    vals.set_size(A.n_rows, 1);
    
    lvecs.set_size(A.n_rows, A.n_rows);
    rvecs.set_size(A.n_rows, A.n_rows);
    
    Mat<T> ltmp(A.n_rows, A.n_rows, arma_nozeros_indicator());
    Mat<T> rtmp(A.n_rows, A.n_rows, arma_nozeros_indicator());
    
    char     jobvl = 'V';
    char     jobvr = 'V';
    blas_int N     = blas_int(A.n_rows);
    blas_int ldvl  = blas_int(ltmp.n_rows);
    blas_int ldvr  = blas_int(rtmp.n_rows);
    blas_int lwork = 64*N;  // lwork_min = (std::max)(blas_int(1), 8*N)
    blas_int info  = 0;
    
    podarray<T> alphar(A.n_rows);
    podarray<T> alphai(A.n_rows);
    podarray<T>   beta(A.n_rows);
    
    podarray<T> work( static_cast<uword>(lwork) );
    
    arma_extra_debug_print("lapack::ggev()");
    lapack::ggev(&jobvl, &jobvr, &N, A.memptr(), &N,  B.memptr(), &N, alphar.memptr(), alphai.memptr(), beta.memptr(), ltmp.memptr(), &ldvl, rtmp.memptr(), &ldvr, work.memptr(), &lwork, &info);
    
    if(info != 0)  { return false; }
    
    arma_extra_debug_print("reformatting eigenvalues and eigenvectors");
    
          eT*   vals_mem =   vals.memptr();
    const  T* alphar_mem = alphar.memptr();
    const  T* alphai_mem = alphai.memptr();
    const  T*   beta_mem =   beta.memptr();
    
    bool beta_has_zero = false;
    
    for(uword j=0; j<A.n_rows; ++j)
      {
      const T alphai_val = alphai_mem[j];
      const T   beta_val =   beta_mem[j];
      
      const T re = alphar_mem[j] / beta_val;
      const T im = alphai_val    / beta_val;
      
      beta_has_zero = (beta_has_zero || (beta_val == T(0)));
      
      vals_mem[j] = std::complex<T>(re, im);
      
      if( (alphai_val > T(0)) && (j < (A.n_rows-1)) )
        {
        ++j;
        vals_mem[j] = std::complex<T>(re,-im);  // force exact conjugate
        }
      }
    
    if(beta_has_zero)  { arma_debug_warn_level(1, "eig_pair(): given matrices appear ill-conditioned"); }
    
    for(uword j=0; j < A.n_rows; ++j)
      {
      if( (j < (A.n_rows-1)) && (vals_mem[j] == std::conj(vals_mem[j+1])) )
        {
        for(uword i=0; i < A.n_rows; ++i)
          {
          lvecs.at(i,j)   = std::complex<T>( ltmp.at(i,j),  ltmp.at(i,j+1) );
          lvecs.at(i,j+1) = std::complex<T>( ltmp.at(i,j), -ltmp.at(i,j+1) );
          rvecs.at(i,j)   = std::complex<T>( rtmp.at(i,j),  rtmp.at(i,j+1) );
          rvecs.at(i,j+1) = std::complex<T>( rtmp.at(i,j), -rtmp.at(i,j+1) );
          }
        ++j;
        }
      else
        {
        for(uword i=0; i<A.n_rows; ++i)
          {
          lvecs.at(i,j) = std::complex<T>(ltmp.at(i,j), T(0));
          rvecs.at(i,j) = std::complex<T>(rtmp.at(i,j), T(0));
          }
        }
      }
    
    return true;
    }
  #else
    {
    arma_ignore(vals);
    arma_ignore(lvecs);
    arma_ignore(rvecs);
    arma_ignore(A_expr);
    arma_ignore(B_expr);
    arma_stop_logic_error("eig_pair(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! two-sided eigendecomposition of general square matrix pair (complex)
template<typename T1, typename T2>
inline
bool
auxlib::eig_pair_twosided
  (
        Mat< std::complex<typename T1::pod_type> >&       vals,
        Mat< std::complex<typename T1::pod_type> >&      lvecs,
        Mat< std::complex<typename T1::pod_type> >&      rvecs,
  const Base< std::complex<typename T1::pod_type>, T1 >& A_expr,
  const Base< std::complex<typename T1::pod_type>, T2 >& B_expr
  )
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type     T;
    typedef typename std::complex<T> eT;
    
    Mat<eT> A(A_expr.get_ref());
    Mat<eT> B(B_expr.get_ref());
    
    arma_debug_check( ((A.is_square() == false) || (B.is_square() == false)), "eig_pair(): given matrices must be square sized" );
    
    arma_debug_check( (A.n_rows != B.n_rows), "eig_pair(): given matrices must have the same size" );
    
    arma_debug_assert_blas_size(A);
    
    if(A.is_empty())  { vals.reset(); lvecs.reset(); rvecs.reset(); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    if(arma_config::check_nonfinite && B.has_nonfinite())  { return false; }
    
    vals.set_size(A.n_rows, 1);
    
    lvecs.set_size(A.n_rows, A.n_rows);
    rvecs.set_size(A.n_rows, A.n_rows);
    
    char     jobvl = 'V';
    char     jobvr = 'V';
    blas_int N     = blas_int(A.n_rows);
    blas_int ldvl  = blas_int(lvecs.n_rows);
    blas_int ldvr  = blas_int(rvecs.n_rows);
    blas_int lwork = 64*N;  // lwork_min = (std::max)(blas_int(1),2*N)
    blas_int info  = 0;
    
    podarray<eT> alpha(A.n_rows);
    podarray<eT>  beta(A.n_rows);
    
    podarray<eT>  work( static_cast<uword>(lwork) );
    podarray<T>  rwork( static_cast<uword>(8*N)   );
    
    arma_extra_debug_print("lapack::cx_ggev()");
    lapack::cx_ggev(&jobvl, &jobvr, &N, A.memptr(), &N, B.memptr(), &N, alpha.memptr(), beta.memptr(), lvecs.memptr(), &ldvl, rvecs.memptr(), &ldvr, work.memptr(), &lwork, rwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
          eT*   vals_mem =  vals.memptr();
    const eT*  alpha_mem = alpha.memptr();
    const eT*   beta_mem =  beta.memptr();
    
    const std::complex<T> zero(T(0), T(0));
    
    bool beta_has_zero = false;
    
    for(uword i=0; i<A.n_rows; ++i)
      {
      const eT& beta_val = beta_mem[i];
      
      vals_mem[i] = alpha_mem[i] / beta_val;
      
      beta_has_zero = (beta_has_zero || (beta_val == zero));
      }
    
    if(beta_has_zero)  { arma_debug_warn_level(1, "eig_pair(): given matrices appear ill-conditioned"); }
    
    return true;
    }
  #else
    {
    arma_ignore(vals);
    arma_ignore(lvecs);
    arma_ignore(rvecs);
    arma_ignore(A_expr);
    arma_ignore(B_expr);
    arma_stop_logic_error("eig_pair(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! eigenvalues of a symmetric real matrix
template<typename eT>
inline
bool
auxlib::eig_sym(Col<eT>& eigval, Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_check( (A.is_square() == false), "eig_sym(): given matrix must be square sized" );
    
    if(A.is_empty())  { eigval.reset(); return true; }
    
    if((arma_config::debug) && (auxlib::rudimentary_sym_check(A) == false))
      {
      arma_debug_warn_level(1, "eig_sym(): given matrix is not symmetric");
      }
    
    if(arma_config::check_nonfinite && trimat_helper::has_nonfinite_triu(A))  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    eigval.set_size(A.n_rows);
    
    char jobz  = 'N';
    char uplo  = 'U';
    
    blas_int N     = blas_int(A.n_rows);
    blas_int lwork = (64+2)*N;  // lwork_min = (std::max)(blas_int(1), 3*N-1)
    blas_int info  = 0;
    
    podarray<eT> work( static_cast<uword>(lwork) );
    
    arma_extra_debug_print("lapack::syev()");
    lapack::syev(&jobz, &uplo, &N, A.memptr(), &N, eigval.memptr(), work.memptr(), &lwork, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(A);
    arma_stop_logic_error("eig_sym(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! eigenvalues of a hermitian complex matrix
template<typename T>
inline
bool
auxlib::eig_sym(Col<T>& eigval, Mat< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename std::complex<T> eT;
    
    arma_debug_check( (A.is_square() == false), "eig_sym(): given matrix must be square sized" );
    
    if(A.is_empty())  { eigval.reset(); return true; }
    
    if((arma_config::debug) && (auxlib::rudimentary_sym_check(A) == false))
      {
      arma_debug_warn_level(1, "eig_sym(): given matrix is not hermitian");
      }
    
    if(arma_config::check_nonfinite && trimat_helper::has_nonfinite_triu(A))  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    eigval.set_size(A.n_rows);
    
    char jobz  = 'N'; 
    char uplo  = 'U';
    
    blas_int N     = blas_int(A.n_rows);
    blas_int lwork = (64+1)*N;  // lwork_min = (std::max)(blas_int(1), 2*N-1)
    blas_int info  = 0;
    
    podarray<eT>  work( static_cast<uword>(lwork) );
    podarray<T>  rwork( static_cast<uword>( (std::max)(blas_int(1), 3*N) ) );
    
    arma_extra_debug_print("lapack::heev()");
    lapack::heev(&jobz, &uplo, &N, A.memptr(), &N, eigval.memptr(), work.memptr(), &lwork, rwork.memptr(), &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(A);
    arma_stop_logic_error("eig_sym(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! eigenvalues and eigenvectors of a symmetric real matrix
template<typename eT>
inline
bool
auxlib::eig_sym(Col<eT>& eigval, Mat<eT>& eigvec, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_check( (X.is_square() == false), "eig_sym(): given matrix must be square sized" );
    
    if(arma_config::check_nonfinite && trimat_helper::has_nonfinite_triu(X))  { return false; }
    
    eigvec = X;
    
    if(eigvec.is_empty())  { eigval.reset(); eigvec.reset(); return true; }
    
    arma_debug_assert_blas_size(eigvec);
    
    eigval.set_size(eigvec.n_rows);
    
    char jobz  = 'V';
    char uplo  = 'U';
    
    blas_int N     = blas_int(eigvec.n_rows);
    blas_int lwork = (64+2)*N;  // lwork_min = (std::max)(blas_int(1), 3*N-1)
    blas_int info  = 0;
    
    podarray<eT> work( static_cast<uword>(lwork) );
    
    arma_extra_debug_print("lapack::syev()");
    lapack::syev(&jobz, &uplo, &N, eigvec.memptr(), &N, eigval.memptr(), work.memptr(), &lwork, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_stop_logic_error("eig_sym(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! eigenvalues and eigenvectors of a hermitian complex matrix
template<typename T>
inline
bool
auxlib::eig_sym(Col<T>& eigval, Mat< std::complex<T> >& eigvec, const Mat< std::complex<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename std::complex<T> eT;
    
    arma_debug_check( (X.is_square() == false), "eig_sym(): given matrix must be square sized" );
    
    if(arma_config::check_nonfinite && trimat_helper::has_nonfinite_triu(X))  { return false; }
    
    eigvec = X;
    
    if(eigvec.is_empty())  { eigval.reset(); eigvec.reset(); return true; }
    
    arma_debug_assert_blas_size(eigvec);
    
    eigval.set_size(eigvec.n_rows);
    
    char jobz  = 'V';
    char uplo  = 'U';
    
    blas_int N     = blas_int(eigvec.n_rows);
    blas_int lwork = (64+1)*N;  // lwork_min = (std::max)(blas_int(1), 2*N-1)
    blas_int info  = 0;
    
    podarray<eT>  work( static_cast<uword>(lwork) );
    podarray<T>  rwork( static_cast<uword>((std::max)(blas_int(1), 3*N)) );
    
    arma_extra_debug_print("lapack::heev()");
    lapack::heev(&jobz, &uplo, &N, eigvec.memptr(), &N, eigval.memptr(), work.memptr(), &lwork, rwork.memptr(), &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_stop_logic_error("eig_sym(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! eigenvalues and eigenvectors of a symmetric real matrix (divide and conquer algorithm)
template<typename eT>
inline
bool
auxlib::eig_sym_dc(Col<eT>& eigval, Mat<eT>& eigvec, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_check( (X.is_square() == false), "eig_sym(): given matrix must be square sized" );
    
    if(arma_config::check_nonfinite && trimat_helper::has_nonfinite_triu(X))  { return false; }
    
    eigvec = X;
    
    if(eigvec.is_empty())  { eigval.reset(); eigvec.reset(); return true; }
    
    arma_debug_assert_blas_size(eigvec);
    
    eigval.set_size(eigvec.n_rows);
    
    char jobz = 'V';
    char uplo = 'U';
    
    blas_int N          = blas_int(eigvec.n_rows);
    blas_int lwork_min  = 1 + 6*N + 2*(N*N);
    blas_int liwork_min = 3 + 5*N;
    blas_int info       = 0;
    
    blas_int  lwork_proposed = 0;
    blas_int liwork_proposed = 0;
    
    if(N >= 32)
      {
      eT        work_query[2] = {};
      blas_int iwork_query[2] = {};
      
      blas_int  lwork_query = -1;
      blas_int liwork_query = -1;
      
      arma_extra_debug_print("lapack::syevd()");
      lapack::syevd(&jobz, &uplo, &N, eigvec.memptr(), &N, eigval.memptr(), &work_query[0], &lwork_query, &iwork_query[0], &liwork_query, &info);
      
      if(info != 0)  { return false; }
      
       lwork_proposed = static_cast<blas_int>( work_query[0] );
      liwork_proposed = iwork_query[0];
      }
    
    blas_int  lwork_final = (std::max)( lwork_proposed,  lwork_min);
    blas_int liwork_final = (std::max)(liwork_proposed, liwork_min);
    
    podarray<eT>        work( static_cast<uword>( lwork_final) );
    podarray<blas_int> iwork( static_cast<uword>(liwork_final) ); 
    
    arma_extra_debug_print("lapack::syevd()");
    lapack::syevd(&jobz, &uplo, &N, eigvec.memptr(), &N, eigval.memptr(), work.memptr(), &lwork_final, iwork.memptr(), &liwork_final, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_stop_logic_error("eig_sym(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! eigenvalues and eigenvectors of a hermitian complex matrix (divide and conquer algorithm)
template<typename T>
inline
bool
auxlib::eig_sym_dc(Col<T>& eigval, Mat< std::complex<T> >& eigvec, const Mat< std::complex<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename std::complex<T> eT;
    
    arma_debug_check( (X.is_square() == false), "eig_sym(): given matrix must be square sized" );
    
    if(arma_config::check_nonfinite && trimat_helper::has_nonfinite_triu(X))  { return false; }
    
    eigvec = X;
    
    if(eigvec.is_empty())  { eigval.reset(); eigvec.reset(); return true; }
    
    arma_debug_assert_blas_size(eigvec);
    
    eigval.set_size(eigvec.n_rows);
    
    char jobz  = 'V';
    char uplo  = 'U';
    
    blas_int N          = blas_int(eigvec.n_rows);
    blas_int lwork_min  = 2*N + N*N;
    blas_int lrwork_min = 1 + 5*N + 2*(N*N);
    blas_int liwork_min = 3 + 5*N;
    blas_int info       = 0;
    
    blas_int  lwork_proposed = 0;
    blas_int lrwork_proposed = 0;
    blas_int liwork_proposed = 0;
    
    if(N >= 32)
      {
      eT        work_query[2] = {};
      T        rwork_query[2] = {};
      blas_int iwork_query[2] = {};
      
      blas_int  lwork_query = -1;
      blas_int lrwork_query = -1;
      blas_int liwork_query = -1;
      
      arma_extra_debug_print("lapack::heevd()");
      lapack::heevd(&jobz, &uplo, &N, eigvec.memptr(), &N, eigval.memptr(), &work_query[0], &lwork_query, &rwork_query[0], &lrwork_query, &iwork_query[0], &liwork_query, &info);
      
      if(info != 0)  { return false; }
      
       lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
      lrwork_proposed = static_cast<blas_int>( rwork_query[0] );
      liwork_proposed = iwork_query[0];
      }
    
    blas_int  lwork_final = (std::max)( lwork_proposed,  lwork_min);
    blas_int lrwork_final = (std::max)(lrwork_proposed, lrwork_min);
    blas_int liwork_final = (std::max)(liwork_proposed, liwork_min);
    
    podarray<eT>        work( static_cast<uword>( lwork_final) );
    podarray< T>       rwork( static_cast<uword>(lrwork_final) );
    podarray<blas_int> iwork( static_cast<uword>(liwork_final) ); 
    
    arma_extra_debug_print("lapack::heevd()");
    lapack::heevd(&jobz, &uplo, &N, eigvec.memptr(), &N, eigval.memptr(), work.memptr(), &lwork_final, rwork.memptr(), &lrwork_final, iwork.memptr(), &liwork_final, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_stop_logic_error("eig_sym(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::chol_simple(Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(X);
    
    char      uplo = 'U';
    blas_int  n    = blas_int(X.n_rows);
    blas_int  info = 0;
    
    arma_extra_debug_print("lapack::potrf()");
    lapack::potrf(&uplo, &n, X.memptr(), &n, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(X);
    
    arma_stop_logic_error("chol(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::chol(Mat<eT>& X, const uword layout)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(X);
    
    char      uplo = (layout == 0) ? 'U' : 'L';
    blas_int  n    = blas_int(X.n_rows);
    blas_int  info = 0;
    
    arma_extra_debug_print("lapack::potrf()");
    lapack::potrf(&uplo, &n, X.memptr(), &n, &info);
    
    if(info != 0)  { return false; }
    
    X = (layout == 0) ? trimatu(X) : trimatl(X);  // trimatu() and trimatl() return the same type
    
    return true;
    }
  #else
    {
    arma_ignore(X);
    arma_ignore(layout);
    
    arma_stop_logic_error("chol(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::chol_band(Mat<eT>& X, const uword KD, const uword layout)
  {
  arma_extra_debug_sigprint();
  
  return auxlib::chol_band_common(X, KD, layout);
  }



template<typename T>
inline
bool
auxlib::chol_band(Mat< std::complex<T> >& X, const uword KD, const uword layout)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_extra_debug_print("auxlib::chol_band(): redirecting to auxlib::chol() due to crippled LAPACK");
    
    arma_ignore(KD);
    
    return auxlib::chol(X, layout);
    }
  #else
    {
    return auxlib::chol_band_common(X, KD, layout);
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::chol_band_common(Mat<eT>& X, const uword KD, const uword layout)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    const uword N = X.n_rows;
    
    const uword KL = (layout == 0) ? uword(0) : KD;
    const uword KU = (layout == 0) ? KD       : uword(0);
    
    Mat<eT> AB;
    band_helper::compress(AB, X, KL, KU, false);
    
    arma_debug_assert_blas_size(AB);
    
    char     uplo = (layout == 0) ? 'U' : 'L';
    blas_int n    = blas_int(N);
    blas_int kd   = blas_int(KD);
    blas_int ldab = blas_int(AB.n_rows);
    blas_int info = 0;
    
    arma_extra_debug_print("lapack::pbtrf()");
    lapack::pbtrf(&uplo, &n, &kd, AB.memptr(), &ldab, &info);
    
    if(info != 0)  { return false; }
    
    band_helper::uncompress(X, AB, KL, KU, false);
    
    return true;
    }
  #else
    {
    arma_ignore(X);
    arma_ignore(KD);
    arma_ignore(layout);
    
    arma_stop_logic_error("chol(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::chol_pivot(Mat<eT>& X, Mat<uword>& P, const uword layout)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename get_pod_type<eT>::result T;
    
    arma_debug_assert_blas_size(X);
    
    char     uplo = (layout == 0) ? 'U' : 'L';
    blas_int n    = blas_int(X.n_rows);
    blas_int rank = 0;
    T        tol  = T(-1);
    blas_int info = 0;
    
    podarray<blas_int> ipiv(  X.n_rows);
    podarray<T>        work(2*X.n_rows);
    
    ipiv.zeros();
    
    arma_extra_debug_print("lapack::pstrf()");
    lapack::pstrf(&uplo, &n, X.memptr(), &n, ipiv.memptr(), &rank, &tol, work.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    X = (layout == 0) ? trimatu(X) : trimatl(X);  // trimatu() and trimatl() return the same type
    
    P.set_size(X.n_rows, 1);
    
    for(uword i=0; i < X.n_rows; ++i)
      {
      P[i] = uword(ipiv[i] - 1);  // take into account that Fortran counts from 1
      }
    
    return true;
    }
  #else
    {
    arma_ignore(X);
    arma_ignore(P);
    arma_ignore(layout);
    
    arma_stop_logic_error("chol(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//
// hessenberg decomposition
template<typename eT, typename T1>
inline
bool
auxlib::hess(Mat<eT>& H, const Base<eT,T1>& X, Col<eT>& tao)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    H = X.get_ref();
    
    arma_debug_check( (H.is_square() == false), "hess(): given matrix must be square sized" );
    
    if(H.is_empty())  { return true; }
    
    arma_debug_assert_blas_size(H);
    
    if(H.n_rows > 2)
      {
      tao.set_size(H.n_rows-1);
      
      blas_int  n      = blas_int(H.n_rows);
      blas_int  ilo    = 1;
      blas_int  ihi    = blas_int(H.n_rows);
      blas_int  lda    = blas_int(H.n_rows);
      blas_int  lwork  = blas_int(H.n_rows) * 64;
      blas_int  info   = 0;
      
      podarray<eT> work(static_cast<uword>(lwork));
      
      arma_extra_debug_print("lapack::gehrd()");
      lapack::gehrd(&n, &ilo, &ihi, H.memptr(), &lda, tao.memptr(), work.memptr(), &lwork, &info);
      
      return (info == 0);
      }
    
    return true;
    }
  #else
    {
    arma_ignore(H);
    arma_ignore(X);
    arma_ignore(tao);
    arma_stop_logic_error("hess(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT, typename T1>
inline
bool
auxlib::qr(Mat<eT>& Q, Mat<eT>& R, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    R = X.get_ref();
    
    const uword R_n_rows = R.n_rows;
    const uword R_n_cols = R.n_cols;
    
    if(R.is_empty())  { Q.eye(R_n_rows, R_n_rows); return true; }
    
    arma_debug_assert_blas_size(R);
    
    blas_int m         = static_cast<blas_int>(R_n_rows);
    blas_int n         = static_cast<blas_int>(R_n_cols);
    blas_int lwork_min = (std::max)(blas_int(1), (std::max)(m,n));  // take into account requirements of geqrf() _and_ orgqr()/ungqr()
    blas_int k         = (std::min)(m,n);
    blas_int info      = 0;
    
    podarray<eT> tau( static_cast<uword>(k) );
    
    eT        work_query[2] = {};
    blas_int lwork_query    = -1;
    
    arma_extra_debug_print("lapack::geqrf()");
    lapack::geqrf(&m, &n, R.memptr(), &m, tau.memptr(), &work_query[0], &lwork_query, &info);
    
    if(info != 0)  { return false; }
    
    blas_int lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
    blas_int lwork_final    = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::geqrf()");
    lapack::geqrf(&m, &n, R.memptr(), &m, tau.memptr(), work.memptr(), &lwork_final, &info);
    
    if(info != 0)  { return false; }
    
    Q.set_size(R_n_rows, R_n_rows);
    
    arrayops::copy( Q.memptr(), R.memptr(), (std::min)(Q.n_elem, R.n_elem) );
    
    //
    // construct R
    
    for(uword col=0; col < R_n_cols; ++col)
      {
      for(uword row=(col+1); row < R_n_rows; ++row)
        {
        R.at(row,col) = eT(0);
        }
      }
    
    
    if( (is_float<eT>::value) || (is_double<eT>::value) )
      {
      arma_extra_debug_print("lapack::orgqr()");
      lapack::orgqr(&m, &m, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &lwork_final, &info);
      }
    else
    if( (is_cx_float<eT>::value) || (is_cx_double<eT>::value) )
      {
      arma_extra_debug_print("lapack::ungqr()");
      lapack::ungqr(&m, &m, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &lwork_final, &info);
      }
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(Q);
    arma_ignore(R);
    arma_ignore(X);
    arma_stop_logic_error("qr(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT, typename T1>
inline
bool 
auxlib::qr_econ(Mat<eT>& Q, Mat<eT>& R, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    if(is_Mat<T1>::value)
      {
      const unwrap<T1>   tmp(X.get_ref());
      const Mat<eT>& M = tmp.M;
      
      if(M.n_rows < M.n_cols)  { return auxlib::qr(Q, R, X); }
      }
    
    Q = X.get_ref();
    
    const uword Q_n_rows = Q.n_rows;
    const uword Q_n_cols = Q.n_cols;
    
    if( Q_n_rows <= Q_n_cols )  { return auxlib::qr(Q, R, Q); }
    
    if(Q.is_empty())  { Q.set_size(Q_n_rows, 0); R.set_size(0, Q_n_cols); return true; }
    
    arma_debug_assert_blas_size(Q);
    
    blas_int m         = static_cast<blas_int>(Q_n_rows);
    blas_int n         = static_cast<blas_int>(Q_n_cols);
    blas_int lwork_min = (std::max)(blas_int(1), (std::max)(m,n));  // take into account requirements of geqrf() _and_ orgqr()/ungqr()
    blas_int k         = (std::min)(m,n);
    blas_int info      = 0;
    
    podarray<eT> tau( static_cast<uword>(k) );
    
    eT        work_query[2] = {};
    blas_int lwork_query    = -1;
    
    arma_extra_debug_print("lapack::geqrf()");
    lapack::geqrf(&m, &n, Q.memptr(), &m, tau.memptr(), &work_query[0], &lwork_query, &info);
    
    if(info != 0)  { return false; }
    
    blas_int lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
    blas_int lwork_final    = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::geqrf()");
    lapack::geqrf(&m, &n, Q.memptr(), &m, tau.memptr(), work.memptr(), &lwork_final, &info);
    
    if(info != 0)  { return false; }
    
    R.zeros(Q_n_cols, Q_n_cols);
    
    //
    // construct R
    
    for(uword col=0; col < Q_n_cols; ++col)
      {
      for(uword row=0; row <= col; ++row)
        {
        R.at(row,col) = Q.at(row,col);
        }
      }
    
    if( (is_float<eT>::value) || (is_double<eT>::value) )
      {
      arma_extra_debug_print("lapack::orgqr()");
      lapack::orgqr(&m, &n, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &lwork_final, &info);
      }
    else
    if( (is_cx_float<eT>::value) || (is_cx_double<eT>::value) )
      {
      arma_extra_debug_print("lapack::ungqr()");
      lapack::ungqr(&m, &n, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &lwork_final, &info);
      }
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(Q);
    arma_ignore(R);
    arma_ignore(X);
    arma_stop_logic_error("qr_econ(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT, typename T1>
inline
bool
auxlib::qr_pivot(Mat<eT>& Q, Mat<eT>& R, Mat<uword>& P, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    R = X.get_ref();
    
    const uword R_n_rows = R.n_rows;
    const uword R_n_cols = R.n_cols;
    
    if(R.is_empty())
      {
      Q.eye(R_n_rows, R_n_rows);
      
      P.set_size(R_n_cols, 1);
      
      for(uword col=0; col < R_n_cols; ++col)  { P.at(col) = col; }
      
      return true;
      }
    
    arma_debug_assert_blas_size(R);
    
    blas_int m         = static_cast<blas_int>(R_n_rows);
    blas_int n         = static_cast<blas_int>(R_n_cols);
    blas_int lwork_min = (std::max)(blas_int(3*n + 1), (std::max)(m,n));  // take into account requirements of geqp3() and orgqr()
    blas_int k         = (std::min)(m,n);
    blas_int info      = 0;
    
    podarray<eT>        tau( static_cast<uword>(k) );
    podarray<blas_int> jpvt( R_n_cols );
    
    jpvt.zeros();
    
    eT        work_query[2] = {};
    blas_int lwork_query    = -1;
    
    arma_extra_debug_print("lapack::geqp3()");
    lapack::geqp3(&m, &n, R.memptr(), &m, jpvt.memptr(), tau.memptr(), &work_query[0], &lwork_query, &info);
    
    if(info != 0)  { return false; }
    
    blas_int lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
    blas_int lwork_final    = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::geqp3()");
    lapack::geqp3(&m, &n, R.memptr(), &m, jpvt.memptr(), tau.memptr(), work.memptr(), &lwork_final, &info);
    
    if(info != 0)  { return false; }
    
    Q.set_size(R_n_rows, R_n_rows);
    
    arrayops::copy( Q.memptr(), R.memptr(), (std::min)(Q.n_elem, R.n_elem) );
    
    //
    // construct R and P
    
    P.set_size(R_n_cols, 1);
    
    for(uword col=0; col < R_n_cols; ++col)
      {
      for(uword row=(col+1); row < R_n_rows; ++row)  { R.at(row,col) = eT(0); }
      
      P.at(col) = jpvt[col] - 1;  // take into account that Fortran counts from 1
      }
    
    arma_extra_debug_print("lapack::orgqr()");
    lapack::orgqr(&m, &m, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &lwork_final, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(Q);
    arma_ignore(R);
    arma_ignore(P);
    arma_ignore(X);
    arma_stop_logic_error("qr(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T, typename T1>
inline
bool
auxlib::qr_pivot(Mat< std::complex<T> >& Q, Mat< std::complex<T> >& R, Mat<uword>& P, const Base<std::complex<T>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename std::complex<T> eT;
    
    R = X.get_ref();
    
    const uword R_n_rows = R.n_rows;
    const uword R_n_cols = R.n_cols;
    
    if(R.is_empty())
      {
      Q.eye(R_n_rows, R_n_rows);
      
      P.set_size(R_n_cols, 1);
      
      for(uword col=0; col < R_n_cols; ++col)  { P.at(col) = col; }
      
      return true;
      }
    
    arma_debug_assert_blas_size(R);
    
    blas_int m         = static_cast<blas_int>(R_n_rows);
    blas_int n         = static_cast<blas_int>(R_n_cols);
    blas_int lwork_min = (std::max)(blas_int(3*n + 1), (std::max)(m,n));  // take into account requirements of geqp3() and ungqr()
    blas_int k         = (std::min)(m,n);
    blas_int info      = 0;
    
    podarray<eT>         tau( static_cast<uword>(k) );
    podarray< T>       rwork( 2*R_n_cols );
    podarray<blas_int>  jpvt( R_n_cols );
    
    jpvt.zeros();
    
    eT        work_query[2] = {};
    blas_int lwork_query    = -1;
    
    arma_extra_debug_print("lapack::geqp3()");
    lapack::cx_geqp3(&m, &n, R.memptr(), &m, jpvt.memptr(), tau.memptr(), &work_query[0], &lwork_query, rwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    blas_int lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
    blas_int lwork_final    = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::geqp3()");
    lapack::cx_geqp3(&m, &n, R.memptr(), &m, jpvt.memptr(), tau.memptr(), work.memptr(), &lwork_final, rwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    Q.set_size(R_n_rows, R_n_rows);
    
    arrayops::copy( Q.memptr(), R.memptr(), (std::min)(Q.n_elem, R.n_elem) );
    
    //
    // construct R and P
    
    P.set_size(R_n_cols, 1);
    
    for(uword col=0; col < R_n_cols; ++col)
      {
      for(uword row=(col+1); row < R_n_rows; ++row)  { R.at(row,col) = eT(0); }
      
      P.at(col) = jpvt[col] - 1;  // take into account that Fortran counts from 1
      }
    
    arma_extra_debug_print("lapack::ungqr()");
    lapack::ungqr(&m, &m, &k, Q.memptr(), &m, tau.memptr(), work.memptr(), &lwork_final, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(Q);
    arma_ignore(R);
    arma_ignore(P);
    arma_ignore(X);
    arma_stop_logic_error("qr(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::svd(Col<eT>& S, Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    if(A.is_empty())  { S.reset(); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    Mat<eT> U(1, 1,        arma_nozeros_indicator());
    Mat<eT> V(1, A.n_cols, arma_nozeros_indicator());
    
    char jobu  = 'N';
    char jobvt = 'N';
    
    blas_int m         = blas_int(A.n_rows);
    blas_int n         = blas_int(A.n_cols);
    blas_int min_mn    = (std::min)(m,n);
    blas_int lda       = blas_int(A.n_rows);
    blas_int ldu       = blas_int(U.n_rows);
    blas_int ldvt      = blas_int(V.n_rows);
    blas_int lwork_min = (std::max)( blas_int(1), (std::max)( (3*min_mn + (std::max)(m,n)), 5*min_mn ) );
    blas_int info      = 0;
    
    S.set_size( static_cast<uword>(min_mn) );
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= 1024)
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = -1;
      
      arma_extra_debug_print("lapack::gesvd()");
      lapack::gesvd<eT>(&jobu, &jobvt, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, &work_query[0], &lwork_query, &info);
      
      if(info != 0)  { return false; }
      
      lwork_proposed = static_cast<blas_int>( work_query[0] );
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::gesvd()");
    lapack::gesvd<eT>(&jobu, &jobvt, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, work.memptr(), &lwork_final, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(S);
    arma_ignore(A);
    arma_stop_logic_error("svd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T>
inline
bool
auxlib::svd(Col<T>& S, Mat< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef std::complex<T> eT;
    
    if(A.is_empty())  { S.reset(); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    Mat<eT> U(1, 1,        arma_nozeros_indicator());
    Mat<eT> V(1, A.n_cols, arma_nozeros_indicator());
    
    char jobu  = 'N';
    char jobvt = 'N';
    
    blas_int  m         = blas_int(A.n_rows);
    blas_int  n         = blas_int(A.n_cols);
    blas_int  min_mn    = (std::min)(m,n);
    blas_int  lda       = blas_int(A.n_rows);
    blas_int  ldu       = blas_int(U.n_rows);
    blas_int  ldvt      = blas_int(V.n_rows);
    blas_int  lwork_min = (std::max)( blas_int(1), 2*min_mn+(std::max)(m,n) );
    blas_int  info      = 0;
    
    S.set_size( static_cast<uword>(min_mn) );
    
    podarray<T> rwork( static_cast<uword>(5*min_mn) );
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= 256)
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = -1;  // query to find optimum size of workspace
      
      arma_extra_debug_print("lapack::cx_gesvd()");
      lapack::cx_gesvd<T>(&jobu, &jobvt, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, &work_query[0], &lwork_query, rwork.memptr(), &info);
      
      if(info != 0)  { return false; }
      
      lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::cx_gesvd()");
    lapack::cx_gesvd<T>(&jobu, &jobvt, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, work.memptr(), &lwork_final, rwork.memptr(), &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(S);
    arma_ignore(A);
    arma_stop_logic_error("svd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::svd(Mat<eT>& U, Col<eT>& S, Mat<eT>& V, Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    if(A.is_empty())  { U.eye(A.n_rows, A.n_rows); S.reset(); V.eye(A.n_cols, A.n_cols); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    U.set_size(A.n_rows, A.n_rows);
    V.set_size(A.n_cols, A.n_cols);
    
    char jobu  = 'A';
    char jobvt = 'A';
    
    blas_int  m         = blas_int(A.n_rows);
    blas_int  n         = blas_int(A.n_cols);
    blas_int  min_mn    = (std::min)(m,n);
    blas_int  lda       = blas_int(A.n_rows);
    blas_int  ldu       = blas_int(U.n_rows);
    blas_int  ldvt      = blas_int(V.n_rows);
    blas_int  lwork_min = (std::max)( blas_int(1), (std::max)( (3*min_mn + (std::max)(m,n)), 5*min_mn ) );
    blas_int  info      = 0;
    
    S.set_size( static_cast<uword>(min_mn) );
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= 1024)
      {
      // query to find optimum size of workspace
      eT        work_query[2] = {};
      blas_int lwork_query    = -1;
      
      arma_extra_debug_print("lapack::gesvd()");
      lapack::gesvd<eT>(&jobu, &jobvt, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, &work_query[0], &lwork_query, &info);
      
      if(info != 0)  { return false; }
      
      lwork_proposed = static_cast<blas_int>( work_query[0] );
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::gesvd()");
    lapack::gesvd<eT>(&jobu, &jobvt, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, work.memptr(), &lwork_final, &info);
    
    if(info != 0)  { return false; }
    
    op_strans::apply_mat_inplace(V);
    
    return true;
    }
  #else
    {
    arma_ignore(U);
    arma_ignore(S);
    arma_ignore(V);
    arma_ignore(A);
    arma_stop_logic_error("svd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T>
inline
bool
auxlib::svd(Mat< std::complex<T> >& U, Col<T>& S, Mat< std::complex<T> >& V, Mat< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef std::complex<T> eT;
    
    if(A.is_empty())  { U.eye(A.n_rows, A.n_rows); S.reset(); V.eye(A.n_cols, A.n_cols); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    U.set_size(A.n_rows, A.n_rows);
    V.set_size(A.n_cols, A.n_cols);
    
    char jobu  = 'A';
    char jobvt = 'A';
    
    blas_int  m         = blas_int(A.n_rows);
    blas_int  n         = blas_int(A.n_cols);
    blas_int  min_mn    = (std::min)(m,n);
    blas_int  lda       = blas_int(A.n_rows);
    blas_int  ldu       = blas_int(U.n_rows);
    blas_int  ldvt      = blas_int(V.n_rows);
    blas_int  lwork_min = (std::max)( blas_int(1), 2*min_mn + (std::max)(m,n) );
    blas_int  info      = 0;
    
    S.set_size( static_cast<uword>(min_mn) );
    
    podarray<T> rwork( static_cast<uword>(5*min_mn) );
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= 256)
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = -1;  // query to find optimum size of workspace
      
      arma_extra_debug_print("lapack::cx_gesvd()");
      lapack::cx_gesvd<T>(&jobu, &jobvt, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, &work_query[0], &lwork_query, rwork.memptr(), &info);
      
      if(info != 0)  { return false; }

      lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::cx_gesvd()");
    lapack::cx_gesvd<T>(&jobu, &jobvt, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, work.memptr(), &lwork_final, rwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    op_htrans::apply_mat_inplace(V);
    
    return true;
    }
  #else
    {
    arma_ignore(U);
    arma_ignore(S);
    arma_ignore(V);
    arma_ignore(A);
    arma_stop_logic_error("svd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::svd_econ(Mat<eT>& U, Col<eT>& S, Mat<eT>& V, Mat<eT>& A, const char mode)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    if(A.is_empty())  { U.eye(); S.reset(); V.eye(); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    blas_int m      = blas_int(A.n_rows);
    blas_int n      = blas_int(A.n_cols);
    blas_int min_mn = (std::min)(m,n);
    blas_int lda    = blas_int(A.n_rows);
    
    S.set_size( static_cast<uword>(min_mn) );
    
    blas_int ldu  = 0;
    blas_int ldvt = 0;
    
    char jobu  = char(0);
    char jobvt = char(0);
    
    if(mode == 'l')
      {
      jobu  = 'S';
      jobvt = 'N';
      
      ldu  = m;
      ldvt = 1;
      
      U.set_size( static_cast<uword>(ldu), static_cast<uword>(min_mn) );
      V.reset();
      }
    
    if(mode == 'r')
      {
      jobu  = 'N';
      jobvt = 'S';
      
      ldu = 1;
      ldvt = (std::min)(m,n);
      
      U.reset();
      V.set_size( static_cast<uword>(ldvt), static_cast<uword>(n) );
      }
    
    if(mode == 'b')
      {
      jobu  = 'S';
      jobvt = 'S';
      
      ldu  = m;
      ldvt = (std::min)(m,n);
      
      U.set_size( static_cast<uword>(ldu),  static_cast<uword>(min_mn) );
      V.set_size( static_cast<uword>(ldvt), static_cast<uword>(n     ) );
      }
    
    
    blas_int lwork_min = (std::max)( blas_int(1), (std::max)( (3*min_mn + (std::max)(m,n)), 5*min_mn ) );
    blas_int info      = 0;
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= 1024)
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = -1;  // query to find optimum size of workspace
      
      arma_extra_debug_print("lapack::gesvd()");
      lapack::gesvd<eT>(&jobu, &jobvt, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, &work_query[0], &lwork_query, &info);
      
      if(info != 0)  { return false; }
      
      lwork_proposed = static_cast<blas_int>(work_query[0]);
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::gesvd()");
    lapack::gesvd<eT>(&jobu, &jobvt, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, work.memptr(), &lwork_final, &info);
    
    if(info != 0)  { return false; }
    
    op_strans::apply_mat_inplace(V);
    
    return true;
    }
  #else
    {
    arma_ignore(U);
    arma_ignore(S);
    arma_ignore(V);
    arma_ignore(A);
    arma_ignore(mode);
    arma_stop_logic_error("svd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T>
inline
bool
auxlib::svd_econ(Mat< std::complex<T> >& U, Col<T>& S, Mat< std::complex<T> >& V, Mat< std::complex<T> >& A, const char mode)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef std::complex<T> eT;
    
    if(A.is_empty())  { U.eye(); S.reset(); V.eye(); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    blas_int m      = blas_int(A.n_rows);
    blas_int n      = blas_int(A.n_cols);
    blas_int min_mn = (std::min)(m,n);
    blas_int lda    = blas_int(A.n_rows);
    
    S.set_size( static_cast<uword>(min_mn) );
    
    blas_int ldu  = 0;
    blas_int ldvt = 0;
    
    char jobu  = char(0);
    char jobvt = char(0);
    
    if(mode == 'l')
      {
      jobu  = 'S';
      jobvt = 'N';
      
      ldu  = m;
      ldvt = 1;
      
      U.set_size( static_cast<uword>(ldu), static_cast<uword>(min_mn) );
      V.reset();
      }
    
    if(mode == 'r')
      {
      jobu  = 'N';
      jobvt = 'S';
      
      ldu  = 1;
      ldvt = (std::min)(m,n);
      
      U.reset();
      V.set_size( static_cast<uword>(ldvt), static_cast<uword>(n) );
      }
    
    if(mode == 'b')
      {
      jobu  = 'S';
      jobvt = 'S';
      
      ldu  = m;
      ldvt = (std::min)(m,n);
      
      U.set_size( static_cast<uword>(ldu),  static_cast<uword>(min_mn) );
      V.set_size( static_cast<uword>(ldvt), static_cast<uword>(n)      );
      }
    
    blas_int lwork_min = (std::max)( blas_int(1), (std::max)( (3*min_mn + (std::max)(m,n)), 5*min_mn ) );
    blas_int info      = 0;
    
    podarray<T> rwork( static_cast<uword>(5*min_mn) );
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= 256)
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = -1;  // query to find optimum size of workspace
      
      arma_extra_debug_print("lapack::cx_gesvd()");
      lapack::cx_gesvd<T>(&jobu, &jobvt, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, &work_query[0], &lwork_query, rwork.memptr(), &info);
      
      if(info != 0)  { return false; }
      
      lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::cx_gesvd()");
    lapack::cx_gesvd<T>(&jobu, &jobvt, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, work.memptr(), &lwork_final, rwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    op_htrans::apply_mat_inplace(V);
    
    return true;
    }
  #else
    {
    arma_ignore(U);
    arma_ignore(S);
    arma_ignore(V);
    arma_ignore(A);
    arma_ignore(mode);
    arma_stop_logic_error("svd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::svd_dc(Col<eT>& S, Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    if(A.is_empty())  { S.reset(); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    Mat<eT> U(1, 1, arma_nozeros_indicator());
    Mat<eT> V(1, 1, arma_nozeros_indicator());
    
    char jobz = 'N';
    
    blas_int  m         = blas_int(A.n_rows);
    blas_int  n         = blas_int(A.n_cols);
    blas_int  min_mn    = (std::min)(m,n);
    blas_int  max_mn    = (std::max)(m,n);
    blas_int  lda       = blas_int(A.n_rows);
    blas_int  ldu       = blas_int(U.n_rows);
    blas_int  ldvt      = blas_int(V.n_rows);
    blas_int  lwork_min = 3*min_mn + (std::max)( max_mn, 7*min_mn );
    blas_int  info      = 0;
    
    S.set_size( static_cast<uword>(min_mn) );
    
    podarray<blas_int> iwork( static_cast<uword>(8*min_mn) );
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= 1024)
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = blas_int(-1);
      
      arma_extra_debug_print("lapack::gesdd()");
      lapack::gesdd<eT>(&jobz, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, &work_query[0], &lwork_query, iwork.memptr(), &info);
      
      if(info != 0)  { return false; }
      
      lwork_proposed = static_cast<blas_int>( work_query[0] );
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::gesdd()");
    lapack::gesdd<eT>(&jobz, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, work.memptr(), &lwork_final, iwork.memptr(), &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(S);
    arma_ignore(A);
    arma_stop_logic_error("svd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T>
inline
bool
auxlib::svd_dc(Col<T>& S, Mat< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef std::complex<T> eT;
    
    if(A.is_empty())  { S.reset(); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    Mat<eT> U(1, 1, arma_nozeros_indicator());
    Mat<eT> V(1, 1, arma_nozeros_indicator());
    
    char jobz = 'N';
    
    blas_int  m         = blas_int(A.n_rows);
    blas_int  n         = blas_int(A.n_cols);
    blas_int  min_mn    = (std::min)(m,n);
    blas_int  max_mn    = (std::max)(m,n);
    blas_int  lda       = blas_int(A.n_rows);
    blas_int  ldu       = blas_int(U.n_rows);
    blas_int  ldvt      = blas_int(V.n_rows);
    blas_int  lwork_min = 2*min_mn + max_mn;
    blas_int  info      = 0;
    
    S.set_size( static_cast<uword>(min_mn) );
    
    podarray<T>        rwork( static_cast<uword>(7*min_mn) );  // from LAPACK 3.8 docs: LAPACK <= v3.6 needs 7*mn
    podarray<blas_int> iwork( static_cast<uword>(8*min_mn) );
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= 256)
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = blas_int(-1);
      
      arma_extra_debug_print("lapack::cx_gesdd()");
      lapack::cx_gesdd<T>(&jobz, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, &work_query[0], &lwork_query, rwork.memptr(), iwork.memptr(), &info);
      
      if(info != 0)  { return false; }
      
      lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::cx_gesdd()");
    lapack::cx_gesdd<T>(&jobz, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, work.memptr(), &lwork_final, rwork.memptr(), iwork.memptr(), &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(S);
    arma_ignore(A);
    arma_stop_logic_error("svd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::svd_dc(Mat<eT>& U, Col<eT>& S, Mat<eT>& V, Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    if(A.is_empty())  { U.eye(A.n_rows, A.n_rows); S.reset(); V.eye(A.n_cols, A.n_cols); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    U.set_size(A.n_rows, A.n_rows);
    V.set_size(A.n_cols, A.n_cols);
    
    char jobz = 'A';
    
    blas_int  m         = blas_int(A.n_rows);
    blas_int  n         = blas_int(A.n_cols);
    blas_int  min_mn    = (std::min)(m,n);
    blas_int  max_mn    = (std::max)(m,n);
    blas_int  lda       = blas_int(A.n_rows);
    blas_int  ldu       = blas_int(U.n_rows);
    blas_int  ldvt      = blas_int(V.n_rows);
    blas_int  lwork1    = 3*min_mn*min_mn + (std::max)(max_mn, 4*min_mn*min_mn + 4*min_mn);  // as per LAPACK 3.2 docs
    blas_int  lwork2    = 4*min_mn*min_mn + 6*min_mn + max_mn;  // as per LAPACK 3.8 docs; consistent with LAPACK 3.4 docs
    blas_int  lwork_min = (std::max)(lwork1, lwork2);  // due to differences between LAPACK 3.2 and 3.8
    blas_int  info      = 0;
    
    S.set_size( static_cast<uword>(min_mn) );
    
    podarray<blas_int> iwork( static_cast<uword>(8*min_mn) );
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= 1024)
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = blas_int(-1);
      
      arma_extra_debug_print("lapack::gesdd()");
      lapack::gesdd<eT>(&jobz, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, &work_query[0], &lwork_query, iwork.memptr(), &info);
      
      if(info != 0)  { return false; }
      
      lwork_proposed = static_cast<blas_int>(work_query[0]);
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::gesdd()");
    lapack::gesdd<eT>(&jobz, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, work.memptr(), &lwork_final, iwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    op_strans::apply_mat_inplace(V);
    
    return true;
    }
  #else
    {
    arma_ignore(U);
    arma_ignore(S);
    arma_ignore(V);
    arma_ignore(A);
    arma_stop_logic_error("svd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T>
inline
bool
auxlib::svd_dc(Mat< std::complex<T> >& U, Col<T>& S, Mat< std::complex<T> >& V, Mat< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef std::complex<T> eT;
    
    if(A.is_empty())  { U.eye(A.n_rows, A.n_rows); S.reset(); V.eye(A.n_cols, A.n_cols); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    U.set_size(A.n_rows, A.n_rows);
    V.set_size(A.n_cols, A.n_cols);
    
    char jobz = 'A';
    
    blas_int m         = blas_int(A.n_rows);
    blas_int n         = blas_int(A.n_cols);
    blas_int min_mn    = (std::min)(m,n);
    blas_int max_mn    = (std::max)(m,n);
    blas_int lda       = blas_int(A.n_rows);
    blas_int ldu       = blas_int(U.n_rows);
    blas_int ldvt      = blas_int(V.n_rows);
    blas_int lwork_min = min_mn*min_mn + 2*min_mn + max_mn;  // as per LAPACK 3.2, 3.4, 3.8 docs
    blas_int lrwork    = min_mn * ((std::max)(5*min_mn+7, 2*max_mn + 2*min_mn+1));   // as per LAPACK 3.4 docs; LAPACK 3.8 uses 5*min_mn+5 instead of 5*min_mn+7
    blas_int info      = 0;
    
    S.set_size( static_cast<uword>(min_mn) );
    
    podarray<T>        rwork( static_cast<uword>(lrwork  ) );
    podarray<blas_int> iwork( static_cast<uword>(8*min_mn) );
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= 256)
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = blas_int(-1);
      
      arma_extra_debug_print("lapack::cx_gesdd()");
      lapack::cx_gesdd<T>(&jobz, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, &work_query[0], &lwork_query, rwork.memptr(), iwork.memptr(), &info);
      
      if(info != 0)  { return false; }
      
      lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::cx_gesdd()");
    lapack::cx_gesdd<T>(&jobz, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, work.memptr(), &lwork_final, rwork.memptr(), iwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    op_htrans::apply_mat_inplace(V);
    
    return true;
    }
  #else
    {
    arma_ignore(U);
    arma_ignore(S);
    arma_ignore(V);
    arma_ignore(A);
    arma_stop_logic_error("svd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::svd_dc_econ(Mat<eT>& U, Col<eT>& S, Mat<eT>& V, Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    char jobz = 'S';
    
    blas_int m         = blas_int(A.n_rows);
    blas_int n         = blas_int(A.n_cols);
    blas_int min_mn    = (std::min)(m,n);
    blas_int max_mn    = (std::max)(m,n);
    blas_int lda       = blas_int(A.n_rows);
    blas_int ldu       = m;
    blas_int ldvt      = min_mn;
    blas_int lwork1    = 3*min_mn*min_mn + (std::max)( max_mn, 4*min_mn*min_mn + 4*min_mn );  // as per LAPACK 3.2 docs
    blas_int lwork2    = 4*min_mn*min_mn + 6*min_mn + max_mn;  // as per LAPACK 3.4 docs; LAPACK 3.8 requires 4*min_mn*min_mn + 7*min_mn
    blas_int lwork_min = (std::max)(lwork1, lwork2);  // due to differences between LAPACK 3.2 and 3.4
    blas_int info      = 0;
    
    if(A.is_empty())
      {
      U.eye();
      S.reset();
      V.eye( static_cast<uword>(n), static_cast<uword>(min_mn) );
      return true;
      }
    
    S.set_size( static_cast<uword>(min_mn) );
    
    U.set_size( static_cast<uword>(m), static_cast<uword>(min_mn) );
    
    V.set_size( static_cast<uword>(min_mn), static_cast<uword>(n) );
    
    podarray<blas_int> iwork( static_cast<uword>(8*min_mn) );
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= 1024)
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = blas_int(-1);
      
      arma_extra_debug_print("lapack::gesdd()");
      lapack::gesdd<eT>(&jobz, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, &work_query[0], &lwork_query, iwork.memptr(), &info);
      
      if(info != 0)  { return false; }
      
      lwork_proposed = static_cast<blas_int>(work_query[0]);
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::gesdd()");
    lapack::gesdd<eT>(&jobz, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, work.memptr(), &lwork_final, iwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    op_strans::apply_mat_inplace(V);
    
    return true;
    }
  #else
    {
    arma_ignore(U);
    arma_ignore(S);
    arma_ignore(V);
    arma_ignore(A);
    arma_stop_logic_error("svd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T>
inline
bool
auxlib::svd_dc_econ(Mat< std::complex<T> >& U, Col<T>& S, Mat< std::complex<T> >& V, Mat< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef std::complex<T> eT;
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    char jobz = 'S';
    
    blas_int m         = blas_int(A.n_rows);
    blas_int n         = blas_int(A.n_cols);
    blas_int min_mn    = (std::min)(m,n);
    blas_int max_mn    = (std::max)(m,n);
    blas_int lda       = blas_int(A.n_rows);
    blas_int ldu       = m;
    blas_int ldvt      = min_mn;
    blas_int lwork_min = min_mn*min_mn + 2*min_mn + max_mn;  // as per LAPACK 3.2 docs
    blas_int lrwork    = min_mn * ((std::max)(5*min_mn+7, 2*max_mn + 2*min_mn+1));  // LAPACK 3.8 uses 5*min_mn+5 instead of 5*min_mn+7
    blas_int info      = 0;
    
    if(A.is_empty())
      {
      U.eye();
      S.reset();
      V.eye( static_cast<uword>(n), static_cast<uword>(min_mn) );
      return true;
      }
    
    S.set_size( static_cast<uword>(min_mn) );
    
    U.set_size( static_cast<uword>(m), static_cast<uword>(min_mn) );
    
    V.set_size( static_cast<uword>(min_mn), static_cast<uword>(n) );
    
    podarray<T>        rwork( static_cast<uword>(lrwork  ) );
    podarray<blas_int> iwork( static_cast<uword>(8*min_mn) );
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= 256)
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = blas_int(-1);
      
      arma_extra_debug_print("lapack::cx_gesdd()");
      lapack::cx_gesdd<T>(&jobz, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, &work_query[0], &lwork_query, rwork.memptr(), iwork.memptr(), &info);
      
      if(info != 0)  { return false; }
      
      lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::cx_gesdd()");
    lapack::cx_gesdd<T>(&jobz, &m, &n, A.memptr(), &lda, S.memptr(), U.memptr(), &ldu, V.memptr(), &ldvt, work.memptr(), &lwork_final, rwork.memptr(), iwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    op_htrans::apply_mat_inplace(V);
    
    return true;
    }
  #else
    {
    arma_ignore(U);
    arma_ignore(S);
    arma_ignore(V);
    arma_ignore(A);
    arma_stop_logic_error("svd(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a system of linear equations via LU decomposition
template<typename T1>
inline
bool
auxlib::solve_square_fast(Mat<typename T1::elem_type>& out, Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  out = B_expr.get_ref();
  
  const uword B_n_rows = out.n_rows;
  const uword B_n_cols = out.n_cols;
  
  arma_debug_check( (A.n_rows != B_n_rows), "solve(): number of rows in given matrices must be the same", [&](){ out.soft_reset(); } );
  
  if(A.is_empty() || out.is_empty())  { out.zeros(A.n_cols, B_n_cols); return true; }
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(A);
    
    blas_int n    = blas_int(A.n_rows);  // assuming A is square
    blas_int lda  = blas_int(A.n_rows);
    blas_int ldb  = blas_int(B_n_rows);
    blas_int nrhs = blas_int(B_n_cols);
    blas_int info = blas_int(0);
    
    podarray<blas_int> ipiv(A.n_rows + 2);  // +2 for paranoia: some versions of Lapack might be trashing memory
    
    arma_extra_debug_print("lapack::gesv()");
    lapack::gesv<eT>(&n, &nrhs, A.memptr(), &lda, ipiv.memptr(), out.memptr(), &ldb, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a system of linear equations via LU decomposition with rcond estimate
template<typename T1>
inline
bool
auxlib::solve_square_rcond(Mat<typename T1::elem_type>& out, typename T1::pod_type& out_rcond, Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::elem_type eT;
    typedef typename T1::pod_type   T;
    
    out_rcond = T(0);
    
    out = B_expr.get_ref();
    
    const uword B_n_rows = out.n_rows;
    const uword B_n_cols = out.n_cols;
    
    arma_debug_check( (A.n_rows != B_n_rows), "solve(): number of rows in given matrices must be the same", [&](){ out.soft_reset(); } );
    
    if(A.is_empty() || out.is_empty())  { out.zeros(A.n_cols, B_n_cols); return true; }
    
    arma_debug_assert_blas_size(A);
    
    char     norm_id  = '1';
    char     trans    = 'N';
    blas_int n        = blas_int(A.n_rows);  // assuming A is square
    blas_int lda      = blas_int(A.n_rows);
    blas_int ldb      = blas_int(B_n_rows);
    blas_int nrhs     = blas_int(B_n_cols);
    blas_int info     = blas_int(0);
    T        norm_val = T(0);
    
    podarray<T>        junk(1);
    podarray<blas_int> ipiv(A.n_rows + 2);  // +2 for paranoia
    
    arma_extra_debug_print("lapack::lange()");
    norm_val = lapack::lange<eT>(&norm_id, &n, &n, A.memptr(), &lda, junk.memptr());
    
    arma_extra_debug_print("lapack::getrf()");
    lapack::getrf<eT>(&n, &n, A.memptr(), &n, ipiv.memptr(), &info);
    
    if(info != blas_int(0))  { return false; }
    
    arma_extra_debug_print("lapack::getrs()");
    lapack::getrs<eT>(&trans, &n, &nrhs, A.memptr(), &lda, ipiv.memptr(), out.memptr(), &ldb, &info);
    
    if(info != blas_int(0))  { return false; }
    
    out_rcond = auxlib::lu_rcond<T>(A, norm_val);
    
    return true;
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a system of linear equations via LU decomposition with refinement (real matrices)
template<typename T1>
inline
bool
auxlib::solve_square_refine(Mat<typename T1::pod_type>& out, typename T1::pod_type& out_rcond, Mat<typename T1::pod_type>& A, const Base<typename T1::pod_type,T1>& B_expr, const bool equilibrate)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type eT;
    
    // Mat<eT> B = B_expr.get_ref();  // B is overwritten by lapack::gesvx() if equilibrate is enabled
    
    quasi_unwrap<T1> UB(B_expr.get_ref());  // deliberately not declaring as const
    
    const Mat<eT>& UB_M_as_Mat = UB.M;  // so we don't confuse the ?: operator below
    
    const bool use_copy = ((equilibrate && UB.is_const) || UB.is_alias(out));
    
    Mat<eT> B_tmp;  if(use_copy)  { B_tmp = UB_M_as_Mat; }
    
    const Mat<eT>& B = (use_copy) ? B_tmp : UB_M_as_Mat;
    
    arma_debug_check( (A.n_rows != B.n_rows), "solve(): number of rows in given matrices must be the same" );
      
    if(A.is_empty() || B.is_empty())  { out.zeros(A.n_rows, B.n_cols); return true; }
    
    arma_debug_assert_blas_size(A,B);
    
    out.set_size(A.n_rows, B.n_cols);
    
    char     fact  = (equilibrate) ? 'E' : 'N'; 
    char     trans = 'N';
    char     equed = char(0);
    blas_int n     = blas_int(A.n_rows);
    blas_int nrhs  = blas_int(B.n_cols);
    blas_int lda   = blas_int(A.n_rows);
    blas_int ldaf  = blas_int(A.n_rows);
    blas_int ldb   = blas_int(A.n_rows);
    blas_int ldx   = blas_int(A.n_rows);
    blas_int info  = blas_int(0);
    eT       rcond = eT(0);
    
    Mat<eT> AF(A.n_rows, A.n_rows, arma_nozeros_indicator());
    
    podarray<blas_int>  IPIV(  A.n_rows);
    podarray<eT>           R(  A.n_rows);
    podarray<eT>           C(  A.n_rows);
    podarray<eT>        FERR(  B.n_cols);
    podarray<eT>        BERR(  B.n_cols);
    podarray<eT>        WORK(4*A.n_rows);
    podarray<blas_int> IWORK(  A.n_rows);
    
    arma_extra_debug_print("lapack::gesvx()");
    lapack::gesvx
      (
      &fact, &trans, &n, &nrhs,
      A.memptr(), &lda,
      AF.memptr(), &ldaf,
      IPIV.memptr(),
      &equed,
      R.memptr(),
      C.memptr(),
      const_cast<eT*>(B.memptr()), &ldb,
      out.memptr(), &ldx,
      &rcond,
      FERR.memptr(),
      BERR.memptr(),
      WORK.memptr(),
      IWORK.memptr(),
      &info
      );
    
    // NOTE: using const_cast<eT*>(B.memptr()) to allow B to be overwritten for equilibration;
    // NOTE: B is created as a copy of B_expr if equilibration is enabled; otherwise B is a reference to B_expr
    
    out_rcond = rcond;
    
    return ((info == 0) || (info == (n+1)));
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_ignore(equilibrate);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a system of linear equations via LU decomposition with refinement (complex matrices)
template<typename T1>
inline
bool
auxlib::solve_square_refine(Mat< std::complex<typename T1::pod_type> >& out, typename T1::pod_type& out_rcond, Mat< std::complex<typename T1::pod_type> >& A, const Base<std::complex<typename T1::pod_type>,T1>& B_expr, const bool equilibrate)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type     T;
    typedef typename std::complex<T> eT;
    
    // Mat<eT> B = B_expr.get_ref();  // B is overwritten by lapack::cx_gesvx() if equilibrate is enabled
    
    quasi_unwrap<T1> UB(B_expr.get_ref());  // deliberately not declaring as const
    
    const Mat<eT>& UB_M_as_Mat = UB.M;  // so we don't confuse the ?: operator below
    
    const bool use_copy = ((equilibrate && UB.is_const) || UB.is_alias(out));
    
    Mat<eT> B_tmp;  if(use_copy)  { B_tmp = UB_M_as_Mat; }
    
    const Mat<eT>& B = (use_copy) ? B_tmp : UB_M_as_Mat;
    
    arma_debug_check( (A.n_rows != B.n_rows), "solve(): number of rows in given matrices must be the same" );
      
    if(A.is_empty() || B.is_empty())  { out.zeros(A.n_rows, B.n_cols); return true; }
    
    arma_debug_assert_blas_size(A,B);
    
    out.set_size(A.n_rows, B.n_cols);
    
    char     fact  = (equilibrate) ? 'E' : 'N'; 
    char     trans = 'N';
    char     equed = char(0);
    blas_int n     = blas_int(A.n_rows);
    blas_int nrhs  = blas_int(B.n_cols);
    blas_int lda   = blas_int(A.n_rows);
    blas_int ldaf  = blas_int(A.n_rows);
    blas_int ldb   = blas_int(A.n_rows);
    blas_int ldx   = blas_int(A.n_rows);
    blas_int info  = blas_int(0);
    T        rcond = T(0);
    
    Mat<eT> AF(A.n_rows, A.n_rows, arma_nozeros_indicator());
    
    podarray<blas_int>  IPIV(  A.n_rows);
    podarray< T>           R(  A.n_rows);
    podarray< T>           C(  A.n_rows);
    podarray< T>        FERR(  B.n_cols);
    podarray< T>        BERR(  B.n_cols);
    podarray<eT>        WORK(2*A.n_rows);
    podarray< T>       RWORK(2*A.n_rows);
    
    arma_extra_debug_print("lapack::cx_gesvx()");
    lapack::cx_gesvx
      (
      &fact, &trans, &n, &nrhs,
      A.memptr(), &lda,
      AF.memptr(), &ldaf,
      IPIV.memptr(),
      &equed,
      R.memptr(),
      C.memptr(),
      const_cast<eT*>(B.memptr()), &ldb,
      out.memptr(), &ldx,
      &rcond,
      FERR.memptr(),
      BERR.memptr(),
      WORK.memptr(),
      RWORK.memptr(),
      &info
      );
    
    // NOTE: using const_cast<eT*>(B.memptr()) to allow B to be overwritten for equilibration;
    // NOTE: B is created as a copy of B_expr if equilibration is enabled; otherwise B is a reference to B_expr
    
    out_rcond = rcond;
    
    return ((info == 0) || (info == (n+1)));
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_ignore(equilibrate);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T1>
inline
bool
auxlib::solve_sympd_fast(Mat<typename T1::elem_type>& out, Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_extra_debug_print("auxlib::solve_sympd_fast(): redirecting to auxlib::solve_square_fast() due to crippled LAPACK");
    
    return auxlib::solve_square_fast(out, A, B_expr);
    }
  #else
    {
    return auxlib::solve_sympd_fast_common(out, A, B_expr);
    }
  #endif
  }



template<typename T1>
inline
bool
auxlib::solve_sympd_fast_common(Mat<typename T1::elem_type>& out, Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  out = B_expr.get_ref();
  
  const uword B_n_rows = out.n_rows;
  const uword B_n_cols = out.n_cols;
  
  arma_debug_check( (A.n_rows != B_n_rows), "solve(): number of rows in given matrices must be the same", [&](){ out.soft_reset(); } );
  
  if(A.is_empty() || out.is_empty())  { out.zeros(A.n_cols, B_n_cols); return true; }
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(A, out);
    
    char     uplo = 'L';
    blas_int n    = blas_int(A.n_rows);  // assuming A is square
    blas_int nrhs = blas_int(B_n_cols);
    blas_int lda  = blas_int(A.n_rows);
    blas_int ldb  = blas_int(B_n_rows);
    blas_int info = blas_int(0);
    
    arma_extra_debug_print("lapack::posv()");
    lapack::posv<eT>(&uplo, &n, &nrhs, A.memptr(), &lda, out.memptr(), &ldb, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a system of linear equations via Cholesky decomposition with rcond estimate (real matrices)
template<typename T1>
inline
bool
auxlib::solve_sympd_rcond(Mat<typename T1::pod_type>& out, bool& out_sympd_state, typename T1::pod_type& out_rcond, Mat<typename T1::pod_type>& A, const Base<typename T1::pod_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::elem_type eT;
    typedef typename T1::pod_type   T;
    
    out_sympd_state = false;
    out_rcond       = T(0);
    
    out = B_expr.get_ref();
    
    const uword B_n_rows = out.n_rows;
    const uword B_n_cols = out.n_cols;
    
    arma_debug_check( (A.n_rows != B_n_rows), "solve(): number of rows in given matrices must be the same", [&](){ out.soft_reset(); } );
    
    if(A.is_empty() || out.is_empty())  { out.zeros(A.n_cols, B_n_cols); return true; }
    
    arma_debug_assert_blas_size(A, out);
    
    char     norm_id  = '1';
    char     uplo     = 'L';
    blas_int n        = blas_int(A.n_rows);  // assuming A is square
    blas_int nrhs     = blas_int(B_n_cols);
    blas_int info     = blas_int(0);
    T        norm_val = T(0);
    
    podarray<T> work(A.n_rows);
    
    arma_extra_debug_print("lapack::lansy()");
    norm_val = lapack::lansy(&norm_id, &uplo, &n, A.memptr(), &n, work.memptr());
    
    arma_extra_debug_print("lapack::potrf()");
    lapack::potrf<eT>(&uplo, &n, A.memptr(), &n, &info);
    
    if(info != 0)  { return false; }
    
    out_sympd_state = true;
    
    arma_extra_debug_print("lapack::potrs()");
    lapack::potrs<eT>(&uplo, &n, &nrhs, A.memptr(), &n, out.memptr(), &n, &info);
    
    if(info != 0)  { return false; }
    
    out_rcond = auxlib::lu_rcond_sympd<T>(A, norm_val);
    
    return true;
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(out_sympd_state);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a system of linear equations via Cholesky decomposition with rcond estimate (complex matrices)
template<typename T1>
inline
bool
auxlib::solve_sympd_rcond(Mat< std::complex<typename T1::pod_type> >& out, bool& out_sympd_state, typename T1::pod_type& out_rcond, Mat< std::complex<typename T1::pod_type> >& A, const Base< std::complex<typename T1::pod_type>,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_extra_debug_print("auxlib::solve_sympd_rcond(): redirecting to auxlib::solve_square_rcond() due to crippled LAPACK");
    
    out_sympd_state = false;
    
    return auxlib::solve_square_rcond(out, out_rcond, A, B_expr);
    }
  #elif defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type     T;
    typedef typename std::complex<T> eT;
    
    out_sympd_state = false;
    out_rcond       = T(0);
    
    out = B_expr.get_ref();
    
    const uword B_n_rows = out.n_rows;
    const uword B_n_cols = out.n_cols;
    
    arma_debug_check( (A.n_rows != B_n_rows), "solve(): number of rows in given matrices must be the same", [&](){ out.soft_reset(); } );
    
    if(A.is_empty() || out.is_empty())  { out.zeros(A.n_cols, B_n_cols); return true; }
    
    arma_debug_assert_blas_size(A, out);
    
    char     norm_id  = '1';
    char     uplo     = 'L';
    blas_int n        = blas_int(A.n_rows);  // assuming A is square
    blas_int nrhs     = blas_int(B_n_cols);
    blas_int info     = blas_int(0);
    T        norm_val = T(0);
    
    podarray<T> work(A.n_rows);
    
    arma_extra_debug_print("lapack::lanhe()");
    norm_val = lapack::lanhe(&norm_id, &uplo, &n, A.memptr(), &n, work.memptr());
    
    arma_extra_debug_print("lapack::potrf()");
    lapack::potrf<eT>(&uplo, &n, A.memptr(), &n, &info);
    
    if(info != 0)  { return false; }
    
    out_sympd_state = true;
    
    arma_extra_debug_print("lapack::potrs()");
    lapack::potrs<eT>(&uplo, &n, &nrhs, A.memptr(), &n, out.memptr(), &n, &info);
    
    if(info != 0)  { return false; }
    
    out_rcond = auxlib::lu_rcond_sympd<T>(A, norm_val);
    
    return true;
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(out_sympd_state);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a system of linear equations via Cholesky decomposition with refinement (real matrices)
template<typename T1>
inline
bool
auxlib::solve_sympd_refine(Mat<typename T1::pod_type>& out, typename T1::pod_type& out_rcond, Mat<typename T1::pod_type>& A, const Base<typename T1::pod_type,T1>& B_expr, const bool equilibrate)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type eT;
    
    // Mat<eT> B = B_expr.get_ref();  // B is overwritten by lapack::posvx() if equilibrate is enabled
    
    quasi_unwrap<T1> UB(B_expr.get_ref());  // deliberately not declaring as const
    
    const Mat<eT>& UB_M_as_Mat = UB.M;  // so we don't confuse the ?: operator below
    
    const bool use_copy = ((equilibrate && UB.is_const) || UB.is_alias(out));
    
    Mat<eT> B_tmp;  if(use_copy)  { B_tmp = UB_M_as_Mat; }
    
    const Mat<eT>& B = (use_copy) ? B_tmp : UB_M_as_Mat;
    
    arma_debug_check( (A.n_rows != B.n_rows), "solve(): number of rows in given matrices must be the same" );
    
    if(A.is_empty() || B.is_empty())  { out.zeros(A.n_rows, B.n_cols); return true; }
    
    arma_debug_assert_blas_size(A,B);
    
    out.set_size(A.n_rows, B.n_cols);
    
    char     fact  = (equilibrate) ? 'E' : 'N'; 
    char     uplo  = 'L';
    char     equed = char(0);
    blas_int n     = blas_int(A.n_rows);
    blas_int nrhs  = blas_int(B.n_cols);
    blas_int lda   = blas_int(A.n_rows);
    blas_int ldaf  = blas_int(A.n_rows);
    blas_int ldb   = blas_int(A.n_rows);
    blas_int ldx   = blas_int(A.n_rows);
    blas_int info  = blas_int(0);
    eT       rcond = eT(0);
    
    Mat<eT> AF(A.n_rows, A.n_rows, arma_nozeros_indicator());
    
    podarray<eT>           S(  A.n_rows);
    podarray<eT>        FERR(  B.n_cols);
    podarray<eT>        BERR(  B.n_cols);
    podarray<eT>        WORK(3*A.n_rows);
    podarray<blas_int> IWORK(  A.n_rows);
    
    arma_extra_debug_print("lapack::posvx()");
    lapack::posvx(&fact, &uplo, &n, &nrhs, A.memptr(), &lda, AF.memptr(), &ldaf, &equed, S.memptr(), const_cast<eT*>(B.memptr()), &ldb, out.memptr(), &ldx, &rcond, FERR.memptr(), BERR.memptr(), WORK.memptr(), IWORK.memptr(), &info);
    
    // NOTE: using const_cast<eT*>(B.memptr()) to allow B to be overwritten for equilibration;
    // NOTE: B is created as a copy of B_expr if equilibration is enabled; otherwise B is a reference to B_expr
    
    // NOTE: lapack::posvx() sets rcond to zero if A is not sympd
    out_rcond = rcond;
    
    return ((info == 0) || (info == (n+1)));
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_ignore(equilibrate);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a system of linear equations via Cholesky decomposition with refinement (complex matrices)
template<typename T1>
inline
bool
auxlib::solve_sympd_refine(Mat< std::complex<typename T1::pod_type> >& out, typename T1::pod_type& out_rcond, Mat< std::complex<typename T1::pod_type> >& A, const Base<std::complex<typename T1::pod_type>,T1>& B_expr, const bool equilibrate)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_extra_debug_print("auxlib::solve_sympd_refine(): redirecting to auxlib::solve_square_refine() due to crippled LAPACK");
    
    return auxlib::solve_square_refine(out, out_rcond, A, B_expr, equilibrate);
    }
  #elif defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type     T;
    typedef typename std::complex<T> eT;
    
    // Mat<eT> B = B_expr.get_ref();  // B is overwritten by lapack::cx_posvx() if equilibrate is enabled
    
    quasi_unwrap<T1> UB(B_expr.get_ref());  // deliberately not declaring as const
    
    const Mat<eT>& UB_M_as_Mat = UB.M;  // so we don't confuse the ?: operator below
    
    const bool use_copy = ((equilibrate && UB.is_const) || UB.is_alias(out));
    
    Mat<eT> B_tmp;  if(use_copy)  { B_tmp = UB_M_as_Mat; }
    
    const Mat<eT>& B = (use_copy) ? B_tmp : UB_M_as_Mat;
    
    arma_debug_check( (A.n_rows != B.n_rows), "solve(): number of rows in given matrices must be the same" );
      
    if(A.is_empty() || B.is_empty())  { out.zeros(A.n_rows, B.n_cols); return true; }
    
    arma_debug_assert_blas_size(A,B);
    
    out.set_size(A.n_rows, B.n_cols);
    
    char     fact  = (equilibrate) ? 'E' : 'N'; 
    char     uplo  = 'L';
    char     equed = char(0);
    blas_int n     = blas_int(A.n_rows);
    blas_int nrhs  = blas_int(B.n_cols);
    blas_int lda   = blas_int(A.n_rows);
    blas_int ldaf  = blas_int(A.n_rows);
    blas_int ldb   = blas_int(A.n_rows);
    blas_int ldx   = blas_int(A.n_rows);
    blas_int info  = blas_int(0);
    T        rcond = T(0);
    
    Mat<eT> AF(A.n_rows, A.n_rows, arma_nozeros_indicator());
    
    podarray< T>           S(  A.n_rows);
    podarray< T>        FERR(  B.n_cols);
    podarray< T>        BERR(  B.n_cols);
    podarray<eT>        WORK(2*A.n_rows);
    podarray< T>       RWORK(  A.n_rows);
    
    arma_extra_debug_print("lapack::cx_posvx()");
    lapack::cx_posvx(&fact, &uplo, &n, &nrhs, A.memptr(), &lda, AF.memptr(), &ldaf, &equed, S.memptr(), const_cast<eT*>(B.memptr()), &ldb, out.memptr(), &ldx, &rcond, FERR.memptr(), BERR.memptr(), WORK.memptr(), RWORK.memptr(), &info);
    
    // NOTE: using const_cast<eT*>(B.memptr()) to allow B to be overwritten for equilibration;
    // NOTE: B is created as a copy of B_expr if equilibration is enabled; otherwise B is a reference to B_expr
    
    // NOTE: lapack::cx_posvx() sets rcond to zero if A is not sympd
    out_rcond = rcond;
    
    return ((info == 0) || (info == (n+1)));
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_ignore(equilibrate);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a non-square full-rank system via QR or LQ decomposition
template<typename T1>
inline
bool
auxlib::solve_rect_fast(Mat<typename T1::elem_type>& out, Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::elem_type eT;
    
    const unwrap<T1>   U(B_expr.get_ref());
    const Mat<eT>& B = U.M;
    
    arma_debug_check( (A.n_rows != B.n_rows), "solve(): number of rows in given matrices must be the same" );
    
    if(A.is_empty() || B.is_empty())  { out.zeros(A.n_cols, B.n_cols); return true; }
    
    arma_debug_assert_blas_size(A,B);
    
    Mat<eT> tmp( (std::max)(A.n_rows, A.n_cols), B.n_cols, arma_nozeros_indicator() );
    
    if(arma::size(tmp) == arma::size(B))
      {
      tmp = B;
      }
    else
      {
      tmp.zeros();
      tmp(0,0, arma::size(B)) = B;
      }
    
    char      trans     = 'N';
    blas_int  m         = blas_int(A.n_rows);
    blas_int  n         = blas_int(A.n_cols);
    blas_int  lda       = blas_int(A.n_rows);
    blas_int  ldb       = blas_int(tmp.n_rows);
    blas_int  nrhs      = blas_int(B.n_cols);
    blas_int  min_mn    = (std::min)(m,n);
    blas_int  lwork_min = (std::max)(blas_int(1), min_mn + (std::max)(min_mn, nrhs));
    blas_int  info      = 0;
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= ((is_cx<eT>::yes) ? uword(256) : uword(1024)))
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = -1;
      
      arma_extra_debug_print("lapack::gels()");
      lapack::gels<eT>( &trans, &m, &n, &nrhs, A.memptr(), &lda, tmp.memptr(), &ldb, &work_query[0], &lwork_query, &info );
      
      if(info != 0)  { return false; }
      
      lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::gels()");
    lapack::gels<eT>( &trans, &m, &n, &nrhs, A.memptr(), &lda, tmp.memptr(), &ldb, work.memptr(), &lwork_final, &info );
    
    if(info != 0)  { return false; }
    
    if(tmp.n_rows == A.n_cols)
      {
      out.steal_mem(tmp);
      }
    else
      {
      out = tmp.head_rows(A.n_cols);
      }
    
    return true;
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a non-square full-rank system via QR or LQ decomposition with rcond estimate (experimental)
template<typename T1>
inline
bool
auxlib::solve_rect_rcond(Mat<typename T1::elem_type>& out, typename T1::pod_type& out_rcond, Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::elem_type eT;
    typedef typename T1::pod_type   T;
    
    out_rcond = T(0);
    
    const unwrap<T1>   U(B_expr.get_ref());
    const Mat<eT>& B = U.M;
    
    arma_debug_check( (A.n_rows != B.n_rows), "solve(): number of rows in given matrices must be the same" );
    
    if(A.is_empty() || B.is_empty())  { out.zeros(A.n_cols, B.n_cols); return true; }
    
    arma_debug_assert_blas_size(A,B);
    
    Mat<eT> tmp( (std::max)(A.n_rows, A.n_cols), B.n_cols, arma_nozeros_indicator() );
    
    if(arma::size(tmp) == arma::size(B))
      {
      tmp = B;
      }
    else
      {
      tmp.zeros();
      tmp(0,0, arma::size(B)) = B;
      }
    
    char      trans     = 'N';
    blas_int  m         = blas_int(A.n_rows);
    blas_int  n         = blas_int(A.n_cols);
    blas_int  lda       = blas_int(A.n_rows);
    blas_int  ldb       = blas_int(tmp.n_rows);
    blas_int  nrhs      = blas_int(B.n_cols);
    blas_int  min_mn    = (std::min)(m,n);
    blas_int  lwork_min = (std::max)(blas_int(1), min_mn + (std::max)(min_mn, nrhs));
    blas_int  info      = 0;
    
    blas_int lwork_proposed = 0;
    
    if(A.n_elem >= ((is_cx<eT>::yes) ? uword(256) : uword(1024)))
      {
      eT        work_query[2] = {};
      blas_int lwork_query    = -1;
      
      arma_extra_debug_print("lapack::gels()");
      lapack::gels<eT>( &trans, &m, &n, &nrhs, A.memptr(), &lda, tmp.memptr(), &ldb, &work_query[0], &lwork_query, &info );
      
      if(info != 0)  { return false; }
      
      lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
      }
    
    blas_int lwork_final = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::gels()");
    lapack::gels<eT>( &trans, &m, &n, &nrhs, A.memptr(), &lda, tmp.memptr(), &ldb, work.memptr(), &lwork_final, &info );
    
    if(info != 0)  { return false; }
    
    if(A.n_rows >= A.n_cols)
      {
      arma_extra_debug_print("estimating rcond via R");
      
      // xGELS  docs: for M >= N, A contains details of its QR decomposition as returned by xGEQRF
      // xGEQRF docs: elements on and above the diagonal contain the min(M,N)-by-N upper trapezoidal matrix R
      
      Mat<eT> R(A.n_cols, A.n_cols, arma_zeros_indicator());
      
      for(uword col=0; col < A.n_cols; ++col)
        {
        for(uword row=0; row <= col; ++row)
          {
          R.at(row,col) = A.at(row,col);
          }
        }
      
      // determine quality of solution
      out_rcond = auxlib::rcond_trimat(R, 0);   // 0: upper triangular; 1: lower triangular
      }
    else
    if(A.n_rows < A.n_cols)
      {
      arma_extra_debug_print("estimating rcond via L");
      
      // xGELS  docs: for M < N, A contains details of its LQ decomposition as returned by xGELQF
      // xGELQF docs: elements on and below the diagonal contain the M-by-min(M,N) lower trapezoidal matrix L
      
      Mat<eT> L(A.n_rows, A.n_rows, arma_zeros_indicator());
      
      for(uword col=0; col < A.n_rows; ++col)
        {
        for(uword row=col; row < A.n_rows; ++row)
          {
          L.at(row,col) = A.at(row,col);
          }
        }
      
      // determine quality of solution
      out_rcond = auxlib::rcond_trimat(L, 1);   // 0: upper triangular; 1: lower triangular
      }
    
    if(tmp.n_rows == A.n_cols)
      {
      out.steal_mem(tmp);
      }
    else
      {
      out = tmp.head_rows(A.n_cols);
      }
    
    return true;
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T1>
inline
bool
auxlib::solve_approx_svd(Mat<typename T1::pod_type>& out, Mat<typename T1::pod_type>& A, const Base<typename T1::pod_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type eT;
    
    const unwrap<T1>   U(B_expr.get_ref());
    const Mat<eT>& B = U.M;
    
    arma_debug_check( (A.n_rows != B.n_rows), "solve(): number of rows in given matrices must be the same" );
    
    if(A.is_empty() || B.is_empty())  { out.zeros(A.n_cols, B.n_cols); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    if(arma_config::check_nonfinite && B.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A,B);
    
    Mat<eT> tmp( (std::max)(A.n_rows, A.n_cols), B.n_cols, arma_nozeros_indicator() );
    
    if(arma::size(tmp) == arma::size(B))
      {
      tmp = B;
      }
    else
      {
      tmp.zeros();
      tmp(0,0, arma::size(B)) = B;
      }
    
    blas_int m      = blas_int(A.n_rows);
    blas_int n      = blas_int(A.n_cols);
    blas_int min_mn = (std::min)(m, n);
    blas_int nrhs   = blas_int(B.n_cols);
    blas_int lda    = blas_int(A.n_rows);
    blas_int ldb    = blas_int(tmp.n_rows);
  //eT       rcond  = eT(-1);  // -1 means "use machine precision"
    eT       rcond  = (std::max)(A.n_rows, A.n_cols) * std::numeric_limits<eT>::epsilon();
    blas_int rank   = blas_int(0);
    blas_int info   = blas_int(0);
    
    podarray<eT> S( static_cast<uword>(min_mn) );
    
    // NOTE: with LAPACK 3.8, can use the workspace query to also obtain liwork,
    // NOTE: which makes the call to lapack::laenv() redundant
    
    blas_int ispec = blas_int(9);
    
    const char* const_name = (is_float<eT>::value) ? "SGELSD" : "DGELSD";
    const char* const_opts = " ";
    
    char* name = const_cast<char*>(const_name);
    char* opts = const_cast<char*>(const_opts);
    
    blas_int n1 = m;
    blas_int n2 = n;
    blas_int n3 = nrhs;
    blas_int n4 = lda;
    
    blas_int laenv_result = (arma_config::hidden_args) ? blas_int(lapack::laenv(&ispec, name, opts, &n1, &n2, &n3, &n4, 6, 1)) : blas_int(0);
    
    blas_int smlsiz    = (std::max)( blas_int(25), laenv_result );
    blas_int smlsiz_p1 = blas_int(1) + smlsiz;
    
    blas_int nlvl   = (std::max)( blas_int(0), blas_int(1) + blas_int( std::log(double(min_mn) / double(smlsiz_p1))/double(0.69314718055994530942) ) );
    blas_int liwork = (std::max)( blas_int(1), (blas_int(3)*min_mn*nlvl + blas_int(11)*min_mn) );
    
    podarray<blas_int> iwork( static_cast<uword>(liwork) );
    
    blas_int lwork_min = blas_int(12)*min_mn + blas_int(2)*min_mn*smlsiz + blas_int(8)*min_mn*nlvl + min_mn*nrhs + smlsiz_p1*smlsiz_p1;
    
    eT        work_query[2] = {};
    blas_int lwork_query    = blas_int(-1);
    
    arma_extra_debug_print("lapack::gelsd()");
    lapack::gelsd(&m, &n, &nrhs, A.memptr(), &lda, tmp.memptr(), &ldb, S.memptr(), &rcond, &rank, &work_query[0], &lwork_query, iwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    // NOTE: in LAPACK 3.8, iwork[0] returns the minimum liwork
    
    blas_int lwork_proposed = static_cast<blas_int>( access::tmp_real(work_query[0]) );
    blas_int lwork_final    = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::gelsd()");
    lapack::gelsd(&m, &n, &nrhs, A.memptr(), &lda, tmp.memptr(), &ldb, S.memptr(), &rcond, &rank, work.memptr(), &lwork_final, iwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    if(tmp.n_rows == A.n_cols)
      {
      out.steal_mem(tmp);
      }
    else
      {
      out = tmp.head_rows(A.n_cols);
      }
    
    return true;
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T1>
inline
bool
auxlib::solve_approx_svd(Mat< std::complex<typename T1::pod_type> >& out, Mat< std::complex<typename T1::pod_type> >& A, const Base<std::complex<typename T1::pod_type>,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type     T;
    typedef typename std::complex<T> eT;
    
    const unwrap<T1>   U(B_expr.get_ref());
    const Mat<eT>& B = U.M;
    
    arma_debug_check( (A.n_rows != B.n_rows), "solve(): number of rows in given matrices must be the same" );
    
    if(A.is_empty() || B.is_empty())  { out.zeros(A.n_cols, B.n_cols); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    if(arma_config::check_nonfinite && B.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A,B);
    
    Mat<eT> tmp( (std::max)(A.n_rows, A.n_cols), B.n_cols, arma_nozeros_indicator() );
    
    if(arma::size(tmp) == arma::size(B))
      {
      tmp = B;
      }
    else
      {
      tmp.zeros();
      tmp(0,0, arma::size(B)) = B;
      }
    
    blas_int m      = blas_int(A.n_rows);
    blas_int n      = blas_int(A.n_cols);
    blas_int min_mn = (std::min)(m, n);
    blas_int nrhs   = blas_int(B.n_cols);
    blas_int lda    = blas_int(A.n_rows);
    blas_int ldb    = blas_int(tmp.n_rows);
  //T        rcond  = T(-1);  // -1 means "use machine precision"
    T        rcond  = (std::max)(A.n_rows, A.n_cols) * std::numeric_limits<T>::epsilon();
    blas_int rank   = blas_int(0);
    blas_int info   = blas_int(0);
    
    podarray<T> S( static_cast<uword>(min_mn) );
    
    blas_int ispec = blas_int(9);
    
    const char* const_name = (is_float<T>::value) ? "CGELSD" : "ZGELSD";
    const char* const_opts = " ";
    
    char* name = const_cast<char*>(const_name);
    char* opts = const_cast<char*>(const_opts);
    
    blas_int n1 = m;
    blas_int n2 = n;
    blas_int n3 = nrhs;
    blas_int n4 = lda;
    
    blas_int laenv_result = (arma_config::hidden_args) ? blas_int(lapack::laenv(&ispec, name, opts, &n1, &n2, &n3, &n4, 6, 1)) : blas_int(0);
    
    blas_int smlsiz    = (std::max)( blas_int(25), laenv_result );
    blas_int smlsiz_p1 = blas_int(1) + smlsiz;
    
    blas_int nlvl = (std::max)( blas_int(0), blas_int(1) + blas_int( std::log(double(min_mn) / double(smlsiz_p1))/double(0.69314718055994530942) ) );
    
    blas_int lrwork = (m >= n)
      ? blas_int(10)*n + blas_int(2)*n*smlsiz + blas_int(8)*n*nlvl + blas_int(3)*smlsiz*nrhs + (std::max)( (smlsiz_p1)*(smlsiz_p1), n*(blas_int(1)+nrhs) + blas_int(2)*nrhs )
      : blas_int(10)*m + blas_int(2)*m*smlsiz + blas_int(8)*m*nlvl + blas_int(3)*smlsiz*nrhs + (std::max)( (smlsiz_p1)*(smlsiz_p1), n*(blas_int(1)+nrhs) + blas_int(2)*nrhs );
    
    blas_int liwork = (std::max)( blas_int(1), (blas_int(3)*blas_int(min_mn)*nlvl + blas_int(11)*blas_int(min_mn)) );
    
    podarray<T>        rwork( static_cast<uword>(lrwork) );
    podarray<blas_int> iwork( static_cast<uword>(liwork) );
    
    blas_int lwork_min = 2*min_mn + min_mn*nrhs;
    
    eT        work_query[2] = {};
    blas_int lwork_query    = blas_int(-1);
    
    arma_extra_debug_print("lapack::cx_gelsd()");
    lapack::cx_gelsd(&m, &n, &nrhs, A.memptr(), &lda, tmp.memptr(), &ldb, S.memptr(), &rcond, &rank, &work_query[0], &lwork_query, rwork.memptr(), iwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    blas_int lwork_proposed = static_cast<blas_int>( access::tmp_real( work_query[0]) );
    blas_int lwork_final    = (std::max)(lwork_proposed, lwork_min);
    
    podarray<eT> work( static_cast<uword>(lwork_final) );
    
    arma_extra_debug_print("lapack::cx_gelsd()");
    lapack::cx_gelsd(&m, &n, &nrhs, A.memptr(), &lda, tmp.memptr(), &ldb, S.memptr(), &rcond, &rank, work.memptr(), &lwork_final, rwork.memptr(), iwork.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    if(tmp.n_rows == A.n_cols)
      {
      out.steal_mem(tmp);
      }
    else
      {
      out = tmp.head_rows(A.n_cols);
      }
    
    return true;
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T1>
inline
bool
auxlib::solve_trimat_fast(Mat<typename T1::elem_type>& out, const Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr, const uword layout)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    out = B_expr.get_ref();
    
    const uword B_n_rows = out.n_rows;
    const uword B_n_cols = out.n_cols;
    
    arma_debug_check( (A.n_rows != B_n_rows), "solve(): number of rows in given matrices must be the same", [&](){ out.soft_reset(); } );
    
    if(A.is_empty() || out.is_empty())  { out.zeros(A.n_cols, B_n_cols); return true; }
    
    arma_debug_assert_blas_size(A,out);
    
    char     uplo  = (layout == 0) ? 'U' : 'L';
    char     trans = 'N';
    char     diag  = 'N';
    blas_int n     = blas_int(A.n_rows);
    blas_int nrhs  = blas_int(B_n_cols);
    blas_int info  = 0;
    
    arma_extra_debug_print("lapack::trtrs()");
    lapack::trtrs(&uplo, &trans, &diag, &n, &nrhs, A.memptr(), &n, out.memptr(), &n, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_ignore(layout);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T1>
inline
bool
auxlib::solve_trimat_rcond(Mat<typename T1::elem_type>& out, typename T1::pod_type& out_rcond, const Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr, const uword layout)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type T;
    
    out_rcond = T(0);
    
    out = B_expr.get_ref();
    
    const uword B_n_rows = out.n_rows;
    const uword B_n_cols = out.n_cols;
    
    arma_debug_check( (A.n_rows != B_n_rows), "solve(): number of rows in given matrices must be the same", [&](){ out.soft_reset(); } );
    
    if(A.is_empty() || out.is_empty())  { out.zeros(A.n_cols, B_n_cols); return true; }
    
    arma_debug_assert_blas_size(A,out);
    
    char     uplo  = (layout == 0) ? 'U' : 'L';
    char     trans = 'N';
    char     diag  = 'N';
    blas_int n     = blas_int(A.n_rows);
    blas_int nrhs  = blas_int(B_n_cols);
    blas_int info  = 0;
    
    arma_extra_debug_print("lapack::trtrs()");
    lapack::trtrs(&uplo, &trans, &diag, &n, &nrhs, A.memptr(), &n, out.memptr(), &n, &info);
    
    if(info != 0)  { return false; }
    
    // determine quality of solution
    out_rcond = auxlib::rcond_trimat(A, layout);
    
    return true;
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_ignore(layout);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a system of linear equations via LU decomposition (real band matrix)
template<typename T1>
inline
bool
auxlib::solve_band_fast(Mat<typename T1::pod_type>& out, Mat<typename T1::pod_type>& A, const uword KL, const uword KU, const Base<typename T1::pod_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  return auxlib::solve_band_fast_common(out, A, KL, KU, B_expr);
  }



//! solve a system of linear equations via LU decomposition (complex band matrix)
template<typename T1>
inline
bool
auxlib::solve_band_fast(Mat< std::complex<typename T1::pod_type> >& out, Mat< std::complex<typename T1::pod_type> >& A, const uword KL, const uword KU, const Base< std::complex<typename T1::pod_type>,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_extra_debug_print("auxlib::solve_band_fast(): redirecting to auxlib::solve_square_fast() due to crippled LAPACK");
    
    arma_ignore(KL);
    arma_ignore(KU);
    
    return auxlib::solve_square_fast(out, A, B_expr);
    }
  #else
    {
    return auxlib::solve_band_fast_common(out, A, KL, KU, B_expr);
    }
  #endif
  }



//! solve a system of linear equations via LU decomposition (band matrix)
template<typename T1>
inline
bool
auxlib::solve_band_fast_common(Mat<typename T1::elem_type>& out, const Mat<typename T1::elem_type>& A, const uword KL, const uword KU, const Base<typename T1::elem_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::elem_type eT;
    
    out = B_expr.get_ref();
    
    const uword B_n_rows = out.n_rows;
    const uword B_n_cols = out.n_cols;
    
    arma_debug_check( (A.n_rows != B_n_rows), "solve(): number of rows in given matrices must be the same", [&](){ out.soft_reset(); } );
    
    if(A.is_empty() || out.is_empty())  { out.zeros(A.n_rows, B_n_cols); return true; }
    
    // for gbsv, matrix AB size: 2*KL+KU+1 x N; band representation of A stored in rows KL+1 to 2*KL+KU+1  (note: fortran counts from 1)
    
    Mat<eT> AB;
    band_helper::compress(AB, A, KL, KU, true);
    
    const uword N = AB.n_cols;  // order of the original square matrix A
    
    arma_debug_assert_blas_size(AB,out);
    
    blas_int n    = blas_int(N);
    blas_int kl   = blas_int(KL);
    blas_int ku   = blas_int(KU);
    blas_int nrhs = blas_int(B_n_cols);
    blas_int ldab = blas_int(AB.n_rows);
    blas_int ldb  = blas_int(B_n_rows);
    blas_int info = blas_int(0);
    
    podarray<blas_int> ipiv(N + 2);  // +2 for paranoia
    
    // NOTE: AB is overwritten
    
    arma_extra_debug_print("lapack::gbsv()");
    lapack::gbsv<eT>(&n, &kl, &ku, &nrhs, AB.memptr(), &ldab, ipiv.memptr(), out.memptr(), &ldb, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(A);
    arma_ignore(KL);
    arma_ignore(KU);
    arma_ignore(B_expr);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a system of linear equations via LU decomposition (real band matrix)
template<typename T1>
inline
bool
auxlib::solve_band_rcond(Mat<typename T1::pod_type>& out, typename T1::pod_type& out_rcond, Mat<typename T1::pod_type>& A, const uword KL, const uword KU, const Base<typename T1::pod_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  return auxlib::solve_band_rcond_common(out, out_rcond, A, KL, KU, B_expr);
  }



//! solve a system of linear equations via LU decomposition (complex band matrix)
template<typename T1>
inline
bool
auxlib::solve_band_rcond(Mat< std::complex<typename T1::pod_type> >& out, typename T1::pod_type& out_rcond, Mat< std::complex<typename T1::pod_type> >& A, const uword KL, const uword KU, const Base< std::complex<typename T1::pod_type>,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_extra_debug_print("auxlib::solve_band_rcond(): redirecting to auxlib::solve_square_rcond() due to crippled LAPACK");
    
    arma_ignore(KL);
    arma_ignore(KU);
    
    return auxlib::solve_square_rcond(out, out_rcond, A, B_expr);
    }
  #else
    {
    return auxlib::solve_band_rcond_common(out, out_rcond, A, KL, KU, B_expr);
    }
  #endif
  }



//! solve a system of linear equations via LU decomposition (band matrix)
template<typename T1>
inline
bool
auxlib::solve_band_rcond_common(Mat<typename T1::elem_type>& out, typename T1::pod_type& out_rcond, const Mat<typename T1::elem_type>& A, const uword KL, const uword KU, const Base<typename T1::elem_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::elem_type eT;
    typedef typename T1::pod_type   T;
    
    out_rcond = T(0);
    
    out = B_expr.get_ref();
    
    const uword B_n_rows = out.n_rows;
    const uword B_n_cols = out.n_cols;
    
    arma_debug_check( (A.n_rows != B_n_rows), "solve(): number of rows in given matrices must be the same", [&](){ out.soft_reset(); } );
    
    if(A.is_empty() || out.is_empty())  { out.zeros(A.n_rows, B_n_cols); return true; }
    
    // for gbtrf, matrix AB size: 2*KL+KU+1 x N; band representation of A stored in rows KL+1 to 2*KL+KU+1  (note: fortran counts from 1)
    
    Mat<eT> AB;
    band_helper::compress(AB, A, KL, KU, true);
    
    const uword N = AB.n_cols;  // order of the original square matrix A
    
    arma_debug_assert_blas_size(AB,out);
    
    char     norm_id  = '1';
    char     trans    = 'N';
    blas_int n        = blas_int(N);  // assuming square matrix
    blas_int kl       = blas_int(KL);
    blas_int ku       = blas_int(KU);
    blas_int nrhs     = blas_int(B_n_cols);
    blas_int ldab     = blas_int(AB.n_rows);
    blas_int ldb      = blas_int(B_n_rows);
    blas_int info     = blas_int(0);
    T        norm_val = T(0);
    
    podarray<T>        junk(1);
    podarray<blas_int> ipiv(N + 2);  // +2 for paranoia
    
    arma_extra_debug_print("lapack::langb()");
    norm_val = lapack::langb<eT>(&norm_id, &n, &kl, &ku, AB.memptr(), &ldab, junk.memptr());
    
    arma_extra_debug_print("lapack::gbtrf()");
    lapack::gbtrf<eT>(&n, &n, &kl, &ku, AB.memptr(), &ldab, ipiv.memptr(), &info);
    
    if(info != 0)  { return false; }
    
    arma_extra_debug_print("lapack::gbtrs()");
    lapack::gbtrs<eT>(&trans, &n, &kl, &ku, &nrhs, AB.memptr(), &ldab, ipiv.memptr(), out.memptr(), &ldb, &info);
    
    if(info != 0)  { return false; }
    
    out_rcond = auxlib::lu_rcond_band<T>(AB, KL, KU, ipiv, norm_val);
    
    return true;
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_ignore(KL);
    arma_ignore(KU);
    arma_ignore(B_expr);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a system of linear equations via LU decomposition with refinement (real band matrices)
template<typename T1>
inline
bool
auxlib::solve_band_refine(Mat<typename T1::pod_type>& out, typename T1::pod_type& out_rcond, Mat<typename T1::pod_type>& A, const uword KL, const uword KU, const Base<typename T1::pod_type,T1>& B_expr, const bool equilibrate)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type eT;
    
    Mat<eT> B = B_expr.get_ref();  // B is overwritten
    
    arma_debug_check( (A.n_rows != B.n_rows), "solve(): number of rows in given matrices must be the same" );
      
    if(A.is_empty() || B.is_empty())  { out.zeros(A.n_rows, B.n_cols); return true; }
    
    // for gbsvx, matrix AB size: KL+KU+1 x N; band representation of A stored in rows 1 to KL+KU+1  (note: fortran counts from 1)
    
    Mat<eT> AB;
    band_helper::compress(AB, A, KL, KU, false);
    
    const uword N = AB.n_cols;
    
    arma_debug_assert_blas_size(AB,B);
    
    out.set_size(N, B.n_cols);
    
    Mat<eT> AFB(2*KL+KU+1, N, arma_nozeros_indicator());
    
    char     fact  = (equilibrate) ? 'E' : 'N'; 
    char     trans = 'N';
    char     equed = char(0);
    blas_int n     = blas_int(N);
    blas_int kl    = blas_int(KL);
    blas_int ku    = blas_int(KU);
    blas_int nrhs  = blas_int(B.n_cols);
    blas_int ldab  = blas_int(AB.n_rows);
    blas_int ldafb = blas_int(AFB.n_rows);
    blas_int ldb   = blas_int(B.n_rows);
    blas_int ldx   = blas_int(N);
    blas_int info  = blas_int(0);
    eT       rcond = eT(0);
    
    podarray<blas_int>  IPIV(  N);
    podarray<eT>           R(  N);
    podarray<eT>           C(  N);
    podarray<eT>        FERR(  B.n_cols);
    podarray<eT>        BERR(  B.n_cols);
    podarray<eT>        WORK(3*N);
    podarray<blas_int> IWORK(  N);
    
    arma_extra_debug_print("lapack::gbsvx()");
    lapack::gbsvx
      (
      &fact, &trans, &n, &kl, &ku, &nrhs, 
      AB.memptr(), &ldab,
      AFB.memptr(), &ldafb,
      IPIV.memptr(),
      &equed,
      R.memptr(),
      C.memptr(),
      B.memptr(), &ldb,
      out.memptr(), &ldx,
      &rcond,
      FERR.memptr(),
      BERR.memptr(),
      WORK.memptr(),
      IWORK.memptr(),
      &info
      );
    
    out_rcond = rcond;
    
    return ((info == 0) || (info == (n+1)));
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_ignore(KL);
    arma_ignore(KU);
    arma_ignore(B_expr);
    arma_ignore(equilibrate);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a system of linear equations via LU decomposition with refinement (complex band matrices)
template<typename T1>
inline
bool
auxlib::solve_band_refine(Mat< std::complex<typename T1::pod_type> >& out, typename T1::pod_type& out_rcond, Mat< std::complex<typename T1::pod_type> >& A, const uword KL, const uword KU, const Base<std::complex<typename T1::pod_type>,T1>& B_expr, const bool equilibrate)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_extra_debug_print("auxlib::solve_band_refine(): redirecting to auxlib::solve_square_refine() due to crippled LAPACK");
    
    arma_ignore(KL);
    arma_ignore(KU);
    
    return auxlib::solve_square_refine(out, out_rcond, A, B_expr, equilibrate);
    }
  #elif defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::pod_type     T;
    typedef typename std::complex<T> eT;
    
    Mat<eT> B = B_expr.get_ref();  // B is overwritten
    
    arma_debug_check( (A.n_rows != B.n_rows), "solve(): number of rows in given matrices must be the same" );
      
    if(A.is_empty() || B.is_empty())  { out.zeros(A.n_rows, B.n_cols); return true; }
    
    // for gbsvx, matrix AB size: KL+KU+1 x N; band representation of A stored in rows 1 to KL+KU+1  (note: fortran counts from 1)
    
    Mat<eT> AB;
    band_helper::compress(AB, A, KL, KU, false);
    
    const uword N = AB.n_cols;
    
    arma_debug_assert_blas_size(AB,B);
    
    out.set_size(N, B.n_cols);
    
    Mat<eT> AFB(2*KL+KU+1, N, arma_nozeros_indicator());
    
    char     fact  = (equilibrate) ? 'E' : 'N'; 
    char     trans = 'N';
    char     equed = char(0);
    blas_int n     = blas_int(N);
    blas_int kl    = blas_int(KL);
    blas_int ku    = blas_int(KU);
    blas_int nrhs  = blas_int(B.n_cols);
    blas_int ldab  = blas_int(AB.n_rows);
    blas_int ldafb = blas_int(AFB.n_rows);
    blas_int ldb   = blas_int(B.n_rows);
    blas_int ldx   = blas_int(N);
    blas_int info  = blas_int(0);
    T        rcond = T(0);
    
    podarray<blas_int>  IPIV(  N);
    podarray< T>           R(  N);
    podarray< T>           C(  N);
    podarray< T>        FERR(  B.n_cols);
    podarray< T>        BERR(  B.n_cols);
    podarray<eT>        WORK(2*N);
    podarray< T>       RWORK(  N);  // NOTE: according to lapack 3.6.1 docs, the size of RWORK in zgbsvx is different to RWORK in dgesvx 
    
    arma_extra_debug_print("lapack::cx_gbsvx()");
    lapack::cx_gbsvx
      (
      &fact, &trans, &n, &kl, &ku, &nrhs,
      AB.memptr(), &ldab,
      AFB.memptr(), &ldafb,
      IPIV.memptr(),
      &equed,
      R.memptr(),
      C.memptr(),
      B.memptr(), &ldb,
      out.memptr(), &ldx,
      &rcond,
      FERR.memptr(),
      BERR.memptr(),
      WORK.memptr(),
      RWORK.memptr(),
      &info
      );
    
    out_rcond = rcond;
    
    return ((info == 0) || (info == (n+1)));
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(out_rcond);
    arma_ignore(A);
    arma_ignore(KL);
    arma_ignore(KU);
    arma_ignore(B_expr);
    arma_ignore(equilibrate);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//! solve a system of linear equations via Gaussian elimination with partial pivoting (real tridiagonal band matrix)
template<typename T1>
inline
bool
auxlib::solve_tridiag_fast(Mat<typename T1::pod_type>& out, Mat<typename T1::pod_type>& A, const Base<typename T1::pod_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  return auxlib::solve_tridiag_fast_common(out, A, B_expr);
  }



//! solve a system of linear equations via Gaussian elimination with partial pivoting (complex tridiagonal band matrix)
template<typename T1>
inline
bool
auxlib::solve_tridiag_fast(Mat< std::complex<typename T1::pod_type> >& out, Mat< std::complex<typename T1::pod_type> >& A, const Base< std::complex<typename T1::pod_type>,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_extra_debug_print("auxlib::solve_tridiag_fast(): redirecting to auxlib::solve_square_fast() due to crippled LAPACK");
    
    return auxlib::solve_square_fast(out, A, B_expr);
    }
  #else
    {
    return auxlib::solve_tridiag_fast_common(out, A, B_expr);
    }
  #endif
  }



//! solve a system of linear equations via Gaussian elimination with partial pivoting (tridiagonal band matrix)
template<typename T1>
inline
bool
auxlib::solve_tridiag_fast_common(Mat<typename T1::elem_type>& out, const Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename T1::elem_type eT;
    
    out = B_expr.get_ref();
    
    const uword B_n_rows = out.n_rows;
    const uword B_n_cols = out.n_cols;
    
    arma_debug_check( (A.n_rows != B_n_rows), "solve(): number of rows in given matrices must be the same", [&](){ out.soft_reset(); } );
    
    if(A.is_empty() || out.is_empty())  { out.zeros(A.n_rows, B_n_cols); return true; }
    
    Mat<eT> tridiag;
    band_helper::extract_tridiag(tridiag, A);
    
    arma_debug_assert_blas_size(tridiag, out);
    
    blas_int n    = blas_int(A.n_rows);
    blas_int nrhs = blas_int(B_n_cols);
    blas_int ldb  = blas_int(B_n_rows);
    blas_int info = blas_int(0);
    
    arma_extra_debug_print("lapack::gtsv()");
    lapack::gtsv<eT>(&n, &nrhs, tridiag.colptr(0), tridiag.colptr(1), tridiag.colptr(2), out.memptr(), &ldb, &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(out);
    arma_ignore(A);
    arma_ignore(B_expr);
    arma_stop_logic_error("solve(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//
// Schur decomposition

template<typename eT, typename T1>
inline
bool
auxlib::schur(Mat<eT>& U, Mat<eT>& S, const Base<eT,T1>& X, const bool calc_U)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    S = X.get_ref();
    
    arma_debug_check( (S.is_square() == false), "schur(): given matrix must be square sized" );
    
    if(S.is_empty())  { U.reset(); S.reset(); return true; }
    
    arma_debug_assert_blas_size(S);
    
    const uword S_n_rows = S.n_rows;
    
    if(calc_U) { U.set_size(S_n_rows, S_n_rows); } else { U.set_size(1,1); }
    
    char      jobvs  = calc_U ? 'V' : 'N';
    char      sort   = 'N';
    void*     select = 0;
    blas_int  n      = blas_int(S_n_rows);
    blas_int  sdim   = 0;
    blas_int  ldvs   = calc_U ? n : blas_int(1);
    blas_int  lwork  = 64*n;  // lwork_min = (std::max)(blas_int(1), 3*n)
    blas_int  info   = 0;
    
    podarray<eT> wr(S_n_rows);
    podarray<eT> wi(S_n_rows);
    
    podarray<eT>        work( static_cast<uword>(lwork) );
    podarray<blas_int> bwork(S_n_rows);
    
    arma_extra_debug_print("lapack::gees()");
    lapack::gees(&jobvs, &sort, select, &n, S.memptr(), &n, &sdim, wr.memptr(), wi.memptr(), U.memptr(), &ldvs, work.memptr(), &lwork, bwork.memptr(), &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(U);
    arma_ignore(S);
    arma_ignore(X);
    arma_ignore(calc_U);
    arma_stop_logic_error("schur(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename T, typename T1>
inline
bool
auxlib::schur(Mat< std::complex<T> >& U, Mat< std::complex<T> >& S, const Base<std::complex<T>,T1>& X, const bool calc_U)
  {
  arma_extra_debug_sigprint();
  
  S = X.get_ref();
  
  arma_debug_check( (S.is_square() == false), "schur(): given matrix must be square sized" );
  
  return auxlib::schur(U,S,calc_U);
  }



template<typename T>
inline
bool
auxlib::schur(Mat< std::complex<T> >& U, Mat< std::complex<T> >& S, const bool calc_U)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef std::complex<T> eT;
    
    if(S.is_empty())  { U.reset(); S.reset(); return true; }
    
    arma_debug_assert_blas_size(S);
    
    const uword S_n_rows = S.n_rows;
    
    if(calc_U) { U.set_size(S_n_rows, S_n_rows); } else { U.set_size(1,1); }
    
    char      jobvs  = calc_U ? 'V' : 'N';
    char      sort   = 'N';
    void*     select = 0;
    blas_int  n      = blas_int(S_n_rows);
    blas_int  sdim   = 0;
    blas_int  ldvs   = calc_U ? n : blas_int(1);
    blas_int  lwork  = 64*n;  // lwork_min = (std::max)(blas_int(1), 2*n)
    blas_int  info   = 0;
    
    podarray<eT>           w(S_n_rows);
    podarray<eT>        work( static_cast<uword>(lwork) );
    podarray< T>       rwork(S_n_rows);
    podarray<blas_int> bwork(S_n_rows);
    
    arma_extra_debug_print("lapack::cx_gees()");
    lapack::cx_gees(&jobvs, &sort, select, &n, S.memptr(), &n, &sdim, w.memptr(), U.memptr(), &ldvs, work.memptr(), &lwork, rwork.memptr(), bwork.memptr(), &info);
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(U);
    arma_ignore(S);
    arma_ignore(calc_U);
    arma_stop_logic_error("schur(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//
// solve the Sylvester equation AX + XB = C

template<typename eT>
inline
bool
auxlib::syl(Mat<eT>& X, const Mat<eT>& A, const Mat<eT>& B, const Mat<eT>& C)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_check( (A.is_square() == false) || (B.is_square() == false), "syl(): given matrices must be square sized" );
      
    arma_debug_check( (C.n_rows != A.n_rows) || (C.n_cols != B.n_cols), "syl(): matrices are not conformant" );
    
    if(A.is_empty() || B.is_empty() || C.is_empty())  { X.reset(); return true; }
    
    Mat<eT> Z1, Z2, T1, T2;
    
    const bool status_sd1 = auxlib::schur(Z1, T1, A);
    const bool status_sd2 = auxlib::schur(Z2, T2, B);
    
    if( (status_sd1 == false) || (status_sd2 == false) )  { return false; }
    
    char     trana = 'N';
    char     tranb = 'N';
    blas_int  isgn = +1;
    blas_int     m = blas_int(T1.n_rows);
    blas_int     n = blas_int(T2.n_cols);
    
    eT       scale = eT(0);
    blas_int  info = 0;
    
    Mat<eT> Y = trans(Z1) * C * Z2;
    
    arma_extra_debug_print("lapack::trsyl()");
    lapack::trsyl<eT>(&trana, &tranb, &isgn, &m, &n, T1.memptr(), &m, T2.memptr(), &n, Y.memptr(), &m, &scale, &info);
    
    if(info < 0)  { return false; }
    
    //Y /= scale;
    Y /= (-scale);
    
    X = Z1 * Y * trans(Z2);
    
    return true;
    }
  #else
    {
    arma_ignore(X);
    arma_ignore(A);
    arma_ignore(B);
    arma_ignore(C);
    arma_stop_logic_error("syl(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }
  
  

//
// QZ decomposition of general square real matrix pair

template<typename T, typename T1, typename T2>
inline
bool
auxlib::qz(Mat<T>& A, Mat<T>& B, Mat<T>& vsl, Mat<T>& vsr, const Base<T,T1>& X_expr, const Base<T,T2>& Y_expr, const char mode)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    A = X_expr.get_ref();
    B = Y_expr.get_ref();
    
    arma_debug_check( ((A.is_square() == false) || (B.is_square() == false)), "qz(): given matrices must be square sized", [&](){ A.soft_reset(); B.soft_reset(); } );
    
    arma_debug_check( (A.n_rows != B.n_rows), "qz(): given matrices must have the same size" );
    
    if(A.is_empty())  { A.reset();  B.reset();  vsl.reset(); vsr.reset(); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    if(arma_config::check_nonfinite && B.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    vsl.set_size(A.n_rows, A.n_rows);
    vsr.set_size(A.n_rows, A.n_rows);
    
    char     jobvsl  = 'V';
    char     jobvsr  = 'V';
    char     eigsort = 'N';
    void*    selctg  = 0;
    blas_int N       = blas_int(A.n_rows);
    blas_int sdim    = 0;
    blas_int lwork   = 64*N+16;  // lwork_min = (std::max)(blas_int(1),8*N+16)
    blas_int info    = 0;
    
         if(mode == 'l')  { eigsort = 'S'; selctg = qz_helper::ptr_cast(&(qz_helper::select_lhp<T>)); }
    else if(mode == 'r')  { eigsort = 'S'; selctg = qz_helper::ptr_cast(&(qz_helper::select_rhp<T>)); }
    else if(mode == 'i')  { eigsort = 'S'; selctg = qz_helper::ptr_cast(&(qz_helper::select_iuc<T>)); }
    else if(mode == 'o')  { eigsort = 'S'; selctg = qz_helper::ptr_cast(&(qz_helper::select_ouc<T>)); }
    
    podarray<T> alphar(A.n_rows);
    podarray<T> alphai(A.n_rows);
    podarray<T>   beta(A.n_rows);
    
    podarray<T>         work( static_cast<uword>(lwork) );
    podarray<blas_int> bwork( static_cast<uword>(N)     );
    
    arma_extra_debug_print("lapack::gges()");
    
    lapack::gges
      (
      &jobvsl, &jobvsr, &eigsort, selctg, &N,
      A.memptr(), &N, B.memptr(), &N, &sdim,
      alphar.memptr(), alphai.memptr(), beta.memptr(),
      vsl.memptr(), &N, vsr.memptr(), &N,
      work.memptr(), &lwork, bwork.memptr(),
      &info
      );
    
    if(info != 0)  { return false; }
    
    op_strans::apply_mat_inplace(vsl);
    
    return true;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(B);
    arma_ignore(vsl);
    arma_ignore(vsr);
    arma_ignore(X_expr);
    arma_ignore(Y_expr);
    arma_ignore(mode);
    arma_stop_logic_error("qz(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



//
// QZ decomposition of general square complex matrix pair

template<typename T, typename T1, typename T2>
inline
bool
auxlib::qz(Mat< std::complex<T> >& A, Mat< std::complex<T> >& B, Mat< std::complex<T> >& vsl, Mat< std::complex<T> >& vsr, const Base< std::complex<T>, T1 >& X_expr, const Base< std::complex<T>, T2 >& Y_expr, const char mode)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename std::complex<T> eT;
    
    A = X_expr.get_ref();
    B = Y_expr.get_ref();
    
    arma_debug_check( ((A.is_square() == false) || (B.is_square() == false)), "qz(): given matrices must be square sized", [&](){ A.soft_reset(); B.soft_reset(); } );
    
    arma_debug_check( (A.n_rows != B.n_rows), "qz(): given matrices must have the same size" );
    
    if(A.is_empty())  { A.reset(); B.reset(); vsl.reset(); vsr.reset(); return true; }
    
    if(arma_config::check_nonfinite && A.has_nonfinite())  { return false; }
    if(arma_config::check_nonfinite && B.has_nonfinite())  { return false; }
    
    arma_debug_assert_blas_size(A);
    
    vsl.set_size(A.n_rows, A.n_rows);
    vsr.set_size(A.n_rows, A.n_rows);
    
    char     jobvsl  = 'V';
    char     jobvsr  = 'V';
    char     eigsort = 'N';
    void*    selctg  = 0;
    blas_int N       = blas_int(A.n_rows);
    blas_int sdim    = 0;
    blas_int lwork   = 64*N;  // lwork_min = (std::max)(blas_int(1),2*N)
    blas_int info    = 0;
    
         if(mode == 'l')  { eigsort = 'S'; selctg = qz_helper::ptr_cast(&(qz_helper::cx_select_lhp<T>)); }
    else if(mode == 'r')  { eigsort = 'S'; selctg = qz_helper::ptr_cast(&(qz_helper::cx_select_rhp<T>)); }
    else if(mode == 'i')  { eigsort = 'S'; selctg = qz_helper::ptr_cast(&(qz_helper::cx_select_iuc<T>)); }
    else if(mode == 'o')  { eigsort = 'S'; selctg = qz_helper::ptr_cast(&(qz_helper::cx_select_ouc<T>)); }
    
    podarray<eT> alpha(A.n_rows);
    podarray<eT>  beta(A.n_rows);
    
    podarray<eT>        work( static_cast<uword>(lwork) );
    podarray< T>       rwork( static_cast<uword>(8*N)   );
    podarray<blas_int> bwork( static_cast<uword>(N)     );
    
    arma_extra_debug_print("lapack::cx_gges()");
    
    lapack::cx_gges
      (
      &jobvsl, &jobvsr, &eigsort, selctg, &N,
      A.memptr(), &N, B.memptr(), &N, &sdim,
      alpha.memptr(), beta.memptr(),
      vsl.memptr(), &N, vsr.memptr(), &N,
      work.memptr(), &lwork, rwork.memptr(), bwork.memptr(),
      &info
      );
    
    if(info != 0)  { return false; }
    
    op_htrans::apply_mat_inplace(vsl);
    
    return true;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(B);
    arma_ignore(vsl);
    arma_ignore(vsr);
    arma_ignore(X_expr);
    arma_ignore(Y_expr);
    arma_ignore(mode);
    arma_stop_logic_error("qz(): use of LAPACK must be enabled");
    return false;
    }
  #endif
  }



template<typename eT>
inline
eT
auxlib::rcond(Mat<eT>& A)
  {
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(A);
    
    char     norm_id  = '1';
    blas_int m        = blas_int(A.n_rows);
    blas_int n        = blas_int(A.n_rows);  // assuming square matrix
    blas_int lda      = blas_int(A.n_rows);
    eT       norm_val = eT(0);
    eT       rcond    = eT(0);
    blas_int info     = blas_int(0);
    
    podarray<eT>        work(4*A.n_rows);
    podarray<blas_int> iwork(  A.n_rows);
    podarray<blas_int> ipiv( (std::min)(A.n_rows, A.n_cols) );
    
    arma_extra_debug_print("lapack::lange()");
    norm_val = lapack::lange(&norm_id, &m, &n, A.memptr(), &lda, work.memptr());
    
    arma_extra_debug_print("lapack::getrf()");
    lapack::getrf(&m, &n, A.memptr(), &lda, ipiv.memptr(), &info);
    
    if(info != blas_int(0))  { return eT(0); }
    
    arma_extra_debug_print("lapack::gecon()");
    lapack::gecon(&norm_id, &n, A.memptr(), &lda, &norm_val, &rcond, work.memptr(), iwork.memptr(), &info);
    
    if(info != blas_int(0))  { return eT(0); }
    
    return rcond;
    }
  #else
    {
    arma_ignore(A);
    arma_stop_logic_error("rcond(): use of LAPACK must be enabled");
    return eT(0);
    }
  #endif
  }



template<typename T>
inline
T
auxlib::rcond(Mat< std::complex<T> >& A)
  {
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename std::complex<T> eT;
    
    arma_debug_assert_blas_size(A);
    
    char     norm_id  = '1';
    blas_int m        = blas_int(A.n_rows);
    blas_int n        = blas_int(A.n_rows);  // assuming square matrix
    blas_int lda      = blas_int(A.n_rows);
    T        norm_val = T(0);
    T        rcond    = T(0);
    blas_int info     = blas_int(0);
    
    podarray< T>        junk(1);
    podarray<eT>        work(2*A.n_rows);
    podarray< T>       rwork(2*A.n_rows);
    podarray<blas_int> ipiv( (std::min)(A.n_rows, A.n_cols) );
    
    arma_extra_debug_print("lapack::lange()");
    norm_val = lapack::lange(&norm_id, &m, &n, A.memptr(), &lda, junk.memptr());
    
    arma_extra_debug_print("lapack::getrf()");
    lapack::getrf(&m, &n, A.memptr(), &lda, ipiv.memptr(), &info);
    
    if(info != blas_int(0))  { return T(0); }
    
    arma_extra_debug_print("lapack::cx_gecon()");
    lapack::cx_gecon(&norm_id, &n, A.memptr(), &lda, &norm_val, &rcond, work.memptr(), rwork.memptr(), &info);
    
    if(info != blas_int(0))  { return T(0); }
    
    return rcond;
    }
  #else
    {
    arma_ignore(A);
    arma_stop_logic_error("rcond(): use of LAPACK must be enabled");
    return T(0);
    }
  #endif
  }



template<typename eT>
inline
eT
auxlib::rcond_sympd(Mat<eT>& A, bool& calc_ok)
  {
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(A);
    
    calc_ok = false;
    
    char     norm_id  = '1';
    char     uplo     = 'L';
    blas_int n        = blas_int(A.n_rows);  // assuming square matrix
    blas_int lda      = blas_int(A.n_rows);
    eT       norm_val = eT(0);
    eT       rcond    = eT(0);
    blas_int info     = blas_int(0);
    
    podarray<eT>        work(3*A.n_rows);
    podarray<blas_int> iwork(  A.n_rows);
    
    arma_extra_debug_print("lapack::lansy()");
    norm_val = lapack::lansy(&norm_id, &uplo, &n, A.memptr(), &lda, work.memptr());
    
    arma_extra_debug_print("lapack::potrf()");
    lapack::potrf(&uplo, &n, A.memptr(), &lda, &info);
    
    if(info != blas_int(0))  { return eT(0); }
    
    arma_extra_debug_print("lapack::pocon()");
    lapack::pocon(&uplo, &n, A.memptr(), &lda, &norm_val, &rcond, work.memptr(), iwork.memptr(), &info);
    
    if(info != blas_int(0))  { return eT(0); }
    
    calc_ok = true;
    
    return rcond;
    }
  #else
    {
    arma_ignore(A);
    calc_ok = false;
    arma_stop_logic_error("rcond(): use of LAPACK must be enabled");
    return eT(0);
    }
  #endif
  }



template<typename T>
inline
T
auxlib::rcond_sympd(Mat< std::complex<T> >& A, bool& calc_ok)
  {
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_extra_debug_print("auxlib::rcond_sympd(): redirecting to auxlib::rcond() due to crippled LAPACK");
    
    calc_ok = true;
    
    return auxlib::rcond(A);
    }
  #elif defined(ARMA_USE_LAPACK)
    {
    typedef typename std::complex<T> eT;
    
    arma_debug_assert_blas_size(A);
    
    calc_ok = false;
    
    char     norm_id  = '1';
    char     uplo     = 'L';
    blas_int n        = blas_int(A.n_rows);  // assuming square matrix
    blas_int lda      = blas_int(A.n_rows);
    T        norm_val = T(0);
    T        rcond    = T(0);
    blas_int info     = blas_int(0);
    
    podarray<eT>  work(2*A.n_rows);
    podarray< T> rwork(  A.n_rows);
    
    arma_extra_debug_print("lapack::lanhe()");
    norm_val = lapack::lanhe(&norm_id, &uplo, &n, A.memptr(), &lda, rwork.memptr());
    
    arma_extra_debug_print("lapack::potrf()");
    lapack::potrf(&uplo, &n, A.memptr(), &lda, &info);
    
    if(info != blas_int(0))  { return T(0); }
    
    arma_extra_debug_print("lapack::cx_pocon()");
    lapack::cx_pocon(&uplo, &n, A.memptr(), &lda, &norm_val, &rcond, work.memptr(), rwork.memptr(), &info);
    
    if(info != blas_int(0))  { return T(0); }
    
    calc_ok = true;
    
    return rcond;
    }
  #else
    {
    arma_ignore(A);
    calc_ok = false;
    arma_stop_logic_error("rcond(): use of LAPACK must be enabled");
    return T(0);
    }
  #endif
  }



template<typename eT>
inline
eT
auxlib::rcond_trimat(const Mat<eT>& A, const uword layout)
  {
  #if defined(ARMA_USE_LAPACK)
    {
    arma_debug_assert_blas_size(A);
    
    char     norm_id = '1';
    char     uplo    = (layout == 0) ? 'U' : 'L';
    char     diag    = 'N';
    blas_int n       = blas_int(A.n_rows);  // assuming square matrix
    eT       rcond   = eT(0);
    blas_int info    = blas_int(0);
    
    podarray<eT>        work(3*A.n_rows);
    podarray<blas_int> iwork(  A.n_rows);
    
    arma_extra_debug_print("lapack::trcon()");
    lapack::trcon(&norm_id, &uplo, &diag, &n, A.memptr(), &n, &rcond, work.memptr(), iwork.memptr(), &info);
    
    if(info != blas_int(0))  { return eT(0); }
    
    return rcond;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(layout);
    arma_stop_logic_error("rcond(): use of LAPACK must be enabled");
    return eT(0);
    }
  #endif
  }



template<typename T>
inline
T
auxlib::rcond_trimat(const Mat< std::complex<T> >& A, const uword layout)
  {
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename std::complex<T> eT;
    
    arma_debug_assert_blas_size(A);
    
    char     norm_id = '1';
    char     uplo    = (layout == 0) ? 'U' : 'L';
    char     diag    = 'N';
    blas_int n       = blas_int(A.n_rows);  // assuming square matrix
    T        rcond   = T(0);
    blas_int info    = blas_int(0);
    
    podarray<eT>  work(2*A.n_rows);
    podarray< T> rwork(  A.n_rows);
    
    arma_extra_debug_print("lapack::cx_trcon()");
    lapack::cx_trcon(&norm_id, &uplo, &diag, &n, A.memptr(), &n, &rcond, work.memptr(), rwork.memptr(), &info);
    
    if(info != blas_int(0))  { return T(0); }
    
    return rcond;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(layout);
    arma_stop_logic_error("rcond(): use of LAPACK must be enabled");
    return T(0);
    }
  #endif
  }



template<typename eT>
inline
eT
auxlib::lu_rcond(const Mat<eT>& A, const eT norm_val)
  {
  #if defined(ARMA_USE_LAPACK)
    {
    char     norm_id = '1';
    blas_int n       = blas_int(A.n_rows);  // assuming square matrix
    blas_int lda     = blas_int(A.n_rows);
    eT       rcond   = eT(0);
    blas_int info    = blas_int(0);
    
    podarray<eT>        work(4*A.n_rows);
    podarray<blas_int> iwork(  A.n_rows);
    
    arma_extra_debug_print("lapack::gecon()");
    lapack::gecon(&norm_id, &n, A.memptr(), &lda, &norm_val, &rcond, work.memptr(), iwork.memptr(), &info);
    
    if(info != blas_int(0))  { return eT(0); }
    
    return rcond;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(norm_val);
    return eT(0);
    }
  #endif
  }



template<typename T>
inline
T
auxlib::lu_rcond(const Mat< std::complex<T> >& A, const T norm_val)
  {
  #if defined(ARMA_USE_LAPACK)
    {
    typedef typename std::complex<T> eT;
    
    char     norm_id = '1';
    blas_int n       = blas_int(A.n_rows);  // assuming square matrix
    blas_int lda     = blas_int(A.n_rows);
    T        rcond   = T(0);
    blas_int info    = blas_int(0);
    
    podarray<eT>  work(2*A.n_rows);
    podarray< T> rwork(2*A.n_rows);
    
    arma_extra_debug_print("lapack::cx_gecon()");
    lapack::cx_gecon(&norm_id, &n, A.memptr(), &lda, &norm_val, &rcond, work.memptr(), rwork.memptr(), &info);
    
    if(info != blas_int(0))  { return T(0); }
    
    return rcond;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(norm_val);
    return T(0);
    }
  #endif
  }



template<typename eT>
inline
eT
auxlib::lu_rcond_sympd(const Mat<eT>& A, const eT norm_val)
  {
  #if defined(ARMA_USE_LAPACK)
    {
    char     uplo  = 'L';
    blas_int n     = blas_int(A.n_rows);  // assuming square matrix
    eT       rcond = eT(0);
    blas_int info  = blas_int(0);
    
    podarray<eT>        work(3*A.n_rows);
    podarray<blas_int> iwork(  A.n_rows);
    
    arma_extra_debug_print("lapack::pocon()");
    lapack::pocon(&uplo, &n, A.memptr(), &n, &norm_val, &rcond, work.memptr(), iwork.memptr(), &info);
    
    if(info != blas_int(0))  { return eT(0); }
    
    return rcond;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(norm_val);
    return eT(0);
    }
  #endif
  }



template<typename T>
inline
T
auxlib::lu_rcond_sympd(const Mat< std::complex<T> >& A, const T norm_val)
  {
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_ignore(A);
    arma_ignore(norm_val);
    return T(0);
    }
  #elif defined(ARMA_USE_LAPACK)
    {
    typedef typename std::complex<T> eT;
    
    char     uplo  = 'L';
    blas_int n     = blas_int(A.n_rows);  // assuming square matrix
    T        rcond = T(0);
    blas_int info  = blas_int(0);
    
    podarray<eT>  work(2*A.n_rows);
    podarray< T> rwork(  A.n_rows);
    
    arma_extra_debug_print("lapack::cx_pocon()");
    lapack::cx_pocon(&uplo, &n, A.memptr(), &n, &norm_val, &rcond, work.memptr(), rwork.memptr(), &info);
    
    if(info != blas_int(0))  { return T(0); }
    
    return rcond;
    }
  #else
    {
    arma_ignore(A);
    arma_ignore(norm_val);
    return T(0);
    }
  #endif
  }



template<typename eT>
inline
eT
auxlib::lu_rcond_band(const Mat<eT>& AB, const uword KL, const uword KU, const podarray<blas_int>& ipiv, const eT norm_val)
  {
  #if defined(ARMA_USE_LAPACK)
    {
    const uword N = AB.n_cols;  // order of the original square matrix A
    
    char     norm_id = '1';
    blas_int n       = blas_int(N);
    blas_int kl      = blas_int(KL);
    blas_int ku      = blas_int(KU);
    blas_int ldab    = blas_int(AB.n_rows);
    eT       rcond   = eT(0);
    blas_int info    = blas_int(0);
    
    podarray<eT>        work(3*N);
    podarray<blas_int> iwork(  N);
    
    arma_extra_debug_print("lapack::gbcon()");
    lapack::gbcon<eT>(&norm_id, &n, &kl, &ku, AB.memptr(), &ldab, ipiv.memptr(), &norm_val, &rcond, work.memptr(), iwork.memptr(), &info);
    
    if(info != blas_int(0))  { return eT(0); }
    
    return rcond;
    }
  #else
    {
    arma_ignore(AB);
    arma_ignore(KL);
    arma_ignore(KU);
    arma_ignore(ipiv);
    arma_ignore(norm_val);
    return eT(0);
    }
  #endif
  }



template<typename T>
inline
T
auxlib::lu_rcond_band(const Mat< std::complex<T> >& AB, const uword KL, const uword KU, const podarray<blas_int>& ipiv, const T norm_val)
  {
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_ignore(AB);
    arma_ignore(KL);
    arma_ignore(KU);
    arma_ignore(ipiv);
    arma_ignore(norm_val);
    return T(0);
    }
  #elif defined(ARMA_USE_LAPACK)
    {
    typedef typename std::complex<T> eT;
    
    const uword N = AB.n_cols;  // order of the original square matrix A
    
    char     norm_id = '1';
    blas_int n       = blas_int(N);
    blas_int kl      = blas_int(KL);
    blas_int ku      = blas_int(KU);
    blas_int ldab    = blas_int(AB.n_rows);
    T        rcond   = T(0);
    blas_int info    = blas_int(0);
    
    podarray<eT>  work(2*N);
    podarray< T> rwork(  N);
    
    arma_extra_debug_print("lapack::cx_gbcon()");
    lapack::cx_gbcon<T>(&norm_id, &n, &kl, &ku, AB.memptr(), &ldab, ipiv.memptr(), &norm_val, &rcond, work.memptr(), rwork.memptr(), &info);
    
    if(info != blas_int(0))  { return T(0); }
    
    return rcond;
    }
  #else
    {
    arma_ignore(AB);
    arma_ignore(KL);
    arma_ignore(KU);
    arma_ignore(ipiv);
    arma_ignore(norm_val);
    return T(0);
    }
  #endif
  }



template<typename T1>
inline
bool
auxlib::crippled_lapack(const Base<typename T1::elem_type, T1>&)
  {
  #if defined(ARMA_CRIPPLED_LAPACK)
    {
    arma_extra_debug_print("auxlib::crippled_lapack(): true");
    
    return (is_cx<typename T1::elem_type>::yes);
    }
  #else
    {
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
auxlib::rudimentary_sym_check(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const uword N   = X.n_rows;
  const uword Nm2 = N-2;
  
  if(N != X.n_cols)  { return false; }
  if(N <= uword(1))  { return true;  }
  
  const eT* X_mem = X.memptr();
  
  const eT* X_offsetA = &(X_mem[Nm2  ]);
  const eT* X_offsetB = &(X_mem[Nm2*N]);
  
  const eT A1 = *(X_offsetA  ); 
  const eT A2 = *(X_offsetA+1);  // bottom-left corner (ie. last value in first column)
  const eT B1 = *(X_offsetB  );  
  const eT B2 = *(X_offsetB+N);  // top-right   corner (ie. first value in last column)
  
  const eT C1 = (std::max)(std::abs(A1), std::abs(B1));
  const eT C2 = (std::max)(std::abs(A2), std::abs(B2));
  
  const eT delta1 = std::abs(A1 - B1);
  const eT delta2 = std::abs(A2 - B2);
  
  const eT tol = eT(10000)*std::numeric_limits<eT>::epsilon();  // allow some leeway
  
  const bool okay1 = ( (delta1 <= tol) || (delta1 <= (C1 * tol)) );
  const bool okay2 = ( (delta2 <= tol) || (delta2 <= (C2 * tol)) );
  
  return (okay1 && okay2);
  }



template<typename T>
inline
bool
auxlib::rudimentary_sym_check(const Mat< std::complex<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: the function name is a misnomer, as it checks for hermitian complex matrices;
  // NOTE: for simplicity of use, the function name is the same as for real matrices
  
  typedef typename std::complex<T> eT;
  
  const uword N   = X.n_rows;
  const uword Nm1 = N-1;
  
  if(N != X.n_cols)  { return false; }
  if(N == uword(0))  { return true;  }
  
  const eT* X_mem = X.memptr();
  
  const T tol = T(10000)*std::numeric_limits<T>::epsilon();  // allow some leeway
  
  if(std::abs(X_mem[0         ].imag()) > tol)  { return false; }  // check top-left
  if(std::abs(X_mem[X.n_elem-1].imag()) > tol)  { return false; }  // check bottom-right
  
  const eT& A = X_mem[Nm1  ];  // bottom-left corner (ie. last value in first column)
  const eT& B = X_mem[Nm1*N];  // top-right   corner (ie. first value in last column)
  
  const T C_real = (std::max)(std::abs(A.real()), std::abs(B.real()));
  const T C_imag = (std::max)(std::abs(A.imag()), std::abs(B.imag()));
  
  const T delta_real = std::abs(A.real() - B.real());
  const T delta_imag = std::abs(A.imag() + B.imag());  // take into account the conjugate
  
  const bool okay_real = ( (delta_real <= tol) || (delta_real <= (C_real * tol)) );
  const bool okay_imag = ( (delta_imag <= tol) || (delta_imag <= (C_imag * tol)) );
  
  return (okay_real && okay_imag);
  }



//



namespace qz_helper
{

// sgges() and dgges() require an external function with three arguments:
// select(alpha_real, alpha_imag, beta)
// where the eigenvalue is defined as complex(alpha_real, alpha_imag) / beta

template<typename T>
inline
blas_int
select_lhp(const T* x_ptr, const T* y_ptr, const T* z_ptr)
  {
  arma_extra_debug_sigprint();
  
  // cout << "select_lhp(): (*x_ptr) = " << (*x_ptr) << endl;
  // cout << "select_lhp(): (*y_ptr) = " << (*y_ptr) << endl;
  // cout << "select_lhp(): (*z_ptr) = " << (*z_ptr) << endl;
  
  arma_ignore(y_ptr);  // ignore imaginary part
  
  const T x = (*x_ptr);
  const T z = (*z_ptr);
  
  if(z == T(0))  { return blas_int(0); }  // consider an infinite eig value not to lie in either lhp or rhp
  
  return ((x/z) < T(0)) ? blas_int(1) : blas_int(0);
  }



template<typename T>
inline
blas_int
select_rhp(const T* x_ptr, const T* y_ptr, const T* z_ptr)
  {
  arma_extra_debug_sigprint();
  
  // cout << "select_rhp(): (*x_ptr) = " << (*x_ptr) << endl;
  // cout << "select_rhp(): (*y_ptr) = " << (*y_ptr) << endl;
  // cout << "select_rhp(): (*z_ptr) = " << (*z_ptr) << endl;
  
  arma_ignore(y_ptr);  // ignore imaginary part
  
  const T x = (*x_ptr);
  const T z = (*z_ptr);
  
  if(z == T(0))  { return blas_int(0); }  // consider an infinite eig value not to lie in either lhp or rhp
  
  return ((x/z) > T(0)) ? blas_int(1) : blas_int(0);
  }



template<typename T>
inline
blas_int
select_iuc(const T* x_ptr, const T* y_ptr, const T* z_ptr)
  {
  arma_extra_debug_sigprint();
  
  // cout << "select_iuc(): (*x_ptr) = " << (*x_ptr) << endl;
  // cout << "select_iuc(): (*y_ptr) = " << (*y_ptr) << endl;
  // cout << "select_iuc(): (*z_ptr) = " << (*z_ptr) << endl;
  
  const T x = (*x_ptr);
  const T y = (*y_ptr);
  const T z = (*z_ptr);
  
  if(z == T(0))  { return blas_int(0); }  // consider an infinite eig value to be outside of the unit circle 
  
  //return (std::abs(std::complex<T>(x,y) / z) < T(1)) ? blas_int(1) : blas_int(0);
  return (std::sqrt(x*x + y*y) < std::abs(z)) ? blas_int(1) : blas_int(0);
  }



template<typename T>
inline
blas_int
select_ouc(const T* x_ptr, const T* y_ptr, const T* z_ptr)
  {
  arma_extra_debug_sigprint();
  
  // cout << "select_ouc(): (*x_ptr) = " << (*x_ptr) << endl;
  // cout << "select_ouc(): (*y_ptr) = " << (*y_ptr) << endl;
  // cout << "select_ouc(): (*z_ptr) = " << (*z_ptr) << endl;
  
  const T x = (*x_ptr);
  const T y = (*y_ptr);
  const T z = (*z_ptr);
  
  if(z == T(0))
    {
    return (x == T(0)) ? blas_int(0) : blas_int(1);  // consider an infinite eig value to be outside of the unit circle 
    }
  
  //return (std::abs(std::complex<T>(x,y) / z) > T(1)) ? blas_int(1) : blas_int(0);
  return (std::sqrt(x*x + y*y) > std::abs(z)) ? blas_int(1) : blas_int(0);
  }



// cgges() and zgges() require an external function with two arguments:
// select(alpha, beta)
// where the complex eigenvalue is defined as (alpha / beta)

template<typename T>
inline
blas_int
cx_select_lhp(const std::complex<T>* x_ptr, const std::complex<T>* y_ptr)
  {
  arma_extra_debug_sigprint();
  
  // cout << "cx_select_lhp(): (*x_ptr) = " << (*x_ptr) << endl;
  // cout << "cx_select_lhp(): (*y_ptr) = " << (*y_ptr) << endl;
  
  const std::complex<T>& x = (*x_ptr);
  const std::complex<T>& y = (*y_ptr);
  
  if( (y.real() == T(0)) && (y.imag() == T(0)) )  { return blas_int(0); }  // consider an infinite eig value not to lie in either lhp or rhp
  
  return (std::real(x / y) < T(0)) ? blas_int(1) : blas_int(0);
  }



template<typename T>
inline
blas_int
cx_select_rhp(const std::complex<T>* x_ptr, const std::complex<T>* y_ptr)
  {
  arma_extra_debug_sigprint();
  
  // cout << "cx_select_rhp(): (*x_ptr) = " << (*x_ptr) << endl;
  // cout << "cx_select_rhp(): (*y_ptr) = " << (*y_ptr) << endl;
  
  const std::complex<T>& x = (*x_ptr);
  const std::complex<T>& y = (*y_ptr);
  
  if( (y.real() == T(0)) && (y.imag() == T(0)) )  { return blas_int(0); }  // consider an infinite eig value not to lie in either lhp or rhp
  
  return (std::real(x / y) > T(0)) ? blas_int(1) : blas_int(0);
  }



template<typename T>
inline
blas_int
cx_select_iuc(const std::complex<T>* x_ptr, const std::complex<T>* y_ptr)
  {
  arma_extra_debug_sigprint();
  
  // cout << "cx_select_iuc(): (*x_ptr) = " << (*x_ptr) << endl;
  // cout << "cx_select_iuc(): (*y_ptr) = " << (*y_ptr) << endl;
  
  const std::complex<T>& x = (*x_ptr);
  const std::complex<T>& y = (*y_ptr);
  
  if( (y.real() == T(0)) && (y.imag() == T(0)) )  { return blas_int(0); }  // consider an infinite eig value to be outside of the unit circle
  
  return (std::abs(x / y) < T(1)) ? blas_int(1) : blas_int(0);
  }



template<typename T>
inline
blas_int
cx_select_ouc(const std::complex<T>* x_ptr, const std::complex<T>* y_ptr)
  {
  arma_extra_debug_sigprint();
  
  // cout << "cx_select_ouc(): (*x_ptr) = " << (*x_ptr) << endl;
  // cout << "cx_select_ouc(): (*y_ptr) = " << (*y_ptr) << endl;
  
  const std::complex<T>& x = (*x_ptr);
  const std::complex<T>& y = (*y_ptr);
  
  if( (y.real() == T(0)) && (y.imag() == T(0)) )
    {
    return ((x.real() == T(0)) && (x.imag() == T(0))) ? blas_int(0) : blas_int(1);  // consider an infinite eig value to be outside of the unit circle
    }
  
  return (std::abs(x / y) > T(1)) ? blas_int(1) : blas_int(0);
  }



// need to do shenanigans with pointers due to:
// - we're using LAPACK ?gges() defined to expect pointer-to-function to be passed as pointer-to-object
// - explicit casting between pointer-to-function and pointer-to-object is a non-standard extension in C
// - the extension is essentially mandatory on POSIX systems
// - some compilers will complain about the extension in pedantic mode

template<typename T>
inline
void_ptr
ptr_cast(blas_int (*function)(const T*, const T*, const T*))
  {
  union converter
    {
    blas_int (*fn)(const T*, const T*, const T*);
    void_ptr obj;
    };
  
  converter tmp;
  
  tmp.obj = 0;
  tmp.fn  = function;
  
  return tmp.obj;
  }



template<typename T>
inline
void_ptr
ptr_cast(blas_int (*function)(const std::complex<T>*, const std::complex<T>*))
  {
  union converter
    {
    blas_int (*fn)(const std::complex<T>*, const std::complex<T>*);
    void_ptr obj;
    };
  
  converter tmp;
  
  tmp.obj = 0;
  tmp.fn  = function;
  
  return tmp.obj;
  }



}  // end of namespace qz_helper


//! @}
