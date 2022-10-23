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


//! low-level interface functions for accessing LAPACK
class auxlib
  {
  public:
  
  //
  // inv
  
  template<typename eT>
  inline static bool inv(Mat<eT>& A);
  
  template<typename eT>
  inline static bool inv(Mat<eT>& out, const Mat<eT>& X);
  
  template<typename eT>
  inline static bool inv_rcond(Mat<eT>& A, typename get_pod_type<eT>::result& out_rcond);
  
  template<typename eT>
  inline static bool inv_tr(Mat<eT>& A, const uword layout);
  
  template<typename eT>
  inline static bool inv_tr_rcond(Mat<eT>& A, typename get_pod_type<eT>::result& out_rcond, const uword layout);
  
  template<typename eT>
  inline static bool inv_sympd(Mat<eT>& A, bool& out_sympd_state);
  
  template<typename eT>
  inline static bool inv_sympd(Mat<eT>& out, const Mat<eT>& X);
  
  template<typename eT>
  inline static bool inv_sympd_rcond(Mat<eT>& A, bool& out_sympd_state, eT& out_rcond, const eT rcond_threshold);
  
  template<typename T>
  inline static bool inv_sympd_rcond(Mat< std::complex<T> >& A, bool& out_sympd_state, T& out_rcond, const T rcond_threshold);
  
  
  //
  // det and log_det
  
  template<typename eT>
  inline static bool det(eT& out_val, Mat<eT>& A);
  
  template<typename eT>
  inline static bool log_det(eT& out_val, typename get_pod_type<eT>::result& out_sign, Mat<eT>& A);
  
  template<typename eT>
  inline static bool log_det_sympd(typename get_pod_type<eT>::result& out_val, Mat<eT>& A);
  
  
  //
  // lu
  
  template<typename eT, typename T1>
  inline static bool lu(Mat<eT>& L, Mat<eT>& U, podarray<blas_int>& ipiv, const Base<eT,T1>& X);
  
  template<typename eT, typename T1>
  inline static bool lu(Mat<eT>& L, Mat<eT>& U, Mat<eT>& P, const Base<eT,T1>& X);
  
  template<typename eT, typename T1>
  inline static bool lu(Mat<eT>& L, Mat<eT>& U, const Base<eT,T1>& X);
  
  
  //
  // eig_gen
  
  template<typename T1>
  inline static bool eig_gen(Mat< std::complex<typename T1::pod_type> >& vals, Mat< std::complex<typename T1::pod_type> >& vecs, const bool vecs_on, const Base<typename T1::pod_type,T1>& expr);
  
  template<typename T1>
  inline static bool eig_gen(Mat< std::complex<typename T1::pod_type> >& vals, Mat< std::complex<typename T1::pod_type> >& vecs, const bool vecs_on, const Base< std::complex<typename T1::pod_type>, T1 >& expr);
  
  
  //
  // eig_gen_balance
  
  template<typename T1>
  inline static bool eig_gen_balance(Mat< std::complex<typename T1::pod_type> >& vals, Mat< std::complex<typename T1::pod_type> >& vecs, const bool vecs_on, const Base<typename T1::pod_type,T1>& expr);
  
  template<typename T1>
  inline static bool eig_gen_balance(Mat< std::complex<typename T1::pod_type> >& vals, Mat< std::complex<typename T1::pod_type> >& vecs, const bool vecs_on, const Base< std::complex<typename T1::pod_type>, T1 >& expr);
  
  
  //
  // eig_gen_twosided
  
  template<typename T1>
  inline static bool eig_gen_twosided(Mat< std::complex<typename T1::pod_type> >& vals, Mat< std::complex<typename T1::pod_type> >& lvecs, Mat< std::complex<typename T1::pod_type> >& rvecs, const Base<typename T1::pod_type,T1>& expr);
  
  template<typename T1>
  inline static bool eig_gen_twosided(Mat< std::complex<typename T1::pod_type> >& vals, Mat< std::complex<typename T1::pod_type> >& lvecs, Mat< std::complex<typename T1::pod_type> >& rvecs, const Base< std::complex<typename T1::pod_type>, T1 >& expr);
  
  
  //
  // eig_gen_twosided_balance
  
  template<typename T1>
  inline static bool eig_gen_twosided_balance(Mat< std::complex<typename T1::pod_type> >& vals, Mat< std::complex<typename T1::pod_type> >& lvecs, Mat< std::complex<typename T1::pod_type> >& rvecs, const Base<typename T1::pod_type,T1>& expr);
  
  template<typename T1>
  inline static bool eig_gen_twosided_balance(Mat< std::complex<typename T1::pod_type> >& vals, Mat< std::complex<typename T1::pod_type> >& lvecs, Mat< std::complex<typename T1::pod_type> >& rvecs, const Base< std::complex<typename T1::pod_type>, T1 >& expr);
  
  
  //
  // eig_pair
  
  template<typename T1, typename T2>
  inline static bool eig_pair(Mat< std::complex<typename T1::pod_type> >& vals, Mat< std::complex<typename T1::pod_type> >& vecs, const bool vecs_on, const Base<typename T1::pod_type,T1>& A_expr, const Base<typename T1::pod_type,T2>& B_expr);
  
  template<typename T1, typename T2>
  inline static bool eig_pair(Mat< std::complex<typename T1::pod_type> >& vals, Mat< std::complex<typename T1::pod_type> >& vecs, const bool vecs_on, const Base< std::complex<typename T1::pod_type>, T1 >& A_expr, const Base< std::complex<typename T1::pod_type>, T2 >& B_expr);
  
  
  //
  // eig_pair_twosided
  
  template<typename T1, typename T2>
  inline static bool eig_pair_twosided(Mat< std::complex<typename T1::pod_type> >& vals, Mat< std::complex<typename T1::pod_type> >& lvecs, Mat< std::complex<typename T1::pod_type> >& rvecs, const Base<typename T1::pod_type,T1>& A_expr, const Base<typename T1::pod_type,T2>& B_expr);
  
  template<typename T1, typename T2>
  inline static bool eig_pair_twosided(Mat< std::complex<typename T1::pod_type> >& vals, Mat< std::complex<typename T1::pod_type> >& lvecs, Mat< std::complex<typename T1::pod_type> >& rvecs, const Base< std::complex<typename T1::pod_type>, T1 >& A_expr, const Base< std::complex<typename T1::pod_type>, T2 >& B_expr);
  
  
  //
  // eig_sym
  
  template<typename eT> 
  inline static bool eig_sym(Col<eT>& eigval, Mat<eT>& A);
  
  template<typename T> 
  inline static bool eig_sym(Col<T>& eigval, Mat< std::complex<T> >& A);
  
  template<typename eT>
  inline static bool eig_sym(Col<eT>& eigval, Mat<eT>& eigvec, const Mat<eT>& X);
  
  template<typename T>
  inline static bool eig_sym(Col<T>& eigval, Mat< std::complex<T> >& eigvec, const Mat< std::complex<T> >& X);
  
  template<typename eT>
  inline static bool eig_sym_dc(Col<eT>& eigval, Mat<eT>& eigvec, const Mat<eT>& X);
  
  template<typename T>
  inline static bool eig_sym_dc(Col<T>& eigval, Mat< std::complex<T> >& eigvec, const Mat< std::complex<T> >& X);
  
  
  //
  // chol
  
  template<typename eT>
  inline static bool chol_simple(Mat<eT>& X);
  
  template<typename eT>
  inline static bool chol(Mat<eT>& X, const uword layout);
  
  template<typename eT>
  inline static bool chol_band(Mat<eT>& X, const uword KD, const uword layout);
  
  template<typename  T>
  inline static bool chol_band(Mat< std::complex<T> >& X, const uword KD, const uword layout);
  
  template<typename eT>
  inline static bool chol_band_common(Mat<eT>& X, const uword KD, const uword layout);
  
  template<typename eT>
  inline static bool chol_pivot(Mat<eT>& X, Mat<uword>& P, const uword layout);
  
  
  //
  // hessenberg decomposition
  
  template<typename eT, typename T1>
  inline static bool hess(Mat<eT>& H, const Base<eT,T1>& X, Col<eT>& tao);
  
  
  //
  // qr
  
  template<typename eT, typename T1>
  inline static bool qr(Mat<eT>& Q, Mat<eT>& R, const Base<eT,T1>& X);
  
  template<typename eT, typename T1>
  inline static bool qr_econ(Mat<eT>& Q, Mat<eT>& R, const Base<eT,T1>& X);
  
  template<typename eT, typename T1>
  inline static bool qr_pivot(Mat<eT>& Q, Mat<eT>& R, Mat<uword>& P, const Base<eT,T1>& X);
  
  template<typename  T, typename T1>
  inline static bool qr_pivot(Mat< std::complex<T> >& Q, Mat< std::complex<T> >& R, Mat<uword>& P, const Base<std::complex<T>,T1>& X);
  
  
  //
  // svd
  
  template<typename eT>
  inline static bool svd(Col<eT>& S, Mat<eT>& A);
  
  template<typename T>
  inline static bool svd(Col<T>& S, Mat< std::complex<T> >& A);
  
  
  template<typename eT>
  inline static bool svd(Mat<eT>& U, Col<eT>& S, Mat<eT>& V, Mat<eT>& A);
  
  template<typename T>
  inline static bool svd(Mat< std::complex<T> >& U, Col<T>& S, Mat< std::complex<T> >& V, Mat< std::complex<T> >& A);
  
  template<typename eT>
  inline static bool svd_econ(Mat<eT>& U, Col<eT>& S, Mat<eT>& V, Mat<eT>& A, const char mode);
  
  template<typename T>
  inline static bool svd_econ(Mat< std::complex<T> >& U, Col<T>& S, Mat< std::complex<T> >& V, Mat< std::complex<T> >& A, const char mode);
  
  
  template<typename eT>
  inline static bool svd_dc(Col<eT>& S, Mat<eT>& A);
  
  template<typename T>
  inline static bool svd_dc(Col<T>& S, Mat< std::complex<T> >& A);
  
  
  template<typename eT>
  inline static bool svd_dc(Mat<eT>& U, Col<eT>& S, Mat<eT>& V, Mat<eT>& A);
  
  template<typename T>
  inline static bool svd_dc(Mat< std::complex<T> >& U, Col<T>& S, Mat< std::complex<T> >& V, Mat< std::complex<T> >& A);
  
  template<typename eT>
  inline static bool svd_dc_econ(Mat<eT>& U, Col<eT>& S, Mat<eT>& V, Mat<eT>& A);
  
  template<typename T>
  inline static bool svd_dc_econ(Mat< std::complex<T> >& U, Col<T>& S, Mat< std::complex<T> >& V, Mat< std::complex<T> >& A);
  
  
  //
  // solve
  
  template<typename T1>
  inline static bool solve_square_fast(Mat<typename T1::elem_type>& out, Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_square_rcond(Mat<typename T1::elem_type>& out, typename T1::pod_type& out_rcond, Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_square_refine(Mat<typename T1::pod_type>& out, typename T1::pod_type& out_rcond, Mat<typename T1::pod_type>& A, const Base<typename T1::pod_type,T1>& B_expr, const bool equilibrate);
  
  template<typename T1>
  inline static bool solve_square_refine(Mat< std::complex<typename T1::pod_type> >& out, typename T1::pod_type& out_rcond, Mat< std::complex<typename T1::pod_type> >& A, const Base<std::complex<typename T1::pod_type>,T1>& B_expr, const bool equilibrate);
  
  //
  
  template<typename T1>
  inline static bool solve_sympd_fast(Mat<typename T1::elem_type>& out, Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_sympd_fast_common(Mat<typename T1::elem_type>& out, Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_sympd_rcond(Mat<typename T1::pod_type>& out, bool& out_sympd_state, typename T1::pod_type& out_rcond, Mat<typename T1::pod_type>& A, const Base<typename T1::pod_type,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_sympd_rcond(Mat< std::complex<typename T1::pod_type> >& out, bool& out_sympd_state, typename T1::pod_type& out_rcond, Mat< std::complex<typename T1::pod_type> >& A, const Base< std::complex<typename T1::pod_type>,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_sympd_refine(Mat<typename T1::pod_type>& out, typename T1::pod_type& out_rcond, Mat<typename T1::pod_type>& A, const Base<typename T1::pod_type,T1>& B_expr, const bool equilibrate);
  
  template<typename T1>
  inline static bool solve_sympd_refine(Mat< std::complex<typename T1::pod_type> >& out, typename T1::pod_type& out_rcond, Mat< std::complex<typename T1::pod_type> >& A, const Base<std::complex<typename T1::pod_type>,T1>& B_expr, const bool equilibrate);
  
  //
  
  template<typename T1>
  inline static bool solve_rect_fast(Mat<typename T1::elem_type>& out, Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_rect_rcond(Mat<typename T1::elem_type>& out, typename T1::pod_type& out_rcond, Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr);
  
  //
  
  template<typename T1>
  inline static bool solve_approx_svd(Mat<typename T1::pod_type>& out, Mat<typename T1::pod_type>& A, const Base<typename T1::pod_type,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_approx_svd(Mat< std::complex<typename T1::pod_type> >& out, Mat< std::complex<typename T1::pod_type> >& A, const Base<std::complex<typename T1::pod_type>,T1>& B_expr);
  
  //
  
  template<typename T1>
  inline static bool solve_trimat_fast(Mat<typename T1::elem_type>& out, const Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr, const uword layout);
  
  template<typename T1>
  inline static bool solve_trimat_rcond(Mat<typename T1::elem_type>& out, typename T1::pod_type& out_rcond, const Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr, const uword layout);
  
  //
  
  template<typename T1>
  inline static bool solve_band_fast(Mat<typename T1::pod_type>& out, Mat<typename T1::pod_type>& A, const uword KL, const uword KU, const Base<typename T1::pod_type,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_band_fast(Mat< std::complex<typename T1::pod_type> >& out, Mat< std::complex<typename T1::pod_type> >& A, const uword KL, const uword KU, const Base< std::complex<typename T1::pod_type>,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_band_fast_common(Mat<typename T1::elem_type>& out, const Mat<typename T1::elem_type>& A, const uword KL, const uword KU, const Base<typename T1::elem_type,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_band_rcond(Mat<typename T1::pod_type>& out, typename T1::pod_type& out_rcond, Mat<typename T1::pod_type>& A, const uword KL, const uword KU, const Base<typename T1::pod_type,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_band_rcond(Mat< std::complex<typename T1::pod_type> >& out, typename T1::pod_type& out_rcond, Mat< std::complex<typename T1::pod_type> >& A, const uword KL, const uword KU, const Base< std::complex<typename T1::pod_type>,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_band_rcond_common(Mat<typename T1::elem_type>& out, typename T1::pod_type& out_rcond, const Mat<typename T1::elem_type>& A, const uword KL, const uword KU, const Base<typename T1::elem_type,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_band_refine(Mat<typename T1::pod_type>& out, typename T1::pod_type& out_rcond, Mat<typename T1::pod_type>& A, const uword KL, const uword KU, const Base<typename T1::pod_type,T1>& B_expr, const bool equilibrate);
  
  template<typename T1>
  inline static bool solve_band_refine(Mat< std::complex<typename T1::pod_type> >& out, typename T1::pod_type& out_rcond, Mat< std::complex<typename T1::pod_type> >& A, const uword KL, const uword KU, const Base<std::complex<typename T1::pod_type>,T1>& B_expr, const bool equilibrate);
  
  //
  
  template<typename T1>
  inline static bool solve_tridiag_fast(Mat<typename T1::pod_type>& out, Mat<typename T1::pod_type>& A, const Base<typename T1::pod_type,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_tridiag_fast(Mat< std::complex<typename T1::pod_type> >& out, Mat< std::complex<typename T1::pod_type> >& A, const Base< std::complex<typename T1::pod_type>,T1>& B_expr);
  
  template<typename T1>
  inline static bool solve_tridiag_fast_common(Mat<typename T1::elem_type>& out, const Mat<typename T1::elem_type>& A, const Base<typename T1::elem_type,T1>& B_expr);
  
  
  //
  // Schur decomposition
  
  template<typename eT, typename T1>
  inline static bool schur(Mat<eT>& U, Mat<eT>& S, const Base<eT,T1>& X, const bool calc_U = true);
  
  template<typename  T, typename T1>
  inline static bool schur(Mat< std::complex<T> >& U, Mat< std::complex<T> >& S, const Base<std::complex<T>,T1>& X, const bool calc_U = true);
  
  template<typename  T>
  inline static bool schur(Mat< std::complex<T> >& U, Mat< std::complex<T> >& S, const bool calc_U = true);
  
  //
  // solve the Sylvester equation AX + XB = C
  
  template<typename eT>
  inline static bool syl(Mat<eT>& X, const Mat<eT>& A, const Mat<eT>& B, const Mat<eT>& C);
  
  
  //
  // QZ decomposition
  
  template<typename T, typename T1, typename T2>
  inline static bool qz(Mat<T>& A, Mat<T>& B, Mat<T>& vsl, Mat<T>& vsr, const Base<T,T1>& X_expr, const Base<T,T2>& Y_expr, const char mode);
  
  template<typename T, typename T1, typename T2>
  inline static bool qz(Mat< std::complex<T> >& A, Mat< std::complex<T> >& B, Mat< std::complex<T> >& vsl, Mat< std::complex<T> >& vsr, const Base< std::complex<T>, T1 >& X_expr, const Base< std::complex<T>, T2 >& Y_expr, const char mode);
  
  
  // 
  // rcond
  
  template<typename eT>
  inline static eT rcond(Mat<eT>& A);
  
  template<typename  T>
  inline static  T rcond(Mat< std::complex<T> >& A);
  
  template<typename eT>
  inline static eT rcond_sympd(Mat<eT>& A, bool& calc_ok);
  
  template<typename T>
  inline static  T rcond_sympd(Mat< std::complex<T> >& A, bool& calc_ok);
  
  template<typename eT>
  inline static eT rcond_trimat(const Mat<eT>& A, const uword layout);
  
  template<typename  T>
  inline static  T rcond_trimat(const Mat< std::complex<T> >& A, const uword layout);
  
  
  //
  // lu_rcond (rcond from pre-computed LU decomposition)
  
  template<typename eT>
  inline static eT lu_rcond(const Mat<eT>& A, const eT norm_val);
  
  template<typename  T>
  inline static  T lu_rcond(const Mat< std::complex<T> >& A, const T norm_val);
  
  template<typename eT>
  inline static eT lu_rcond_sympd(const Mat<eT>& A, const eT norm_val);
  
  template<typename  T>
  inline static  T lu_rcond_sympd(const Mat< std::complex<T> >& A, const T norm_val);
  
  template<typename eT>
  inline static eT lu_rcond_band(const Mat<eT>& AB, const uword KL, const uword KU, const podarray<blas_int>& ipiv, const eT norm_val);
  
  template<typename  T>
  inline static  T lu_rcond_band(const Mat< std::complex<T> >& AB, const uword KL, const uword KU, const podarray<blas_int>& ipiv, const T norm_val);
  
  
  //
  // misc
  
  template<typename T1>
  inline static bool crippled_lapack(const Base<typename T1::elem_type, T1>&);
  
  template<typename eT>
  inline static bool rudimentary_sym_check(const Mat<eT>& X);
  
  template<typename T>
  inline static bool rudimentary_sym_check(const Mat< std::complex<T> >& X);
  };



namespace qz_helper
  {
  template<typename T> inline blas_int select_lhp(const T* x_ptr, const T* y_ptr, const T* z_ptr);
  template<typename T> inline blas_int select_rhp(const T* x_ptr, const T* y_ptr, const T* z_ptr);
  template<typename T> inline blas_int select_iuc(const T* x_ptr, const T* y_ptr, const T* z_ptr);
  template<typename T> inline blas_int select_ouc(const T* x_ptr, const T* y_ptr, const T* z_ptr);
  
  template<typename T> inline blas_int cx_select_lhp(const std::complex<T>* x_ptr, const std::complex<T>* y_ptr);
  template<typename T> inline blas_int cx_select_rhp(const std::complex<T>* x_ptr, const std::complex<T>* y_ptr);
  template<typename T> inline blas_int cx_select_iuc(const std::complex<T>* x_ptr, const std::complex<T>* y_ptr);
  template<typename T> inline blas_int cx_select_ouc(const std::complex<T>* x_ptr, const std::complex<T>* y_ptr);
  
  template<typename T> inline void_ptr ptr_cast(blas_int (*function)(const T*, const T*, const T*));
  template<typename T> inline void_ptr ptr_cast(blas_int (*function)(const std::complex<T>*, const std::complex<T>*));
  }



//! @}
