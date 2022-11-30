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



//! \addtogroup op_pinv
//! @{



template<typename T1>
inline
void
op_pinv_default::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_pinv_default>& in)
  {
  arma_extra_debug_sigprint();
  
  const bool status = op_pinv_default::apply_direct(out, in.m);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("pinv(): svd failed");
    }
  }



template<typename T1>
inline
bool
op_pinv_default::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  constexpr T     tol       = T(0);
  constexpr uword method_id = uword(0);
  
  return op_pinv::apply_direct(out, expr, tol, method_id);
  }



//



template<typename T1>
inline
void
op_pinv::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_pinv>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const T     tol       = access::tmp_real(in.aux);
  const uword method_id = in.aux_uword_a;
  
  const bool status = op_pinv::apply_direct(out, in.m, tol, method_id);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("pinv(): svd failed");
    }
  }



template<typename T1>
inline
bool
op_pinv::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, typename T1::pod_type tol, const uword method_id)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  arma_debug_check((tol < T(0)), "pinv(): tolerance must be >= 0");
  
  // method_id = 0 -> default setting
  // method_id = 1 -> use standard algorithm
  // method_id = 2 -> use divide and conquer algorithm
  
  Mat<eT> A(expr.get_ref());
  
  if(A.is_empty())  { out.set_size(A.n_cols,A.n_rows); return true; }
  
  if(is_op_diagmat<T1>::value || A.is_diagmat())
    {
    arma_extra_debug_print("op_pinv: detected diagonal matrix");
    
    return op_pinv::apply_diag(out, A, tol);
    }
  
  bool do_sym   = false;
  bool do_sympd = false;
  
  const bool is_sym_size_ok = (A.n_rows > (is_cx<eT>::yes ? uword(20) : uword(40)));
  const bool is_arg_default = ((tol == T(0)) && (method_id == uword(0)));
  
  if( (arma_config::optimise_sympd) && (auxlib::crippled_lapack(A) == false) && (is_arg_default || is_sym_size_ok) )
    {
    bool is_approx_sym   = false;
    bool is_approx_sympd = false;
    
    sympd_helper::analyse_matrix(is_approx_sym, is_approx_sympd, A);
    
    do_sym   = is_sym_size_ok && ((is_cx<eT>::no) ? (is_approx_sym) : (is_approx_sym && is_approx_sympd));
    do_sympd = is_arg_default && is_approx_sympd;
    }
  
  if(do_sympd)
    {
    arma_extra_debug_print("op_pinv: attempting sympd optimisation");
    
    out = A;
    
          bool is_sympd_junk   = false;
          T    rcond_calc      = T(0);
    const T    rcond_threshold = T((std::max)(uword(100), uword(A.n_rows))) * std::numeric_limits<T>::epsilon();
    
    const bool status = auxlib::inv_sympd_rcond(out, is_sympd_junk, rcond_calc, rcond_threshold);
    
    if(status && arma_isfinite(rcond_calc))  { return true; }
    
    arma_extra_debug_print("op_pinv: sympd optimisation failed");
    // auxlib::inv_sympd_rcond() will fail if A isn't really positive definite or its rcond is below rcond_threshold
    }
  
  if(do_sym)
    {
    arma_extra_debug_print("op_pinv: symmetric/hermitian optimisation");
    
    return op_pinv::apply_sym(out, A, tol, method_id);
    }
  
  return op_pinv::apply_gen(out, A, tol, method_id);
  }



template<typename eT>
inline
bool
op_pinv::apply_diag(Mat<eT>& out, const Mat<eT>& A, typename get_pod_type<eT>::result tol)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  out.zeros(A.n_cols, A.n_rows);
  
  const uword N = (std::min)(A.n_rows, A.n_cols);
  
  podarray<T> diag_abs_vals(N);
  
  T max_abs_Aii = T(0);
  
  for(uword i=0; i<N; ++i)
    {
    const eT     Aii = A.at(i,i);
    const  T abs_Aii = std::abs(Aii);
    
    if(arma_isnan(Aii))  { return false; }
    
    diag_abs_vals[i] = abs_Aii;
    
    max_abs_Aii = (abs_Aii > max_abs_Aii) ? abs_Aii : max_abs_Aii;
    }
  
  if(tol == T(0))  { tol = (std::max)(A.n_rows, A.n_cols) * max_abs_Aii * std::numeric_limits<T>::epsilon(); }
  
  for(uword i=0; i<N; ++i)
    {
    if(diag_abs_vals[i] >= tol)
      {
      const eT Aii = A.at(i,i);
      
      if(Aii != eT(0))  { out.at(i,i) = eT(eT(1) / Aii); }
      }
    }
  
  return true;
  }



template<typename eT>
inline
bool
op_pinv::apply_sym(Mat<eT>& out, const Mat<eT>& A, typename get_pod_type<eT>::result tol, const uword method_id)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  Col< T> eigval;
  Mat<eT> eigvec;
  
  const bool status = ((method_id == uword(0)) || (method_id == uword(2))) ? auxlib::eig_sym_dc(eigval, eigvec, A) : auxlib::eig_sym(eigval, eigvec, A);
  
  if(status == false)  { return false; }
  
  if(eigval.n_elem == 0)  { out.zeros(A.n_cols, A.n_rows); return true; }
  
  Col<T> abs_eigval = arma::abs(eigval);
  
  const uvec indices = sort_index(abs_eigval, "descend");
  
  abs_eigval = abs_eigval.elem(indices);
      eigval =     eigval.elem(indices);
      eigvec =     eigvec.cols(indices);
  
  // set tolerance to default if it hasn't been specified
  if(tol == T(0))  { tol = (std::max)(A.n_rows, A.n_cols) * abs_eigval[0] * std::numeric_limits<T>::epsilon(); }
  
  uword count = 0;
  
  for(uword i=0; i < abs_eigval.n_elem; ++i)  { count += (abs_eigval[i] >= tol) ? uword(1) : uword(0); }
  
  if(count == 0)  { out.zeros(A.n_cols, A.n_rows); return true; }
  
  Col<T> eigval2(count, arma_nozeros_indicator());
  
  uword count2 = 0;
  
  for(uword i=0; i < eigval.n_elem; ++i)
    {
    const T abs_val = abs_eigval[i];
    const T     val =     eigval[i];
    
    if(abs_val >= tol)  { eigval2[count2] = (val != T(0)) ? T(T(1) / val) : T(0); ++count2; }
    }
  
  const Mat<eT> eigvec_use(eigvec.memptr(), eigvec.n_rows, count, false);
  
  out = (eigvec_use * diagmat(eigval2)).eval() * eigvec_use.t();
  
  return true;
  }




template<typename eT>
inline
bool
op_pinv::apply_gen(Mat<eT>& out, Mat<eT>& A, typename get_pod_type<eT>::result tol, const uword method_id)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const uword n_rows = A.n_rows;
  const uword n_cols = A.n_cols;
  
  // economical SVD decomposition 
  Mat<eT> U;
  Col< T> s;
  Mat<eT> V;
  
  if(n_cols > n_rows)  { A = trans(A); }
  
  const bool status = ((method_id == uword(0)) || (method_id == uword(2))) ? auxlib::svd_dc_econ(U, s, V, A) : auxlib::svd_econ(U, s, V, A, 'b');
  
  if(status == false)  { return false; }
  
  // set tolerance to default if it hasn't been specified
  if( (tol == T(0)) && (s.n_elem > 0) )  { tol = (std::max)(n_rows, n_cols) * s[0] * std::numeric_limits<T>::epsilon(); }
  
  uword count = 0;
  
  for(uword i=0; i < s.n_elem; ++i)  { count += (s[i] >= tol) ? uword(1) : uword(0); }
  
  if(count == 0)  { out.zeros(n_cols, n_rows); return true; }
  
  Col<T> s2(count, arma_nozeros_indicator());
  
  uword count2 = 0;
  
  for(uword i=0; i < s.n_elem; ++i)
    {
    const T val = s[i];
    
    if(val >= tol)  { s2[count2] = (val > T(0)) ? T(T(1) / val) : T(0); ++count2; }
    }
  
  const Mat<eT> U_use(U.memptr(), U.n_rows, count, false);
  const Mat<eT> V_use(V.memptr(), V.n_rows, count, false);
  
  Mat<eT> tmp;
  
  if(n_rows >= n_cols)
    {
    // out = ( (V.n_cols > count) ? V.cols(0,count-1) : V ) * diagmat(s2) * trans( (U.n_cols > count) ? U.cols(0,count-1) : U );
    
    tmp = V_use * diagmat(s2);
    
    out = tmp * trans(U_use);
    }
  else
    {
    // out = ( (U.n_cols > count) ? U.cols(0,count-1) : U ) * diagmat(s2) * trans( (V.n_cols > count) ? V.cols(0,count-1) : V );
    
    tmp = U_use * diagmat(s2);
    
    out = tmp * trans(V_use);
    }
  
  return true;
  }



//! @}
