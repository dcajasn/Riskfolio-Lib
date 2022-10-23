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



//! \addtogroup op_rank
//! @{



template<typename T1>
inline
bool
op_rank::apply(uword& out, const Base<typename T1::elem_type,T1>& expr, const typename T1::pod_type tol)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> A(expr.get_ref());
  
  if(A.is_empty())  { out = uword(0); return true; }
  
  if(is_op_diagmat<T1>::value || A.is_diagmat())
    {
    arma_extra_debug_print("op_rank::apply(): detected diagonal matrix");
    
    return op_rank::apply_diag(out, A, tol);
    }
  
  bool do_sym = false;
  
  if((arma_config::optimise_sympd) && (auxlib::crippled_lapack(A) == false) && (A.n_rows >= (is_cx<eT>::yes ? uword(64) : uword(128))))
    {
    bool is_approx_sym   = false;
    bool is_approx_sympd = false;
    
    sympd_helper::analyse_matrix(is_approx_sym, is_approx_sympd, A);
    
    do_sym = (is_cx<eT>::no) ? (is_approx_sym) : (is_approx_sym && is_approx_sympd);
    }
  
  if(do_sym)
    {
    arma_extra_debug_print("op_rank::apply(): symmetric/hermitian optimisation");
    
    return op_rank::apply_sym(out, A, tol);
    }
  
  return op_rank::apply_gen(out, A, tol);
  }



template<typename eT>
inline
bool
op_rank::apply_diag(uword& out, Mat<eT>& A, typename get_pod_type<eT>::result tol)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const uword N = (std::min)(A.n_rows, A.n_cols);
  
  podarray<T> diag_abs_vals(N);
  
  T max_abs_Aii = T(0);
  
  for(uword i=0; i<N; ++i)
    {
    const eT     Aii = A.at(i,i);
    const  T abs_Aii = std::abs(Aii);
    
    if(arma_isnan(Aii))  { out = uword(0); return false; }
    
    diag_abs_vals[i] = abs_Aii;
    
    max_abs_Aii = (abs_Aii > max_abs_Aii) ? abs_Aii : max_abs_Aii;
    }
  
  // set tolerance to default if it hasn't been specified
  if(tol == T(0))  { tol = (std::max)(A.n_rows, A.n_cols) * max_abs_Aii * std::numeric_limits<T>::epsilon(); }
  
  uword count = 0;
  
  for(uword i=0; i<N; ++i)  { count += (diag_abs_vals[i] > tol) ? uword(1) : uword(0); }
  
  out = count;
  
  return true;
  }



template<typename eT>
inline
bool
op_rank::apply_sym(uword& out, Mat<eT>& A, typename get_pod_type<eT>::result tol)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  if(A.is_square() == false)  { out = uword(0); return false; }
  
  Col<T> v;
  
  const bool status = auxlib::eig_sym(v, A);
  
  if(status == false)  { out = uword(0); return false; }
  
  const uword v_n_elem = v.n_elem;
        T*    v_mem    = v.memptr();
  
  if(v_n_elem == 0)  { out = uword(0); return true; }
  
  T max_abs_v = T(0);
  
  for(uword i=0; i < v_n_elem; ++i)  { const T val = std::abs(v_mem[i]); v_mem[i] = val; if(val > max_abs_v) { max_abs_v = val; } }
  
  // set tolerance to default if it hasn't been specified
  if(tol == T(0))  { tol = (std::max)(A.n_rows, A.n_cols) * max_abs_v * std::numeric_limits<T>::epsilon(); }
  
  uword count = 0;
  
  for(uword i=0; i < v_n_elem; ++i)  { count += (v_mem[i] > tol) ? uword(1) : uword(0); }
  
  out = count;
  
  return true;
  }



template<typename eT>
inline
bool
op_rank::apply_gen(uword& out, Mat<eT>& A, typename get_pod_type<eT>::result tol)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  Col<T> s;
  
  const bool status = auxlib::svd_dc(s, A);
  
  if(status == false)  { out = uword(0); return false; }
  
  const uword s_n_elem = s.n_elem;
  const T*    s_mem    = s.memptr();
  
  if(s_n_elem == 0)  { out = uword(0); return true; }
  
  // set tolerance to default if it hasn't been specified
  if(tol == T(0))  { tol = (std::max)(A.n_rows, A.n_cols) * s_mem[0] * std::numeric_limits<T>::epsilon(); }
  
  uword count = 0;
  
  for(uword i=0; i < s_n_elem; ++i)  { count += (s_mem[i] > tol) ? uword(1) : uword(0); }
  
  out = count;
  
  return true;
  }



//! @}
