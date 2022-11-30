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


//! \addtogroup op_inv_gen
//! @{



template<typename T1>
inline
void
op_inv_gen_default::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_gen_default>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = op_inv_gen_default::apply_direct(out, X.m, "inv()");
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("inv(): matrix is singular");
    }
  }



template<typename T1>
inline
bool
op_inv_gen_default::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, const char* caller_sig)
  {
  arma_extra_debug_sigprint();
  
  return op_inv_gen_full::apply_direct<T1,false>(out, expr, caller_sig, uword(0));
  }



//



template<typename T1>
inline
void
op_inv_gen_full::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_gen_full>& X)
  {
  arma_extra_debug_sigprint();
  
  const uword flags = X.aux_uword_a;
  
  const bool status = op_inv_gen_full::apply_direct(out, X.m, "inv()", flags);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("inv(): matrix is singular");
    }
  }



template<typename T1, const bool has_user_flags>
inline
bool
op_inv_gen_full::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, const char* caller_sig, const uword flags)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  if(has_user_flags == true )  { arma_extra_debug_print("op_inv_gen_full: has_user_flags = true");  }
  if(has_user_flags == false)  { arma_extra_debug_print("op_inv_gen_full: has_user_flags = false"); }
  
  const bool tiny         = has_user_flags && bool(flags & inv_opts::flag_tiny        );
  const bool allow_approx = has_user_flags && bool(flags & inv_opts::flag_allow_approx);
  const bool likely_sympd = has_user_flags && bool(flags & inv_opts::flag_likely_sympd);
  const bool no_sympd     = has_user_flags && bool(flags & inv_opts::flag_no_sympd    );
  const bool no_ugly      = has_user_flags && bool(flags & inv_opts::flag_no_ugly     );
  
  if(has_user_flags)
    {
    arma_extra_debug_print("op_inv_gen_full: enabled flags:");
    
    if(tiny        )  { arma_extra_debug_print("tiny");         }
    if(allow_approx)  { arma_extra_debug_print("allow_approx"); }
    if(likely_sympd)  { arma_extra_debug_print("likely_sympd"); }
    if(no_sympd    )  { arma_extra_debug_print("no_sympd");     }
    if(no_ugly     )  { arma_extra_debug_print("no_ugly");      }
    
    arma_debug_check( (no_sympd && likely_sympd), "inv(): options 'no_sympd' and 'likely_sympd' are mutually exclusive" );
    arma_debug_check( (no_ugly  && allow_approx), "inv(): options 'no_ugly' and 'allow_approx' are mutually exclusive"  );
    }
  
  if(no_ugly)
    {
    op_inv_gen_state<T> inv_state;
    
    const bool status = op_inv_gen_rcond::apply_direct(out, inv_state, expr);
    
    const T local_rcond = inv_state.rcond;  // workaround for bug in gcc 4.8
    
    if((status == false) || (local_rcond < std::numeric_limits<T>::epsilon()) || arma_isnan(local_rcond))  { return false; }
    
    return true;
    }
  
  if(allow_approx)
    {
    op_inv_gen_state<T> inv_state;
    
    Mat<eT> tmp;
    
    const bool status = op_inv_gen_rcond::apply_direct(tmp, inv_state, expr);
    
    const T local_rcond = inv_state.rcond;  // workaround for bug in gcc 4.8

    if((status == false) || (local_rcond < std::numeric_limits<T>::epsilon()) || arma_isnan(local_rcond))
      {
      Mat<eT> A = expr.get_ref();
      
      if(inv_state.is_diag)  { return op_pinv::apply_diag(out, A, T(0)          ); }
      if(inv_state.is_sym )  { return op_pinv::apply_sym (out, A, T(0), uword(0)); }
      
      return op_pinv::apply_gen(out, A, T(0), uword(0));
      }
    
    out.steal_mem(tmp);
    
    return true;
    }
  
  out = expr.get_ref();
  
  arma_debug_check( (out.is_square() == false), caller_sig, ": given matrix must be square sized", [&](){ out.soft_reset(); } );
  
  const uword N = out.n_rows;
  
  if(N == 0)  { return true; }
  
  if(is_cx<eT>::no)
    {
    if(N == 1)
      {
      const eT a = out[0];
      
      out[0] = eT(1) / a;
      
      return (a != eT(0));
      }
    else
    if(N == 2)
      {
      const bool status = op_inv_gen_full::apply_tiny_2x2(out);
      
      if(status)  { return true; }
      }
    else
    if((N == 3) && tiny)
      {
      const bool status = op_inv_gen_full::apply_tiny_3x3(out);
      
      if(status)  { return true; }
      }
    else
    if((N == 4) && tiny)
      {
      const bool status = op_inv_gen_full::apply_tiny_4x4(out);
      
      if(status)  { return true; }
      }
    
    // fallthrough if optimisation failed
    }
  
  if(is_op_diagmat<T1>::value || out.is_diagmat())
    {
    arma_extra_debug_print("op_inv_gen_full: detected diagonal matrix");
    
    eT* colmem = out.memptr();
    
    for(uword i=0; i<N; ++i)
      {
      eT& out_ii = colmem[i];
      
      const eT src_val = out_ii;
      const eT inv_val = eT(1) / src_val;
      
      if(src_val == eT(0))  { return false; }
      
      out_ii = inv_val;
      
      colmem += N;
      }
    
    return true;
    }
  
  const strip_trimat<T1> strip(expr.get_ref());
  
  const bool is_triu_expr = strip.do_triu;
  const bool is_tril_expr = strip.do_tril;
  
  const bool is_triu_mat = (is_triu_expr || is_tril_expr) ? false : (                        trimat_helper::is_triu(out));
  const bool is_tril_mat = (is_triu_expr || is_tril_expr) ? false : ((is_triu_mat) ? false : trimat_helper::is_tril(out));
  
  if(is_triu_expr || is_tril_expr || is_triu_mat || is_tril_mat)
    {
    return auxlib::inv_tr(out, ((is_triu_expr || is_triu_mat) ? uword(0) : uword(1)));
    }
  
  const bool try_sympd = arma_config::optimise_sympd && ((no_sympd) ? false : (likely_sympd ? true : sympd_helper::guess_sympd(out)));
  
  if(try_sympd)
    {
    arma_extra_debug_print("op_inv_gen_full: attempting sympd optimisation");
    
    Mat<eT> tmp = out;
    
    bool sympd_state = false;
    
    const bool status = auxlib::inv_sympd(tmp, sympd_state);
    
    if(status)  { out.steal_mem(tmp); return true; }
    
    if((status == false) && (sympd_state == true))  { return false; }
    
    arma_extra_debug_print("op_inv_gen_full: sympd optimisation failed");
    
    // fallthrough if optimisation failed
    }
  
  return auxlib::inv(out);
  }



template<typename eT>
arma_cold
inline
bool
op_inv_gen_full::apply_tiny_2x2(Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  // NOTE: assuming matrix X is square sized
  
  constexpr T det_min =        std::numeric_limits<T>::epsilon();
  constexpr T det_max = T(1) / std::numeric_limits<T>::epsilon();
  
  eT* Xm = X.memptr();
  
  const eT a = Xm[pos<0,0>::n2];
  const eT b = Xm[pos<0,1>::n2];
  const eT c = Xm[pos<1,0>::n2];
  const eT d = Xm[pos<1,1>::n2];
  
  const eT     det_val = (a*d - b*c);
  const  T abs_det_val = std::abs(det_val);
  
  if((abs_det_val < det_min) || (abs_det_val > det_max))  { return false; }
  
  Xm[pos<0,0>::n2] =  d / det_val;
  Xm[pos<0,1>::n2] = -b / det_val;
  Xm[pos<1,0>::n2] = -c / det_val;
  Xm[pos<1,1>::n2] =  a / det_val;
  
  return true;
  }



template<typename eT>
arma_cold
inline
bool
op_inv_gen_full::apply_tiny_3x3(Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  // NOTE: assuming matrix X is square sized
  
  constexpr T det_min =        std::numeric_limits<T>::epsilon();
  constexpr T det_max = T(1) / std::numeric_limits<T>::epsilon();
  
  Mat<eT> Y(3, 3, arma_nozeros_indicator());
  
  eT* Xm = X.memptr();
  eT* Ym = Y.memptr();
  
  const eT     det_val = op_det::apply_tiny_3x3(X);
  const  T abs_det_val = std::abs(det_val);
  
  if((abs_det_val < det_min) || (abs_det_val > det_max))  { return false; }
  
  Ym[pos<0,0>::n3] =  (Xm[pos<2,2>::n3]*Xm[pos<1,1>::n3] - Xm[pos<2,1>::n3]*Xm[pos<1,2>::n3]) / det_val;
  Ym[pos<1,0>::n3] = -(Xm[pos<2,2>::n3]*Xm[pos<1,0>::n3] - Xm[pos<2,0>::n3]*Xm[pos<1,2>::n3]) / det_val;
  Ym[pos<2,0>::n3] =  (Xm[pos<2,1>::n3]*Xm[pos<1,0>::n3] - Xm[pos<2,0>::n3]*Xm[pos<1,1>::n3]) / det_val;
  
  Ym[pos<0,1>::n3] = -(Xm[pos<2,2>::n3]*Xm[pos<0,1>::n3] - Xm[pos<2,1>::n3]*Xm[pos<0,2>::n3]) / det_val;
  Ym[pos<1,1>::n3] =  (Xm[pos<2,2>::n3]*Xm[pos<0,0>::n3] - Xm[pos<2,0>::n3]*Xm[pos<0,2>::n3]) / det_val;
  Ym[pos<2,1>::n3] = -(Xm[pos<2,1>::n3]*Xm[pos<0,0>::n3] - Xm[pos<2,0>::n3]*Xm[pos<0,1>::n3]) / det_val;
  
  Ym[pos<0,2>::n3] =  (Xm[pos<1,2>::n3]*Xm[pos<0,1>::n3] - Xm[pos<1,1>::n3]*Xm[pos<0,2>::n3]) / det_val;
  Ym[pos<1,2>::n3] = -(Xm[pos<1,2>::n3]*Xm[pos<0,0>::n3] - Xm[pos<1,0>::n3]*Xm[pos<0,2>::n3]) / det_val;
  Ym[pos<2,2>::n3] =  (Xm[pos<1,1>::n3]*Xm[pos<0,0>::n3] - Xm[pos<1,0>::n3]*Xm[pos<0,1>::n3]) / det_val;
  
  const eT check_val = Xm[pos<0,0>::n3]*Ym[pos<0,0>::n3] + Xm[pos<0,1>::n3]*Ym[pos<1,0>::n3] + Xm[pos<0,2>::n3]*Ym[pos<2,0>::n3];
  
  const  T max_diff  = (is_float<T>::value) ? T(1e-4) : T(1e-10);  // empirically determined; may need tuning
  
  if(std::abs(T(1) - check_val) >= max_diff)  { return false; }
  
  arrayops::copy(Xm, Ym, uword(3*3));
  
  return true;
  }



template<typename eT>
arma_cold
inline
bool
op_inv_gen_full::apply_tiny_4x4(Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  // NOTE: assuming matrix X is square sized
  
  constexpr T det_min =        std::numeric_limits<T>::epsilon();
  constexpr T det_max = T(1) / std::numeric_limits<T>::epsilon();
  
  Mat<eT> Y(4, 4, arma_nozeros_indicator());
  
  eT* Xm = X.memptr();
  eT* Ym = Y.memptr();
  
  const eT     det_val = op_det::apply_tiny_4x4(X);
  const  T abs_det_val = std::abs(det_val);
  
  if((abs_det_val < det_min) || (abs_det_val > det_max))  { return false; }
  
  Ym[pos<0,0>::n4] = ( Xm[pos<1,2>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,1>::n4] - Xm[pos<1,3>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,1>::n4] + Xm[pos<1,3>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,2>::n4] - Xm[pos<1,1>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,2>::n4] - Xm[pos<1,2>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,3>::n4] + Xm[pos<1,1>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,3>::n4] ) / det_val;
  Ym[pos<1,0>::n4] = ( Xm[pos<1,3>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,0>::n4] - Xm[pos<1,2>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,0>::n4] - Xm[pos<1,3>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,2>::n4] + Xm[pos<1,0>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,2>::n4] + Xm[pos<1,2>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,3>::n4] - Xm[pos<1,0>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,3>::n4] ) / det_val;
  Ym[pos<2,0>::n4] = ( Xm[pos<1,1>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,0>::n4] - Xm[pos<1,3>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,0>::n4] + Xm[pos<1,3>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,1>::n4] - Xm[pos<1,0>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,1>::n4] - Xm[pos<1,1>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,3>::n4] + Xm[pos<1,0>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,3>::n4] ) / det_val;
  Ym[pos<3,0>::n4] = ( Xm[pos<1,2>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,0>::n4] - Xm[pos<1,1>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,0>::n4] - Xm[pos<1,2>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,1>::n4] + Xm[pos<1,0>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,1>::n4] + Xm[pos<1,1>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,2>::n4] - Xm[pos<1,0>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,2>::n4] ) / det_val;
  
  Ym[pos<0,1>::n4] = ( Xm[pos<0,3>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,2>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,3>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,2>::n4] + Xm[pos<0,1>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,2>::n4] + Xm[pos<0,2>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,3>::n4] - Xm[pos<0,1>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,3>::n4] ) / det_val;
  Ym[pos<1,1>::n4] = ( Xm[pos<0,2>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,3>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,0>::n4] + Xm[pos<0,3>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,2>::n4] - Xm[pos<0,0>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,2>::n4] - Xm[pos<0,2>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,3>::n4] + Xm[pos<0,0>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,3>::n4] ) / det_val;
  Ym[pos<2,1>::n4] = ( Xm[pos<0,3>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,1>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,3>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,1>::n4] + Xm[pos<0,0>::n4]*Xm[pos<2,3>::n4]*Xm[pos<3,1>::n4] + Xm[pos<0,1>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,3>::n4] - Xm[pos<0,0>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,3>::n4] ) / det_val;
  Ym[pos<3,1>::n4] = ( Xm[pos<0,1>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,2>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,0>::n4] + Xm[pos<0,2>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,0>::n4]*Xm[pos<2,2>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,1>::n4]*Xm[pos<2,0>::n4]*Xm[pos<3,2>::n4] + Xm[pos<0,0>::n4]*Xm[pos<2,1>::n4]*Xm[pos<3,2>::n4] ) / det_val;
  
  Ym[pos<0,2>::n4] = ( Xm[pos<0,2>::n4]*Xm[pos<1,3>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,3>::n4]*Xm[pos<1,2>::n4]*Xm[pos<3,1>::n4] + Xm[pos<0,3>::n4]*Xm[pos<1,1>::n4]*Xm[pos<3,2>::n4] - Xm[pos<0,1>::n4]*Xm[pos<1,3>::n4]*Xm[pos<3,2>::n4] - Xm[pos<0,2>::n4]*Xm[pos<1,1>::n4]*Xm[pos<3,3>::n4] + Xm[pos<0,1>::n4]*Xm[pos<1,2>::n4]*Xm[pos<3,3>::n4] ) / det_val;
  Ym[pos<1,2>::n4] = ( Xm[pos<0,3>::n4]*Xm[pos<1,2>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,2>::n4]*Xm[pos<1,3>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,3>::n4]*Xm[pos<1,0>::n4]*Xm[pos<3,2>::n4] + Xm[pos<0,0>::n4]*Xm[pos<1,3>::n4]*Xm[pos<3,2>::n4] + Xm[pos<0,2>::n4]*Xm[pos<1,0>::n4]*Xm[pos<3,3>::n4] - Xm[pos<0,0>::n4]*Xm[pos<1,2>::n4]*Xm[pos<3,3>::n4] ) / det_val;
  Ym[pos<2,2>::n4] = ( Xm[pos<0,1>::n4]*Xm[pos<1,3>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,3>::n4]*Xm[pos<1,1>::n4]*Xm[pos<3,0>::n4] + Xm[pos<0,3>::n4]*Xm[pos<1,0>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,0>::n4]*Xm[pos<1,3>::n4]*Xm[pos<3,1>::n4] - Xm[pos<0,1>::n4]*Xm[pos<1,0>::n4]*Xm[pos<3,3>::n4] + Xm[pos<0,0>::n4]*Xm[pos<1,1>::n4]*Xm[pos<3,3>::n4] ) / det_val;
  Ym[pos<3,2>::n4] = ( Xm[pos<0,2>::n4]*Xm[pos<1,1>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,1>::n4]*Xm[pos<1,2>::n4]*Xm[pos<3,0>::n4] - Xm[pos<0,2>::n4]*Xm[pos<1,0>::n4]*Xm[pos<3,1>::n4] + Xm[pos<0,0>::n4]*Xm[pos<1,2>::n4]*Xm[pos<3,1>::n4] + Xm[pos<0,1>::n4]*Xm[pos<1,0>::n4]*Xm[pos<3,2>::n4] - Xm[pos<0,0>::n4]*Xm[pos<1,1>::n4]*Xm[pos<3,2>::n4] ) / det_val;
  
  Ym[pos<0,3>::n4] = ( Xm[pos<0,3>::n4]*Xm[pos<1,2>::n4]*Xm[pos<2,1>::n4] - Xm[pos<0,2>::n4]*Xm[pos<1,3>::n4]*Xm[pos<2,1>::n4] - Xm[pos<0,3>::n4]*Xm[pos<1,1>::n4]*Xm[pos<2,2>::n4] + Xm[pos<0,1>::n4]*Xm[pos<1,3>::n4]*Xm[pos<2,2>::n4] + Xm[pos<0,2>::n4]*Xm[pos<1,1>::n4]*Xm[pos<2,3>::n4] - Xm[pos<0,1>::n4]*Xm[pos<1,2>::n4]*Xm[pos<2,3>::n4] ) / det_val;
  Ym[pos<1,3>::n4] = ( Xm[pos<0,2>::n4]*Xm[pos<1,3>::n4]*Xm[pos<2,0>::n4] - Xm[pos<0,3>::n4]*Xm[pos<1,2>::n4]*Xm[pos<2,0>::n4] + Xm[pos<0,3>::n4]*Xm[pos<1,0>::n4]*Xm[pos<2,2>::n4] - Xm[pos<0,0>::n4]*Xm[pos<1,3>::n4]*Xm[pos<2,2>::n4] - Xm[pos<0,2>::n4]*Xm[pos<1,0>::n4]*Xm[pos<2,3>::n4] + Xm[pos<0,0>::n4]*Xm[pos<1,2>::n4]*Xm[pos<2,3>::n4] ) / det_val;
  Ym[pos<2,3>::n4] = ( Xm[pos<0,3>::n4]*Xm[pos<1,1>::n4]*Xm[pos<2,0>::n4] - Xm[pos<0,1>::n4]*Xm[pos<1,3>::n4]*Xm[pos<2,0>::n4] - Xm[pos<0,3>::n4]*Xm[pos<1,0>::n4]*Xm[pos<2,1>::n4] + Xm[pos<0,0>::n4]*Xm[pos<1,3>::n4]*Xm[pos<2,1>::n4] + Xm[pos<0,1>::n4]*Xm[pos<1,0>::n4]*Xm[pos<2,3>::n4] - Xm[pos<0,0>::n4]*Xm[pos<1,1>::n4]*Xm[pos<2,3>::n4] ) / det_val;
  Ym[pos<3,3>::n4] = ( Xm[pos<0,1>::n4]*Xm[pos<1,2>::n4]*Xm[pos<2,0>::n4] - Xm[pos<0,2>::n4]*Xm[pos<1,1>::n4]*Xm[pos<2,0>::n4] + Xm[pos<0,2>::n4]*Xm[pos<1,0>::n4]*Xm[pos<2,1>::n4] - Xm[pos<0,0>::n4]*Xm[pos<1,2>::n4]*Xm[pos<2,1>::n4] - Xm[pos<0,1>::n4]*Xm[pos<1,0>::n4]*Xm[pos<2,2>::n4] + Xm[pos<0,0>::n4]*Xm[pos<1,1>::n4]*Xm[pos<2,2>::n4] ) / det_val;
  
  const eT check_val = Xm[pos<0,0>::n4]*Ym[pos<0,0>::n4] + Xm[pos<0,1>::n4]*Ym[pos<1,0>::n4] + Xm[pos<0,2>::n4]*Ym[pos<2,0>::n4] + Xm[pos<0,3>::n4]*Ym[pos<3,0>::n4];
  
  const  T max_diff  = (is_float<T>::value) ? T(1e-4) : T(1e-10);  // empirically determined; may need tuning
  
  if(std::abs(T(1) - check_val) >= max_diff)  { return false; }
  
  arrayops::copy(Xm, Ym, uword(4*4));
  
  return true;
  }



template<typename T1>
inline
bool
op_inv_gen_rcond::apply_direct(Mat<typename T1::elem_type>& out, op_inv_gen_state<typename T1::pod_type>& out_state, const Base<typename T1::elem_type,T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  out             = expr.get_ref();
  out_state.rcond = T(0);
  
  arma_debug_check( (out.is_square() == false), "inv(): given matrix must be square sized", [&](){ out.soft_reset(); } );
  
  if(is_op_diagmat<T1>::value || out.is_diagmat())
    {
    arma_extra_debug_print("op_inv_gen_rcond: detected diagonal matrix");
    
    out_state.is_diag = true;
    
    eT* colmem = out.memptr();
    
    T max_abs_src_val = T(0);
    T max_abs_inv_val = T(0);
    
    const uword N = out.n_rows;
    
    for(uword i=0; i<N; ++i)
      {
      eT& out_ii = colmem[i];
      
      const eT src_val = out_ii;
      const eT inv_val = eT(1) / src_val;
      
      if(src_val == eT(0))  { return false; }
      
      out_ii = inv_val;
      
      const T abs_src_val = std::abs(src_val);
      const T abs_inv_val = std::abs(inv_val);
      
      max_abs_src_val = (abs_src_val > max_abs_src_val) ? abs_src_val : max_abs_src_val;
      max_abs_inv_val = (abs_inv_val > max_abs_inv_val) ? abs_inv_val : max_abs_inv_val;
      
      colmem += N;
      }
    
    out_state.rcond = T(1) / (max_abs_src_val * max_abs_inv_val);
    
    return true;
    }
  
  const strip_trimat<T1> strip(expr.get_ref());
  
  const bool is_triu_expr = strip.do_triu;
  const bool is_tril_expr = strip.do_tril;
  
  const bool is_triu_mat = (is_triu_expr || is_tril_expr) ? false : (                        trimat_helper::is_triu(out));
  const bool is_tril_mat = (is_triu_expr || is_tril_expr) ? false : ((is_triu_mat) ? false : trimat_helper::is_tril(out));
  
  if(is_triu_expr || is_tril_expr || is_triu_mat || is_tril_mat)
    {
    return auxlib::inv_tr_rcond(out, out_state.rcond, ((is_triu_expr || is_triu_mat) ? uword(0) : uword(1)));
    }
  
  const bool try_sympd = arma_config::optimise_sympd && ((auxlib::crippled_lapack(out)) ? false : sympd_helper::guess_sympd(out));
  
  if(try_sympd)
    {
    arma_extra_debug_print("op_inv_gen_rcond: attempting sympd optimisation");
    
    out_state.is_sym = true;
    
    Mat<eT> tmp = out;
    
    bool sympd_state = false;
    
    const bool status = auxlib::inv_sympd_rcond(tmp, sympd_state, out_state.rcond, T(-1));
    
    if(status)  { out.steal_mem(tmp); return true; }
    
    if((status == false) && (sympd_state == true))  { return false; }
    
    arma_extra_debug_print("op_inv_gen_rcond: sympd optimisation failed");
    
    // fallthrough if optimisation failed
    }
  
  return auxlib::inv_rcond(out, out_state.rcond);
  }



//! @}
