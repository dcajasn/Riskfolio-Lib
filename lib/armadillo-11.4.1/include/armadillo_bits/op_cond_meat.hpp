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


//! \addtogroup op_cond
//! @{



template<typename T1>
inline
typename T1::pod_type
op_cond::apply(const Base<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  Mat<eT> A(X.get_ref());
  
  if(A.n_elem == 0)  { return T(0); }
  
  if(is_op_diagmat<T1>::value || A.is_diagmat())
    {
    arma_extra_debug_print("op_cond::apply(): detected diagonal matrix");
    
    return op_cond::apply_diag(A);
    }
  
  bool is_approx_sym   = false;
  bool is_approx_sympd = false;
  
  sympd_helper::analyse_matrix(is_approx_sym, is_approx_sympd, A);
  
  const bool do_sym = (is_cx<eT>::no) ? (is_approx_sym) : (is_approx_sym && is_approx_sympd);
  
  if(do_sym)
    {
    arma_extra_debug_print("op_cond: symmetric/hermitian optimisation");
    
    return op_cond::apply_sym(A);
    }
  
  return op_cond::apply_gen(A);
  }



template<typename eT>
inline
typename get_pod_type<eT>::result
op_cond::apply_diag(const Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const uword N = (std::min)(A.n_rows, A.n_cols);
  
  T abs_min = Datum<T>::inf;
  T abs_max = T(0);
  
  for(uword i=0; i < N; ++i)
    {
    const T abs_val = std::abs(A.at(i,i));
    
    if(arma_isnan(abs_val))
      {
      arma_debug_warn_level(3, "cond(): failed");
      
      return Datum<T>::nan;
      }
    
    abs_min = (abs_val < abs_min) ? abs_val : abs_min;
    abs_max = (abs_val > abs_max) ? abs_val : abs_max;
    }
  
  if((abs_min == T(0)) || (abs_max == T(0)))  { return Datum<T>::inf; }
  
  return T(abs_max / abs_min);
  }



template<typename eT>
inline
typename get_pod_type<eT>::result
op_cond::apply_sym(Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  Col<T> eigval;
  
  const bool status = auxlib::eig_sym(eigval, A);
  
  if(status == false)
    {
    arma_debug_warn_level(3, "cond(): failed");
    
    return Datum<T>::nan;
    }
  
  if(eigval.n_elem == 0)  { return T(0); }
  
  const T* eigval_mem = eigval.memptr();
  
  T abs_min = std::abs(eigval_mem[0]);
  T abs_max = abs_min;
  
  for(uword i=1; i < eigval.n_elem; ++i)
    {
    const T abs_val = std::abs(eigval_mem[i]);
    
    abs_min = (abs_val < abs_min) ? abs_val : abs_min;
    abs_max = (abs_val > abs_max) ? abs_val : abs_max;
    }
  
  if((abs_min == T(0)) || (abs_max == T(0)))  { return Datum<T>::inf; }
  
  return T(abs_max / abs_min);
  }



template<typename eT>
inline
typename get_pod_type<eT>::result
op_cond::apply_gen(Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  Col<T> S;
  
  const bool status = auxlib::svd_dc(S, A);
  
  if(status == false)
    {
    arma_debug_warn_level(3, "cond(): failed");
    
    return Datum<T>::nan;
    }
  
  if(S.n_elem == 0)  { return T(0); }
  
  const T S_max = S[0];
  const T S_min = S[S.n_elem-1];
  
  if((S_max == T(0)) || (S_min == T(0)))  { return Datum<T>::inf; }
  
  return T(S_max / S_min);
  }



//! @}
