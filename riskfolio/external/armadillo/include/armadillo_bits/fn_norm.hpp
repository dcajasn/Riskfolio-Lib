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


//! \addtogroup fn_norm
//! @{



template<typename T1>
inline
arma_warn_unused
typename enable_if2< is_arma_type<T1>::value, typename T1::pod_type >::result
norm
  (
  const T1&   X,
  const uword k = uword(2),
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type T;
  
  const Proxy<T1> P(X);
  
  if(P.get_n_elem() == 0)  { return T(0); }
  
  const bool is_vec = (T1::is_xvec) || (T1::is_row) || (T1::is_col) || (P.get_n_rows() == 1) || (P.get_n_cols() == 1);
  
  if(is_vec)
    {
    if(k == uword(1))  { return op_norm::vec_norm_1(P); }
    if(k == uword(2))  { return op_norm::vec_norm_2(P); }
    
    arma_debug_check( (k == 0), "norm(): k must be greater than zero" );
    
    return op_norm::vec_norm_k(P, int(k));
    }
  else
    {
    const quasi_unwrap<typename Proxy<T1>::stored_type> U(P.Q);
    
    if(k == uword(1))  { return op_norm::mat_norm_1(U.M); }
    if(k == uword(2))  { return op_norm::mat_norm_2(U.M); }
      
    arma_stop_logic_error("norm(): unsupported matrix norm type");
    }
  
  return T(0);
  }



template<typename T1>
inline
arma_warn_unused
typename enable_if2< is_arma_type<T1>::value, typename T1::pod_type >::result
norm
  (
  const T1&   X,
  const char* method,
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type T;
  
  const Proxy<T1> P(X);
  
  if(P.get_n_elem() == 0)  { return T(0); }
  
  const char sig    = (method != nullptr) ? method[0] : char(0);
  const bool is_vec = (T1::is_xvec) || (T1::is_row) || (T1::is_col) || (P.get_n_rows() == 1) || (P.get_n_cols() == 1);
  
  if(is_vec)
    {
    if( (sig == 'i') || (sig == 'I') || (sig == '+') )  { return op_norm::vec_norm_max(P); }
    if( (sig == '-')                                 )  { return op_norm::vec_norm_min(P); }
    if( (sig == 'f') || (sig == 'F')                 )  { return op_norm::vec_norm_2(P);   }
    
    arma_stop_logic_error("norm(): unsupported vector norm type");
    }
  else
    {
    if( (sig == 'i') || (sig == 'I') || (sig == '+') )   // inf norm
      {
      const quasi_unwrap<typename Proxy<T1>::stored_type> U(P.Q);
      
      return op_norm::mat_norm_inf(U.M);
      }
    else
    if( (sig == 'f') || (sig == 'F') )
      {
      return op_norm::vec_norm_2(P);
      }
    
    arma_stop_logic_error("norm(): unsupported matrix norm type");
    }
  
  return T(0);
  }



template<typename T1>
inline
arma_warn_unused
typename enable_if2< is_arma_type<T1>::value, double >::result
norm
  (
  const T1&   X,
  const uword k = uword(2),
  const typename arma_integral_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  if(resolves_to_colvector<T1>::value)  { return norm(conv_to< Col<double> >::from(X), k); }
  if(resolves_to_rowvector<T1>::value)  { return norm(conv_to< Row<double> >::from(X), k); }
  
  return norm(conv_to< Mat<double> >::from(X), k);
  }



template<typename T1>
inline
arma_warn_unused
typename enable_if2< is_arma_type<T1>::value, double >::result
norm
  (
  const T1&   X,
  const char* method,
  const typename arma_integral_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  if(resolves_to_colvector<T1>::value)  { return norm(conv_to< Col<double> >::from(X), method); }
  if(resolves_to_rowvector<T1>::value)  { return norm(conv_to< Row<double> >::from(X), method); }
  
  return norm(conv_to< Mat<double> >::from(X), method);
  }



//
// norms for sparse matrices


template<typename T1>
inline
arma_warn_unused
typename enable_if2< is_arma_sparse_type<T1>::value, typename T1::pod_type >::result
norm
  (
  const T1&   expr,
  const uword k = uword(2),
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  if(is_SpSubview_col<T1>::value)
    {
    const SpSubview_col<eT>& sv = reinterpret_cast< const SpSubview_col<eT>& >(expr);
    
    if(sv.n_rows == sv.m.n_rows)
      {
      const SpMat<eT>& m   = sv.m;
      const uword      col = sv.aux_col1;
      const eT*        mem = &(m.values[ m.col_ptrs[col] ]);
      
      return spop_norm::vec_norm_k(mem, sv.n_nonzero, k);
      }
    }
  
  const unwrap_spmat<T1> U(expr);
  const SpMat<eT>& X   = U.M;
  
  if(X.n_nonzero == 0)  { return T(0); }
  
  const bool is_vec = (T1::is_xvec) || (T1::is_row) || (T1::is_col) || (X.n_rows == 1) || (X.n_cols == 1);
  
  if(is_vec)
    {
    return spop_norm::vec_norm_k(X.values, X.n_nonzero, k);
    }
  else
    {
    if(k == uword(1))  { return spop_norm::mat_norm_1(X); }
    if(k == uword(2))  { return spop_norm::mat_norm_2(X); }
    
    arma_stop_logic_error("norm(): unsupported or unimplemented norm type for sparse matrices");
    }
  
  return T(0);
  }



template<typename T1>
inline
arma_warn_unused
typename enable_if2< is_arma_sparse_type<T1>::value, typename T1::pod_type >::result
norm
  (
  const T1&   expr,
  const char* method,
  const typename arma_real_or_cx_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const unwrap_spmat<T1> U(expr);
  const SpMat<eT>& X   = U.M;
  
  if(X.n_nonzero == 0)  { return T(0); }
  
  // create a fake dense vector to allow reuse of code for dense vectors
  Col<eT> fake_vector( access::rwp(X.values), X.n_nonzero, false );
  
  const Proxy< Col<eT> > P_fake_vector(fake_vector);
  
  
  const char sig    = (method != nullptr) ? method[0] : char(0);
  const bool is_vec = (T1::is_xvec) || (T1::is_row) || (T1::is_col) || (X.n_rows == 1) || (X.n_cols == 1);
  
  if(is_vec)
    {
    if( (sig == 'i') || (sig == 'I') || (sig == '+') )   // max norm
      {
      return op_norm::vec_norm_max(P_fake_vector);
      }
    else
    if(sig == '-')   // min norm
      {
      const T val = op_norm::vec_norm_min(P_fake_vector);
      
      return (X.n_nonzero < X.n_elem) ? T((std::min)(T(0), val)) : T(val);
      }
    else
    if( (sig == 'f') || (sig == 'F') )
      {
      return op_norm::vec_norm_2(P_fake_vector);
      }
    
    arma_stop_logic_error("norm(): unsupported vector norm type");
    }
  else
    {
    if( (sig == 'i') || (sig == 'I') || (sig == '+') )   // inf norm
      {
      return spop_norm::mat_norm_inf(X);
      }
    else
    if( (sig == 'f') || (sig == 'F') )
      {
      return op_norm::vec_norm_2(P_fake_vector);
      }
    
    arma_stop_logic_error("norm(): unsupported matrix norm type");
    }
  
  return T(0);
  }



//! @}
