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


//! \addtogroup op_norm
//! @{



template<typename eT>
inline
typename get_pod_type<eT>::result
spop_norm::mat_norm_1(const SpMat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  // TODO: this can be sped up with a dedicated implementation
  return as_scalar( max( sum(abs(X), 0), 1) );
  }



template<typename eT>
inline
typename get_pod_type<eT>::result
spop_norm::mat_norm_2(const SpMat<eT>& X, const typename arma_real_only<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  // norm = sqrt( largest eigenvalue of (A^H)*A ), where ^H is the conjugate transpose
  // http://math.stackexchange.com/questions/4368/computing-the-largest-eigenvalue-of-a-very-large-sparse-matrix
  
  typedef typename get_pod_type<eT>::result T;
  
  const SpMat<eT>& A = X;
  const SpMat<eT>  B = trans(A);
  
  const SpMat<eT>  C = (A.n_rows <= A.n_cols) ? (A*B) : (B*A);
  
  Col<T> eigval;
  eigs_sym(eigval, C, 1);
  
  return (eigval.n_elem > 0) ? T(std::sqrt(eigval[0])) : T(0);
  }



template<typename eT>
inline
typename get_pod_type<eT>::result
spop_norm::mat_norm_2(const SpMat<eT>& X, const typename arma_cx_only<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename get_pod_type<eT>::result T;
  
  // we're calling eigs_gen(), which currently requires ARPACK
  #if !defined(ARMA_USE_ARPACK)
    {
    arma_stop_logic_error("norm(): use of ARPACK must be enabled for norm of complex matrices");
    return T(0);
    }
  #endif
  
  const SpMat<eT>& A = X;
  const SpMat<eT>  B = trans(A);
  
  const SpMat<eT>  C = (A.n_rows <= A.n_cols) ? (A*B) : (B*A);
  
  Col<eT> eigval;
  eigs_gen(eigval, C, 1);
  
  return (eigval.n_elem > 0) ? T(std::sqrt(std::real(eigval[0]))) : T(0);
  }



template<typename eT>
inline
typename get_pod_type<eT>::result
spop_norm::mat_norm_inf(const SpMat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  // TODO: this can be sped up with a dedicated implementation
  return as_scalar( max( sum(abs(X), 1), 0) );
  }



template<typename eT>
inline
typename get_pod_type<eT>::result
spop_norm::vec_norm_k(const eT* mem, const uword N, const uword k)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (k == 0), "norm(): k must be greater than zero" );
  
  // create a fake dense vector to allow reuse of code for dense vectors
  Col<eT> fake_vector( access::rwp(mem), N, false );
  
  const Proxy< Col<eT> > P_fake_vector(fake_vector);
  
  if(k == uword(1))  { return op_norm::vec_norm_1(P_fake_vector); }
  if(k == uword(2))  { return op_norm::vec_norm_2(P_fake_vector); }
  
  return op_norm::vec_norm_k(P_fake_vector, int(k));
  }



//! @}
