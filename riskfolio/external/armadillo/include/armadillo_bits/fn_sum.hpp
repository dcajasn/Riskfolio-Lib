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


//! \addtogroup fn_sum
//! @{


template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_arma_type<T1>::value && resolves_to_vector<T1>::yes, typename T1::elem_type >::result
sum(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return accu(X);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value && resolves_to_vector<T1>::no, const Op<T1, op_sum> >::result
sum(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_sum>(X, 0, 0);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const Op<T1, op_sum> >::result
sum(const T1& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_sum>(X, dim, 0);
  }



template<typename T>
arma_warn_unused
arma_inline
typename arma_scalar_only<T>::result
sum(const T& x)
  {
  return x;
  }



//! sum of cube
template<typename T1>
arma_warn_unused
arma_inline
const OpCube<T1, op_sum>
sum
  (
  const BaseCube<typename T1::elem_type,T1>& X,
  const uword dim = 0
  )
  {
  arma_extra_debug_sigprint();
  
  return OpCube<T1, op_sum>(X.get_ref(), dim, 0);
  }



//! sum of sparse object
template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value && resolves_to_sparse_vector<T1>::yes,
  typename T1::elem_type
  >::result
sum(const T1& x)
  {
  arma_extra_debug_sigprint();
  
  // sum elements
  return accu(x);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value && resolves_to_sparse_vector<T1>::no,
  const SpOp<T1,spop_sum>
  >::result
sum(const T1& x)
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1,spop_sum>(x, 0, 0);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value,
  const SpOp<T1,spop_sum>
  >::result
sum(const T1& x, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1,spop_sum>(x, dim, 0);
  }



//! @}
