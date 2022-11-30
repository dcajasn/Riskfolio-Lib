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


//! \addtogroup fn_min
//! @{


template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_arma_type<T1>::value && resolves_to_vector<T1>::yes, typename T1::elem_type >::result
min(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return op_min::min(X);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value && resolves_to_vector<T1>::no, const Op<T1, op_min> >::result
min(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_min>(X, 0, 0);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_arma_type<T1>::value, const Op<T1, op_min> >::result
min(const T1& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_min>(X, dim, 0);
  }



template<typename T>
arma_warn_unused
arma_inline
typename arma_scalar_only<T>::result
min(const T& x)
  {
  return x;
  }



//! element-wise minimum
template<typename T1, typename T2>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  ( is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value ),
  const Glue<T1, T2, glue_min>
  >::result
min
  (
  const T1& X,
  const T2& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_min>(X, Y);
  }



template<typename T1>
arma_warn_unused
arma_inline
const OpCube<T1, op_min>
min
  (
  const BaseCube<typename T1::elem_type, T1>& X,
  const uword dim = 0
  )
  {
  arma_extra_debug_sigprint();
  
  return OpCube<T1, op_min>(X.get_ref(), dim, 0);
  }



template<typename T1, typename T2>
arma_warn_unused
arma_inline
const GlueCube<T1, T2, glue_min>
min
  (
  const BaseCube<typename T1::elem_type, T1>& X,
  const BaseCube<typename T1::elem_type, T2>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_min>(X.get_ref(), Y.get_ref());
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value && resolves_to_sparse_vector<T1>::yes,
  typename T1::elem_type
  >::result
min(const T1& x)
  {
  arma_extra_debug_sigprint();
  
  return spop_min::vector_min(x);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value && resolves_to_sparse_vector<T1>::no,
  const SpOp<T1, spop_min>
  >::result
min(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1, spop_min>(X, 0, 0);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value,
  const SpOp<T1, spop_min>
  >::result
min(const T1& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1, spop_min>(X, dim, 0);
  }



// elementwise sparse min
template<typename T1, typename T2>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const SpGlue<T1, T2, spglue_min>
  >::result
min(const T1& x, const T2& y)
  {
  arma_extra_debug_sigprint();

  return SpGlue<T1, T2, spglue_min>(x, y);
  }



//! elementwise min of dense and sparse objects with the same element type
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  Mat<typename T1::elem_type>
  >::result
min
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  Mat<typename T1::elem_type> out;
  
  spglue_min::dense_sparse_min(out, x, y);
  
  return out;
  }



//! elementwise min of sparse and dense objects with the same element type
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  Mat<typename T1::elem_type>
  >::result
min
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  Mat<typename T1::elem_type> out;
  
  // Just call the other order (these operations are commutative)
  // TODO: if there is a matrix size mismatch, the debug assert will print the matrix sizes in wrong order
  spglue_min::dense_sparse_min(out, y, x);
  
  return out;
  }



arma_warn_unused
inline
uword
min(const SizeMat& s)
  {
  return (std::min)(s.n_rows, s.n_cols);
  }



arma_warn_unused
inline
uword
min(const SizeCube& s)
  {
  return (std::min)( (std::min)(s.n_rows, s.n_cols), s.n_slices );
  }



//! @}
