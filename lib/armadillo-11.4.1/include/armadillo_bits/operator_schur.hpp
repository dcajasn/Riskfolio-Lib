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


//! \addtogroup operator_schur
//! @{


// operator %, which we define it to do a schur product (element-wise multiplication)


//! element-wise multiplication of user-accessible Armadillo objects with same element type
template<typename T1, typename T2>
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value,
  const eGlue<T1, T2, eglue_schur>
  >::result
operator%
  (
  const T1& X,
  const T2& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return eGlue<T1, T2, eglue_schur>(X, Y);
  }



//! element-wise multiplication of user-accessible Armadillo objects with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && (is_same_type<typename T1::elem_type, typename T2::elem_type>::no)),
  const mtGlue<typename promote_type<typename T1::elem_type, typename T2::elem_type>::result, T1, T2, glue_mixed_schur>
  >::result
operator%
  (
  const T1& X,
  const T2& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  return mtGlue<out_eT, T1, T2, glue_mixed_schur>( X, Y );
  }



//! element-wise multiplication of two sparse matrices
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  SpGlue<T1,T2,spglue_schur>
  >::result
operator%
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  return SpGlue<T1,T2,spglue_schur>(x, y);
  }



//! element-wise multiplication of one dense and one sparse object
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  SpMat<typename T1::elem_type>
  >::result
operator%
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  SpMat<typename T1::elem_type> out;
  
  spglue_schur_misc::dense_schur_sparse(out, x, y);
  
  return out;
  }



//! element-wise multiplication of one sparse and one dense object
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  SpMat<typename T1::elem_type>
  >::result
operator%
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  SpMat<typename T1::elem_type> out;
  
  // Just call the other order (these operations are commutative)
  // TODO: if there is a matrix size mismatch, the debug assert will print the matrix sizes in wrong order
  spglue_schur_misc::dense_schur_sparse(out, y, x);
  
  return out;
  }



//! optimization: sparse % (sparse +/- scalar) can be done without forming the dense result of the (sparse +/- scalar) term
template<typename T1, typename T2, typename op_type>
inline
typename
enable_if2
  <
  (
  is_arma_sparse_type<T1>::value && is_arma_sparse_type<T2>::value &&
  is_same_type<typename T1::elem_type, typename T2::elem_type>::yes &&
      (is_same_type<op_type, op_sp_plus>::value ||
       is_same_type<op_type, op_sp_minus_pre>::value ||
       is_same_type<op_type, op_sp_minus_post>::value)
  ),
  SpMat<typename T1::elem_type>
  >::result
operator%
  (
  const T1& x,
  const SpToDOp<T2, op_type>& y
  )
  {
  arma_extra_debug_sigprint();
  
  SpMat<typename T1::elem_type> out;
  
  op_type::apply_inside_schur(out, x, y);
  
  return out;
  }



//! optimization: (sparse +/- scalar) % sparse can be done without forming the dense result of the (sparse +/- scalar) term
template<typename T1, typename T2, typename op_type>
inline
typename
enable_if2
  <
  (
  is_arma_sparse_type<T1>::value && is_arma_sparse_type<T2>::value &&
  is_same_type<typename T1::elem_type, typename T2::elem_type>::yes &&
      (is_same_type<op_type, op_sp_plus>::value ||
       is_same_type<op_type, op_sp_minus_pre>::value ||
       is_same_type<op_type, op_sp_minus_post>::value)
  ),
  SpMat<typename T1::elem_type>
  >::result
operator%
  (
  const SpToDOp<T1, op_type>& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  SpMat<typename T1::elem_type> out;
  
  // Just call the other order (these operations are commutative)
  // TODO: if there is a matrix size mismatch, the debug assert will print the matrix sizes in wrong order
  op_type::apply_inside_schur(out, y, x);
  
  return out;
  }



//! element-wise multiplication of two sparse objects with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::no),
  const mtSpGlue< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result, T1, T2, spglue_schur_mixed >
  >::result
operator%
  (
  const T1& X,
  const T2& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  return mtSpGlue<out_eT, T1, T2, spglue_schur_mixed>( X, Y );
  }



//! element-wise multiplication of one dense and one sparse object with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::no),
  SpMat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >
  >::result
operator%
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  SpMat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result > out;
  
  spglue_schur_mixed::dense_schur_sparse(out, x, y);
  
  return out;
  }



//! element-wise multiplication of one sparse and one dense object with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::no),
  SpMat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >
  >::result
operator%
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  SpMat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result > out;
  
  // Just call the other order (these operations are commutative)
  // TODO: if there is a matrix size mismatch, the debug assert will print the matrix sizes in wrong order
  spglue_schur_mixed::dense_schur_sparse(out, y, x);
  
  return out;
  }



template<typename parent, unsigned int mode, typename T2>
inline
Mat<typename parent::elem_type>
operator%
  (
  const subview_each1<parent,mode>&          X,
  const Base<typename parent::elem_type,T2>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return subview_each1_aux::operator_schur(X, Y.get_ref());
  }



template<typename T1, typename parent, unsigned int mode>
arma_inline
Mat<typename parent::elem_type>
operator%
  (
  const Base<typename parent::elem_type,T1>& X,
  const subview_each1<parent,mode>&          Y
  )
  {
  arma_extra_debug_sigprint();
  
  return subview_each1_aux::operator_schur(Y, X.get_ref());  // NOTE: swapped order
  }



template<typename parent, unsigned int mode, typename TB, typename T2>
inline
Mat<typename parent::elem_type>
operator%
  (
  const subview_each2<parent,mode,TB>&       X,
  const Base<typename parent::elem_type,T2>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return subview_each2_aux::operator_schur(X, Y.get_ref());
  }



template<typename T1, typename parent, unsigned int mode, typename TB>
arma_inline
Mat<typename parent::elem_type>
operator%
  (
  const Base<typename parent::elem_type,T1>& X,
  const subview_each2<parent,mode,TB>&       Y
  )
  {
  arma_extra_debug_sigprint();
  
  return subview_each2_aux::operator_schur(Y, X.get_ref());  // NOTE: swapped order
  }



//! @}
