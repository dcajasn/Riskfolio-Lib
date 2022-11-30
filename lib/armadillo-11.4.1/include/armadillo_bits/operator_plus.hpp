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


//! \addtogroup operator_plus
//! @{



//! unary plus operation (does nothing, but is required for completeness)
template<typename T1>
arma_inline
typename enable_if2< is_arma_type<T1>::value, const T1& >::result
operator+
(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return X;
  }



//! Base + scalar
template<typename T1>
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_scalar_plus> >::result
operator+
(const T1& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_plus>(X, k);
  }



//! scalar + Base
template<typename T1>
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_scalar_plus> >::result
operator+
(const typename T1::elem_type k, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_plus>(X, k);  // NOTE: order is swapped
  }



//! non-complex Base + complex scalar
template<typename T1>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no),
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_plus>
  >::result
operator+
  (
  const T1&                                  X,
  const std::complex<typename T1::pod_type>& k
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_plus>('j', X, k);
  }



//! complex scalar + non-complex Base
template<typename T1>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no),
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_plus>
  >::result
operator+
  (
  const std::complex<typename T1::pod_type>& k,
  const T1&                                  X
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_plus>('j', X, k);  // NOTE: order is swapped
  }



//! addition of user-accessible Armadillo objects with same element type
template<typename T1, typename T2>
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value,
  const eGlue<T1, T2, eglue_plus>
  >::result
operator+
  (
  const T1& X,
  const T2& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return eGlue<T1, T2, eglue_plus>(X, Y);
  }



//! addition of user-accessible Armadillo objects with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && (is_same_type<typename T1::elem_type, typename T2::elem_type>::no)),
  const mtGlue<typename promote_type<typename T1::elem_type, typename T2::elem_type>::result, T1, T2, glue_mixed_plus>
  >::result
operator+
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
  
  return mtGlue<out_eT, T1, T2, glue_mixed_plus>( X, Y );
  }



//! addition of two sparse objects
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  SpGlue<T1,T2,spglue_plus>
  >::result
operator+
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  return SpGlue<T1,T2,spglue_plus>(x, y);
  }



//! addition of one dense and one sparse object
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  Mat<typename T1::elem_type>
  >::result
operator+
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  Mat<typename T1::elem_type> result(x);
  
  const SpProxy<T2> pb(y);
  
  arma_debug_assert_same_size( result.n_rows, result.n_cols, pb.get_n_rows(), pb.get_n_cols(), "addition" );
  
  typename SpProxy<T2>::const_iterator_type it     = pb.begin();
  typename SpProxy<T2>::const_iterator_type it_end = pb.end();
  
  while(it != it_end)
    {
    result.at(it.row(), it.col()) += (*it);
    ++it;
    }
  
  return result;
  }



//! addition of one sparse and one dense object
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  Mat<typename T1::elem_type>
  >::result
operator+
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<T1> pa(x);
  
  Mat<typename T1::elem_type> result(y);
  
  arma_debug_assert_same_size( pa.get_n_rows(), pa.get_n_cols(), result.n_rows, result.n_cols, "addition" );
  
  typename SpProxy<T1>::const_iterator_type it     = pa.begin();
  typename SpProxy<T1>::const_iterator_type it_end = pa.end();
  
  while(it != it_end)
    {
    result.at(it.row(), it.col()) += (*it);
    ++it;
    }
  
  return result;
  }



//! addition of two sparse objects with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::no),
  const mtSpGlue< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result, T1, T2, spglue_plus_mixed >
  >::result
operator+
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
  
  return mtSpGlue<out_eT, T1, T2, spglue_plus_mixed>( X, Y );
  }



//! addition of sparse and non-sparse objects with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::no),
  Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >
  >::result
operator+
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result > out;
  
  spglue_plus_mixed::dense_plus_sparse(out, x, y);
  
  return out;
  }



//! addition of sparse and non-sparse objects with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::no),
  Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >
  >::result
operator+
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result > out;
  
  // Just call the other order (these operations are commutative)
  // TODO: if there is a matrix size mismatch, the debug assert will print the matrix sizes in wrong order
  spglue_plus_mixed::dense_plus_sparse(out, y, x);
  
  return out;
  }



//! addition of sparse object with scalar
template<typename T1>
inline
typename enable_if2< is_arma_sparse_type<T1>::value, const SpToDOp<T1, op_sp_plus> >::result
operator+
  (
  const T1&                    X,
  const typename T1::elem_type k
  )
  {
  arma_extra_debug_sigprint();
  
  return SpToDOp<T1, op_sp_plus>(X, k);
  }



template<typename T1>
inline
typename enable_if2< is_arma_sparse_type<T1>::value, const SpToDOp<T1, op_sp_plus> >::result
operator+
  (
  const typename T1::elem_type k,
  const T1&                    X
  )
  {
  arma_extra_debug_sigprint();

  return SpToDOp<T1, op_sp_plus>(X, k);  // NOTE: swapped order
  }



// TODO: this is an uncommon use case; remove?
//! multiple applications of add/subtract scalars can be condensed
template<typename T1, typename op_type>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value &&
      (is_same_type<op_type, op_sp_plus>::value ||
       is_same_type<op_type, op_sp_minus_post>::value)),
  const SpToDOp<T1, op_sp_plus>
  >::result
operator+
  (
  const SpToDOp<T1, op_type>&  x,
  const typename T1::elem_type k
  )
  {
  arma_extra_debug_sigprint();

  const typename T1::elem_type aux = (is_same_type<op_type, op_sp_plus>::value) ? x.aux : -x.aux;

  return SpToDOp<T1, op_sp_plus>(x.m, aux + k);
  }



// TODO: this is an uncommon use case; remove?
//! multiple applications of add/subtract scalars can be condensed
template<typename T1, typename op_type>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value &&
       is_same_type<op_type, op_sp_minus_pre>::value),
  const SpToDOp<T1, op_sp_minus_pre>
  >::result
operator+
  (
  const SpToDOp<T1, op_type>&  x,
  const typename T1::elem_type k
  )
  {
  arma_extra_debug_sigprint();

  return SpToDOp<T1, op_sp_minus_pre>(x.m, x.aux + k);
  }



// TODO: this is an uncommon use case; remove?
//! multiple applications of add/subtract scalars can be condensed
template<typename T1, typename op_type>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value &&
      (is_same_type<op_type, op_sp_plus>::value ||
       is_same_type<op_type, op_sp_minus_post>::value)),
  const SpToDOp<T1, op_sp_plus>
  >::result
operator+
  (
  const typename T1::elem_type k,
  const SpToDOp<T1, op_type>&  x
  )
  {
  arma_extra_debug_sigprint();

  const typename T1::elem_type aux = (is_same_type<op_type, op_sp_plus>::value) ? x.aux : -x.aux;

  return SpToDOp<T1, op_sp_plus>(x.m, aux + k);
  }



// TODO: this is an uncommon use case; remove?
//! multiple applications of add/subtract scalars can be condensed
template<typename T1, typename op_type>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value &&
       is_same_type<op_type, op_sp_minus_pre>::value),
  const SpToDOp<T1, op_sp_minus_pre>
  >::result
operator+
  (
  const typename T1::elem_type k,
  const SpToDOp<T1, op_type>&  x
  )
  {
  arma_extra_debug_sigprint();

  return SpToDOp<T1, op_sp_minus_pre>(x.m, x.aux + k);
  }




template<typename parent, unsigned int mode, typename T2>
arma_inline
Mat<typename parent::elem_type>
operator+
  (
  const subview_each1<parent,mode>&          X,
  const Base<typename parent::elem_type,T2>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return subview_each1_aux::operator_plus(X, Y.get_ref());
  }



template<typename T1, typename parent, unsigned int mode>
arma_inline
Mat<typename parent::elem_type>
operator+
  (
  const Base<typename parent::elem_type,T1>& X,
  const subview_each1<parent,mode>&          Y
  )
  {
  arma_extra_debug_sigprint();
  
  return subview_each1_aux::operator_plus(Y, X.get_ref());  // NOTE: swapped order
  }



template<typename parent, unsigned int mode, typename TB, typename T2>
arma_inline
Mat<typename parent::elem_type>
operator+
  (
  const subview_each2<parent,mode,TB>&       X,
  const Base<typename parent::elem_type,T2>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return subview_each2_aux::operator_plus(X, Y.get_ref());
  }



template<typename T1, typename parent, unsigned int mode, typename TB>
arma_inline
Mat<typename parent::elem_type>
operator+
  (
  const Base<typename parent::elem_type,T1>& X,
  const subview_each2<parent,mode,TB>&       Y
  )
  {
  arma_extra_debug_sigprint();
  
  return subview_each2_aux::operator_plus(Y, X.get_ref());  // NOTE: swapped order
  }



//! @}
