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


//! \addtogroup operator_times
//! @{



//! Base * scalar
template<typename T1>
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_scalar_times> >::result
operator*
(const T1& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_times>(X,k);
  }



//! scalar * Base
template<typename T1>
arma_inline
typename enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_scalar_times> >::result
operator*
(const typename T1::elem_type k, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_times>(X,k);  // NOTE: order is swapped
  }



//! non-complex Base * complex scalar
template<typename T1>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no),
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>
  >::result
operator*
  (
  const T1&                                  X,
  const std::complex<typename T1::pod_type>& k
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>('j', X, k);
  }



//! complex scalar * non-complex Base
template<typename T1>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no),
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>
  >::result
operator*
  (
  const std::complex<typename T1::pod_type>& k,
  const T1&                                  X
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>('j', X, k);
  }



//! scalar * trans(T1)
template<typename T1>
arma_inline
const Op<T1, op_htrans2>
operator*
(const typename T1::elem_type k, const Op<T1, op_htrans>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_htrans2>(X.m, k);
  }



//! trans(T1) * scalar
template<typename T1>
arma_inline
const Op<T1, op_htrans2>
operator*
(const Op<T1, op_htrans>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_htrans2>(X.m, k);
  }



//! Base * diagmat
template<typename T1, typename T2>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const Glue<T1, Op<T2, op_diagmat>, glue_times_diag>
  >::result
operator*
(const T1& X, const Op<T2, op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, Op<T2, op_diagmat>, glue_times_diag>(X, Y);
  }



//! diagmat * Base
template<typename T1, typename T2>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const Glue<Op<T1, op_diagmat>, T2, glue_times_diag>
  >::result
operator*
(const Op<T1, op_diagmat>& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Op<T1, op_diagmat>, T2, glue_times_diag>(X, Y);
  }



//! diagmat * diagmat
template<typename T1, typename T2>
inline
Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >
operator*
(const Op<T1, op_diagmat>& X, const Op<T2, op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const diagmat_proxy<T1> A(X.m);
  const diagmat_proxy<T2> B(Y.m);
  
  arma_debug_assert_mul_size(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  Mat<out_eT> out(A.n_rows, B.n_cols, arma_zeros_indicator());
  
  const uword A_length = (std::min)(A.n_rows, A.n_cols);
  const uword B_length = (std::min)(B.n_rows, B.n_cols);
  
  const uword N = (std::min)(A_length, B_length);
  
  for(uword i=0; i<N; ++i)
    {
    out.at(i,i) = upgrade_val<eT1,eT2>::apply( A[i] ) * upgrade_val<eT1,eT2>::apply( B[i] );
    }
  
  return out;
  }



//! multiplication of Base objects with same element type
template<typename T1, typename T2>
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value,
  const Glue<T1, T2, glue_times>
  >::result
operator*
(const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X, Y);
  }



//! multiplication of Base objects with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && (is_same_type<typename T1::elem_type, typename T2::elem_type>::no)),
  const mtGlue< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result, T1, T2, glue_mixed_times >
  >::result
operator*
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
  
  return mtGlue<out_eT, T1, T2, glue_mixed_times>( X, Y );
  }



//! sparse multiplied by scalar
template<typename T1>
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value,
  SpOp<T1,spop_scalar_times>
  >::result
operator*
  (
  const T1& X,
  const typename T1::elem_type k
  )
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1,spop_scalar_times>(X, k);
  }



template<typename T1>
inline
typename
enable_if2
  <
  is_arma_sparse_type<T1>::value,
  SpOp<T1,spop_scalar_times>
  >::result
operator*
  (
  const typename T1::elem_type k,
  const T1& X
  )
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1,spop_scalar_times>(X, k);
  }



//! non-complex sparse * complex scalar
template<typename T1>
arma_inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_cx<typename T1::elem_type>::no),
  const mtSpOp<typename std::complex<typename T1::pod_type>, T1, spop_cx_scalar_times>
  >::result
operator*
  (
  const T1&                                  X,
  const std::complex<typename T1::pod_type>& k
  )
  {
  arma_extra_debug_sigprint();
  
  return mtSpOp<typename std::complex<typename T1::pod_type>, T1, spop_cx_scalar_times>('j', X, k);
  }



//! complex scalar * non-complex sparse
template<typename T1>
arma_inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_cx<typename T1::elem_type>::no),
  const mtSpOp<typename std::complex<typename T1::pod_type>, T1, spop_cx_scalar_times>
  >::result
operator*
  (
  const std::complex<typename T1::pod_type>& k,
  const T1&                                  X
  )
  {
  arma_extra_debug_sigprint();
  
  return mtSpOp<typename std::complex<typename T1::pod_type>, T1, spop_cx_scalar_times>('j', X, k);
  }



//! multiplication of two sparse objects
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const SpGlue<T1,T2,spglue_times>
  >::result
operator*
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();

  return SpGlue<T1,T2,spglue_times>(x, y);
  }



//! multiplication of one sparse and one dense object
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  Mat<typename T1::elem_type>
  >::result
operator*
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> result;
  
  spglue_times_misc::sparse_times_dense(result, x, y);
  
  return result;
  }



//! multiplication of one dense and one sparse object
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  Mat<typename T1::elem_type>
  >::result
operator*
  (
  const T1& x,
  const T2& y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> result;
  
  spglue_times_misc::dense_times_sparse(result, x, y);
  
  return result;
  }



//! multiplication of two sparse objects with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_sparse_type<T2>::value && (is_same_type<typename T1::elem_type, typename T2::elem_type>::no)),
  const mtSpGlue< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result, T1, T2, spglue_times_mixed >
  >::result
operator*
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
  
  return mtSpGlue<out_eT, T1, T2, spglue_times_mixed>( X, Y );
  }



//! multiplication of one sparse and one dense object with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::no),
  Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >
  >::result
operator*
  (
  const T1& X,
  const T2& Y
  )
  {
  arma_extra_debug_sigprint();
  
  Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result > out;
  
  spglue_times_mixed::sparse_times_dense(out, X, Y);
  
  return out;
  }



//! multiplication of one dense and one sparse object with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::no),
  Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >
  >::result
operator*
  (
  const T1& X,
  const T2& Y
  )
  {
  arma_extra_debug_sigprint();
  
  Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result > out;
  
  spglue_times_mixed::dense_times_sparse(out, X, Y);
  
  return out;
  }



//! @}
