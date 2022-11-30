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


//! \addtogroup op_sp_minus
//! @{


// scalar - SpBase
template<typename T1>
inline
void
op_sp_minus_pre::apply(Mat<typename T1::elem_type>& out, const SpToDOp<T1,op_sp_minus_pre>& in)
  {
  arma_extra_debug_sigprint();
  
  // Note that T1 will be a sparse type, so we use SpProxy.
  const SpProxy<T1> proxy(in.m);
  
  out.set_size(proxy.get_n_rows(), proxy.get_n_cols());
  out.fill(in.aux);
  
  typename SpProxy<T1>::const_iterator_type it     = proxy.begin();
  typename SpProxy<T1>::const_iterator_type it_end = proxy.end();
  
  for(; it != it_end; ++it)
    {
    out.at(it.row(), it.col()) -= (*it);
    }
  }



// force apply into SpMat
template<typename T1>
inline
void
op_sp_minus_pre::apply(SpMat<typename T1::elem_type>& out, const SpToDOp<T1,op_sp_minus_pre>& in)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;

  // Note that T1 will be a sparse type, so we use SpProxy.
  const SpProxy<T1> proxy(in.m);
  
  const uword n_rows = proxy.get_n_rows();
  const uword n_cols = proxy.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
  const eT k = in.aux;
  
  for(uword c = 0; c < n_cols; ++c)
  for(uword r = 0; r < n_rows; ++r)
    {
    out.at(r, c) = k - proxy.at(r, c);
    }
  }



// used for the optimization of sparse % (scalar - sparse)
template<typename eT, typename T2, typename T3>
inline
void
op_sp_minus_pre::apply_inside_schur(SpMat<eT>& out, const T2& x, const SpToDOp<T3, op_sp_minus_pre>& y)
  {
  arma_extra_debug_sigprint();

  const SpProxy<T2> proxy2(x);
  const SpProxy<T3> proxy3(y.m);

  arma_debug_assert_same_size(proxy2.get_n_rows(), proxy2.get_n_cols(), proxy3.get_n_rows(), proxy3.get_n_cols(), "element-wise multiplication");

  out.zeros(proxy2.get_n_rows(), proxy2.get_n_cols());
  
  typename SpProxy<T2>::const_iterator_type it     = proxy2.begin();
  typename SpProxy<T2>::const_iterator_type it_end = proxy2.end();
  
  const eT k = y.aux;
  
  for(; it != it_end; ++it)
    {
    const uword it_row = it.row();
    const uword it_col = it.col();
    
    out.at(it_row, it_col) = (*it) * (k - proxy3.at(it_row, it_col));
    }
  }



// used for the optimization of sparse / (scalar - sparse)
template<typename eT, typename T2, typename T3>
inline
void
op_sp_minus_pre::apply_inside_div(SpMat<eT>& out, const T2& x, const SpToDOp<T3, op_sp_minus_pre>& y)
  {
  arma_extra_debug_sigprint();

  const SpProxy<T2> proxy2(x);
  const SpProxy<T3> proxy3(y.m);

  arma_debug_assert_same_size(proxy2.get_n_rows(), proxy2.get_n_cols(), proxy3.get_n_rows(), proxy3.get_n_cols(), "element-wise multiplication");

  out.zeros(proxy2.get_n_rows(), proxy2.get_n_cols());
  
  typename SpProxy<T2>::const_iterator_type it     = proxy2.begin();
  typename SpProxy<T2>::const_iterator_type it_end = proxy2.end();
  
  const eT k = y.aux;
  
  for(; it != it_end; ++it)
    {
    const uword it_row = it.row();
    const uword it_col = it.col();
    
    out.at(it_row, it_col) = (*it) / (k - proxy3.at(it_row, it_col));
    }
  }



// SpBase - scalar
template<typename T1>
inline
void
op_sp_minus_post::apply(Mat<typename T1::elem_type>& out, const SpToDOp<T1,op_sp_minus_post>& in)
  {
  arma_extra_debug_sigprint();
  
  // Note that T1 will be a sparse type, so we use SpProxy.
  const SpProxy<T1> proxy(in.m);
  
  out.set_size(proxy.get_n_rows(), proxy.get_n_cols());
  out.fill(-in.aux);
  
  typename SpProxy<T1>::const_iterator_type it     = proxy.begin();
  typename SpProxy<T1>::const_iterator_type it_end = proxy.end();
  
  for(; it != it_end; ++it)
    {
    out.at(it.row(), it.col()) += (*it);
    }
  }



// force apply into sparse matrix
template<typename T1>
inline
void
op_sp_minus_post::apply(SpMat<typename T1::elem_type>& out, const SpToDOp<T1,op_sp_minus_post>& in)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;

  // Note that T1 will be a sparse type, so we use SpProxy.
  const SpProxy<T1> proxy(in.m);
  
  const uword n_rows = proxy.get_n_rows();
  const uword n_cols = proxy.get_n_cols();

  out.set_size(n_rows, n_cols);
  
  const eT k = in.aux;

  for(uword c = 0; c < n_cols; ++c)
  for(uword r = 0; r < n_rows; ++r)
    {
    out.at(r, c) = proxy.at(r, c) - k;
    }
  }



// used for the optimization of sparse % (sparse - scalar)
template<typename eT, typename T2, typename T3>
inline
void
op_sp_minus_post::apply_inside_schur(SpMat<eT>& out, const T2& x, const SpToDOp<T3, op_sp_minus_post>& y)
  {
  arma_extra_debug_sigprint();

  const SpProxy<T2> proxy2(x);
  const SpProxy<T3> proxy3(y.m);

  arma_debug_assert_same_size(proxy2.get_n_rows(), proxy2.get_n_cols(), proxy3.get_n_rows(), proxy3.get_n_cols(), "element-wise multiplication");

  out.zeros(proxy2.get_n_rows(), proxy2.get_n_cols());
  
  typename SpProxy<T2>::const_iterator_type it     = proxy2.begin();
  typename SpProxy<T2>::const_iterator_type it_end = proxy2.end();
  
  const eT k = y.aux;
  
  for(; it != it_end; ++it)
    {
    const uword it_row = it.row();
    const uword it_col = it.col();
    
    out.at(it_row, it_col) = (*it) * (proxy3.at(it_row, it_col) - k);
    }
  }



// used for the optimization of sparse / (sparse - scalar)
template<typename eT, typename T2, typename T3>
inline
void
op_sp_minus_post::apply_inside_div(SpMat<eT>& out, const T2& x, const SpToDOp<T3, op_sp_minus_post>& y)
  {
  arma_extra_debug_sigprint();

  const SpProxy<T2> proxy2(x);
  const SpProxy<T3> proxy3(y.m);

  arma_debug_assert_same_size(proxy2.get_n_rows(), proxy2.get_n_cols(), proxy3.get_n_rows(), proxy3.get_n_cols(), "element-wise multiplication");

  out.zeros(proxy2.get_n_rows(), proxy2.get_n_cols());
  
  typename SpProxy<T2>::const_iterator_type it     = proxy2.begin();
  typename SpProxy<T2>::const_iterator_type it_end = proxy2.end();
  
  const eT k = y.aux;
  
  for(; it != it_end; ++it)
    {
    const uword it_row = it.row();
    const uword it_col = it.col();
    
    out.at(it_row, it_col) = (*it) / (proxy3.at(it_row, it_col) - k);
    }
  }



//! @}
