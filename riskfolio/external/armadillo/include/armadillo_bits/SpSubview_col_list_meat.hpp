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


//! \addtogroup SpSubview_col_list
//! @{



template<typename eT, typename T1>
inline
SpSubview_col_list<eT,T1>::~SpSubview_col_list()
  {
  arma_extra_debug_sigprint();
  }



template<typename eT, typename T1>
arma_inline
SpSubview_col_list<eT,T1>::SpSubview_col_list
  (
  const SpMat<eT>&      in_m,
  const Base<uword,T1>& in_ci
  )
  : m   (in_m           )
  , U_ci(in_ci.get_ref())
  {
  arma_extra_debug_sigprint();
  
  const umat&  ci        = U_ci.M;
  const uword* ci_mem    = ci.memptr();
  const uword  ci_n_elem = ci.n_elem;
  
  arma_debug_check
    (
    ( (ci.is_vec() == false) && (ci.is_empty() == false) ),
    "SpMat::cols(): given object must be a vector"
    );
  
  for(uword ci_count=0; ci_count < ci_n_elem; ++ci_count)
    {
    const uword i = ci_mem[ci_count];
    
    arma_debug_check_bounds( (i >= in_m.n_cols), "SpMat::cols(): index out of bounds" );
    }
  }



//! apply a functor to each element
template<typename eT, typename T1>
template<typename functor>
inline
void
SpSubview_col_list<eT,T1>::for_each(functor F)
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT> tmp(*this);
  
  tmp.for_each(F);
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
template<typename functor>
inline
void
SpSubview_col_list<eT,T1>::for_each(functor F) const
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(*this);
  
  tmp.for_each(F);
  }



//! transform each element using a functor
template<typename eT, typename T1>
template<typename functor>
inline
void
SpSubview_col_list<eT,T1>::transform(functor F)
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT> tmp(*this);
  
  tmp.transform(F);
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::replace(const eT old_val, const eT new_val)
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT> tmp(*this);
  
  tmp.replace(old_val, new_val);
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::clean(const typename get_pod_type<eT>::result threshold)
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT> tmp(*this);
  
  tmp.clean(threshold);
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp(m.n_rows, U_ci.M.n_elem, arma_nozeros_indicator());  tmp.fill(val);
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::zeros()
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT>& m_local = const_cast< SpMat<eT>& >(m);
  
  const umat&  ci        = U_ci.M;
  const uword* ci_mem    = ci.memptr();
  const uword  ci_n_elem = ci.n_elem;
  
  m_local.sync_csc();
  m_local.invalidate_cache();
  
  for(uword ci_count=0; ci_count < ci_n_elem; ++ci_count)
    {
    const uword i = ci_mem[ci_count];
    
    const uword col_n_nonzero = m_local.col_ptrs[i+1] - m_local.col_ptrs[i];
    
    uword offset = m_local.col_ptrs[i];
    
    for(uword j=0; j < col_n_nonzero; ++j)
      {
      access::rw(m_local.values[offset]) = eT(0);
      
      ++offset;
      }
    }
  
  m_local.remove_zeros();
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::ones()
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> tmp(m.n_rows, U_ci.M.n_elem, fill::ones);
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::operator+= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp1(*this);
  
  Mat<eT> tmp2(tmp1.n_rows, tmp1.n_cols, arma_nozeros_indicator());  tmp2.fill(val);
  
  const Mat<eT> tmp3 = tmp1 + tmp2;
  
  (*this).operator=(tmp3);
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::operator-= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp1(*this);
  
  Mat<eT> tmp2(tmp1.n_rows, tmp1.n_cols, arma_nozeros_indicator());  tmp2.fill(val);
  
  const Mat<eT> tmp3 = tmp1 - tmp2;
  
  (*this).operator=(tmp3);
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::operator*= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(val == eT(0))  { (*this).zeros(); return; }
  
  SpMat<eT>& m_local = const_cast< SpMat<eT>& >(m);
  
  const umat&  ci        = U_ci.M;
  const uword* ci_mem    = ci.memptr();
  const uword  ci_n_elem = ci.n_elem;
  
  m_local.sync_csc();
  m_local.invalidate_cache();
  
  bool has_zero = false;
  
  for(uword ci_count=0; ci_count < ci_n_elem; ++ci_count)
    {
    const uword i = ci_mem[ci_count];
    
    const uword col_n_nonzero = m_local.col_ptrs[i+1] - m_local.col_ptrs[i];
    
    uword offset = m_local.col_ptrs[i];
    
    for(uword j=0; j < col_n_nonzero; ++j)
      {
      eT& m_local_val = access::rw(m_local.values[offset]);
      
      m_local_val *= val;
      
      if(m_local_val == eT(0))  { has_zero = true; }
      
      ++offset;
      }
    }
  
  if(has_zero)  { m_local.remove_zeros(); }
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::operator/= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp1(*this);
  
  Mat<eT> tmp2(tmp1.n_rows, tmp1.n_cols, arma_nozeros_indicator());  tmp2.fill(val);
  
  const SpMat<eT> tmp3 = tmp1 / tmp2;
  
  (*this).operator=(tmp3);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
SpSubview_col_list<eT,T1>::operator= (const Base<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<expr> U(x.get_ref());
  const Mat<eT>&       X = U.M;
  
  SpMat<eT>& m_local = const_cast< SpMat<eT>& >(m);
  
  const umat&  ci        = U_ci.M;
  const uword* ci_mem    = ci.memptr();
  const uword  ci_n_elem = ci.n_elem;
  
  arma_debug_assert_same_size( m_local.n_rows, ci_n_elem, X.n_rows, X.n_cols, "SpMat::cols()" );
  
  const uword X_n_elem = X.n_elem;
  const eT*   X_mem    = X.memptr();
  
  uword X_n_nonzero = 0;
  
  for(uword i=0; i < X_n_elem; ++i)  { X_n_nonzero += (X_mem[i] != eT(0)) ? uword(1) : uword(0); }
  
  SpMat<eT> Y(arma_reserve_indicator(), X.n_rows, m_local.n_cols, X_n_nonzero);
  
  uword count = 0;
  
  for(uword ci_count=0; ci_count < ci_n_elem; ++ci_count)
    {
    const uword i = ci_mem[ci_count];
    
    for(uword row=0; row < X.n_rows; ++row)
      {
      const eT X_val = (*X_mem);  ++X_mem;
      
      if(X_val != eT(0))
        {
        access::rw(Y.row_indices[count]) = row;
        access::rw(Y.values     [count]) = X_val;
        ++count;
        ++access::rw(Y.col_ptrs[i + 1]);
        }
      }
    }
  
  // fix the column pointers
  for(uword i = 0; i < Y.n_cols; ++i)
    {
    access::rw(Y.col_ptrs[i+1]) += Y.col_ptrs[i];
    }
  
  (*this).zeros();
  
  SpMat<eT> tmp = m_local + Y;
  
  m_local.steal_mem(tmp);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
SpSubview_col_list<eT,T1>::operator+= (const Base<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> tmp = SpMat<eT>(*this) + x.get_ref();
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
SpSubview_col_list<eT,T1>::operator-= (const Base<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT> tmp = SpMat<eT>(*this) - x.get_ref();
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
SpSubview_col_list<eT,T1>::operator%= (const Base<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp = SpMat<eT>(*this) % x.get_ref();
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
SpSubview_col_list<eT,T1>::operator/= (const Base<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp = SpMat<eT>(*this) / x.get_ref();
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::operator= (const SpSubview_col_list<eT,T1>& x)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(x);
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
SpSubview_col_list<eT,T1>::operator= (const SpSubview_col_list<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(x);
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
SpSubview_col_list<eT,T1>::operator= (const SpBase<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<expr> U(x.get_ref());
  const SpMat<eT>&     X = U.M;
  
  if(U.is_alias(m))
    {
    const SpMat<eT> tmp(X);
    
    (*this).operator=(tmp);
    
    return;
    }
  
  SpMat<eT>& m_local = const_cast< SpMat<eT>& >(m);
  
  const umat&  ci        = U_ci.M;
  const uword* ci_mem    = ci.memptr();
  const uword  ci_n_elem = ci.n_elem;
  
  arma_debug_assert_same_size( m_local.n_rows, ci_n_elem, X.n_rows, X.n_cols, "SpMat::cols()" );
  
  SpMat<eT> Y(arma_reserve_indicator(), X.n_rows, m_local.n_cols, X.n_nonzero);
  
  uword count = 0;
  
  for(uword ci_count=0; ci_count < ci_n_elem; ++ci_count)
    {
    const uword i = ci_mem[ci_count];
    
    typename SpMat<eT>::const_col_iterator X_col_it     = X.begin_col(ci_count);
    typename SpMat<eT>::const_col_iterator X_col_it_end = X.end_col(ci_count);
    
    while(X_col_it != X_col_it_end)
      {
      access::rw(Y.row_indices[count]) = X_col_it.row();
      access::rw(Y.values     [count]) = (*X_col_it);
      ++count;
      ++access::rw(Y.col_ptrs[i + 1]);
      ++X_col_it;
      }
    }
  
  // fix the column pointers
  for(uword i = 0; i < Y.n_cols; ++i)
    {
    access::rw(Y.col_ptrs[i+1]) += Y.col_ptrs[i];
    }
  
  (*this).zeros();
  
  SpMat<eT> tmp = m_local + Y;
  
  m_local.steal_mem(tmp);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
SpSubview_col_list<eT,T1>::operator+= (const SpBase<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp = SpMat<eT>(*this) + x.get_ref();
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
SpSubview_col_list<eT,T1>::operator-= (const SpBase<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp = SpMat<eT>(*this) - x.get_ref();
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
SpSubview_col_list<eT,T1>::operator%= (const SpBase<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp = SpMat<eT>(*this) % x.get_ref();
  
  (*this).operator=(tmp);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
SpSubview_col_list<eT,T1>::operator/= (const SpBase<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT> tmp(*this);
  
  tmp /= x.get_ref();
  
  (*this).operator=(tmp);
  }



//
//



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::extract(SpMat<eT>& out, const SpSubview_col_list<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: aliasing is handled by SpMat<eT>::operator=(const SpSubview_col_list<eT,T1>& in)
  
  const umat&  ci        = in.U_ci.M;
  const uword* ci_mem    = ci.memptr();
  const uword  ci_n_elem = ci.n_elem;
  
  const SpMat<eT>& in_m = in.m;
  
  in_m.sync_csc();
  
  uword total_n_nonzero = 0;
  
  for(uword ci_count=0; ci_count < ci_n_elem; ++ci_count)
    {
    const uword i = ci_mem[ci_count];
    
    const uword col_n_nonzero = in_m.col_ptrs[i+1] - in_m.col_ptrs[i];
    
    total_n_nonzero += col_n_nonzero;
    }
  
  out.reserve(in.m.n_rows, ci_n_elem, total_n_nonzero);
  
  uword out_n_nonzero = 0;
  uword out_col_count = 0;
  
  for(uword ci_count=0; ci_count < ci_n_elem; ++ci_count)
    {
    const uword i = ci_mem[ci_count];
    
    const uword col_n_nonzero = in_m.col_ptrs[i+1] - in_m.col_ptrs[i];
    
    uword offset = in_m.col_ptrs[i];
    
    for(uword j=0; j < col_n_nonzero; ++j)
      {
      const eT    val = in_m.values     [ offset ];
      const uword row = in_m.row_indices[ offset ];
      
      ++offset;
      
      access::rw(out.values     [out_n_nonzero]) = val;
      access::rw(out.row_indices[out_n_nonzero]) = row;
      
      access::rw(out.col_ptrs[out_col_count+1])++;
      
      ++out_n_nonzero;
      }
    
    ++out_col_count;
    }
  
  // fix the column pointers
  for(uword i = 0; i < out.n_cols; ++i)
    {
    access::rw(out.col_ptrs[i+1]) += out.col_ptrs[i];
    }
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::plus_inplace(SpMat<eT>& out, const SpSubview_col_list& in)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(in);
  
  out += tmp;
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::minus_inplace(SpMat<eT>& out, const SpSubview_col_list& in)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(in);
  
  out -= tmp;
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::schur_inplace(SpMat<eT>& out, const SpSubview_col_list& in)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(in);
  
  out %= tmp;
  }



template<typename eT, typename T1>
inline
void
SpSubview_col_list<eT,T1>::div_inplace(SpMat<eT>& out, const SpSubview_col_list& in)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(in);
  
  out /= tmp;
  }



//! @}
