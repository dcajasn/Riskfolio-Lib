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


//! \addtogroup spdiagview
//! @{


template<typename eT>
inline
spdiagview<eT>::~spdiagview()
  {
  arma_extra_debug_sigprint();
  }


template<typename eT>
arma_inline
spdiagview<eT>::spdiagview(const SpMat<eT>& in_m, const uword in_row_offset, const uword in_col_offset, const uword in_len)
  : m(in_m)
  , row_offset(in_row_offset)
  , col_offset(in_col_offset)
  , n_rows(in_len)
  , n_elem(in_len)
  {
  arma_extra_debug_sigprint();
  }



//! set a diagonal of our matrix using a diagonal from a foreign matrix
template<typename eT>
inline
void
spdiagview<eT>::operator= (const spdiagview<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  spdiagview<eT>& d = *this;
  
  arma_debug_check( (d.n_elem != x.n_elem), "spdiagview: diagonals have incompatible lengths" );
  
        SpMat<eT>& d_m = const_cast< SpMat<eT>& >(d.m);
  const SpMat<eT>& x_m = x.m;
  
  if( (&d_m == &x_m) || ((d.row_offset == 0) && (d.col_offset == 0)) )
    {
    const Mat<eT> tmp(x);
    
    (*this).operator=(tmp);
    }
  else
    {
    const uword d_n_elem     = d.n_elem;
    const uword d_row_offset = d.row_offset;
    const uword d_col_offset = d.col_offset;
    
    const uword x_row_offset = x.row_offset;
    const uword x_col_offset = x.col_offset;
    
    for(uword i=0; i < d_n_elem; ++i)
      {
      d_m.at(i + d_row_offset, i + d_col_offset) = x_m.at(i + x_row_offset, i + x_col_offset);
      }
    }
  }



template<typename eT>
inline
void
spdiagview<eT>::operator+=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(val == eT(0))  { return; }
  
  SpMat<eT>& t_m = const_cast< SpMat<eT>& >(m);
  
  const uword t_n_elem     = n_elem;
  const uword t_row_offset = row_offset;
  const uword t_col_offset = col_offset;
  
  for(uword i=0; i < t_n_elem; ++i)
    {
    t_m.at(i + t_row_offset, i + t_col_offset) += val;
    }
  }



template<typename eT>
inline
void
spdiagview<eT>::operator-=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(val == eT(0))  { return; }
  
  SpMat<eT>& t_m = const_cast< SpMat<eT>& >(m);
  
  const uword t_n_elem     = n_elem;
  const uword t_row_offset = row_offset;
  const uword t_col_offset = col_offset;
  
  for(uword i=0; i < t_n_elem; ++i)
    {
    t_m.at(i + t_row_offset, i + t_col_offset) -= val;
    }
  }



template<typename eT>
inline
void
spdiagview<eT>::operator*=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(val == eT(0))  { (*this).zeros(); return; }
  
  SpMat<eT>& t_m = const_cast< SpMat<eT>& >(m);
  
  const uword t_n_elem     = n_elem;
  const uword t_row_offset = row_offset;
  const uword t_col_offset = col_offset;
  
  for(uword i=0; i < t_n_elem; ++i)
    {
    t_m.at(i + t_row_offset, i + t_col_offset) *= val;
    }
  }



template<typename eT>
inline
void
spdiagview<eT>::operator/=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT>& t_m = const_cast< SpMat<eT>& >(m);
  
  const uword t_n_elem     = n_elem;
  const uword t_row_offset = row_offset;
  const uword t_col_offset = col_offset;
  
  for(uword i=0; i < t_n_elem; ++i)
    {
    t_m.at(i + t_row_offset, i + t_col_offset) /= val;
    }
  }



//! set a diagonal of our matrix using data from a foreign object
template<typename eT>
template<typename T1>
inline
void
spdiagview<eT>::operator= (const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  spdiagview<eT>& d = *this;
  
  SpMat<eT>& d_m = const_cast< SpMat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
  
  if(is_same_type< T1, Gen<Col<eT>, gen_zeros> >::yes)
    {
    const Proxy<T1> P(o.get_ref());
    
    arma_debug_check( (d_n_elem != P.get_n_elem()), "spdiagview: given object has incompatible size" );
    
    (*this).zeros();
    
    return;
    }
  
  if(is_same_type< T1, Gen<Col<eT>, gen_ones> >::yes)
    {
    const Proxy<T1> P(o.get_ref());
    
    arma_debug_check( (d_n_elem != P.get_n_elem()), "spdiagview: given object has incompatible size" );
    
    (*this).ones();
    
    return;
    }
  
  const quasi_unwrap<T1> U(o.get_ref());
  const Mat<eT>& x     = U.M;
  
  const eT* x_mem = x.memptr();
  
  arma_debug_check
    (
    ( (d_n_elem != x.n_elem) || ((x.n_rows != 1) && (x.n_cols != 1)) ),
    "spdiagview: given object has incompatible size"
    );
  
  if( (d_row_offset == 0) && (d_col_offset == 0) )
    {
    SpMat<eT> tmp1;
    
    tmp1.eye(d_m.n_rows, d_m.n_cols);
    
    bool has_zero = false;
    
    for(uword i=0; i < d_n_elem; ++i)
      {
      const eT val = x_mem[i];
      
      access::rw(tmp1.values[i]) = val;
      
      if(val == eT(0))  { has_zero = true; }
      }
    
    if(has_zero)  { tmp1.remove_zeros(); }
    
    if(tmp1.n_nonzero == 0)  { (*this).zeros(); return; }
    
    SpMat<eT> tmp2;
    
    spglue_merge::diagview_merge(tmp2, d_m, tmp1);
    
    d_m.steal_mem(tmp2);
    }
  else
    {
    for(uword i=0; i < d_n_elem; ++i)
      {
      d_m.at(i + d_row_offset, i + d_col_offset) = x_mem[i];
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
spdiagview<eT>::operator+=(const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  spdiagview<eT>& d = *this;
  
  SpMat<eT>& d_m = const_cast< SpMat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
    
  const Proxy<T1> P( o.get_ref() );
  
  arma_debug_check
    (
    ( (d_n_elem != P.get_n_elem()) || ((P.get_n_rows() != 1) && (P.get_n_cols() != 1)) ),
    "spdiagview: given object has incompatible size"
    );
  
  if( (is_Mat<typename Proxy<T1>::stored_type>::value) || (Proxy<T1>::use_at) )
    {
    const unwrap<typename Proxy<T1>::stored_type> tmp(P.Q);
    const Mat<eT>& x = tmp.M;
    
    const eT* x_mem = x.memptr();

    for(uword i=0; i < d_n_elem; ++i)
      {
      d_m.at(i + d_row_offset, i + d_col_offset) += x_mem[i];
      }
    }
  else
    {
    typename Proxy<T1>::ea_type Pea = P.get_ea();
      
    for(uword i=0; i < d_n_elem; ++i)
      {
      d_m.at(i + d_row_offset, i + d_col_offset) += Pea[i];
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
spdiagview<eT>::operator-=(const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  spdiagview<eT>& d = *this;
  
  SpMat<eT>& d_m = const_cast< SpMat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
    
  const Proxy<T1> P( o.get_ref() );
  
  arma_debug_check
    (
    ( (d_n_elem != P.get_n_elem()) || ((P.get_n_rows() != 1) && (P.get_n_cols() != 1)) ),
    "spdiagview: given object has incompatible size"
    );
  
  if( (is_Mat<typename Proxy<T1>::stored_type>::value) || (Proxy<T1>::use_at) )
    {
    const unwrap<typename Proxy<T1>::stored_type> tmp(P.Q);
    const Mat<eT>& x = tmp.M;
    
    const eT* x_mem = x.memptr();

    for(uword i=0; i < d_n_elem; ++i)
      {
      d_m.at(i + d_row_offset, i + d_col_offset) -= x_mem[i];
      }
    }
  else
    {
    typename Proxy<T1>::ea_type Pea = P.get_ea();
      
    for(uword i=0; i < d_n_elem; ++i)
      {
      d_m.at(i + d_row_offset, i + d_col_offset) -= Pea[i];
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
spdiagview<eT>::operator%=(const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  spdiagview<eT>& d = *this;
  
  SpMat<eT>& d_m = const_cast< SpMat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
    
  const Proxy<T1> P( o.get_ref() );
  
  arma_debug_check
    (
    ( (d_n_elem != P.get_n_elem()) || ((P.get_n_rows() != 1) && (P.get_n_cols() != 1)) ),
    "spdiagview: given object has incompatible size"
    );
  
  if( (is_Mat<typename Proxy<T1>::stored_type>::value) || (Proxy<T1>::use_at) )
    {
    const unwrap<typename Proxy<T1>::stored_type> tmp(P.Q);
    const Mat<eT>& x = tmp.M;
    
    const eT* x_mem = x.memptr();

    for(uword i=0; i < d_n_elem; ++i)
      {
      d_m.at(i + d_row_offset, i + d_col_offset) *= x_mem[i];
      }
    }
  else
    {
    typename Proxy<T1>::ea_type Pea = P.get_ea();
      
    for(uword i=0; i < d_n_elem; ++i)
      {
      d_m.at(i + d_row_offset, i + d_col_offset) *= Pea[i];
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
spdiagview<eT>::operator/=(const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  spdiagview<eT>& d = *this;
  
  SpMat<eT>& d_m = const_cast< SpMat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
    
  const Proxy<T1> P( o.get_ref() );
  
  arma_debug_check
    (
    ( (d_n_elem != P.get_n_elem()) || ((P.get_n_rows() != 1) && (P.get_n_cols() != 1)) ),
    "spdiagview: given object has incompatible size"
    );
  
  if( (is_Mat<typename Proxy<T1>::stored_type>::value) || (Proxy<T1>::use_at) )
    {
    const unwrap<typename Proxy<T1>::stored_type> tmp(P.Q);
    const Mat<eT>& x = tmp.M;
    
    const eT* x_mem = x.memptr();

    for(uword i=0; i < d_n_elem; ++i)
      {
      d_m.at(i + d_row_offset, i + d_col_offset) /= x_mem[i];
      }
    }
  else
    {
    typename Proxy<T1>::ea_type Pea = P.get_ea();
      
    for(uword i=0; i < d_n_elem; ++i)
      {
      d_m.at(i + d_row_offset, i + d_col_offset) /= Pea[i];
      }
    }
  }



//! set a diagonal of our matrix using data from a foreign object
template<typename eT>
template<typename T1>
inline
void
spdiagview<eT>::operator= (const SpBase<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> U( o.get_ref() );
  const SpMat<eT>& x   = U.M;
  
  arma_debug_check
    (
    ( (n_elem != x.n_elem) || ((x.n_rows != 1) && (x.n_cols != 1)) ),
    "spdiagview: given object has incompatible size"
    );
  
  const Mat<eT> tmp(x);
  
  (*this).operator=(tmp);
  }



template<typename eT>
template<typename T1>
inline
void
spdiagview<eT>::operator+=(const SpBase<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  spdiagview<eT>& d = *this;
  
  SpMat<eT>& d_m = const_cast< SpMat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
  
  const SpProxy<T1> P( o.get_ref() );
  
  arma_debug_check
    (
    ( (d_n_elem != P.get_n_elem()) || ((P.get_n_rows() != 1) && (P.get_n_cols() != 1)) ),
    "spdiagview: given object has incompatible size"
    );
  
  if( SpProxy<T1>::use_iterator || P.is_alias(d_m) )
    {
    const SpMat<eT> tmp(P.Q);
    
    if(tmp.n_cols == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) += tmp.at(i,0); }
      }
    else
    if(tmp.n_rows == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) += tmp.at(0,i); }
      }
    }
  else
    {
    if(P.get_n_cols() == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) += P.at(i,0); }
      }
    else
    if(P.get_n_rows() == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) += P.at(0,i); }
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
spdiagview<eT>::operator-=(const SpBase<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  spdiagview<eT>& d = *this;
  
  SpMat<eT>& d_m = const_cast< SpMat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
  
  const SpProxy<T1> P( o.get_ref() );
  
  arma_debug_check
    (
    ( (d_n_elem != P.get_n_elem()) || ((P.get_n_rows() != 1) && (P.get_n_cols() != 1)) ),
    "spdiagview: given object has incompatible size"
    );
  
  if( SpProxy<T1>::use_iterator || P.is_alias(d_m) )
    {
    const SpMat<eT> tmp(P.Q);
    
    if(tmp.n_cols == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) -= tmp.at(i,0); }
      }
    else
    if(tmp.n_rows == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) -= tmp.at(0,i); }
      }
    }
  else
    {
    if(P.get_n_cols() == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) -= P.at(i,0); }
      }
    else
    if(P.get_n_rows() == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) -= P.at(0,i); }
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
spdiagview<eT>::operator%=(const SpBase<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  spdiagview<eT>& d = *this;
  
  SpMat<eT>& d_m = const_cast< SpMat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
  
  const SpProxy<T1> P( o.get_ref() );
  
  arma_debug_check
    (
    ( (d_n_elem != P.get_n_elem()) || ((P.get_n_rows() != 1) && (P.get_n_cols() != 1)) ),
    "spdiagview: given object has incompatible size"
    );
  
  if( SpProxy<T1>::use_iterator || P.is_alias(d_m) )
    {
    const SpMat<eT> tmp(P.Q);
    
    if(tmp.n_cols == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) *= tmp.at(i,0); }
      }
    else
    if(tmp.n_rows == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) *= tmp.at(0,i); }
      }
    }
  else
    {
    if(P.get_n_cols() == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) *= P.at(i,0); }
      }
    else
    if(P.get_n_rows() == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) *= P.at(0,i); }
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
spdiagview<eT>::operator/=(const SpBase<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  spdiagview<eT>& d = *this;
  
  SpMat<eT>& d_m = const_cast< SpMat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
  
  const SpProxy<T1> P( o.get_ref() );
  
  arma_debug_check
    (
    ( (d_n_elem != P.get_n_elem()) || ((P.get_n_rows() != 1) && (P.get_n_cols() != 1)) ),
    "spdiagview: given object has incompatible size"
    );
  
  if( SpProxy<T1>::use_iterator || P.is_alias(d_m) )
    {
    const SpMat<eT> tmp(P.Q);
    
    if(tmp.n_cols == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) /= tmp.at(i,0); }
      }
    else
    if(tmp.n_rows == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) /= tmp.at(0,i); }
      }
    }
  else
    {
    if(P.get_n_cols() == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) /= P.at(i,0); }
      }
    else
    if(P.get_n_rows() == 1)
      {
      for(uword i=0; i < d_n_elem; ++i)  { d_m.at(i + d_row_offset, i + d_col_offset) /= P.at(0,i); }
      }
    }
  }



template<typename eT>
inline
void
spdiagview<eT>::extract(SpMat<eT>& out, const spdiagview<eT>& d)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT>& d_m = d.m;
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
  
  Col<eT> cache(d_n_elem, arma_nozeros_indicator());
  eT* cache_mem = cache.memptr();
  
  uword d_n_nonzero = 0;
  
  for(uword i=0; i < d_n_elem; ++i)
    {
    const eT val = d_m.at(i + d_row_offset, i + d_col_offset);
    
    cache_mem[i] = val;
    
    d_n_nonzero += (val != eT(0)) ? uword(1) : uword(0);
    }
  
  out.reserve(d_n_elem, 1, d_n_nonzero);
  
  uword count = 0;
  for(uword i=0; i < d_n_elem; ++i)
    {
    const eT val = cache_mem[i];
    
    if(val != eT(0))
      {
      access::rw(out.row_indices[count]) = i;
      access::rw(out.values[count])      = val;
      ++count;
      }
    }
  
  access::rw(out.col_ptrs[0]) = 0;
  access::rw(out.col_ptrs[1]) = d_n_nonzero;
  }



//! extract a diagonal and store it as a dense column vector
template<typename eT>
inline
void
spdiagview<eT>::extract(Mat<eT>& out, const spdiagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: we're assuming that the 'out' matrix has already been set to the correct size;
  // size setting is done by either the Mat contructor or Mat::operator=()
  
  const SpMat<eT>& in_m = in.m;
  
  const uword in_n_elem     = in.n_elem;
  const uword in_row_offset = in.row_offset;
  const uword in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  for(uword i=0; i < in_n_elem; ++i)
    {
    out_mem[i] = in_m.at(i + in_row_offset, i + in_col_offset);
    }
  }



template<typename eT>
inline
SpMat_MapMat_val<eT>
spdiagview<eT>::operator[](const uword i)
  {
  return (const_cast< SpMat<eT>& >(m)).at(i+row_offset, i+col_offset);
  }



template<typename eT>
inline
eT
spdiagview<eT>::operator[](const uword i) const
  {
  return m.at(i+row_offset, i+col_offset);
  }



template<typename eT>
inline
SpMat_MapMat_val<eT>
spdiagview<eT>::at(const uword i)
  {
  return (const_cast< SpMat<eT>& >(m)).at(i+row_offset, i+col_offset);
  }



template<typename eT>
inline
eT
spdiagview<eT>::at(const uword i) const
  {
  return m.at(i+row_offset, i+col_offset);
  }



template<typename eT>
inline
SpMat_MapMat_val<eT>
spdiagview<eT>::operator()(const uword i)
  {
  arma_debug_check_bounds( (i >= n_elem), "spdiagview::operator(): out of bounds" );
  
  return (const_cast< SpMat<eT>& >(m)).at(i+row_offset, i+col_offset);
  }



template<typename eT>
inline
eT
spdiagview<eT>::operator()(const uword i) const
  {
  arma_debug_check_bounds( (i >= n_elem), "spdiagview::operator(): out of bounds" );
  
  return m.at(i+row_offset, i+col_offset);
  }



template<typename eT>
inline
SpMat_MapMat_val<eT>
spdiagview<eT>::at(const uword row, const uword)
  {
  return (const_cast< SpMat<eT>& >(m)).at(row+row_offset, row+col_offset);
  }



template<typename eT>
inline
eT
spdiagview<eT>::at(const uword row, const uword) const
  {
  return m.at(row+row_offset, row+col_offset);
  }



template<typename eT>
inline
SpMat_MapMat_val<eT>
spdiagview<eT>::operator()(const uword row, const uword col)
  {
  arma_debug_check_bounds( ((row >= n_elem) || (col > 0)), "spdiagview::operator(): out of bounds" );
  
  return (const_cast< SpMat<eT>& >(m)).at(row+row_offset, row+col_offset);
  }



template<typename eT>
inline
eT
spdiagview<eT>::operator()(const uword row, const uword col) const
  {
  arma_debug_check_bounds( ((row >= n_elem) || (col > 0)), "spdiagview::operator(): out of bounds" );
  
  return m.at(row+row_offset, row+col_offset);
  }



template<typename eT>
inline
void
spdiagview<eT>::replace(const eT old_val, const eT new_val)
  {
  arma_extra_debug_sigprint();
  
  if(old_val == eT(0))
    {
    arma_debug_warn_level(1, "spdiagview::replace(): replacement not done, as old_val = 0");
    }
  else
    {
    Mat<eT> tmp(*this);
    
    tmp.replace(old_val, new_val);
    
    (*this).operator=(tmp);
    }
  }



template<typename eT>
inline
void
spdiagview<eT>::clean(const typename get_pod_type<eT>::result threshold)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp(*this);
  
  tmp.clean(threshold);
  
  (*this).operator=(tmp);
  }



template<typename eT>
inline
void
spdiagview<eT>::clamp(const eT min_val, const eT max_val)
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT> tmp(*this);
  
  tmp.clamp(min_val, max_val);
  
  (*this).operator=(tmp);
  }



template<typename eT>
inline
void
spdiagview<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  if( (row_offset == 0) && (col_offset == 0) && (m.sync_state != 1) )
    {
    if(val == eT(0))
      {
      SpMat<eT> tmp(arma_reserve_indicator(), m.n_rows, m.n_cols, m.n_nonzero);  // worst case scenario
      
      typename SpMat<eT>::const_iterator it     = m.begin();
      typename SpMat<eT>::const_iterator it_end = m.end();
      
      uword count = 0;
        
      for(; it != it_end; ++it)
        {
        const uword row = it.row();
        const uword col = it.col();
        
        if(row != col)
          {
          access::rw(tmp.values[count])      = (*it);
          access::rw(tmp.row_indices[count]) = row;
          access::rw(tmp.col_ptrs[col + 1])++;
          ++count;
          }
        }
      
      for(uword i=0; i < tmp.n_cols; ++i)
        {
        access::rw(tmp.col_ptrs[i + 1]) += tmp.col_ptrs[i];
        }
      
      // quick resize without reallocating memory and copying data
      access::rw(         tmp.n_nonzero) = count;
      access::rw(     tmp.values[count]) = eT(0);
      access::rw(tmp.row_indices[count]) = uword(0);
      
      access::rw(m).steal_mem(tmp);
      }
    else  // val != eT(0)
      {
      SpMat<eT> tmp1;
      
      tmp1.eye(m.n_rows, m.n_cols);
      
      if(val != eT(1))  { tmp1 *= val; }
      
      SpMat<eT> tmp2;
      
      spglue_merge::diagview_merge(tmp2, m, tmp1);
      
      access::rw(m).steal_mem(tmp2);
      }
    }
  else
    {
    SpMat<eT>& x = const_cast< SpMat<eT>& >(m);
    
    const uword local_n_elem = n_elem;
    
    for(uword i=0; i < local_n_elem; ++i)
      {
      x.at(i+row_offset, i+col_offset) = val;
      }
    }
  }



template<typename eT>
inline
void
spdiagview<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  (*this).fill(eT(0));
  }



template<typename eT>
inline
void
spdiagview<eT>::ones()
  {
  arma_extra_debug_sigprint();
  
  (*this).fill(eT(1));
  }



template<typename eT>
inline
void
spdiagview<eT>::randu()
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT>& x = const_cast< SpMat<eT>& >(m);
  
  const uword local_n_elem = n_elem;
  
  for(uword i=0; i < local_n_elem; ++i)
    {
    x.at(i+row_offset, i+col_offset) = eT(arma_rng::randu<eT>());
    }
  }



template<typename eT>
inline
void
spdiagview<eT>::randn()
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT>& x = const_cast< SpMat<eT>& >(m);
  
  const uword local_n_elem = n_elem;
  
  for(uword i=0; i < local_n_elem; ++i)
    {
    x.at(i+row_offset, i+col_offset) = eT(arma_rng::randn<eT>());
    }
  }



//! @}
