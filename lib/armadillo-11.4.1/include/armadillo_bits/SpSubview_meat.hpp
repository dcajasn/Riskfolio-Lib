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


//! \addtogroup SpSubview
//! @{


template<typename eT>
inline
SpSubview<eT>::~SpSubview()
  {
  arma_extra_debug_sigprint_this(this);
  }



template<typename eT>
inline
SpSubview<eT>::SpSubview(const SpMat<eT>& in_m, const uword in_row1, const uword in_col1, const uword in_n_rows, const uword in_n_cols)
  : m(in_m)
  , aux_row1(in_row1)
  , aux_col1(in_col1)
  , n_rows(in_n_rows)
  , n_cols(in_n_cols)
  , n_elem(in_n_rows * in_n_cols)
  , n_nonzero(0)
  {
  arma_extra_debug_sigprint_this(this);
  
  m.sync_csc();
  
  // There must be a O(1) way to do this
  uword lend     = m.col_ptrs[in_col1 + in_n_cols];
  uword lend_row = in_row1 + in_n_rows;
  uword count   = 0;
  
  for(uword i = m.col_ptrs[in_col1]; i < lend; ++i)
    {
    const uword m_row_indices_i = m.row_indices[i];
    
    const bool condition = (m_row_indices_i >= in_row1) && (m_row_indices_i < lend_row);
    
    count += condition ? uword(1) : uword(0);
    }
  
  access::rw(n_nonzero) = count;
  }



template<typename eT>
inline
SpSubview<eT>::SpSubview(const SpSubview<eT>& in)
  : m        (in.m        )
  , aux_row1 (in.aux_row1 )
  , aux_col1 (in.aux_col1 )
  , n_rows   (in.n_rows   )
  , n_cols   (in.n_cols   )
  , n_elem   (in.n_elem   )
  , n_nonzero(in.n_nonzero)
  {
  arma_extra_debug_sigprint(arma_str::format("this = %x   in = %x") % this % &in);
  }



template<typename eT>
inline
SpSubview<eT>::SpSubview(SpSubview<eT>&& in)
  : m        (in.m        )
  , aux_row1 (in.aux_row1 )
  , aux_col1 (in.aux_col1 )
  , n_rows   (in.n_rows   )
  , n_cols   (in.n_cols   )
  , n_elem   (in.n_elem   )
  , n_nonzero(in.n_nonzero)
  {
  arma_extra_debug_sigprint(arma_str::format("this = %x   in = %x") % this % &in);
  
  // for paranoia
  
  access::rw(in.aux_row1 ) = 0;
  access::rw(in.aux_col1 ) = 0;
  access::rw(in.n_rows   ) = 0;
  access::rw(in.n_cols   ) = 0;
  access::rw(in.n_elem   ) = 0;
  access::rw(in.n_nonzero) = 0;
  }



template<typename eT>
inline
const SpSubview<eT>&
SpSubview<eT>::operator+=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(val == eT(0))  { return *this; }
  
  Mat<eT> tmp( (*this).n_rows, (*this).n_cols, arma_nozeros_indicator() );
  
  tmp.fill(val);
  
  return (*this).operator=( (*this) + tmp );
  }



template<typename eT>
inline
const SpSubview<eT>&
SpSubview<eT>::operator-=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(val == eT(0))  { return *this; }
  
  Mat<eT> tmp( (*this).n_rows, (*this).n_cols, arma_nozeros_indicator() );
  
  tmp.fill(val);
  
  return (*this).operator=( (*this) - tmp );
  }



template<typename eT>
inline
const SpSubview<eT>&
SpSubview<eT>::operator*=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(val == eT(0))  { (*this).zeros(); return *this; }
  
  if((n_elem == 0) || (n_nonzero == 0))  { return *this; }
  
  m.sync_csc();
  m.invalidate_cache();
  
  const uword lstart_row = aux_row1;
  const uword lend_row   = aux_row1 + n_rows;
  
  const uword lstart_col = aux_col1;
  const uword lend_col   = aux_col1 + n_cols;
  
  const uword* m_row_indices = m.row_indices;
        eT*    m_values      = access::rwp(m.values);
  
  bool has_zero = false;
  
  for(uword c = lstart_col; c < lend_col; ++c)
    {
    const uword r_start = m.col_ptrs[c    ];
    const uword r_end   = m.col_ptrs[c + 1];
    
    for(uword r = r_start; r < r_end; ++r)
      {
      const uword m_row_indices_r = m_row_indices[r];
      
      if( (m_row_indices_r >= lstart_row) && (m_row_indices_r < lend_row) )
        {
        eT& m_values_r = m_values[r];
        
        m_values_r *= val;
        
        if(m_values_r == eT(0))  { has_zero = true; }
        }
      }
    }
  
  if(has_zero)
    {
    const uword old_m_n_nonzero = m.n_nonzero;
    
    access::rw(m).remove_zeros();
    
    if(m.n_nonzero != old_m_n_nonzero)
      {
      access::rw(n_nonzero) = n_nonzero - (old_m_n_nonzero - m.n_nonzero); 
      }
    }
  
  return *this;
  }



template<typename eT>
inline
const SpSubview<eT>&
SpSubview<eT>::operator/=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (val == eT(0)), "element-wise division: division by zero" );
  
  m.sync_csc();
  m.invalidate_cache();
  
  const uword lstart_row = aux_row1;
  const uword lend_row   = aux_row1 + n_rows;
  
  const uword lstart_col = aux_col1;
  const uword lend_col   = aux_col1 + n_cols;
  
  const uword* m_row_indices = m.row_indices;
        eT*    m_values      = access::rwp(m.values);
  
  bool has_zero = false;
  
  for(uword c = lstart_col; c < lend_col; ++c)
    {
    const uword r_start = m.col_ptrs[c    ];
    const uword r_end   = m.col_ptrs[c + 1];
    
    for(uword r = r_start; r < r_end; ++r)
      {
      const uword m_row_indices_r = m_row_indices[r];
      
      if( (m_row_indices_r >= lstart_row) && (m_row_indices_r < lend_row) )
        {
        eT& m_values_r = m_values[r];
        
        m_values_r /= val;
        
        if(m_values_r == eT(0))  { has_zero = true; }
        }
      }
    }
  
  if(has_zero)
    {
    const uword old_m_n_nonzero = m.n_nonzero;
    
    access::rw(m).remove_zeros();
    
    if(m.n_nonzero != old_m_n_nonzero)
      {
      access::rw(n_nonzero) = n_nonzero - (old_m_n_nonzero - m.n_nonzero); 
      }
    }
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator=(const Base<eT, T1>& in)
  {
  arma_extra_debug_sigprint();
  
  if(is_same_type< T1, Gen<Mat<eT>, gen_zeros> >::yes)
    {
    const Proxy<T1> P(in.get_ref());
    
    arma_debug_assert_same_size(n_rows, n_cols, P.get_n_rows(), P.get_n_cols(), "insertion into sparse submatrix");
    
    (*this).zeros();
    
    return *this;
    }
  
  if(is_same_type< T1, Gen<Mat<eT>, gen_eye> >::yes)
    {
    const Proxy<T1> P(in.get_ref());
    
    arma_debug_assert_same_size(n_rows, n_cols, P.get_n_rows(), P.get_n_cols(), "insertion into sparse submatrix");
    
    (*this).eye();
    
    return *this;
    }
  
  const quasi_unwrap<T1> U(in.get_ref());
  
  arma_debug_assert_same_size(n_rows, n_cols, U.M.n_rows, U.M.n_cols, "insertion into sparse submatrix");
  
  spglue_merge::subview_merge(*this, U.M);
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator+=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  return (*this).operator=( (*this) + x.get_ref() );
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator-=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  return (*this).operator=( (*this) - x.get_ref() );
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator*=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT> tmp(*this);
  
  tmp *= x.get_ref();
  
  return (*this).operator=(tmp);
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator%=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  return (*this).operator=( (*this) % x.get_ref() );
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator/=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  return (*this).operator=( (*this) / x.get_ref() );
  }



template<typename eT>
inline
const SpSubview<eT>&
SpSubview<eT>::operator=(const SpSubview<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  return (*this).operator_equ_common(x);
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator=(const SpBase<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  return (*this).operator_equ_common( x.get_ref() );
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator_equ_common(const SpBase<eT, T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> U(in.get_ref());
  
  arma_debug_assert_same_size(n_rows, n_cols, U.M.n_rows, U.M.n_cols, "insertion into sparse submatrix");
  
  if(U.is_alias(m))
    {
    const SpMat<eT> tmp(U.M);
    
    spglue_merge::subview_merge(*this, tmp);
    }
  else
    {
    spglue_merge::subview_merge(*this, U.M);
    }
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator+=(const SpBase<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  // TODO: implement dedicated machinery
  return (*this).operator=( (*this) + x.get_ref() );
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator-=(const SpBase<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  // TODO: implement dedicated machinery
  return (*this).operator=( (*this) - x.get_ref() );
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator*=(const SpBase<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  return (*this).operator=( (*this) * x.get_ref() );
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator%=(const SpBase<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  // TODO: implement dedicated machinery
  return (*this).operator=( (*this) % x.get_ref() );
  }



template<typename eT>
template<typename T1>
inline
const SpSubview<eT>&
SpSubview<eT>::operator/=(const SpBase<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: use of this function is not advised; it is implemented only for completeness
  
  SpProxy<T1> p(x.get_ref());
  
  arma_debug_assert_same_size(n_rows, n_cols, p.get_n_rows(), p.get_n_cols(), "element-wise division");
  
  if(p.is_alias(m) == false)
    {
    for(uword lcol = 0; lcol < n_cols; ++lcol)
    for(uword lrow = 0; lrow < n_rows; ++lrow)
      {
      at(lrow,lcol) /= p.at(lrow,lcol);
      }
    }
  else
    {
    const SpMat<eT> tmp(p.Q);
    
    (*this).operator/=(tmp);
    }
  
  return *this;
  }



//! apply a functor to each element
template<typename eT>
template<typename functor>
inline
void
SpSubview<eT>::for_each(functor F)
  {
  arma_extra_debug_sigprint();
  
  m.sync_csc();
  m.invalidate_cache();
  
  const uword lstart_row = aux_row1;
  const uword lend_row   = aux_row1 + n_rows;
  
  const uword lstart_col = aux_col1;
  const uword lend_col   = aux_col1 + n_cols;
  
  const uword* m_row_indices = m.row_indices;
        eT*    m_values      = access::rwp(m.values);
  
  bool has_zero = false;
  
  for(uword c = lstart_col; c < lend_col; ++c)
    {
    const uword r_start = m.col_ptrs[c    ];
    const uword r_end   = m.col_ptrs[c + 1];
    
    for(uword r = r_start; r < r_end; ++r)
      {
      const uword m_row_indices_r = m_row_indices[r];
      
      if( (m_row_indices_r >= lstart_row) && (m_row_indices_r < lend_row) )
        {
        eT& m_values_r = m_values[r];
        
        F(m_values_r);
        
        if(m_values_r == eT(0))  { has_zero = true; }
        }
      }
    }
  
  if(has_zero)
    {
    const uword old_m_n_nonzero = m.n_nonzero;
    
    access::rw(m).remove_zeros();
    
    if(m.n_nonzero != old_m_n_nonzero)
      {
      access::rw(n_nonzero) = n_nonzero - (old_m_n_nonzero - m.n_nonzero); 
      }
    }
  }



template<typename eT>
template<typename functor>
inline
void
SpSubview<eT>::for_each(functor F) const
  {
  arma_extra_debug_sigprint();
  
  m.sync_csc();
  
  const uword lstart_row = aux_row1;
  const uword lend_row   = aux_row1 + n_rows;
  
  const uword lstart_col = aux_col1;
  const uword lend_col   = aux_col1 + n_cols;
  
  const uword* m_row_indices = m.row_indices;
  
  for(uword c = lstart_col; c < lend_col; ++c)
    {
    const uword r_start = m.col_ptrs[c    ];
    const uword r_end   = m.col_ptrs[c + 1];
    
    for(uword r = r_start; r < r_end; ++r)
      {
      const uword m_row_indices_r = m_row_indices[r];
      
      if( (m_row_indices_r >= lstart_row) && (m_row_indices_r < lend_row) )
        {
        F(m.values[r]);
        }
      }
    }
  }



//! transform each element using a functor
template<typename eT>
template<typename functor>
inline
void
SpSubview<eT>::transform(functor F)
  {
  arma_extra_debug_sigprint();
  
  m.sync_csc();
  m.invalidate_cache();
  
  const uword lstart_row = aux_row1;
  const uword lend_row   = aux_row1 + n_rows;
  
  const uword lstart_col = aux_col1;
  const uword lend_col   = aux_col1 + n_cols;
  
  const uword* m_row_indices = m.row_indices;
        eT*    m_values      = access::rwp(m.values);
  
  bool has_zero = false;
  
  for(uword c = lstart_col; c < lend_col; ++c)
    {
    const uword r_start = m.col_ptrs[c    ];
    const uword r_end   = m.col_ptrs[c + 1];
    
    for(uword r = r_start; r < r_end; ++r)
      {
      const uword m_row_indices_r = m_row_indices[r];
      
      if( (m_row_indices_r >= lstart_row) && (m_row_indices_r < lend_row) )
        {
        eT& m_values_r = m_values[r];
        
        m_values_r = eT( F(m_values_r) );
        
        if(m_values_r == eT(0))  { has_zero = true; }
        }
      }
    }
  
  if(has_zero)
    {
    const uword old_m_n_nonzero = m.n_nonzero;
    
    access::rw(m).remove_zeros();
    
    if(m.n_nonzero != old_m_n_nonzero)
      {
      access::rw(n_nonzero) = n_nonzero - (old_m_n_nonzero - m.n_nonzero); 
      }
    }
  }



template<typename eT>
inline
void
SpSubview<eT>::replace(const eT old_val, const eT new_val)
  {
  arma_extra_debug_sigprint();
  
  if(old_val == eT(0))
    {
    if(new_val != eT(0))
      {
      Mat<eT> tmp(*this);
      
      tmp.replace(old_val, new_val);
      
      (*this).operator=(tmp);
      }
    
    return;
    }
  
  m.sync_csc();
  m.invalidate_cache();
  
  const uword lstart_row = aux_row1;
  const uword lend_row   = aux_row1 + n_rows;
  
  const uword lstart_col = aux_col1;
  const uword lend_col   = aux_col1 + n_cols;
  
  const uword* m_row_indices = m.row_indices;
        eT*    m_values      = access::rwp(m.values);
  
  if(arma_isnan(old_val))
    {
    for(uword c = lstart_col; c < lend_col; ++c)
      {
      const uword r_start = m.col_ptrs[c    ];
      const uword r_end   = m.col_ptrs[c + 1];
      
      for(uword r = r_start; r < r_end; ++r)
        {
        const uword m_row_indices_r = m_row_indices[r];
        
        if( (m_row_indices_r >= lstart_row) && (m_row_indices_r < lend_row) )
          {
          eT& val = m_values[r];
          
          val = (arma_isnan(val)) ? new_val : val;
          }
        }
      }
    }
  else
    {
    for(uword c = lstart_col; c < lend_col; ++c)
      {
      const uword r_start = m.col_ptrs[c    ];
      const uword r_end   = m.col_ptrs[c + 1];
      
      for(uword r = r_start; r < r_end; ++r)
        {
        const uword m_row_indices_r = m_row_indices[r];
        
        if( (m_row_indices_r >= lstart_row) && (m_row_indices_r < lend_row) )
          {
          eT& val = m_values[r];
          
          val = (val == old_val) ? new_val : val;
          }
        }
      }
    }
  
  if(new_val == eT(0))  { access::rw(m).remove_zeros(); }
  }



template<typename eT>
inline
void
SpSubview<eT>::clean(const typename get_pod_type<eT>::result threshold)
  {
  arma_extra_debug_sigprint();
  
  if((n_elem == 0) || (n_nonzero == 0))  { return; }
  
  // TODO: replace with a more efficient implementation
  
  SpMat<eT> tmp(*this);
  
  tmp.clean(threshold);
  
  if(is_cx<eT>::yes)
    {
    (*this).operator=(tmp);
    }
  else
  if(tmp.n_nonzero != n_nonzero)
    {
    (*this).operator=(tmp);
    }
  }



template<typename eT>
inline
void
SpSubview<eT>::clamp(const eT min_val, const eT max_val)
  {
  arma_extra_debug_sigprint();
  
  if(is_cx<eT>::no)
    {
    arma_debug_check( (access::tmp_real(min_val) > access::tmp_real(max_val)), "SpSubview::clamp(): min_val must be less than max_val" );
    }
  else
    {
    arma_debug_check( (access::tmp_real(min_val) > access::tmp_real(max_val)), "SpSubview::clamp(): real(min_val) must be less than real(max_val)" );
    arma_debug_check( (access::tmp_imag(min_val) > access::tmp_imag(max_val)), "SpSubview::clamp(): imag(min_val) must be less than imag(max_val)" );
    }
  
  if((n_elem == 0) || (n_nonzero == 0))  { return; }
  
  // TODO: replace with a more efficient implementation
  
  SpMat<eT> tmp(*this);
  
  tmp.clamp(min_val, max_val);
  
  (*this).operator=(tmp);
  }



template<typename eT>
inline
void
SpSubview<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(val != eT(0))
    {
    Mat<eT> tmp( (*this).n_rows, (*this).n_cols, arma_nozeros_indicator() );
    
    tmp.fill(val);
    
    (*this).operator=(tmp);
    }
  else
    {
    (*this).zeros();
    }
  }



template<typename eT>
inline
void
SpSubview<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  if((n_elem == 0) || (n_nonzero == 0))  { return; }
  
  if((m.n_nonzero - n_nonzero) == 0)
    {
    access::rw(m).zeros();
    access::rw(n_nonzero) = 0;
    return;
    }
  
  SpMat<eT> tmp(arma_reserve_indicator(), m.n_rows, m.n_cols, m.n_nonzero - n_nonzero);
  
  const uword sv_row_start = aux_row1;
  const uword sv_col_start = aux_col1;
  
  const uword sv_row_end   = aux_row1 + n_rows - 1;
  const uword sv_col_end   = aux_col1 + n_cols - 1;
  
  typename SpMat<eT>::const_iterator m_it     = m.begin();
  typename SpMat<eT>::const_iterator m_it_end = m.end();
  
  uword tmp_count = 0;
  
  for(; m_it != m_it_end; ++m_it)
    {
    const uword m_it_row = m_it.row();
    const uword m_it_col = m_it.col();
    
    const bool inside_box = ((m_it_row >= sv_row_start) && (m_it_row <= sv_row_end)) && ((m_it_col >= sv_col_start) && (m_it_col <= sv_col_end));
    
    if(inside_box == false)
      {
      access::rw(tmp.values[tmp_count])      = (*m_it);
      access::rw(tmp.row_indices[tmp_count]) = m_it_row;
      access::rw(tmp.col_ptrs[m_it_col + 1])++;
      ++tmp_count;
      }
    }
  
  for(uword i=0; i < tmp.n_cols; ++i)
    {
    access::rw(tmp.col_ptrs[i + 1]) += tmp.col_ptrs[i];
    }
  
  access::rw(m).steal_mem(tmp);
  
  access::rw(n_nonzero) = 0;
  }



template<typename eT>
inline
void
SpSubview<eT>::ones()
  {
  arma_extra_debug_sigprint();
  
  (*this).fill(eT(1));
  }



template<typename eT>
inline
void
SpSubview<eT>::eye()
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT> tmp;
  
  tmp.eye( (*this).n_rows, (*this).n_cols );
  
  (*this).operator=(tmp);
  }



template<typename eT>
inline
void
SpSubview<eT>::randu()
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp( (*this).n_rows, (*this).n_cols, fill::randu );
  
  (*this).operator=(tmp);
  }



template<typename eT>
inline
void
SpSubview<eT>::randn()
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp( (*this).n_rows, (*this).n_cols, fill::randn );
  
  (*this).operator=(tmp);
  }



template<typename eT>
arma_hot
inline
SpSubview_MapMat_val<eT>
SpSubview<eT>::operator[](const uword i)
  {
  const uword lrow = i % n_rows;
  const uword lcol = i / n_rows;
  
  return (*this).at(lrow, lcol);
  }



template<typename eT>
arma_hot
inline
eT
SpSubview<eT>::operator[](const uword i) const
  {
  const uword lrow = i % n_rows;
  const uword lcol = i / n_rows;
  
  return (*this).at(lrow, lcol);
  }



template<typename eT>
arma_hot
inline
SpSubview_MapMat_val<eT>
SpSubview<eT>::operator()(const uword i)
  {
  arma_debug_check_bounds( (i >= n_elem), "SpSubview::operator(): index out of bounds" );
  
  const uword lrow = i % n_rows;
  const uword lcol = i / n_rows;
  
  return (*this).at(lrow, lcol);
  }



template<typename eT>
arma_hot
inline
eT
SpSubview<eT>::operator()(const uword i) const
  {
  arma_debug_check_bounds( (i >= n_elem), "SpSubview::operator(): index out of bounds" );
  
  const uword lrow = i % n_rows;
  const uword lcol = i / n_rows;
  
  return (*this).at(lrow, lcol);
  }



template<typename eT>
arma_hot
inline
SpSubview_MapMat_val<eT>
SpSubview<eT>::operator()(const uword in_row, const uword in_col)
  {
  arma_debug_check_bounds( (in_row >= n_rows) || (in_col >= n_cols), "SpSubview::operator(): index out of bounds" );
  
  return (*this).at(in_row, in_col);
  }



template<typename eT>
arma_hot
inline
eT
SpSubview<eT>::operator()(const uword in_row, const uword in_col) const
  {
  arma_debug_check_bounds( (in_row >= n_rows) || (in_col >= n_cols), "SpSubview::operator(): index out of bounds" );
  
  return (*this).at(in_row, in_col);
  }



template<typename eT>
arma_hot
inline
SpSubview_MapMat_val<eT>
SpSubview<eT>::at(const uword i)
  {
  const uword lrow = i % n_rows;
  const uword lcol = i / n_cols;
  
  return (*this).at(lrow, lcol);
  }



template<typename eT>
arma_hot
inline
eT
SpSubview<eT>::at(const uword i) const
  {
  const uword lrow = i % n_rows;
  const uword lcol = i / n_cols;
  
  return (*this).at(lrow, lcol);
  }



template<typename eT>
arma_hot
inline
SpSubview_MapMat_val<eT>
SpSubview<eT>::at(const uword in_row, const uword in_col)
  {
  return SpSubview_MapMat_val<eT>((*this), m.cache, aux_row1 + in_row, aux_col1 + in_col);
  }



template<typename eT>
arma_hot
inline
eT
SpSubview<eT>::at(const uword in_row, const uword in_col) const
  {
  return m.at(aux_row1 + in_row, aux_col1 + in_col);
  }



template<typename eT>
inline
bool
SpSubview<eT>::check_overlap(const SpSubview<eT>& x) const
  {
  const SpSubview<eT>& t = *this;
  
  if(&t.m != &x.m)
    {
    return false;
    }
  else
    {
    if( (t.n_elem == 0) || (x.n_elem == 0) )
      {
      return false;
      }
    else
      {
      const uword t_row_start  = t.aux_row1;
      const uword t_row_end_p1 = t_row_start + t.n_rows;
      
      const uword t_col_start  = t.aux_col1;
      const uword t_col_end_p1 = t_col_start + t.n_cols;
      
      const uword x_row_start  = x.aux_row1;
      const uword x_row_end_p1 = x_row_start + x.n_rows;
      
      const uword x_col_start  = x.aux_col1;
      const uword x_col_end_p1 = x_col_start + x.n_cols;
      
      const bool outside_rows = ( (x_row_start >= t_row_end_p1) || (t_row_start >= x_row_end_p1) );
      const bool outside_cols = ( (x_col_start >= t_col_end_p1) || (t_col_start >= x_col_end_p1) );
      
      return ( (outside_rows == false) && (outside_cols == false) );
      }
    }
  }



template<typename eT>
inline
bool
SpSubview<eT>::is_vec() const
  {
  return ( (n_rows == 1) || (n_cols == 1) );
  }



template<typename eT>
inline
SpSubview_row<eT>
SpSubview<eT>::row(const uword row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds(row_num >= n_rows, "SpSubview::row(): out of bounds");
  
  return SpSubview_row<eT>(const_cast< SpMat<eT>& >(m), row_num + aux_row1, aux_col1, n_cols);
  }



template<typename eT>
inline
const SpSubview_row<eT>
SpSubview<eT>::row(const uword row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds(row_num >= n_rows, "SpSubview::row(): out of bounds");
  
  return SpSubview_row<eT>(m, row_num + aux_row1, aux_col1, n_cols);
  }



template<typename eT>
inline
SpSubview_col<eT>
SpSubview<eT>::col(const uword col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds(col_num >= n_cols, "SpSubview::col(): out of bounds");
  
  return SpSubview_col<eT>(const_cast< SpMat<eT>& >(m), col_num + aux_col1, aux_row1, n_rows);
  }



template<typename eT>
inline
const SpSubview_col<eT>
SpSubview<eT>::col(const uword col_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds(col_num >= n_cols, "SpSubview::col(): out of bounds");
  
  return SpSubview_col<eT>(m, col_num + aux_col1, aux_row1, n_rows);
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "SpSubview::rows(): indices out of bounds or incorrectly used"
    );
  
  return submat(in_row1, 0, in_row2, n_cols - 1);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::rows(const uword in_row1, const uword in_row2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "SpSubview::rows(): indices out of bounds or incorrectly used"
    );

  return submat(in_row1, 0, in_row2, n_cols - 1);
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "SpSubview::cols(): indices out of bounds or incorrectly used"
    );
  
  return submat(0, in_col1, n_rows - 1, in_col2);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::cols(const uword in_col1, const uword in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "SpSubview::cols(): indices out of bounds or incorrectly used"
    );
  
  return submat(0, in_col1, n_rows - 1, in_col2);
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_row1 > in_row2) || (in_col1 > in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "SpSubview::submat(): indices out of bounds or incorrectly used"
    );
  
  return access::rw(m).submat(in_row1 + aux_row1, in_col1 + aux_col1, in_row2 + aux_row1, in_col2 + aux_col1);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_row1 > in_row2) || (in_col1 > in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "SpSubview::submat(): indices out of bounds or incorrectly used"
    );
  
  return m.submat(in_row1 + aux_row1, in_col1 + aux_col1, in_row2 + aux_row1, in_col2 + aux_col1);
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::submat(const span& row_span, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;
  const bool col_all = row_span.whole;
  
  const uword in_row1 = row_all ? 0      : row_span.a;
  const uword in_row2 = row_all ? n_rows : row_span.b;
  
  const uword in_col1 = col_all ? 0      : col_span.a;
  const uword in_col2 = col_all ? n_cols : col_span.b;
  
  arma_debug_check_bounds
    (
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= n_rows)))
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= n_cols))),
    "SpSubview::submat(): indices out of bounds or incorrectly used"
    );
  
  return submat(in_row1, in_col1, in_row2, in_col2);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::submat(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;
  const bool col_all = row_span.whole;
  
  const uword in_row1 = row_all ? 0          : row_span.a;
  const uword in_row2 = row_all ? n_rows - 1 : row_span.b;
  
  const uword in_col1 = col_all ? 0          : col_span.a;
  const uword in_col2 = col_all ? n_cols - 1 : col_span.b;
  
  arma_debug_check_bounds
    (
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= n_rows)))
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= n_cols))),
    "SpSubview::submat(): indices out of bounds or incorrectly used"
    );
  
  return submat(in_row1, in_col1, in_row2, in_col2);
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::operator()(const uword row_num, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  return submat(span(row_num, row_num), col_span);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::operator()(const uword row_num, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  return submat(span(row_num, row_num), col_span);
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::operator()(const span& row_span, const uword col_num)
  {
  arma_extra_debug_sigprint();
  
  return submat(row_span, span(col_num, col_num));
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::operator()(const span& row_span, const uword col_num) const
  {
  arma_extra_debug_sigprint();
  
  return submat(row_span, span(col_num, col_num));
  }



template<typename eT>
inline
SpSubview<eT>
SpSubview<eT>::operator()(const span& row_span, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  return submat(row_span, col_span);
  }



template<typename eT>
inline
const SpSubview<eT>
SpSubview<eT>::operator()(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  return submat(row_span, col_span);
  }



template<typename eT>
inline
void
SpSubview<eT>::swap_rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check((in_row1 >= n_rows) || (in_row2 >= n_rows), "SpSubview::swap_rows(): invalid row index");
  
  const uword lstart_col = aux_col1;
  const uword lend_col   = aux_col1 + n_cols;
  
  for(uword c = lstart_col; c < lend_col; ++c)
    {
    const eT val = access::rw(m).at(in_row1 + aux_row1, c);
    access::rw(m).at(in_row2 + aux_row1, c) = eT( access::rw(m).at(in_row1 + aux_row1, c) );
    access::rw(m).at(in_row1 + aux_row1, c) = val;
    }
  }



template<typename eT>
inline
void
SpSubview<eT>::swap_cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check((in_col1 >= n_cols) || (in_col2 >= n_cols), "SpSubview::swap_cols(): invalid column index");
  
  const uword lstart_row = aux_row1;
  const uword lend_row   = aux_row1 + n_rows;
  
  for(uword r = lstart_row; r < lend_row; ++r)
    {
    const eT val = access::rw(m).at(r, in_col1 + aux_col1);
    access::rw(m).at(r, in_col1 + aux_col1) = eT( access::rw(m).at(r, in_col2 + aux_col1) );
    access::rw(m).at(r, in_col2 + aux_col1) = val;
    }
  }



template<typename eT>
inline
typename SpSubview<eT>::iterator
SpSubview<eT>::begin()
  {
  return iterator(*this);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_iterator
SpSubview<eT>::begin() const
  {
  m.sync_csc();
  
  return const_iterator(*this);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_iterator
SpSubview<eT>::cbegin() const
  {
  m.sync_csc();
  
  return const_iterator(*this);
  }



template<typename eT>
inline
typename SpSubview<eT>::iterator
SpSubview<eT>::begin_col(const uword col_num)
  {
  m.sync_csc();
  
  return iterator(*this, 0, col_num);
  }


template<typename eT>
inline
typename SpSubview<eT>::const_iterator
SpSubview<eT>::begin_col(const uword col_num) const
  {
  m.sync_csc();
  
  return const_iterator(*this, 0, col_num);
  }



template<typename eT>
inline
typename SpSubview<eT>::row_iterator
SpSubview<eT>::begin_row(const uword row_num)
  {
  m.sync_csc();
  
  return row_iterator(*this, row_num, 0);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_row_iterator
SpSubview<eT>::begin_row(const uword row_num) const
  {
  m.sync_csc();
  
  return const_row_iterator(*this, row_num, 0);
  }



template<typename eT>
inline
typename SpSubview<eT>::iterator
SpSubview<eT>::end()
  {
  m.sync_csc();
  
  return iterator(*this, 0, n_cols, n_nonzero, m.n_nonzero - n_nonzero);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_iterator
SpSubview<eT>::end() const
  {
  m.sync_csc();
  
  return const_iterator(*this, 0, n_cols, n_nonzero, m.n_nonzero - n_nonzero);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_iterator
SpSubview<eT>::cend() const
  {
  m.sync_csc();
  
  return const_iterator(*this, 0, n_cols, n_nonzero, m.n_nonzero - n_nonzero);
  }



template<typename eT>
inline
typename SpSubview<eT>::row_iterator
SpSubview<eT>::end_row()
  {
  m.sync_csc();
  
  return row_iterator(*this, n_nonzero);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_row_iterator
SpSubview<eT>::end_row() const
  {
  m.sync_csc();
  
  return const_row_iterator(*this, n_nonzero);
  }



template<typename eT>
inline
typename SpSubview<eT>::row_iterator
SpSubview<eT>::end_row(const uword row_num)
  {
  m.sync_csc();
  
  return row_iterator(*this, row_num + 1, 0);
  }



template<typename eT>
inline
typename SpSubview<eT>::const_row_iterator
SpSubview<eT>::end_row(const uword row_num) const
  {
  m.sync_csc();
  
  return const_row_iterator(*this, row_num + 1, 0);
  }



template<typename eT>
arma_inline
bool
SpSubview<eT>::is_alias(const SpMat<eT>& X) const
  {
  return m.is_alias(X);
  }



template<typename eT>
inline
arma_warn_unused
eT&
SpSubview<eT>::insert_element(const uword in_row, const uword in_col, const eT in_val)
  {
  arma_extra_debug_sigprint();
  
  // This may not actually insert an element.
  const uword old_n_nonzero = m.n_nonzero;
  eT& retval = access::rw(m).insert_element(in_row + aux_row1, in_col + aux_col1, in_val);
  // Update n_nonzero (if necessary).
  access::rw(n_nonzero) += (m.n_nonzero - old_n_nonzero);
  
  return retval;
  }



template<typename eT>
inline
void
SpSubview<eT>::delete_element(const uword in_row, const uword in_col)
  {
  arma_extra_debug_sigprint();
  
  // This may not actually delete an element.
  const uword old_n_nonzero = m.n_nonzero;
  access::rw(m).delete_element(in_row + aux_row1, in_col + aux_col1);
  access::rw(n_nonzero) -= (old_n_nonzero - m.n_nonzero);
  }



template<typename eT>
inline
void
SpSubview<eT>::invalidate_cache() const
  {
  arma_extra_debug_sigprint();
  
  m.invalidate_cache();
  }



//
//
//



template<typename eT>
inline
SpSubview_col<eT>::SpSubview_col(const SpMat<eT>& in_m, const uword in_col)
  : SpSubview<eT>(in_m, 0, in_col, in_m.n_rows, 1)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
SpSubview_col<eT>::SpSubview_col(const SpMat<eT>& in_m, const uword in_col, const uword in_row1, const uword in_n_rows)
  : SpSubview<eT>(in_m, in_row1, in_col, in_n_rows, 1)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
void
SpSubview_col<eT>::operator=(const SpSubview<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  SpSubview<eT>::operator=(x);
  }



template<typename eT>
inline
void
SpSubview_col<eT>::operator=(const SpSubview_col<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  SpSubview<eT>::operator=(x); // interprets 'SpSubview_col' as 'SpSubview'
  }



template<typename eT>
template<typename T1>
inline
void
SpSubview_col<eT>::operator=(const SpBase<eT,T1>& x)
  {
  arma_extra_debug_sigprint();
  
  SpSubview<eT>::operator=(x);
  }



template<typename eT>
template<typename T1>
inline
void
SpSubview_col<eT>::operator=(const Base<eT,T1>& x)
  {
  arma_extra_debug_sigprint();
  
  SpSubview<eT>::operator=(x);
  }



template<typename eT>
inline
arma_warn_unused
const SpOp<SpSubview_col<eT>,spop_htrans>
SpSubview_col<eT>::t() const
  {
  return SpOp<SpSubview_col<eT>,spop_htrans>(*this);
  }



template<typename eT>
inline
arma_warn_unused
const SpOp<SpSubview_col<eT>,spop_htrans>
SpSubview_col<eT>::ht() const
  {
  return SpOp<SpSubview_col<eT>,spop_htrans>(*this);
  }



template<typename eT>
inline
arma_warn_unused
const SpOp<SpSubview_col<eT>,spop_strans>
SpSubview_col<eT>::st() const
  {
  return SpOp<SpSubview_col<eT>,spop_strans>(*this);
  }



//
//
//



template<typename eT>
inline
SpSubview_row<eT>::SpSubview_row(const SpMat<eT>& in_m, const uword in_row)
  : SpSubview<eT>(in_m, in_row, 0, 1, in_m.n_cols)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
SpSubview_row<eT>::SpSubview_row(const SpMat<eT>& in_m, const uword in_row, const uword in_col1, const uword in_n_cols)
  : SpSubview<eT>(in_m, in_row, in_col1, 1, in_n_cols)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
void
SpSubview_row<eT>::operator=(const SpSubview<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  SpSubview<eT>::operator=(x);
  }



template<typename eT>
inline
void
SpSubview_row<eT>::operator=(const SpSubview_row<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  SpSubview<eT>::operator=(x); // interprets 'SpSubview_row' as 'SpSubview'
  }



template<typename eT>
template<typename T1>
inline
void
SpSubview_row<eT>::operator=(const SpBase<eT,T1>& x)
  {
  arma_extra_debug_sigprint();
  
  SpSubview<eT>::operator=(x);
  }



template<typename eT>
template<typename T1>
inline
void
SpSubview_row<eT>::operator=(const Base<eT,T1>& x)
  {
  arma_extra_debug_sigprint();
  
  SpSubview<eT>::operator=(x);
  }



template<typename eT>
inline
arma_warn_unused
const SpOp<SpSubview_row<eT>,spop_htrans>
SpSubview_row<eT>::t() const
  {
  return SpOp<SpSubview_row<eT>,spop_htrans>(*this);
  }



template<typename eT>
inline
arma_warn_unused
const SpOp<SpSubview_row<eT>,spop_htrans>
SpSubview_row<eT>::ht() const
  {
  return SpOp<SpSubview_row<eT>,spop_htrans>(*this);
  }



template<typename eT>
inline
arma_warn_unused
const SpOp<SpSubview_row<eT>,spop_strans>
SpSubview_row<eT>::st() const
  {
  return SpOp<SpSubview_row<eT>,spop_strans>(*this);
  }



//! @}
