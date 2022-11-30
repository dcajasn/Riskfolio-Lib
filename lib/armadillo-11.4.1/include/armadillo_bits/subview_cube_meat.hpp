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


//! \addtogroup subview_cube
//! @{


template<typename eT>
inline
subview_cube<eT>::~subview_cube()
  {
  arma_extra_debug_sigprint_this(this);
  }



template<typename eT>
arma_inline
subview_cube<eT>::subview_cube
  (
  const Cube<eT>& in_m,
  const uword     in_row1,
  const uword     in_col1,
  const uword     in_slice1,
  const uword     in_n_rows,
  const uword     in_n_cols,
  const uword     in_n_slices
  )
  : m           (in_m)
  , aux_row1    (in_row1)
  , aux_col1    (in_col1)
  , aux_slice1  (in_slice1)
  , n_rows      (in_n_rows)
  , n_cols      (in_n_cols)
  , n_elem_slice(in_n_rows * in_n_cols)
  , n_slices    (in_n_slices)
  , n_elem      (n_elem_slice * in_n_slices)
  {
  arma_extra_debug_sigprint_this(this);
  }



template<typename eT>
inline
subview_cube<eT>::subview_cube(const subview_cube<eT>& in)
  : m           (in.m           )
  , aux_row1    (in.aux_row1    )
  , aux_col1    (in.aux_col1    )
  , aux_slice1  (in.aux_slice1  )
  , n_rows      (in.n_rows      )
  , n_cols      (in.n_cols      )
  , n_elem_slice(in.n_elem_slice)
  , n_slices    (in.n_slices    )
  , n_elem      (in.n_elem      )
  {
  arma_extra_debug_sigprint(arma_str::format("this = %x   in = %x") % this % &in);
  }



template<typename eT>
inline
subview_cube<eT>::subview_cube(subview_cube<eT>&& in)
  : m           (in.m           )
  , aux_row1    (in.aux_row1    )
  , aux_col1    (in.aux_col1    )
  , aux_slice1  (in.aux_slice1  )
  , n_rows      (in.n_rows      )
  , n_cols      (in.n_cols      )
  , n_elem_slice(in.n_elem_slice)
  , n_slices    (in.n_slices    )
  , n_elem      (in.n_elem      )
  {
  arma_extra_debug_sigprint(arma_str::format("this = %x   in = %x") % this % &in);
  
  // for paranoia
  
  access::rw(in.aux_row1    ) = 0;
  access::rw(in.aux_col1    ) = 0;
  access::rw(in.aux_slice1  ) = 0;
  access::rw(in.n_rows      ) = 0;
  access::rw(in.n_cols      ) = 0;
  access::rw(in.n_elem_slice) = 0;
  access::rw(in.n_slices    ) = 0;
  access::rw(in.n_elem      ) = 0;
  }



template<typename eT>
template<typename op_type>
inline
void
subview_cube<eT>::inplace_op(const eT val)
  {
  arma_extra_debug_sigprint();
  
  subview_cube<eT>& t = *this;
  
  const uword t_n_rows   = t.n_rows;
  const uword t_n_cols   = t.n_cols;
  const uword t_n_slices = t.n_slices;
  
  for(uword s=0; s < t_n_slices; ++s)
  for(uword c=0; c < t_n_cols;   ++c)
    {
    if(is_same_type<op_type, op_internal_plus >::yes)  { arrayops::inplace_plus ( slice_colptr(s,c), val, t_n_rows ); }
    if(is_same_type<op_type, op_internal_minus>::yes)  { arrayops::inplace_minus( slice_colptr(s,c), val, t_n_rows ); }
    if(is_same_type<op_type, op_internal_schur>::yes)  { arrayops::inplace_mul  ( slice_colptr(s,c), val, t_n_rows ); }
    if(is_same_type<op_type, op_internal_div  >::yes)  { arrayops::inplace_div  ( slice_colptr(s,c), val, t_n_rows ); }
    }
  }






template<typename eT>
template<typename op_type, typename T1>
inline
void
subview_cube<eT>::inplace_op(const BaseCube<eT,T1>& in, const char* identifier)
  {
  arma_extra_debug_sigprint();
  
  const ProxyCube<T1> P(in.get_ref());
  
  subview_cube<eT>& t = *this;
  
  const uword t_n_rows   = t.n_rows;
  const uword t_n_cols   = t.n_cols;
  const uword t_n_slices = t.n_slices;
  
  arma_debug_assert_same_size(t, P, identifier);
  
  const bool use_mp      = arma_config::openmp && ProxyCube<T1>::use_mp && mp_gate<eT>::eval(t.n_elem);
  const bool has_overlap = P.has_overlap(t);
  
  if(has_overlap)  { arma_extra_debug_print("aliasing or overlap detected"); }
  
  if( (is_Cube<typename ProxyCube<T1>::stored_type>::value) || (use_mp) || (has_overlap) )
    {
    const unwrap_cube_check<typename ProxyCube<T1>::stored_type> tmp(P.Q, has_overlap);
    const Cube<eT>& B = tmp.M;
    
    if( (is_same_type<op_type, op_internal_equ>::yes) && (t.aux_row1 == 0) && (t_n_rows == t.m.n_rows) )
      {
      for(uword s=0; s < t_n_slices; ++s)
        {
        arrayops::copy( t.slice_colptr(s,0), B.slice_colptr(s,0), t.n_elem_slice );
        }
      }
    else
      {
      for(uword s=0; s < t_n_slices; ++s)
      for(uword c=0; c < t_n_cols;   ++c)
        {
        if(is_same_type<op_type, op_internal_equ  >::yes)  { arrayops::copy         ( t.slice_colptr(s,c), B.slice_colptr(s,c), t_n_rows ); }
        if(is_same_type<op_type, op_internal_plus >::yes)  { arrayops::inplace_plus ( t.slice_colptr(s,c), B.slice_colptr(s,c), t_n_rows ); }
        if(is_same_type<op_type, op_internal_minus>::yes)  { arrayops::inplace_minus( t.slice_colptr(s,c), B.slice_colptr(s,c), t_n_rows ); }
        if(is_same_type<op_type, op_internal_schur>::yes)  { arrayops::inplace_mul  ( t.slice_colptr(s,c), B.slice_colptr(s,c), t_n_rows ); }
        if(is_same_type<op_type, op_internal_div  >::yes)  { arrayops::inplace_div  ( t.slice_colptr(s,c), B.slice_colptr(s,c), t_n_rows ); }
        }
      }
    }
  else  // use the Proxy
    {
    if(ProxyCube<T1>::use_at)
      {
      for(uword s=0; s < t_n_slices; ++s)
      for(uword c=0; c < t_n_cols;   ++c)
        {
        eT* t_col_data = t.slice_colptr(s,c);
        
        for(uword r=0; r < t_n_rows; ++r)
          {
          const eT tmp = P.at(r,c,s);
          
          if(is_same_type<op_type, op_internal_equ  >::yes)  { (*t_col_data) =  tmp; t_col_data++; }
          if(is_same_type<op_type, op_internal_plus >::yes)  { (*t_col_data) += tmp; t_col_data++; }
          if(is_same_type<op_type, op_internal_minus>::yes)  { (*t_col_data) -= tmp; t_col_data++; }
          if(is_same_type<op_type, op_internal_schur>::yes)  { (*t_col_data) *= tmp; t_col_data++; }
          if(is_same_type<op_type, op_internal_div  >::yes)  { (*t_col_data) /= tmp; t_col_data++; }
          }
        }
      }
    else
      {
      typename ProxyCube<T1>::ea_type Pea = P.get_ea();
      
      uword count = 0;
      
      for(uword s=0; s < t_n_slices; ++s)
      for(uword c=0; c < t_n_cols;   ++c)
        {
        eT* t_col_data = t.slice_colptr(s,c);
        
        for(uword r=0; r < t_n_rows; ++r)
          {
          const eT tmp = Pea[count];  count++;
          
          if(is_same_type<op_type, op_internal_equ  >::yes)  { (*t_col_data) =  tmp; t_col_data++; }
          if(is_same_type<op_type, op_internal_plus >::yes)  { (*t_col_data) += tmp; t_col_data++; }
          if(is_same_type<op_type, op_internal_minus>::yes)  { (*t_col_data) -= tmp; t_col_data++; }
          if(is_same_type<op_type, op_internal_schur>::yes)  { (*t_col_data) *= tmp; t_col_data++; }
          if(is_same_type<op_type, op_internal_div  >::yes)  { (*t_col_data) /= tmp; t_col_data++; }
          }
        }
      }
    }
  }



template<typename eT>
template<typename op_type>
inline
void
subview_cube<eT>::inplace_op(const subview_cube<eT>& x, const char* identifier)
  {
  arma_extra_debug_sigprint();
  
  if(check_overlap(x))
    {
    const Cube<eT> tmp(x);
    
    if(is_same_type<op_type, op_internal_equ  >::yes)  { (*this).operator= (tmp); }
    if(is_same_type<op_type, op_internal_plus >::yes)  { (*this).operator+=(tmp); }
    if(is_same_type<op_type, op_internal_minus>::yes)  { (*this).operator-=(tmp); }
    if(is_same_type<op_type, op_internal_schur>::yes)  { (*this).operator%=(tmp); }
    if(is_same_type<op_type, op_internal_div  >::yes)  { (*this).operator/=(tmp); }
    
    return;
    }
  
  subview_cube<eT>& t = *this;
  
  arma_debug_assert_same_size(t, x, identifier);
  
  const uword t_n_rows   = t.n_rows;
  const uword t_n_cols   = t.n_cols;
  const uword t_n_slices = t.n_slices;
  
  for(uword s=0; s < t_n_slices; ++s)
  for(uword c=0; c < t_n_cols;   ++c)
    {
    if(is_same_type<op_type, op_internal_equ  >::yes)  { arrayops::copy         ( t.slice_colptr(s,c), x.slice_colptr(s,c), t_n_rows ); }
    if(is_same_type<op_type, op_internal_plus >::yes)  { arrayops::inplace_plus ( t.slice_colptr(s,c), x.slice_colptr(s,c), t_n_rows ); }
    if(is_same_type<op_type, op_internal_minus>::yes)  { arrayops::inplace_minus( t.slice_colptr(s,c), x.slice_colptr(s,c), t_n_rows ); }
    if(is_same_type<op_type, op_internal_schur>::yes)  { arrayops::inplace_mul  ( t.slice_colptr(s,c), x.slice_colptr(s,c), t_n_rows ); }
    if(is_same_type<op_type, op_internal_div  >::yes)  { arrayops::inplace_div  ( t.slice_colptr(s,c), x.slice_colptr(s,c), t_n_rows ); }
    }
  }



template<typename eT>
inline
void
subview_cube<eT>::operator= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(n_elem != 1)
    {
    arma_debug_assert_same_size(n_rows, n_cols, n_slices, 1, 1, 1, "copy into subcube");
    }
  
  Cube<eT>& Q = const_cast< Cube<eT>& >(m);
  
  Q.at(aux_row1, aux_col1, aux_slice1) = val;
  }



template<typename eT>
inline
void
subview_cube<eT>::operator+= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_plus>(val);
  }



template<typename eT>
inline
void
subview_cube<eT>::operator-= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_minus>(val);
  }



template<typename eT>
inline
void
subview_cube<eT>::operator*= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_schur>(val);
  }



template<typename eT>
inline
void
subview_cube<eT>::operator/= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_div>(val);
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator= (const BaseCube<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_equ>(in, "copy into subcube");
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator+= (const BaseCube<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_plus>(in, "addition");
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator-= (const BaseCube<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_minus>(in, "subtraction");
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator%= (const BaseCube<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_schur>(in, "element-wise multiplication");
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator/= (const BaseCube<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_div>(in, "element-wise division");
  }



//! x.subcube(...) = y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::operator= (const subview_cube<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_equ>(x, "copy into subcube");
  }



template<typename eT>
inline
void
subview_cube<eT>::operator+= (const subview_cube<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_plus>(x, "addition");
  }



template<typename eT>
inline
void
subview_cube<eT>::operator-= (const subview_cube<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_minus>(x, "subtraction");
  }



template<typename eT>
inline
void
subview_cube<eT>::operator%= (const subview_cube<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_schur>(x, "element-wise multiplication");
  }



template<typename eT>
inline
void
subview_cube<eT>::operator/= (const subview_cube<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_div>(x, "element-wise division");
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<T1> tmp(in.get_ref());
  
  const Mat<eT>&          x = tmp.M;
        subview_cube<eT>& t = *this;
  
  const uword t_n_rows   = t.n_rows;
  const uword t_n_cols   = t.n_cols;
  const uword t_n_slices = t.n_slices;
  
  const uword x_n_rows   = x.n_rows;
  const uword x_n_cols   = x.n_cols;
  
  if( ((x_n_rows == 1) || (x_n_cols == 1)) && (t_n_rows == 1) && (t_n_cols == 1) && (x.n_elem == t_n_slices) )
    {
    Cube<eT>& Q = const_cast< Cube<eT>& >(t.m);
    
    const uword t_aux_row1   = t.aux_row1;
    const uword t_aux_col1   = t.aux_col1;
    const uword t_aux_slice1 = t.aux_slice1;
    
    const eT* x_mem = x.memptr();
    
    uword i,j;
    for(i=0, j=1; j < t_n_slices; i+=2, j+=2)
      {
      const eT tmp_i = x_mem[i];
      const eT tmp_j = x_mem[j];
      
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + i) = tmp_i;
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + j) = tmp_j;
      }
    
    if(i < t_n_slices)
      {
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + i) = x_mem[i];
      }
    }
  else
  if( (t_n_rows == x_n_rows) && (t_n_cols == x_n_cols) && (t_n_slices == 1) )
    {
    // interpret the matrix as a cube with one slice
    
    for(uword col = 0; col < t_n_cols; ++col)
      {
      arrayops::copy( t.slice_colptr(0, col), x.colptr(col), t_n_rows );
      }
    }
  else
  if( (t_n_rows == x_n_rows) && (t_n_cols == 1) && (t_n_slices == x_n_cols) )
    {
    for(uword i=0; i < t_n_slices; ++i)
      {
      arrayops::copy( t.slice_colptr(i, 0), x.colptr(i), t_n_rows );
      }
    }
  else
  if( (t_n_rows == 1) && (t_n_cols == x_n_rows) && (t_n_slices == x_n_cols) )
    {
    Cube<eT>& Q = const_cast< Cube<eT>& >(t.m);
    
    const uword t_aux_row1   = t.aux_row1;
    const uword t_aux_col1   = t.aux_col1;
    const uword t_aux_slice1 = t.aux_slice1;
    
    for(uword slice=0; slice < t_n_slices; ++slice)
      {
      const uword mod_slice = t_aux_slice1 + slice;
      
      const eT* x_colptr = x.colptr(slice);
      
      uword i,j;
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        const eT tmp_i = x_colptr[i];
        const eT tmp_j = x_colptr[j];
        
        Q.at(t_aux_row1, t_aux_col1 + i, mod_slice) = tmp_i;
        Q.at(t_aux_row1, t_aux_col1 + j, mod_slice) = tmp_j;
        }
      
      if(i < t_n_cols)
        {
        Q.at(t_aux_row1, t_aux_col1 + i, mod_slice) = x_colptr[i];
        }
      }
    }
  else
    {
    if(arma_config::debug)
      {
      arma_stop_logic_error( arma_incompat_size_string(t, x, "copy into subcube") );
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator+= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<T1> tmp(in.get_ref());
  
  const Mat<eT>&          x = tmp.M;
        subview_cube<eT>& t = *this;
  
  const uword t_n_rows   = t.n_rows;
  const uword t_n_cols   = t.n_cols;
  const uword t_n_slices = t.n_slices;
  
  const uword x_n_rows   = x.n_rows;
  const uword x_n_cols   = x.n_cols;
  
  if( ((x_n_rows == 1) || (x_n_cols == 1)) && (t_n_rows == 1) && (t_n_cols == 1) && (x.n_elem == t_n_slices) )
    {
    Cube<eT>& Q = const_cast< Cube<eT>& >(t.m);
    
    const uword t_aux_row1   = t.aux_row1;
    const uword t_aux_col1   = t.aux_col1;
    const uword t_aux_slice1 = t.aux_slice1;
    
    const eT* x_mem = x.memptr();
    
    uword i,j;
    for(i=0, j=1; j < t_n_slices; i+=2, j+=2)
      {
      const eT tmp_i = x_mem[i];
      const eT tmp_j = x_mem[j];
      
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + i) += tmp_i;
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + j) += tmp_j;
      }
    
    if(i < t_n_slices)
      {
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + i) += x_mem[i];
      }
    }
  else
  if( (t_n_rows == x_n_rows) && (t_n_cols == x_n_cols) && (t_n_slices == 1) )
    {
    for(uword col = 0; col < t_n_cols; ++col)
      {
      arrayops::inplace_plus( t.slice_colptr(0, col), x.colptr(col), t_n_rows );
      }
    }
  else
  if( (t_n_rows == x_n_rows) && (t_n_cols == 1) && (t_n_slices == x_n_cols) )
    {
    for(uword i=0; i < t_n_slices; ++i)
      {
      arrayops::inplace_plus( t.slice_colptr(i, 0), x.colptr(i), t_n_rows );
      }
    }
  else
  if( (t_n_rows == 1) && (t_n_cols == x_n_rows) && (t_n_slices == x_n_cols) )
    {
    Cube<eT>& Q = const_cast< Cube<eT>& >(t.m);
    
    const uword t_aux_row1   = t.aux_row1;
    const uword t_aux_col1   = t.aux_col1;
    const uword t_aux_slice1 = t.aux_slice1;
    
    for(uword slice=0; slice < t_n_slices; ++slice)
      {
      const uword mod_slice = t_aux_slice1 + slice;
      
      const eT* x_colptr = x.colptr(slice);
      
      uword i,j;
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        const eT tmp_i = x_colptr[i];
        const eT tmp_j = x_colptr[j];
        
        Q.at(t_aux_row1, t_aux_col1 + i, mod_slice) += tmp_i;
        Q.at(t_aux_row1, t_aux_col1 + j, mod_slice) += tmp_j;
        }
      
      if(i < t_n_cols)
        {
        Q.at(t_aux_row1, t_aux_col1 + i, mod_slice) += x_colptr[i];
        }
      }
    }
  else
    {
    if(arma_config::debug)
      {
      arma_stop_logic_error( arma_incompat_size_string(t, x, "addition") );
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator-= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<T1> tmp(in.get_ref());
  
  const Mat<eT>&          x = tmp.M;
        subview_cube<eT>& t = *this;
  
  const uword t_n_rows   = t.n_rows;
  const uword t_n_cols   = t.n_cols;
  const uword t_n_slices = t.n_slices;
  
  const uword x_n_rows   = x.n_rows;
  const uword x_n_cols   = x.n_cols;
  
  if( ((x_n_rows == 1) || (x_n_cols == 1)) && (t_n_rows == 1) && (t_n_cols == 1) && (x.n_elem == t_n_slices) )
    {
    Cube<eT>& Q = const_cast< Cube<eT>& >(t.m);
    
    const uword t_aux_row1   = t.aux_row1;
    const uword t_aux_col1   = t.aux_col1;
    const uword t_aux_slice1 = t.aux_slice1;
    
    const eT* x_mem = x.memptr();
    
    uword i,j;
    for(i=0, j=1; j < t_n_slices; i+=2, j+=2)
      {
      const eT tmp_i = x_mem[i];
      const eT tmp_j = x_mem[j];
      
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + i) -= tmp_i;
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + j) -= tmp_j;
      }
    
    if(i < t_n_slices)
      {
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + i) -= x_mem[i];
      }
    }
  else
  if( (t_n_rows == x_n_rows) && (t_n_cols == x_n_cols) && (t_n_slices == 1) )
    {
    for(uword col = 0; col < t_n_cols; ++col)
      {
      arrayops::inplace_minus( t.slice_colptr(0, col), x.colptr(col), t_n_rows );
      }
    }
  else
  if( (t_n_rows == x_n_rows) && (t_n_cols == 1) && (t_n_slices == x_n_cols) )
    {
    for(uword i=0; i < t_n_slices; ++i)
      {
      arrayops::inplace_minus( t.slice_colptr(i, 0), x.colptr(i), t_n_rows );
      }
    }
  else
  if( (t_n_rows == 1) && (t_n_cols == x_n_rows) && (t_n_slices == x_n_cols) )
    {
    Cube<eT>& Q = const_cast< Cube<eT>& >(t.m);
    
    const uword t_aux_row1   = t.aux_row1;
    const uword t_aux_col1   = t.aux_col1;
    const uword t_aux_slice1 = t.aux_slice1;
    
    for(uword slice=0; slice < t_n_slices; ++slice)
      {
      const uword mod_slice = t_aux_slice1 + slice;
      
      const eT* x_colptr = x.colptr(slice);
      
      uword i,j;
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        const eT tmp_i = x_colptr[i];
        const eT tmp_j = x_colptr[j];
        
        Q.at(t_aux_row1, t_aux_col1 + i, mod_slice) -= tmp_i;
        Q.at(t_aux_row1, t_aux_col1 + j, mod_slice) -= tmp_j;
        }
      
      if(i < t_n_cols)
        {
        Q.at(t_aux_row1, t_aux_col1 + i, mod_slice) -= x_colptr[i];
        }
      }
    }
  else
    {
    if(arma_config::debug)
      {
      arma_stop_logic_error( arma_incompat_size_string(t, x, "subtraction") );
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator%= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<T1> tmp(in.get_ref());
  
  const Mat<eT>&          x = tmp.M;
        subview_cube<eT>& t = *this;
  
  const uword t_n_rows   = t.n_rows;
  const uword t_n_cols   = t.n_cols;
  const uword t_n_slices = t.n_slices;
  
  const uword x_n_rows   = x.n_rows;
  const uword x_n_cols   = x.n_cols;
  
  if( ((x_n_rows == 1) || (x_n_cols == 1)) && (t_n_rows == 1) && (t_n_cols == 1) && (x.n_elem == t_n_slices) )
    {
    Cube<eT>& Q = const_cast< Cube<eT>& >(t.m);
    
    const uword t_aux_row1   = t.aux_row1;
    const uword t_aux_col1   = t.aux_col1;
    const uword t_aux_slice1 = t.aux_slice1;
    
    const eT* x_mem = x.memptr();
    
    uword i,j;
    for(i=0, j=1; j < t_n_slices; i+=2, j+=2)
      {
      const eT tmp_i = x_mem[i];
      const eT tmp_j = x_mem[j];
      
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + i) *= tmp_i;
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + j) *= tmp_j;
      }
    
    if(i < t_n_slices)
      {
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + i) *= x_mem[i];
      }
    }
  else
  if( (t_n_rows == x_n_rows) && (t_n_cols == x_n_cols) && (t_n_slices == 1) )
    {
    for(uword col = 0; col < t_n_cols; ++col)
      {
      arrayops::inplace_mul( t.slice_colptr(0, col), x.colptr(col), t_n_rows );
      }
    }
  else
  if( (t_n_rows == x_n_rows) && (t_n_cols == 1) && (t_n_slices == x_n_cols) )
    {
    for(uword i=0; i < t_n_slices; ++i)
      {
      arrayops::inplace_mul( t.slice_colptr(i, 0), x.colptr(i), t_n_rows );
      }
    }
  else
  if( (t_n_rows == 1) && (t_n_cols == x_n_rows) && (t_n_slices == x_n_cols) )
    {
    Cube<eT>& Q = const_cast< Cube<eT>& >(t.m);
    
    const uword t_aux_row1   = t.aux_row1;
    const uword t_aux_col1   = t.aux_col1;
    const uword t_aux_slice1 = t.aux_slice1;
    
    for(uword slice=0; slice < t_n_slices; ++slice)
      {
      const uword mod_slice = t_aux_slice1 + slice;
      
      const eT* x_colptr = x.colptr(slice);
      
      uword i,j;
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        const eT tmp_i = x_colptr[i];
        const eT tmp_j = x_colptr[j];
        
        Q.at(t_aux_row1, t_aux_col1 + i, mod_slice) *= tmp_i;
        Q.at(t_aux_row1, t_aux_col1 + j, mod_slice) *= tmp_j;
        }
      
      if(i < t_n_cols)
        {
        Q.at(t_aux_row1, t_aux_col1 + i, mod_slice) *= x_colptr[i];
        }
      }
    }
  else
    {
    if(arma_config::debug)
      {
      arma_stop_logic_error( arma_incompat_size_string(t, x, "element-wise multiplication") );
      }
    }
  }



template<typename eT>
template<typename T1>
inline
void
subview_cube<eT>::operator/= (const Base<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<T1> tmp(in.get_ref());
  
  const Mat<eT>&          x = tmp.M;
        subview_cube<eT>& t = *this;
  
  const uword t_n_rows   = t.n_rows;
  const uword t_n_cols   = t.n_cols;
  const uword t_n_slices = t.n_slices;
  
  const uword x_n_rows   = x.n_rows;
  const uword x_n_cols   = x.n_cols;
  
  if( ((x_n_rows == 1) || (x_n_cols == 1)) && (t_n_rows == 1) && (t_n_cols == 1) && (x.n_elem == t_n_slices) )
    {
    Cube<eT>& Q = const_cast< Cube<eT>& >(t.m);
    
    const uword t_aux_row1   = t.aux_row1;
    const uword t_aux_col1   = t.aux_col1;
    const uword t_aux_slice1 = t.aux_slice1;
    
    const eT* x_mem = x.memptr();
    
    uword i,j;
    for(i=0, j=1; j < t_n_slices; i+=2, j+=2)
      {
      const eT tmp_i = x_mem[i];
      const eT tmp_j = x_mem[j];
      
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + i) /= tmp_i;
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + j) /= tmp_j;
      }
    
    if(i < t_n_slices)
      {
      Q.at(t_aux_row1, t_aux_col1, t_aux_slice1 + i) /= x_mem[i];
      }
    }
  else
  if( (t_n_rows == x_n_rows) && (t_n_cols == x_n_cols) && (t_n_slices == 1) )
    {
    for(uword col = 0; col < t_n_cols; ++col)
      {
      arrayops::inplace_div( t.slice_colptr(0, col), x.colptr(col), t_n_rows );
      }
    }
  else
  if( (t_n_rows == x_n_rows) && (t_n_cols == 1) && (t_n_slices == x_n_cols) )
    {
    for(uword i=0; i < t_n_slices; ++i)
      {
      arrayops::inplace_div( t.slice_colptr(i, 0), x.colptr(i), t_n_rows );
      }
    }
  else
  if( (t_n_rows == 1) && (t_n_cols == x_n_rows) && (t_n_slices == x_n_cols) )
    {
    Cube<eT>& Q = const_cast< Cube<eT>& >(t.m);
    
    const uword t_aux_row1   = t.aux_row1;
    const uword t_aux_col1   = t.aux_col1;
    const uword t_aux_slice1 = t.aux_slice1;
    
    for(uword slice=0; slice < t_n_slices; ++slice)
      {
      const uword mod_slice = t_aux_slice1 + slice;
      
      const eT* x_colptr = x.colptr(slice);
      
      uword i,j;
      for(i=0, j=1; j < t_n_cols; i+=2, j+=2)
        {
        const eT tmp_i = x_colptr[i];
        const eT tmp_j = x_colptr[j];
        
        Q.at(t_aux_row1, t_aux_col1 + i, mod_slice) /= tmp_i;
        Q.at(t_aux_row1, t_aux_col1 + j, mod_slice) /= tmp_j;
        }
      
      if(i < t_n_cols)
        {
        Q.at(t_aux_row1, t_aux_col1 + i, mod_slice) /= x_colptr[i];
        }
      }
    }
  else
    {
    if(arma_config::debug)
      {
      arma_stop_logic_error( arma_incompat_size_string(t, x, "element-wise division") );
      }
    }
  }



template<typename eT>
template<typename gen_type>
inline
void
subview_cube<eT>::operator= (const GenCube<eT,gen_type>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(n_rows, n_cols, n_slices, in.n_rows, in.n_cols, in.n_slices, "copy into subcube");
  
  in.apply(*this);
  }



//! apply a functor to each element
template<typename eT>
template<typename functor>
inline
void
subview_cube<eT>::for_each(functor F)
  {
  arma_extra_debug_sigprint();
  
  Cube<eT>& Q = const_cast< Cube<eT>& >(m);
  
  const uword start_col   = aux_col1;
  const uword start_row   = aux_row1;
  const uword start_slice = aux_slice1;
  
  const uword end_col_plus1   = start_col   + n_cols;
  const uword end_row_plus1   = start_row   + n_rows;
  const uword end_slice_plus1 = start_slice + n_slices;
  
  for(uword uslice = start_slice; uslice < end_slice_plus1; ++uslice)
  for(uword ucol   = start_col;     ucol < end_col_plus1;   ++ucol  )
  for(uword urow   = start_row;     urow < end_row_plus1;   ++urow  )
    {
    F( Q.at(urow, ucol, uslice) );
    }
  }



template<typename eT>
template<typename functor>
inline
void
subview_cube<eT>::for_each(functor F) const
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT>& Q = m;
  
  const uword start_col   = aux_col1;
  const uword start_row   = aux_row1;
  const uword start_slice = aux_slice1;
  
  const uword end_col_plus1   = start_col   + n_cols;
  const uword end_row_plus1   = start_row   + n_rows;
  const uword end_slice_plus1 = start_slice + n_slices;
  
  for(uword uslice = start_slice; uslice < end_slice_plus1; ++uslice)
  for(uword ucol   = start_col;     ucol < end_col_plus1;   ++ucol  )
  for(uword urow   = start_row;     urow < end_row_plus1;   ++urow  )
    {
    F( Q.at(urow, ucol, uslice) );
    }
  }



//! transform each element in the subview using a functor
template<typename eT>
template<typename functor>
inline
void
subview_cube<eT>::transform(functor F)
  {
  arma_extra_debug_sigprint();
  
  Cube<eT>& Q = const_cast< Cube<eT>& >(m);
  
  const uword start_col   = aux_col1;
  const uword start_row   = aux_row1;
  const uword start_slice = aux_slice1;
  
  const uword end_col_plus1   = start_col   + n_cols;
  const uword end_row_plus1   = start_row   + n_rows;
  const uword end_slice_plus1 = start_slice + n_slices;
  
  for(uword uslice = start_slice; uslice < end_slice_plus1; ++uslice)
  for(uword ucol   = start_col;     ucol < end_col_plus1;   ++ucol  )
  for(uword urow   = start_row;     urow < end_row_plus1;   ++urow  )
    {
    Q.at(urow, ucol, uslice) = eT( F( Q.at(urow, ucol, uslice) ) );
    }
  }



//! imbue (fill) the subview with values provided by a functor
template<typename eT>
template<typename functor>
inline
void
subview_cube<eT>::imbue(functor F)
  {
  arma_extra_debug_sigprint();
  
  Cube<eT>& Q = const_cast< Cube<eT>& >(m);
  
  const uword start_col   = aux_col1;
  const uword start_row   = aux_row1;
  const uword start_slice = aux_slice1;
  
  const uword end_col_plus1   = start_col   + n_cols;
  const uword end_row_plus1   = start_row   + n_rows;
  const uword end_slice_plus1 = start_slice + n_slices;
  
  for(uword uslice = start_slice; uslice < end_slice_plus1; ++uslice)
  for(uword ucol   = start_col;     ucol < end_col_plus1;   ++ucol  )
  for(uword urow   = start_row;     urow < end_row_plus1;   ++urow  )
    {
    Q.at(urow, ucol, uslice) = eT( F() );
    }
  }



//! apply a lambda function to each slice, where each slice is interpreted as a matrix
template<typename eT>
inline
void
subview_cube<eT>::each_slice(const std::function< void(Mat<eT>&) >& F)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp1(n_rows, n_cols, arma_nozeros_indicator());
  Mat<eT> tmp2('j', tmp1.memptr(), n_rows, n_cols);
  
  for(uword slice_id=0; slice_id < n_slices; ++slice_id)
    {
    for(uword col_id=0; col_id < n_cols; ++col_id)
      {
      arrayops::copy( tmp1.colptr(col_id), slice_colptr(slice_id, col_id), n_rows );
      }
    
    F(tmp2);
    
    for(uword col_id=0; col_id < n_cols; ++col_id)
      {
      arrayops::copy( slice_colptr(slice_id, col_id), tmp1.colptr(col_id), n_rows );
      }
    }
  }



template<typename eT>
inline
void
subview_cube<eT>::each_slice(const std::function< void(const Mat<eT>&) >& F) const
  {
  arma_extra_debug_sigprint();
  
        Mat<eT> tmp1(n_rows, n_cols, arma_nozeros_indicator());
  const Mat<eT> tmp2('j', tmp1.memptr(), n_rows, n_cols);
  
  for(uword slice_id=0; slice_id < n_slices; ++slice_id)
    {
    for(uword col_id=0; col_id < n_cols; ++col_id)
      {
      arrayops::copy( tmp1.colptr(col_id), slice_colptr(slice_id, col_id), n_rows );
      }
    
    F(tmp2);
    }
  }



template<typename eT>
inline
void
subview_cube<eT>::replace(const eT old_val, const eT new_val)
  {
  arma_extra_debug_sigprint();
  
  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  for(uword slice = 0; slice < local_n_slices; ++slice)
    {
    for(uword col = 0; col < local_n_cols; ++col)
      {
      arrayops::replace(slice_colptr(slice,col), local_n_rows, old_val, new_val);
      }
    }
  }



template<typename eT>
inline
void
subview_cube<eT>::clean(const typename get_pod_type<eT>::result threshold)
  {
  arma_extra_debug_sigprint();

  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  for(uword slice = 0; slice < local_n_slices; ++slice)
    {
    for(uword col = 0; col < local_n_cols; ++col)
      {
      arrayops::clean( slice_colptr(slice,col), local_n_rows, threshold );
      }
    }
  }



template<typename eT>
inline
void
subview_cube<eT>::clamp(const eT min_val, const eT max_val)
  {
  arma_extra_debug_sigprint();
  
  if(is_cx<eT>::no)
    {
    arma_debug_check( (access::tmp_real(min_val) > access::tmp_real(max_val)), "subview_cube::clamp(): min_val must be less than max_val" );
    }
  else
    {
    arma_debug_check( (access::tmp_real(min_val) > access::tmp_real(max_val)), "subview_cube::clamp(): real(min_val) must be less than real(max_val)" );
    arma_debug_check( (access::tmp_imag(min_val) > access::tmp_imag(max_val)), "subview_cube::clamp(): imag(min_val) must be less than imag(max_val)" );
    }
  
  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  for(uword slice = 0; slice < local_n_slices; ++slice)
    {
    for(uword col = 0; col < local_n_cols; ++col)
      {
      arrayops::clamp( slice_colptr(slice,col), local_n_rows, min_val, max_val );
      }
    }
  }



template<typename eT>
inline
void
subview_cube<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();

  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  for(uword slice = 0; slice < local_n_slices; ++slice)
    {
    for(uword col = 0; col < local_n_cols; ++col)
      {
      arrayops::inplace_set( slice_colptr(slice,col), val, local_n_rows );
      }
    }
  }



template<typename eT>
inline
void
subview_cube<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  for(uword slice = 0; slice < local_n_slices; ++slice)
    {
    for(uword col = 0; col < local_n_cols; ++col)
      {
      arrayops::fill_zeros( slice_colptr(slice,col), local_n_rows );
      }
    }
  }



template<typename eT>
inline
void
subview_cube<eT>::ones()
  {
  arma_extra_debug_sigprint();
  
  fill(eT(1));
  }



template<typename eT>
inline
void
subview_cube<eT>::randu()
  {
  arma_extra_debug_sigprint();
  
  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  for(uword slice = 0; slice < local_n_slices; ++slice)
    {
    for(uword col = 0; col < local_n_cols; ++col)
      {
      arma_rng::randu<eT>::fill( slice_colptr(slice,col), local_n_rows );
      }
    }
  }



template<typename eT>
inline
void
subview_cube<eT>::randn()
  {
  arma_extra_debug_sigprint();
  
  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  for(uword slice = 0; slice < local_n_slices; ++slice)
    {
    for(uword col = 0; col < local_n_cols; ++col)
      {
      arma_rng::randn<eT>::fill( slice_colptr(slice,col), local_n_rows );
      }
    }
  }



template<typename eT>
inline
arma_warn_unused
bool
subview_cube<eT>::is_finite() const
  {
  arma_extra_debug_sigprint();
  
  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  for(uword slice = 0; slice < local_n_slices; ++slice)
    {
    for(uword col = 0; col < local_n_cols; ++col)
      {
      if(arrayops::is_finite(slice_colptr(slice,col), local_n_rows) == false)  { return false; }
      }
    }
  
  return true;
  }



template<typename eT>
inline
arma_warn_unused
bool
subview_cube<eT>::is_zero(const typename get_pod_type<eT>::result tol) const
  {
  arma_extra_debug_sigprint();
  
  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  for(uword slice = 0; slice < local_n_slices; ++slice)
    {
    for(uword col = 0; col < local_n_cols; ++col)
      {
      if(arrayops::is_zero(slice_colptr(slice,col), local_n_rows, tol) == false)  { return false; }
      }
    }
  
  return true;
  }



template<typename eT>
inline
arma_warn_unused
bool
subview_cube<eT>::has_inf() const
  {
  arma_extra_debug_sigprint();
  
  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  for(uword slice = 0; slice < local_n_slices; ++slice)
    {
    for(uword col = 0; col < local_n_cols; ++col)
      {
      if(arrayops::has_inf(slice_colptr(slice,col), local_n_rows))  { return true; }
      }
    }
  
  return false;
  }



template<typename eT>
inline
arma_warn_unused
bool
subview_cube<eT>::has_nan() const
  {
  arma_extra_debug_sigprint();
  
  const uword local_n_rows   = n_rows;
  const uword local_n_cols   = n_cols;
  const uword local_n_slices = n_slices;
  
  for(uword slice = 0; slice < local_n_slices; ++slice)
    {
    for(uword col = 0; col < local_n_cols; ++col)
      {
      if(arrayops::has_nan(slice_colptr(slice,col), local_n_rows))  { return true; }
      }
    }
  
  return false;
  }



template<typename eT>
inline
eT
subview_cube<eT>::at_alt(const uword i) const
  {
  return operator[](i);
  }



template<typename eT>
inline
eT&
subview_cube<eT>::operator[](const uword i)
  {
  const uword in_slice = i / n_elem_slice;
  const uword offset   = in_slice * n_elem_slice;
  const uword j        = i - offset;
  
  const uword in_col   = j / n_rows;
  const uword in_row   = j % n_rows;

  const uword index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  
  return access::rw( (const_cast< Cube<eT>& >(m)).mem[index] );
  }



template<typename eT>
inline
eT
subview_cube<eT>::operator[](const uword i) const
  {
  const uword in_slice = i / n_elem_slice;
  const uword offset   = in_slice * n_elem_slice;
  const uword j        = i - offset;
  
  const uword in_col   = j / n_rows;
  const uword in_row   = j % n_rows;

  const uword index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  
  return m.mem[index];
  }



template<typename eT>
inline
eT&
subview_cube<eT>::operator()(const uword i)
  {
  arma_debug_check_bounds( (i >= n_elem), "subview_cube::operator(): index out of bounds" );
  
  const uword in_slice = i / n_elem_slice;
  const uword offset   = in_slice * n_elem_slice;
  const uword j        = i - offset;
  
  const uword in_col   = j / n_rows;
  const uword in_row   = j % n_rows;

  const uword index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  
  return access::rw( (const_cast< Cube<eT>& >(m)).mem[index] );
  }



template<typename eT>
inline
eT
subview_cube<eT>::operator()(const uword i) const
  {
  arma_debug_check_bounds( (i >= n_elem), "subview_cube::operator(): index out of bounds" );
  
  const uword in_slice = i / n_elem_slice;
  const uword offset   = in_slice * n_elem_slice;
  const uword j        = i - offset;
  
  const uword in_col   = j / n_rows;
  const uword in_row   = j % n_rows;

  const uword index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  
  return m.mem[index];
  }



template<typename eT>
arma_inline
eT&
subview_cube<eT>::operator()(const uword in_row, const uword in_col, const uword in_slice)
  {
  arma_debug_check_bounds( ( (in_row >= n_rows) || (in_col >= n_cols) || (in_slice >= n_slices) ), "subview_cube::operator(): location out of bounds" );
  
  const uword index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  
  return access::rw( (const_cast< Cube<eT>& >(m)).mem[index] );
  }



template<typename eT>
arma_inline
eT
subview_cube<eT>::operator()(const uword in_row, const uword in_col, const uword in_slice) const
  {
  arma_debug_check_bounds( ( (in_row >= n_rows) || (in_col >= n_cols) || (in_slice >= n_slices) ), "subview_cube::operator(): location out of bounds" );
  
  const uword index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  
  return m.mem[index];
  }



template<typename eT>
arma_inline
eT&
subview_cube<eT>::at(const uword in_row, const uword in_col, const uword in_slice)
  {
  const uword index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  
  return access::rw( (const_cast< Cube<eT>& >(m)).mem[index] );
  }



template<typename eT>
arma_inline
eT
subview_cube<eT>::at(const uword in_row, const uword in_col, const uword in_slice) const
  {
  const uword index = (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 + in_row;
  
  return m.mem[index];
  }



template<typename eT>
arma_inline
eT*
subview_cube<eT>::slice_colptr(const uword in_slice, const uword in_col)
  {
  return & access::rw((const_cast< Cube<eT>& >(m)).mem[  (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1  ]);
  }



template<typename eT>
arma_inline
const eT*
subview_cube<eT>::slice_colptr(const uword in_slice, const uword in_col) const
  {
  return & m.mem[ (in_slice + aux_slice1)*m.n_elem_slice + (in_col + aux_col1)*m.n_rows + aux_row1 ];
  }



template<typename eT>
template<typename eT2>
inline
bool
subview_cube<eT>::check_overlap(const subview_cube<eT2>& x) const
  {
  if(is_same_type<eT,eT2>::value == false)  { return false; }
  
  const subview_cube<eT>& t = (*this);
  
  if(void_ptr(&(t.m)) != void_ptr(&(x.m)))  { return false; }
  
  if( (t.n_elem == 0) || (x.n_elem == 0) )  { return false; }
  
  const uword t_row_start  = t.aux_row1;
  const uword t_row_end_p1 = t_row_start + t.n_rows;
  
  const uword t_col_start  = t.aux_col1;
  const uword t_col_end_p1 = t_col_start + t.n_cols;
  
  const uword t_slice_start  = t.aux_slice1;
  const uword t_slice_end_p1 = t_slice_start + t.n_slices;
  
  
  const uword x_row_start  = x.aux_row1;
  const uword x_row_end_p1 = x_row_start + x.n_rows;
  
  const uword x_col_start  = x.aux_col1;
  const uword x_col_end_p1 = x_col_start + x.n_cols;
  
  const uword x_slice_start  = x.aux_slice1;
  const uword x_slice_end_p1 = x_slice_start + x.n_slices;
  
  
  const bool outside_rows   = ( (x_row_start   >= t_row_end_p1  ) || (t_row_start   >= x_row_end_p1  ) );
  const bool outside_cols   = ( (x_col_start   >= t_col_end_p1  ) || (t_col_start   >= x_col_end_p1  ) );
  const bool outside_slices = ( (x_slice_start >= t_slice_end_p1) || (t_slice_start >= x_slice_end_p1) );
  
  return ( (outside_rows == false) && (outside_cols == false) && (outside_slices == false) );
  }



template<typename eT>
inline
bool
subview_cube<eT>::check_overlap(const Mat<eT>& x) const
  {
  const subview_cube<eT>& t = *this;
  
  const uword t_aux_slice1        = t.aux_slice1;
  const uword t_aux_slice2_plus_1 = t_aux_slice1 + t.n_slices;
  
  for(uword slice = t_aux_slice1; slice < t_aux_slice2_plus_1; ++slice)
    {
    if(t.m.mat_ptrs[slice] != nullptr)
      {
      const Mat<eT>& y = *(t.m.mat_ptrs[slice]);
      
      if( x.memptr() == y.memptr() )  { return true; }
      }
    }
  
  return false;
  }



//! cube X = Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::extract(Cube<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();

  // NOTE: we're assuming that the cube has already been set to the correct size and there is no aliasing;
  // size setting and alias checking is done by either the Cube contructor or operator=()
  
  const uword n_rows   = in.n_rows;
  const uword n_cols   = in.n_cols;
  const uword n_slices = in.n_slices;
  
  arma_extra_debug_print(arma_str::format("out.n_rows = %d   out.n_cols = %d    out.n_slices = %d    in.m.n_rows = %d   in.m.n_cols = %d   in.m.n_slices = %d") % out.n_rows % out.n_cols % out.n_slices % in.m.n_rows % in.m.n_cols % in.m.n_slices);
  
  if( (in.aux_row1 == 0) && (n_rows == in.m.n_rows) )
    {
    for(uword s=0; s < n_slices; ++s)
      {
      arrayops::copy( out.slice_colptr(s,0), in.slice_colptr(s,0), in.n_elem_slice );
      }
    
    return;
    }
  
  for(uword s=0; s < n_slices; ++s)
  for(uword c=0; c < n_cols;   ++c)
    {
    arrayops::copy( out.slice_colptr(s,c), in.slice_colptr(s,c), n_rows );
    }
  }



//! cube X += Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::plus_inplace(Cube<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "addition");
  
  const uword n_rows   = out.n_rows;
  const uword n_cols   = out.n_cols;
  const uword n_slices = out.n_slices;
  
  for(uword slice = 0; slice<n_slices; ++slice)
    {
    for(uword col = 0; col<n_cols; ++col)
      {
      arrayops::inplace_plus( out.slice_colptr(slice,col), in.slice_colptr(slice,col), n_rows );
      }
    }
  }



//! cube X -= Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::minus_inplace(Cube<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "subtraction");
  
  const uword n_rows   = out.n_rows;
  const uword n_cols   = out.n_cols;
  const uword n_slices = out.n_slices;
  
  for(uword slice = 0; slice<n_slices; ++slice)
    {
    for(uword col = 0; col<n_cols; ++col)
      {
      arrayops::inplace_minus( out.slice_colptr(slice,col), in.slice_colptr(slice,col), n_rows );
      }
    }
  }



//! cube X %= Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::schur_inplace(Cube<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "element-wise multiplication");
  
  const uword n_rows   = out.n_rows;
  const uword n_cols   = out.n_cols;
  const uword n_slices = out.n_slices;
  
  for(uword slice = 0; slice<n_slices; ++slice)
    {
    for(uword col = 0; col<n_cols; ++col)
      {
      arrayops::inplace_mul( out.slice_colptr(slice,col), in.slice_colptr(slice,col), n_rows );
      }
    }
  }



//! cube X /= Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::div_inplace(Cube<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, in, "element-wise division");
  
  const uword n_rows   = out.n_rows;
  const uword n_cols   = out.n_cols;
  const uword n_slices = out.n_slices;
  
  for(uword slice = 0; slice<n_slices; ++slice)
    {
    for(uword col = 0; col<n_cols; ++col)
      {
      arrayops::inplace_div( out.slice_colptr(slice,col), in.slice_colptr(slice,col), n_rows );
      }
    }
  }



//! mat X = Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::extract(Mat<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_cube_as_mat(out, in, "copy into matrix", false);
  
  const uword in_n_rows   = in.n_rows;
  const uword in_n_cols   = in.n_cols;
  const uword in_n_slices = in.n_slices;
  
  const uword out_vec_state = out.vec_state;
  
  if(in_n_slices == 1)
    {
    out.set_size(in_n_rows, in_n_cols);
    
    for(uword col=0; col < in_n_cols; ++col)
      {
      arrayops::copy( out.colptr(col), in.slice_colptr(0, col), in_n_rows );
      }
    }
  else
    {
    if(out_vec_state == 0)
      {
      if(in_n_cols == 1)
        {
        out.set_size(in_n_rows, in_n_slices);
        
        for(uword i=0; i < in_n_slices; ++i)
          {
          arrayops::copy( out.colptr(i), in.slice_colptr(i, 0), in_n_rows );
          }
        }
      else
      if(in_n_rows == 1)
        {
        const Cube<eT>& Q = in.m;
        
        const uword in_aux_row1   = in.aux_row1;
        const uword in_aux_col1   = in.aux_col1;
        const uword in_aux_slice1 = in.aux_slice1;
        
        out.set_size(in_n_cols, in_n_slices);
        
        for(uword slice=0; slice < in_n_slices; ++slice)
          {
          const uword mod_slice = in_aux_slice1 + slice;
          
          eT* out_colptr = out.colptr(slice);
          
          uword i,j;
          for(i=0, j=1; j < in_n_cols; i+=2, j+=2)
            {
            const eT tmp_i = Q.at(in_aux_row1, in_aux_col1 + i, mod_slice);
            const eT tmp_j = Q.at(in_aux_row1, in_aux_col1 + j, mod_slice);
            
            out_colptr[i] = tmp_i;
            out_colptr[j] = tmp_j;
            }
          
          if(i < in_n_cols)
            {
            out_colptr[i] = Q.at(in_aux_row1, in_aux_col1 + i, mod_slice);
            }
          }
        }
      }
    else
      {
      out.set_size(in_n_slices);
      
      eT* out_mem = out.memptr();
      
      const Cube<eT>& Q = in.m;
      
      const uword in_aux_row1   = in.aux_row1;
      const uword in_aux_col1   = in.aux_col1;
      const uword in_aux_slice1 = in.aux_slice1;
      
      for(uword i=0; i<in_n_slices; ++i)
        {
        out_mem[i] = Q.at(in_aux_row1, in_aux_col1, in_aux_slice1 + i);
        }
      }
    }
  }



//! mat X += Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::plus_inplace(Mat<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_cube_as_mat(out, in, "addition", true);
  
  const uword in_n_rows   = in.n_rows;
  const uword in_n_cols   = in.n_cols;
  const uword in_n_slices = in.n_slices;
  
  const uword out_n_rows    = out.n_rows;
  const uword out_n_cols    = out.n_cols;
  const uword out_vec_state = out.vec_state;
  
  if(in_n_slices == 1)
    {
    if( (arma_config::debug) && ((out_n_rows != in_n_rows) || (out_n_cols != in_n_cols)) )
      {
      std::ostringstream tmp;
      
      tmp
        << "in-place addition: "
        << out_n_rows << 'x' << out_n_cols << " output matrix is incompatible with "
        <<  in_n_rows << 'x' <<  in_n_cols << 'x' << in_n_slices << " cube interpreted as "
        <<  in_n_rows << 'x' <<  in_n_cols << " matrix";
      
      arma_stop_logic_error(tmp.str());
      }
    
    for(uword col=0; col < in_n_cols; ++col)
      {
      arrayops::inplace_plus( out.colptr(col), in.slice_colptr(0, col), in_n_rows );
      }
    }
  else
    {
    if(out_vec_state == 0)
      {
      if( (in_n_rows == out_n_rows) && (in_n_cols == 1) && (in_n_slices == out_n_cols) )
        {
        for(uword i=0; i < in_n_slices; ++i)
          {
          arrayops::inplace_plus( out.colptr(i), in.slice_colptr(i, 0), in_n_rows );
          }
        }
      else
      if( (in_n_rows == 1) && (in_n_cols == out_n_rows) && (in_n_slices == out_n_cols) )
        {
        const Cube<eT>& Q = in.m;
        
        const uword in_aux_row1   = in.aux_row1;
        const uword in_aux_col1   = in.aux_col1;
        const uword in_aux_slice1 = in.aux_slice1;
        
        for(uword slice=0; slice < in_n_slices; ++slice)
          {
          const uword mod_slice = in_aux_slice1 + slice;
          
          eT* out_colptr = out.colptr(slice);
          
          uword i,j;
          for(i=0, j=1; j < in_n_cols; i+=2, j+=2)
            {
            const eT tmp_i = Q.at(in_aux_row1, in_aux_col1 + i, mod_slice);
            const eT tmp_j = Q.at(in_aux_row1, in_aux_col1 + j, mod_slice);
            
            out_colptr[i] += tmp_i;
            out_colptr[j] += tmp_j;
            }
          
          if(i < in_n_cols)
            {
            out_colptr[i] += Q.at(in_aux_row1, in_aux_col1 + i, mod_slice);
            }
          }
        }
      }
    else
      {
      eT* out_mem = out.memptr();
      
      const Cube<eT>& Q = in.m;
      
      const uword in_aux_row1   = in.aux_row1;
      const uword in_aux_col1   = in.aux_col1;
      const uword in_aux_slice1 = in.aux_slice1;
      
      for(uword i=0; i<in_n_slices; ++i)
        {
        out_mem[i] += Q.at(in_aux_row1, in_aux_col1, in_aux_slice1 + i);
        }
      }
    }
  }



//! mat X -= Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::minus_inplace(Mat<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_cube_as_mat(out, in, "subtraction", true);
  
  const uword in_n_rows   = in.n_rows;
  const uword in_n_cols   = in.n_cols;
  const uword in_n_slices = in.n_slices;
  
  const uword out_n_rows    = out.n_rows;
  const uword out_n_cols    = out.n_cols;
  const uword out_vec_state = out.vec_state;
  
  if(in_n_slices == 1)
    {
    if( (arma_config::debug) && ((out_n_rows != in_n_rows) || (out_n_cols != in_n_cols)) )
      {
      std::ostringstream tmp;
      
      tmp
        << "in-place subtraction: "
        << out_n_rows << 'x' << out_n_cols << " output matrix is incompatible with "
        <<  in_n_rows << 'x' <<  in_n_cols << 'x' << in_n_slices << " cube interpreted as "
        <<  in_n_rows << 'x' <<  in_n_cols << " matrix";
      
      arma_stop_logic_error(tmp.str());
      }
    
    for(uword col=0; col < in_n_cols; ++col)
      {
      arrayops::inplace_minus( out.colptr(col), in.slice_colptr(0, col), in_n_rows );
      }
    }
  else
    {
    if(out_vec_state == 0)
      {
      if( (in_n_rows == out_n_rows) && (in_n_cols == 1) && (in_n_slices == out_n_cols) )
        {
        for(uword i=0; i < in_n_slices; ++i)
          {
          arrayops::inplace_minus( out.colptr(i), in.slice_colptr(i, 0), in_n_rows );
          }
        }
      else
      if( (in_n_rows == 1) && (in_n_cols == out_n_rows) && (in_n_slices == out_n_cols) )
        {
        const Cube<eT>& Q = in.m;
        
        const uword in_aux_row1   = in.aux_row1;
        const uword in_aux_col1   = in.aux_col1;
        const uword in_aux_slice1 = in.aux_slice1;
        
        for(uword slice=0; slice < in_n_slices; ++slice)
          {
          const uword mod_slice = in_aux_slice1 + slice;
          
          eT* out_colptr = out.colptr(slice);
          
          uword i,j;
          for(i=0, j=1; j < in_n_cols; i+=2, j+=2)
            {
            const eT tmp_i = Q.at(in_aux_row1, in_aux_col1 + i, mod_slice);
            const eT tmp_j = Q.at(in_aux_row1, in_aux_col1 + j, mod_slice);
            
            out_colptr[i] -= tmp_i;
            out_colptr[j] -= tmp_j;
            }
          
          if(i < in_n_cols)
            {
            out_colptr[i] -= Q.at(in_aux_row1, in_aux_col1 + i, mod_slice);
            }
          }
        }
      }
    else
      {
      eT* out_mem = out.memptr();
      
      const Cube<eT>& Q = in.m;
      
      const uword in_aux_row1   = in.aux_row1;
      const uword in_aux_col1   = in.aux_col1;
      const uword in_aux_slice1 = in.aux_slice1;
      
      for(uword i=0; i<in_n_slices; ++i)
        {
        out_mem[i] -= Q.at(in_aux_row1, in_aux_col1, in_aux_slice1 + i);
        }
      }
    }
  }



//! mat X %= Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::schur_inplace(Mat<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_cube_as_mat(out, in, "element-wise multiplication", true);
  
  const uword in_n_rows   = in.n_rows;
  const uword in_n_cols   = in.n_cols;
  const uword in_n_slices = in.n_slices;
  
  const uword out_n_rows    = out.n_rows;
  const uword out_n_cols    = out.n_cols;
  const uword out_vec_state = out.vec_state;
  
  if(in_n_slices == 1)
    {
    if( (arma_config::debug) && ((out_n_rows != in_n_rows) || (out_n_cols != in_n_cols)) )
      {
      std::ostringstream tmp;
      
      tmp
        << "in-place element-wise multiplication: "
        << out_n_rows << 'x' << out_n_cols << " output matrix is incompatible with "
        <<  in_n_rows << 'x' <<  in_n_cols << 'x' << in_n_slices << " cube interpreted as "
        <<  in_n_rows << 'x' <<  in_n_cols << " matrix";
      
      arma_stop_logic_error(tmp.str());
      }
    
    for(uword col=0; col < in_n_cols; ++col)
      {
      arrayops::inplace_mul( out.colptr(col), in.slice_colptr(0, col), in_n_rows );
      }
    }
  else
    {
    if(out_vec_state == 0)
      {
      if( (in_n_rows == out_n_rows) && (in_n_cols == 1) && (in_n_slices == out_n_cols) )
        {
        for(uword i=0; i < in_n_slices; ++i)
          {
          arrayops::inplace_mul( out.colptr(i), in.slice_colptr(i, 0), in_n_rows );
          }
        }
      else
      if( (in_n_rows == 1) && (in_n_cols == out_n_rows) && (in_n_slices == out_n_cols) )
        {
        const Cube<eT>& Q = in.m;
        
        const uword in_aux_row1   = in.aux_row1;
        const uword in_aux_col1   = in.aux_col1;
        const uword in_aux_slice1 = in.aux_slice1;
        
        for(uword slice=0; slice < in_n_slices; ++slice)
          {
          const uword mod_slice = in_aux_slice1 + slice;
          
          eT* out_colptr = out.colptr(slice);
          
          uword i,j;
          for(i=0, j=1; j < in_n_cols; i+=2, j+=2)
            {
            const eT tmp_i = Q.at(in_aux_row1, in_aux_col1 + i, mod_slice);
            const eT tmp_j = Q.at(in_aux_row1, in_aux_col1 + j, mod_slice);
            
            out_colptr[i] *= tmp_i;
            out_colptr[j] *= tmp_j;
            }
          
          if(i < in_n_cols)
            {
            out_colptr[i] *= Q.at(in_aux_row1, in_aux_col1 + i, mod_slice);
            }
          }
        }
      }
    else
      {
      eT* out_mem = out.memptr();
      
      const Cube<eT>& Q = in.m;
      
      const uword in_aux_row1   = in.aux_row1;
      const uword in_aux_col1   = in.aux_col1;
      const uword in_aux_slice1 = in.aux_slice1;
      
      for(uword i=0; i<in_n_slices; ++i)
        {
        out_mem[i] *= Q.at(in_aux_row1, in_aux_col1, in_aux_slice1 + i);
        }
      }
    }
  }



//! mat X /= Y.subcube(...)
template<typename eT>
inline
void
subview_cube<eT>::div_inplace(Mat<eT>& out, const subview_cube<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_cube_as_mat(out, in, "element-wise division", true);
  
  const uword in_n_rows   = in.n_rows;
  const uword in_n_cols   = in.n_cols;
  const uword in_n_slices = in.n_slices;
  
  const uword out_n_rows    = out.n_rows;
  const uword out_n_cols    = out.n_cols;
  const uword out_vec_state = out.vec_state;
  
  if(in_n_slices == 1)
    {
    if( (arma_config::debug) && ((out_n_rows != in_n_rows) || (out_n_cols != in_n_cols)) )
      {
      std::ostringstream tmp;
      
      tmp
        << "in-place element-wise division: "
        << out_n_rows << 'x' << out_n_cols << " output matrix is incompatible with "
        <<  in_n_rows << 'x' <<  in_n_cols << 'x' << in_n_slices << " cube interpreted as "
        <<  in_n_rows << 'x' <<  in_n_cols << " matrix";
      
      arma_stop_logic_error(tmp.str());
      }
    
    for(uword col=0; col < in_n_cols; ++col)
      {
      arrayops::inplace_div( out.colptr(col), in.slice_colptr(0, col), in_n_rows );
      }
    }
  else
    {
    if(out_vec_state == 0)
      {
      if( (in_n_rows == out_n_rows) && (in_n_cols == 1) && (in_n_slices == out_n_cols) )
        {
        for(uword i=0; i < in_n_slices; ++i)
          {
          arrayops::inplace_div( out.colptr(i), in.slice_colptr(i, 0), in_n_rows );
          }
        }
      else
      if( (in_n_rows == 1) && (in_n_cols == out_n_rows) && (in_n_slices == out_n_cols) )
        {
        const Cube<eT>& Q = in.m;
        
        const uword in_aux_row1   = in.aux_row1;
        const uword in_aux_col1   = in.aux_col1;
        const uword in_aux_slice1 = in.aux_slice1;
        
        for(uword slice=0; slice < in_n_slices; ++slice)
          {
          const uword mod_slice = in_aux_slice1 + slice;
          
          eT* out_colptr = out.colptr(slice);
          
          uword i,j;
          for(i=0, j=1; j < in_n_cols; i+=2, j+=2)
            {
            const eT tmp_i = Q.at(in_aux_row1, in_aux_col1 + i, mod_slice);
            const eT tmp_j = Q.at(in_aux_row1, in_aux_col1 + j, mod_slice);
            
            out_colptr[i] /= tmp_i;
            out_colptr[j] /= tmp_j;
            }
          
          if(i < in_n_cols)
            {
            out_colptr[i] /= Q.at(in_aux_row1, in_aux_col1 + i, mod_slice);
            }
          }
        }
      }
    else
      {
      eT* out_mem = out.memptr();
      
      const Cube<eT>& Q = in.m;
      
      const uword in_aux_row1   = in.aux_row1;
      const uword in_aux_col1   = in.aux_col1;
      const uword in_aux_slice1 = in.aux_slice1;
      
      for(uword i=0; i<in_n_slices; ++i)
        {
        out_mem[i] /= Q.at(in_aux_row1, in_aux_col1, in_aux_slice1 + i);
        }
      }
    }
  }



template<typename eT>
inline
typename subview_cube<eT>::iterator
subview_cube<eT>::begin()
  {
  return iterator(*this, aux_row1, aux_col1, aux_slice1);
  }



template<typename eT>
inline
typename subview_cube<eT>::const_iterator
subview_cube<eT>::begin() const
  {
  return const_iterator(*this, aux_row1, aux_col1, aux_slice1);
  }



template<typename eT>
inline
typename subview_cube<eT>::const_iterator
subview_cube<eT>::cbegin() const
  {
  return const_iterator(*this, aux_row1, aux_col1, aux_slice1);
  }



template<typename eT>
inline
typename subview_cube<eT>::iterator
subview_cube<eT>::end()
  {
  return iterator(*this, aux_row1, aux_col1, aux_slice1 + n_slices);
  }



template<typename eT>
inline
typename subview_cube<eT>::const_iterator
subview_cube<eT>::end() const
  {
  return const_iterator(*this, aux_row1, aux_col1, aux_slice1 + n_slices);
  }



template<typename eT>
inline
typename subview_cube<eT>::const_iterator
subview_cube<eT>::cend() const
  {
  return const_iterator(*this, aux_row1, aux_col1, aux_slice1 + n_slices);
  }



//
//
//



template<typename eT>
inline
subview_cube<eT>::iterator::iterator()
  : M            (nullptr)
  , current_ptr  (nullptr)
  , current_row  (0      )
  , current_col  (0      )
  , current_slice(0      )
  , aux_row1     (0      )
  , aux_col1     (0      )
  , aux_row2_p1  (0      )
  , aux_col2_p1  (0      )
  {
  arma_extra_debug_sigprint();
  // Technically this iterator is invalid (it does not point to a valid element)
  }



template<typename eT>
inline
subview_cube<eT>::iterator::iterator(const iterator& X)
  : M            (X.M            )
  , current_ptr  (X.current_ptr  )
  , current_row  (X.current_row  )
  , current_col  (X.current_col  )
  , current_slice(X.current_slice)
  , aux_row1     (X.aux_row1     )
  , aux_col1     (X.aux_col1     )
  , aux_row2_p1  (X.aux_row2_p1  )
  , aux_col2_p1  (X.aux_col2_p1  )
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
subview_cube<eT>::iterator::iterator(subview_cube<eT>& in_sv, const uword in_row, const uword in_col, const uword in_slice)
  : M            (&(const_cast< Cube<eT>& >(in_sv.m)))
  , current_ptr  (&(M->at(in_row,in_col,in_slice))   )
  , current_row  (in_row                             )
  , current_col  (in_col                             )
  , current_slice(in_slice                           )
  , aux_row1     (in_sv.aux_row1                     )
  , aux_col1     (in_sv.aux_col1                     )
  , aux_row2_p1  (in_sv.aux_row1 + in_sv.n_rows      )
  , aux_col2_p1  (in_sv.aux_col1 + in_sv.n_cols      )
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
arma_warn_unused
eT&
subview_cube<eT>::iterator::operator*()
  {
  return (*current_ptr);
  }



template<typename eT>
inline
typename subview_cube<eT>::iterator&
subview_cube<eT>::iterator::operator++()
  {
  current_row++;
  
  if(current_row == aux_row2_p1)
    {
    current_row = aux_row1;
    current_col++;
    
    if(current_col == aux_col2_p1)
      {
      current_col = aux_col1;
      current_slice++;
      }
    
    current_ptr = &( (*M).at(current_row,current_col,current_slice) );
    }
  else
    {
    current_ptr++;
    }
  
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename subview_cube<eT>::iterator
subview_cube<eT>::iterator::operator++(int)
  {
  typename subview_cube<eT>::iterator temp(*this);
  
  ++(*this);
  
  return temp;
  }



template<typename eT>
inline
arma_warn_unused
bool
subview_cube<eT>::iterator::operator==(const iterator& rhs) const
  {
  return (current_ptr == rhs.current_ptr);
  }



template<typename eT>
inline
arma_warn_unused
bool
subview_cube<eT>::iterator::operator!=(const iterator& rhs) const
  {
  return (current_ptr != rhs.current_ptr);
  }



template<typename eT>
inline
arma_warn_unused
bool
subview_cube<eT>::iterator::operator==(const const_iterator& rhs) const
  {
  return (current_ptr == rhs.current_ptr);
  }



template<typename eT>
inline
arma_warn_unused
bool
subview_cube<eT>::iterator::operator!=(const const_iterator& rhs) const
  {
  return (current_ptr != rhs.current_ptr);
  }



//
//
//



template<typename eT>
inline
subview_cube<eT>::const_iterator::const_iterator()
  : M            (nullptr)
  , current_ptr  (nullptr)
  , current_row  (0   )
  , current_col  (0   )
  , current_slice(0   )
  , aux_row1     (0   )
  , aux_col1     (0   )
  , aux_row2_p1  (0   )
  , aux_col2_p1  (0   )
  {
  arma_extra_debug_sigprint();
  // Technically this iterator is invalid (it does not point to a valid element)
  }



template<typename eT>
inline
subview_cube<eT>::const_iterator::const_iterator(const iterator& X)
  : M            (X.M            )
  , current_ptr  (X.current_ptr  )
  , current_row  (X.current_row  )
  , current_col  (X.current_col  )
  , current_slice(X.current_slice)
  , aux_row1     (X.aux_row1     )
  , aux_col1     (X.aux_col1     )
  , aux_row2_p1  (X.aux_row2_p1  )
  , aux_col2_p1  (X.aux_col2_p1  )
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
subview_cube<eT>::const_iterator::const_iterator(const const_iterator& X)
  : M            (X.M            )
  , current_ptr  (X.current_ptr  )
  , current_row  (X.current_row  )
  , current_col  (X.current_col  )
  , current_slice(X.current_slice)
  , aux_row1     (X.aux_row1     )
  , aux_col1     (X.aux_col1     )
  , aux_row2_p1  (X.aux_row2_p1  )
  , aux_col2_p1  (X.aux_col2_p1  )
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
subview_cube<eT>::const_iterator::const_iterator(const subview_cube<eT>& in_sv, const uword in_row, const uword in_col, const uword in_slice)
  : M            (&(in_sv.m)                      )
  , current_ptr  (&(M->at(in_row,in_col,in_slice)))
  , current_row  (in_row                          )
  , current_col  (in_col                          )
  , current_slice(in_slice                        )
  , aux_row1     (in_sv.aux_row1                  )
  , aux_col1     (in_sv.aux_col1                  )
  , aux_row2_p1  (in_sv.aux_row1 + in_sv.n_rows   )
  , aux_col2_p1  (in_sv.aux_col1 + in_sv.n_cols   )
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
arma_warn_unused
const eT&
subview_cube<eT>::const_iterator::operator*()
  {
  return (*current_ptr);
  }



template<typename eT>
inline
typename subview_cube<eT>::const_iterator&
subview_cube<eT>::const_iterator::operator++()
  {
  current_row++;
  
  if(current_row == aux_row2_p1)
    {
    current_row = aux_row1;
    current_col++;
    
    if(current_col == aux_col2_p1)
      {
      current_col = aux_col1;
      current_slice++;
      }
    
    current_ptr = &( (*M).at(current_row,current_col,current_slice) );
    }
  else
    {
    current_ptr++;
    }
  
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename subview_cube<eT>::const_iterator
subview_cube<eT>::const_iterator::operator++(int)
  {
  typename subview_cube<eT>::const_iterator temp(*this);
  
  ++(*this);
  
  return temp;
  }



template<typename eT>
inline
arma_warn_unused
bool
subview_cube<eT>::const_iterator::operator==(const iterator& rhs) const
  {
  return (current_ptr == rhs.current_ptr);
  }



template<typename eT>
inline
arma_warn_unused
bool
subview_cube<eT>::const_iterator::operator!=(const iterator& rhs) const
  {
  return (current_ptr != rhs.current_ptr);
  }



template<typename eT>
inline
arma_warn_unused
bool
subview_cube<eT>::const_iterator::operator==(const const_iterator& rhs) const
  {
  return (current_ptr == rhs.current_ptr);
  }



template<typename eT>
inline
arma_warn_unused
bool
subview_cube<eT>::const_iterator::operator!=(const const_iterator& rhs) const
  {
  return (current_ptr != rhs.current_ptr);
  }



//! @}
