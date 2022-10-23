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


//! \addtogroup subview_cube_slices
//! @{


template<typename eT, typename T1>
inline
subview_cube_slices<eT,T1>::~subview_cube_slices()
  {
  arma_extra_debug_sigprint();
  }


template<typename eT, typename T1>
arma_inline
subview_cube_slices<eT,T1>::subview_cube_slices
  (
  const Cube<eT>&       in_m,
  const Base<uword,T1>& in_si
  )
  : m      (in_m )
  , base_si(in_si)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::inplace_rand(const uword rand_mode)
  {
  arma_extra_debug_sigprint();
  
  Cube<eT>& m_local = const_cast< Cube<eT>& >(m);
  
  const uword m_n_slices     = m_local.n_slices;
  const uword m_n_elem_slice = m_local.n_elem_slice;
  
  const quasi_unwrap<T1> U(base_si.get_ref());
  const umat& si       = U.M;
  
  arma_debug_check
    (
    ( (si.is_vec() == false) && (si.is_empty() == false) ),
    "Cube::slices(): given object must be a vector"
    );
  
  const uword* si_mem    = si.memptr();
  const uword  si_n_elem = si.n_elem;
  
  for(uword si_count=0; si_count < si_n_elem; ++si_count)
    {
    const uword i = si_mem[si_count];
    
    arma_debug_check_bounds( (i >= m_n_slices), "Cube::slices(): index out of bounds" );
    
    eT* m_slice_ptr = m_local.slice_memptr(i);
    
    if(rand_mode == 0)  { arma_rng::randu<eT>::fill(m_slice_ptr, m_n_elem_slice); }
    if(rand_mode == 1)  { arma_rng::randn<eT>::fill(m_slice_ptr, m_n_elem_slice); }
    }
  }



template<typename eT, typename T1>
template<typename op_type>
inline
void
subview_cube_slices<eT,T1>::inplace_op(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Cube<eT>& m_local = const_cast< Cube<eT>& >(m);
  
  const uword m_n_slices     = m_local.n_slices;
  const uword m_n_elem_slice = m_local.n_elem_slice;
  
  const quasi_unwrap<T1> U(base_si.get_ref());
  const umat& si       = U.M;
  
  arma_debug_check
    (
    ( (si.is_vec() == false) && (si.is_empty() == false) ),
    "Cube::slices(): given object must be a vector"
    );
  
  const uword* si_mem    = si.memptr();
  const uword  si_n_elem = si.n_elem;
  
  for(uword si_count=0; si_count < si_n_elem; ++si_count)
    {
    const uword i = si_mem[si_count];
    
    arma_debug_check_bounds( (i >= m_n_slices), "Cube::slices(): index out of bounds" );
    
    eT* m_slice_ptr = m_local.slice_memptr(i);
    
    if(is_same_type<op_type, op_internal_equ  >::yes) { arrayops::inplace_set  (m_slice_ptr, val, m_n_elem_slice); }
    if(is_same_type<op_type, op_internal_plus >::yes) { arrayops::inplace_plus (m_slice_ptr, val, m_n_elem_slice); }
    if(is_same_type<op_type, op_internal_minus>::yes) { arrayops::inplace_minus(m_slice_ptr, val, m_n_elem_slice); }
    if(is_same_type<op_type, op_internal_schur>::yes) { arrayops::inplace_mul  (m_slice_ptr, val, m_n_elem_slice); }
    if(is_same_type<op_type, op_internal_div  >::yes) { arrayops::inplace_div  (m_slice_ptr, val, m_n_elem_slice); }
    }
  }



template<typename eT, typename T1>
template<typename op_type, typename expr>
inline
void
subview_cube_slices<eT,T1>::inplace_op(const BaseCube<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  Cube<eT>& m_local = const_cast< Cube<eT>& >(m);
  
  const uword m_n_slices     = m_local.n_slices;
  const uword m_n_elem_slice = m_local.n_elem_slice;
  
  const quasi_unwrap<T1> U(base_si.get_ref());
  const umat& si       = U.M;
  
  arma_debug_check
    (
    ( (si.is_vec() == false) && (si.is_empty() == false) ),
    "Cube::slices(): given object must be a vector"
    );
  
  const uword* si_mem    = si.memptr();
  const uword  si_n_elem = si.n_elem;
  
  const unwrap_cube_check<expr> tmp(x.get_ref(), m_local);
  const Cube<eT>& X           = tmp.M;
  
  arma_debug_assert_same_size( m_local.n_rows, m_local.n_cols, si_n_elem, X.n_rows, X.n_cols, X.n_slices, "Cube::slices()" );
  
  for(uword si_count=0; si_count < si_n_elem; ++si_count)
    {
    const uword i = si_mem[si_count];
    
    arma_debug_check_bounds( (i >= m_n_slices), "Cube::slices(): index out of bounds" );
    
          eT* m_slice_ptr = m_local.slice_memptr(i);
    const eT* X_slice_ptr =       X.slice_memptr(si_count);
    
    if(is_same_type<op_type, op_internal_equ  >::yes) { arrayops::copy         (m_slice_ptr, X_slice_ptr, m_n_elem_slice); }
    if(is_same_type<op_type, op_internal_plus >::yes) { arrayops::inplace_plus (m_slice_ptr, X_slice_ptr, m_n_elem_slice); }
    if(is_same_type<op_type, op_internal_minus>::yes) { arrayops::inplace_minus(m_slice_ptr, X_slice_ptr, m_n_elem_slice); }
    if(is_same_type<op_type, op_internal_schur>::yes) { arrayops::inplace_mul  (m_slice_ptr, X_slice_ptr, m_n_elem_slice); }
    if(is_same_type<op_type, op_internal_div  >::yes) { arrayops::inplace_div  (m_slice_ptr, X_slice_ptr, m_n_elem_slice); }
    }
  }



//
//



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_equ>(val);
  }



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::zeros()
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_equ>(eT(0));
  }



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::ones()
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_equ>(eT(1));
  }



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::randu()
  {
  arma_extra_debug_sigprint();
  
  inplace_rand(0);
  }



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::randn()
  {
  arma_extra_debug_sigprint();
  
  inplace_rand(1);
  }



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::operator+= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_plus>(val);
  }



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::operator-= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_minus>(val);
  }



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::operator*= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_schur>(val);
  }



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::operator/= (const eT val)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_div>(val);
  }



//
//



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_cube_slices<eT,T1>::operator_equ(const subview_cube_slices<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_equ>(x);
  }




template<typename eT, typename T1>
template<typename T2>
inline
void
subview_cube_slices<eT,T1>::operator= (const subview_cube_slices<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  (*this).operator_equ(x);
  }



//! work around compiler bugs
template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::operator= (const subview_cube_slices<eT,T1>& x)
  {
  arma_extra_debug_sigprint();
  
  (*this).operator_equ(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_cube_slices<eT,T1>::operator+= (const subview_cube_slices<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_plus>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_cube_slices<eT,T1>::operator-= (const subview_cube_slices<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_minus>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_cube_slices<eT,T1>::operator%= (const subview_cube_slices<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_schur>(x);
  }



template<typename eT, typename T1>
template<typename T2>
inline
void
subview_cube_slices<eT,T1>::operator/= (const subview_cube_slices<eT,T2>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_div>(x);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
subview_cube_slices<eT,T1>::operator= (const BaseCube<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_equ>(x);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
subview_cube_slices<eT,T1>::operator+= (const BaseCube<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_plus>(x);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
subview_cube_slices<eT,T1>::operator-= (const BaseCube<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_minus>(x);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
subview_cube_slices<eT,T1>::operator%= (const BaseCube<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_schur>(x);
  }



template<typename eT, typename T1>
template<typename expr>
inline
void
subview_cube_slices<eT,T1>::operator/= (const BaseCube<eT,expr>& x)
  {
  arma_extra_debug_sigprint();
  
  inplace_op<op_internal_div>(x);
  }



//
//



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::extract(Cube<eT>& out, const subview_cube_slices<eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT>& m_local = in.m;
  
  const uword m_n_slices     = m_local.n_slices;
  const uword m_n_elem_slice = m_local.n_elem_slice;
  
  const quasi_unwrap<T1> U(in.base_si.get_ref());
  const umat& si = U.M;
  
  arma_debug_check
    (
    ( (si.is_vec() == false) && (si.is_empty() == false) ),
    "Cube::slices(): given object must be a vector"
    );
  
  const uword* si_mem    = si.memptr();
  const uword  si_n_elem = si.n_elem;
  
  out.set_size(m_local.n_rows, m_local.n_cols, si_n_elem);
  
  for(uword si_count=0; si_count < si_n_elem; ++si_count)
    {
    const uword i = si_mem[si_count];
    
    arma_debug_check_bounds( (i >= m_n_slices), "Cube::slices(): index out of bounds" );
    
          eT* out_slice_ptr =     out.slice_memptr(si_count);
    const eT*   m_slice_ptr = m_local.slice_memptr(i);
    
    arrayops::copy(out_slice_ptr, m_slice_ptr, m_n_elem_slice);
    }
  }



// TODO: implement a dedicated function instead of creating a temporary
template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::plus_inplace(Cube<eT>& out, const subview_cube_slices& in)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT> tmp(in);
  
  out += tmp;
  }



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::minus_inplace(Cube<eT>& out, const subview_cube_slices& in)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT> tmp(in);
  
  out -= tmp;
  }



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::schur_inplace(Cube<eT>& out, const subview_cube_slices& in)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT> tmp(in);
  
  out %= tmp;
  }



template<typename eT, typename T1>
inline
void
subview_cube_slices<eT,T1>::div_inplace(Cube<eT>& out, const subview_cube_slices& in)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT> tmp(in);
  
  out /= tmp;
  }



//! @}
