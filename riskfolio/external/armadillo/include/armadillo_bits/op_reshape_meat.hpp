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



//! \addtogroup op_reshape
//! @{



template<typename T1>
inline
void
op_reshape::apply(Mat<typename T1::elem_type>& actual_out, const Op<T1,op_reshape>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword new_n_rows = in.aux_uword_a;
  const uword new_n_cols = in.aux_uword_b;
  
  if(is_Mat<T1>::value || (arma_config::openmp && Proxy<T1>::use_mp))
    {
    const unwrap<T1>   U(in.m);
    const Mat<eT>& A = U.M;
    
    if(&actual_out == &A)
      {
      op_reshape::apply_mat_inplace(actual_out, new_n_rows, new_n_cols);
      }
    else
      {
      op_reshape::apply_mat_noalias(actual_out, A, new_n_rows, new_n_cols);
      }
    }
  else
    {
    const Proxy<T1> P(in.m);
    
    const bool is_alias = P.is_alias(actual_out);
    
    Mat<eT>  tmp;
    Mat<eT>& out = (is_alias) ? tmp : actual_out;
    
    if(is_Mat<typename Proxy<T1>::stored_type>::value)
      {
      const quasi_unwrap<typename Proxy<T1>::stored_type> U(P.Q);
      
      op_reshape::apply_mat_noalias(out, U.M, new_n_rows, new_n_cols);
      }
    else
      {
      op_reshape::apply_proxy_noalias(out, P, new_n_rows, new_n_cols);
      }
    
    if(is_alias)  { actual_out.steal_mem(tmp); }
    }
  }



template<typename eT>
inline
void
op_reshape::apply_mat_inplace(Mat<eT>& A, const uword new_n_rows, const uword new_n_cols)
  {
  arma_extra_debug_sigprint();
  
  const uword new_n_elem = new_n_rows * new_n_cols;
  
  if(A.n_elem == new_n_elem)  { A.set_size(new_n_rows, new_n_cols); return; }
  
  Mat<eT> B;
  
  op_reshape::apply_mat_noalias(B, A, new_n_rows, new_n_cols);
  
  A.steal_mem(B);
  }



template<typename eT>
inline
void
op_reshape::apply_mat_noalias(Mat<eT>& out, const Mat<eT>& A, const uword new_n_rows, const uword new_n_cols)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(new_n_rows, new_n_cols);
  
  const uword n_elem_to_copy = (std::min)(A.n_elem, out.n_elem);
  
  eT* out_mem = out.memptr();
  
  arrayops::copy( out_mem, A.memptr(), n_elem_to_copy );
  
  if(n_elem_to_copy < out.n_elem)
    {
    const uword n_elem_leftover = out.n_elem - n_elem_to_copy;
    
    arrayops::fill_zeros(&(out_mem[n_elem_to_copy]), n_elem_leftover);
    }
  }



template<typename T1>
inline
void
op_reshape::apply_proxy_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword new_n_rows, const uword new_n_cols)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  out.set_size(new_n_rows, new_n_cols);
  
  const uword n_elem_to_copy = (std::min)(P.get_n_elem(), out.n_elem);
  
  eT* out_mem = out.memptr();
  
  if(Proxy<T1>::use_at == false)
    {
    typename Proxy<T1>::ea_type Pea = P.get_ea();
    
    for(uword i=0; i < n_elem_to_copy; ++i)  { out_mem[i] = Pea[i]; }
    }
  else
    {
    uword i = 0;
    
    const uword P_n_rows = P.get_n_rows();
    const uword P_n_cols = P.get_n_cols();
    
    for(uword col=0; col < P_n_cols; ++col)
    for(uword row=0; row < P_n_rows; ++row)
      {
      if(i >= n_elem_to_copy)  { goto nested_loop_end; }
      
      out_mem[i] = P.at(row,col);
      
      ++i;
      }
    
    nested_loop_end: ;
    }
  
  if(n_elem_to_copy < out.n_elem)
    {
    const uword n_elem_leftover = out.n_elem - n_elem_to_copy;
    
    arrayops::fill_zeros(&(out_mem[n_elem_to_copy]), n_elem_leftover);
    }
  }



template<typename T1>
inline
void
op_reshape::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_reshape>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube<T1> U(in.m);
  const Cube<eT>& A   = U.M;
  
  const uword new_n_rows   = in.aux_uword_a;
  const uword new_n_cols   = in.aux_uword_b;
  const uword new_n_slices = in.aux_uword_c;
  
  if(&out == &A)
    {
    op_reshape::apply_cube_inplace(out, new_n_rows, new_n_cols, new_n_slices);
    }
  else
    {
    op_reshape::apply_cube_noalias(out, A, new_n_rows, new_n_cols, new_n_slices);
    }
  }
  


template<typename eT>
inline
void
op_reshape::apply_cube_inplace(Cube<eT>& A, const uword new_n_rows, const uword new_n_cols, const uword new_n_slices)
  {
  arma_extra_debug_sigprint();
  
  const uword new_n_elem = new_n_rows * new_n_cols * new_n_slices;
  
  if(A.n_elem == new_n_elem)  { A.set_size(new_n_rows, new_n_cols, new_n_slices); return; }
  
  Cube<eT> B;
  
  op_reshape::apply_cube_noalias(B, A, new_n_rows, new_n_cols, new_n_slices);
  
  A.steal_mem(B);
  }



template<typename eT>
inline
void
op_reshape::apply_cube_noalias(Cube<eT>& out, const Cube<eT>& A, const uword new_n_rows, const uword new_n_cols, const uword new_n_slices)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(new_n_rows, new_n_cols, new_n_slices);
  
  const uword n_elem_to_copy = (std::min)(A.n_elem, out.n_elem);
  
  eT* out_mem = out.memptr();
  
  arrayops::copy( out_mem, A.memptr(), n_elem_to_copy );
  
  if(n_elem_to_copy < out.n_elem)
    {
    const uword n_elem_leftover = out.n_elem - n_elem_to_copy;
    
    arrayops::fill_zeros(&(out_mem[n_elem_to_copy]), n_elem_leftover);
    }
  }



//



template<typename T1>
arma_cold
inline
void
op_reshape_old::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_reshape_old>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword new_n_rows = in.aux_uword_a;
  const uword new_n_cols = in.aux_uword_b;
  const uword dim        = in.aux_uword_c;
  
  const unwrap<T1>   U(in.m);
  const Mat<eT>& A = U.M;
  
  if(&out == &A)
    {
    op_reshape_old::apply_mat_inplace(out, new_n_rows, new_n_cols, dim);
    }
  else
    {
    op_reshape_old::apply_mat_noalias(out, A, new_n_rows, new_n_cols, dim);
    }
  }



template<typename eT>
arma_cold
inline
void
op_reshape_old::apply_mat_inplace(Mat<eT>& A, const uword new_n_rows, const uword new_n_cols, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  if(dim == 0)
    {
    op_reshape::apply_mat_inplace(A, new_n_rows, new_n_cols);
    }
  else
  if(dim == 1)
    {
    Mat<eT> tmp;
    
    op_strans::apply_mat_noalias(tmp, A);
    
    op_reshape::apply_mat_noalias(A, tmp, new_n_rows, new_n_cols);
    }
  }



template<typename eT>
arma_cold
inline
void
op_reshape_old::apply_mat_noalias(Mat<eT>& out, const Mat<eT>& A, const uword new_n_rows, const uword new_n_cols, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  if(dim == 0)
    {
    op_reshape::apply_mat_noalias(out, A, new_n_rows, new_n_cols);
    }
  else
  if(dim == 1)
    {
    Mat<eT> tmp;
    
    op_strans::apply_mat_noalias(tmp, A);
    
    op_reshape::apply_mat_noalias(out, tmp, new_n_rows, new_n_cols);
    }
  }



//! @}
