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



//! \addtogroup op_vectorise
//! @{



template<typename T1>
inline
void
op_vectorise_col::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_vectorise_col>& in)
  {
  arma_extra_debug_sigprint();
  
  op_vectorise_col::apply_direct(out, in.m);
  }



template<typename T1>
inline
void
op_vectorise_col::apply_direct(Mat<typename T1::elem_type>& out, const T1& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // allow detection of in-place operation
  if(is_Mat<T1>::value || (arma_config::openmp && Proxy<T1>::use_mp))
    {
    const unwrap<T1> U(expr);
    
    if(&out == &(U.M))
      {
      // output matrix is the same as the input matrix
      
      out.set_size(out.n_elem, 1);  // set_size() doesn't destroy data as long as the number of elements in the matrix remains the same
      }
    else
      {
      out.set_size(U.M.n_elem, 1);
      
      arrayops::copy(out.memptr(), U.M.memptr(), U.M.n_elem);
      }
    }
  else
  if(is_subview<T1>::value)
    {
    const subview<eT>& sv = reinterpret_cast< const subview<eT>& >(expr);
    
    if(&out == &(sv.m))
      {
      Mat<eT> tmp;
      
      op_vectorise_col::apply_subview(tmp, sv);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_vectorise_col::apply_subview(out, sv);
      }
    }
  else
    {
    const Proxy<T1> P(expr);
    
    const bool is_alias = P.is_alias(out);
    
    if(is_Mat<typename Proxy<T1>::stored_type>::value)
      {
      const quasi_unwrap<typename Proxy<T1>::stored_type> U(P.Q);
      
      if(is_alias)
        {
        Mat<eT> tmp(U.M.memptr(), U.M.n_elem, 1);
        
        out.steal_mem(tmp);
        }
      else
        {
        out.set_size(U.M.n_elem, 1);
        
        arrayops::copy(out.memptr(), U.M.memptr(), U.M.n_elem);
        }
      }
    else
      {
      if(is_alias)
        {
        Mat<eT> tmp;
        
        op_vectorise_col::apply_proxy(tmp, P);
        
        out.steal_mem(tmp);
        }
      else
        {
        op_vectorise_col::apply_proxy(out, P);
        }
      }
    }
  }



template<typename eT>
inline
void
op_vectorise_col::apply_subview(Mat<eT>& out, const subview<eT>& sv)
  {
  arma_extra_debug_sigprint();
  
  const uword sv_n_rows = sv.n_rows;
  const uword sv_n_cols = sv.n_cols;
  
  out.set_size(sv.n_elem, 1);
  
  eT* out_mem = out.memptr();
  
  for(uword col=0; col < sv_n_cols; ++col)
    {
    arrayops::copy(out_mem, sv.colptr(col), sv_n_rows);
    
    out_mem += sv_n_rows;
    }
  }



template<typename T1>
inline
void
op_vectorise_col::apply_proxy(Mat<typename T1::elem_type>& out, const Proxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword N = P.get_n_elem();
  
  out.set_size(N, 1);
  
  eT* outmem = out.memptr();
  
  if(Proxy<T1>::use_at == false)
    {
    // TODO: add handling of aligned access ?
    
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    uword i,j;
    
    for(i=0, j=1; j < N; i+=2, j+=2)
      {
      const eT tmp_i = A[i];
      const eT tmp_j = A[j];
      
      outmem[i] = tmp_i;
      outmem[j] = tmp_j;
      }
    
    if(i < N)
      {
      outmem[i] = A[i];
      }
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    if(n_rows == 1)
      {
      for(uword i=0; i < n_cols; ++i)
        {
        outmem[i] = P.at(0,i);
        }
      }
    else
      {
      for(uword col=0; col < n_cols; ++col)
      for(uword row=0; row < n_rows; ++row)
        {
        *outmem = P.at(row,col);
        outmem++;
        }
      }
    }
  }



template<typename T1>
inline
void
op_vectorise_row::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_vectorise_row>& in)
  {
  arma_extra_debug_sigprint();
  
  op_vectorise_row::apply_direct(out, in.m);
  }



template<typename T1>
inline
void
op_vectorise_row::apply_direct(Mat<typename T1::elem_type>& out, const T1& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(expr);
  
  if(P.is_alias(out))
    {
    Mat<eT> tmp;
    
    op_vectorise_row::apply_proxy(tmp, P);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_vectorise_row::apply_proxy(out, P);
    }
  }



template<typename T1>
inline
void
op_vectorise_row::apply_proxy(Mat<typename T1::elem_type>& out, const Proxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  const uword n_elem = P.get_n_elem();
  
  out.set_size(1, n_elem);
  
  eT* outmem = out.memptr();
  
  if(n_cols == 1)
    {
    if(is_Mat<typename Proxy<T1>::stored_type>::value)
      {
      const unwrap<typename Proxy<T1>::stored_type> tmp(P.Q);
      
      arrayops::copy(out.memptr(), tmp.M.memptr(), n_elem);
      }
    else
      {
      for(uword i=0; i < n_elem; ++i)  { outmem[i] = P.at(i,0); }
      }
    }
  else
    {
    for(uword row=0; row < n_rows; ++row)
      {
      uword i,j;
      
      for(i=0, j=1; j < n_cols; i+=2, j+=2)
        {
        const eT tmp_i = P.at(row,i);
        const eT tmp_j = P.at(row,j);
        
        *outmem = tmp_i; outmem++;
        *outmem = tmp_j; outmem++;
        }
      
      if(i < n_cols)
        {
        *outmem = P.at(row,i); outmem++;
        }
      }
    }
  }



template<typename T1>
inline
void
op_vectorise_all::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_vectorise_all>& in)
  {
  arma_extra_debug_sigprint();
  
  const uword dim = in.aux_uword_a;
  
  if(dim == 0)
    {
    op_vectorise_col::apply_direct(out, in.m);
    }
  else
    {
    op_vectorise_row::apply_direct(out, in.m);
    }
  }



//



template<typename T1>
inline
void
op_vectorise_cube_col::apply(Mat<typename T1::elem_type>& out, const CubeToMatOp<T1, op_vectorise_cube_col>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_same_type< T1, subview_cube<eT> >::yes)
    {
    op_vectorise_cube_col::apply_subview(out, reinterpret_cast< const subview_cube<eT>& >(in.m));
    }
  else
    {
    if(is_Cube<T1>::value || (arma_config::openmp && ProxyCube<T1>::use_mp))
      {
      op_vectorise_cube_col::apply_unwrap(out, in.m);
      }
    else
      {
      op_vectorise_cube_col::apply_proxy(out, in.m);
      }
    }
  }



template<typename eT>
inline
void
op_vectorise_cube_col::apply_subview(Mat<eT>& out, const subview_cube<eT>& sv)
  {
  arma_extra_debug_sigprint();
  
  const uword sv_nr = sv.n_rows;
  const uword sv_nc = sv.n_cols;
  const uword sv_ns = sv.n_slices;
  
  out.set_size(sv.n_elem, 1);
  
  eT* out_mem = out.memptr();
  
  for(uword s=0; s < sv_ns; ++s)
  for(uword c=0; c < sv_nc; ++c)
    {
    arrayops::copy(out_mem, sv.slice_colptr(s,c), sv_nr);
    
    out_mem += sv_nr;
    }
  }



template<typename T1>
inline
void
op_vectorise_cube_col::apply_unwrap(Mat<typename T1::elem_type>& out, const T1& expr)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<T1> U(expr);
  
  out.set_size(U.M.n_elem, 1);
  
  arrayops::copy(out.memptr(), U.M.memptr(), U.M.n_elem);
  }



template<typename T1>
inline
void
op_vectorise_cube_col::apply_proxy(Mat<typename T1::elem_type>& out, const T1& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const ProxyCube<T1> P(expr);
  
  if(is_Cube<typename ProxyCube<T1>::stored_type>::value)
    {
    op_vectorise_cube_col::apply_unwrap(out, P.Q);
    
    return;
    }
  
  const uword N = P.get_n_elem();
  
  out.set_size(N, 1);
  
  eT* outmem = out.memptr();
  
  if(ProxyCube<T1>::use_at == false)
    {
    typename ProxyCube<T1>::ea_type A = P.get_ea();
    
    uword i,j;
    
    for(i=0, j=1; j < N; i+=2, j+=2)
      {
      const eT tmp_i = A[i];
      const eT tmp_j = A[j];
      
      outmem[i] = tmp_i;
      outmem[j] = tmp_j;
      }
    
    if(i < N)
      {
      outmem[i] = A[i];
      }
    }
  else
    {
    const uword nr = P.get_n_rows();
    const uword nc = P.get_n_cols();
    const uword ns = P.get_n_slices();
    
    for(uword s=0; s < ns; ++s)
    for(uword c=0; c < nc; ++c)
    for(uword r=0; r < nr; ++r)
      {
      *outmem = P.at(r,c,s);
      outmem++;
      }
    }
  }



//! @}
