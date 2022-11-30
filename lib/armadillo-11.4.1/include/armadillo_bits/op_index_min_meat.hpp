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


//! \addtogroup op_index_min
//! @{



template<typename T1>
inline
void
op_index_min::apply(Mat<uword>& out, const mtOp<uword,T1,op_index_min>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  arma_debug_check( (dim > 1), "index_min(): parameter 'dim' must be 0 or 1" );
  
  const quasi_unwrap<T1> U(in.m);
  const Mat<eT>& X = U.M;
  
  if(U.is_alias(out) == false)
    {
    op_index_min::apply_noalias(out, X, dim);
    }
  else
    {
    Mat<uword> tmp;
    
    op_index_min::apply_noalias(tmp, X, dim);
    
    out.steal_mem(tmp);
    }
  }



template<typename eT>
inline
void
op_index_min::apply_noalias(Mat<uword>& out, const Mat<eT>& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_index_min::apply(): dim = 0");
    
    out.set_size((X_n_rows > 0) ? 1 : 0, X_n_cols);
    
    if(X_n_rows == 0)  { return; }
    
    uword* out_mem = out.memptr();
    
    for(uword col=0; col < X_n_cols; ++col)
      {
      op_min::direct_min( X.colptr(col), X_n_rows, out_mem[col] );
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_index_min::apply(): dim = 1");
    
    out.zeros(X_n_rows, (X_n_cols > 0) ? 1 : 0);
    
    if(X_n_cols == 0)  { return; }
    
    uword* out_mem = out.memptr();
    
    Col<T> tmp(X_n_rows, arma_nozeros_indicator());
    
    T* tmp_mem = tmp.memptr();
    
    if(is_cx<eT>::yes)
      {
      const eT* col_mem = X.colptr(0);
      
      for(uword row=0; row < X_n_rows; ++row)
        {
        tmp_mem[row] = eop_aux::arma_abs(col_mem[row]);
        }
      }
    else
      {
      arrayops::copy(tmp_mem, (T*)(X.colptr(0)), X_n_rows);
      }
    
    for(uword col=1; col < X_n_cols; ++col)
      {
      const eT* col_mem = X.colptr(col);
      
      for(uword row=0; row < X_n_rows; ++row)
        {
        T& min_val = tmp_mem[row];
        T  col_val = (is_cx<eT>::yes) ? T(eop_aux::arma_abs(col_mem[row])) : T(access::tmp_real(col_mem[row]));
        
        if(min_val > col_val)
          {
          min_val = col_val;
          
          out_mem[row] = col;
          }
        }
      }
    }
  }



template<typename T1>
inline
void
op_index_min::apply(Cube<uword>& out, const mtOpCube<uword, T1, op_index_min>& in)
  {
  arma_extra_debug_sigprint();
  
  const uword dim = in.aux_uword_a;
  arma_debug_check( (dim > 2), "index_min(): parameter 'dim' must be 0 or 1 or 2" );
  
  const unwrap_cube<T1> U(in.m);
  
  if(U.is_alias(out) == false)
    {
    op_index_min::apply_noalias(out, U.M, dim);
    }
  else
    {
    Cube<uword> tmp;
    
    op_index_min::apply_noalias(tmp, U.M, dim);
    
    out.steal_mem(tmp);
    }
  }



template<typename eT>
inline
void
op_index_min::apply_noalias(Cube<uword>& out, const Cube<eT>& X, const uword dim, const typename arma_not_cx<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword X_n_rows   = X.n_rows;
  const uword X_n_cols   = X.n_cols;
  const uword X_n_slices = X.n_slices;
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_index_min::apply(): dim = 0");
    
    out.set_size((X_n_rows > 0) ? 1 : 0, X_n_cols, X_n_slices);
    
    if(out.is_empty() || X.is_empty())  { return; }
    
    for(uword slice=0; slice < X_n_slices; ++slice)
      {
      uword* out_mem = out.slice_memptr(slice);
      
      for(uword col=0; col < X_n_cols; ++col)
        {
        op_min::direct_min( X.slice_colptr(slice,col), X_n_rows, out_mem[col] );
        }
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_index_min::apply(): dim = 1");
    
    out.zeros(X_n_rows, (X_n_cols > 0) ? 1 : 0, X_n_slices);
    
    if(out.is_empty() || X.is_empty())  { return; }
    
    Col<eT> tmp(X_n_rows, arma_nozeros_indicator());
    
    eT* tmp_mem = tmp.memptr();
    
    for(uword slice=0; slice < X_n_slices; ++slice)
      {
      uword* out_mem = out.slice_memptr(slice);
      
      arrayops::copy(tmp_mem, X.slice_colptr(slice,0), X_n_rows);
      
      for(uword col=1; col < X_n_cols; ++col)
        {
        const eT* col_mem = X.slice_colptr(slice,col);
        
        for(uword row=0; row < X_n_rows; ++row)
          {
          const eT val = col_mem[row];
          
          if(val < tmp_mem[row])
            {
            tmp_mem[row] = val;
            out_mem[row] = col;
            }
          }
        }
      }
    }
  else
  if(dim == 2)
    {
    arma_extra_debug_print("op_index_min::apply(): dim = 2");
    
    out.zeros(X_n_rows, X_n_cols, (X_n_slices > 0) ? 1 : 0);
    
    if(out.is_empty() || X.is_empty())  { return; }
    
    Mat<eT> tmp(X.slice_memptr(0), X_n_rows, X_n_cols);  // copy slice 0
    
    eT*    tmp_mem = tmp.memptr();
    uword* out_mem = out.memptr();
    
    const uword N = X.n_elem_slice;
    
    for(uword slice=1; slice < X_n_slices; ++slice)
      {
      const eT* X_slice_mem = X.slice_memptr(slice);
      
      for(uword i=0; i < N; ++i)
        {
        const eT val = X_slice_mem[i];
        
        if(val < tmp_mem[i])
          {
          tmp_mem[i] = val;
          out_mem[i] = slice;
          }
        }
      }
    }
  }



template<typename eT>
inline
void
op_index_min::apply_noalias(Cube<uword>& out, const Cube<eT>& X, const uword dim, const typename arma_cx_only<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename get_pod_type<eT>::result T;
  
  const uword X_n_rows   = X.n_rows;
  const uword X_n_cols   = X.n_cols;
  const uword X_n_slices = X.n_slices;
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_index_min::apply(): dim = 0");
    
    out.set_size((X_n_rows > 0) ? 1 : 0, X_n_cols, X_n_slices);
    
    if(out.is_empty() || X.is_empty())  { return; }
    
    for(uword slice=0; slice < X_n_slices; ++slice)
      {
      uword* out_mem = out.slice_memptr(slice);
      
      for(uword col=0; col < X_n_cols; ++col)
        {
        op_min::direct_min( X.slice_colptr(slice,col), X_n_rows, out_mem[col] );
        }
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_index_min::apply(): dim = 1");
    
    out.zeros(X_n_rows, (X_n_cols > 0) ? 1 : 0, X_n_slices);
    
    if(out.is_empty() || X.is_empty())  { return; }
    
    Col<T> tmp(X_n_rows, arma_nozeros_indicator());
    
    T* tmp_mem = tmp.memptr();
    
    for(uword slice=0; slice < X_n_slices; ++slice)
      {
      uword* out_mem = out.slice_memptr(slice);
      
      const eT* col0_mem = X.slice_colptr(slice,0);
      
      for(uword row=0; row < X_n_rows; ++row)
        {
        tmp_mem[row] = std::abs( col0_mem[row] );
        }
      
      for(uword col=1; col < X_n_cols; ++col)
        {
        const eT* col_mem = X.slice_colptr(slice,col);
        
        for(uword row=0; row < X_n_rows; ++row)
          {
          const T val = std::abs( col_mem[row] );
          
          if(val < tmp_mem[row])
            {
            tmp_mem[row] = val;
            out_mem[row] = col;
            }
          }
        }
      }
    }
  else
  if(dim == 2)
    {
    arma_extra_debug_print("op_index_min::apply(): dim = 2");
    
    out.zeros(X_n_rows, X_n_cols, (X_n_slices > 0) ? 1 : 0);
    
    if(out.is_empty() || X.is_empty())  { return; }
    
    uword* out_mem = out.memptr();
    
    Mat<T> tmp(X_n_rows, X_n_cols, arma_nozeros_indicator());
    
           T*      tmp_mem = tmp.memptr();
    const eT* X_slice0_mem = X.slice_memptr(0);
    
    const uword N = X.n_elem_slice;
    
    for(uword i=0; i<N; ++i)
      {
      tmp_mem[i] = std::abs( X_slice0_mem[i] );
      }
    
    for(uword slice=1; slice < X_n_slices; ++slice)
      {
      const eT* X_slice_mem = X.slice_memptr(slice);
      
      for(uword i=0; i < N; ++i)
        {
        const T val = std::abs( X_slice_mem[i] );
        
        if(val < tmp_mem[i])
          {
          tmp_mem[i] = val;
          out_mem[i] = slice;
          }
        }
      }
    }
  }



template<typename T1>
inline
void
op_index_min::apply(Mat<uword>& out, const SpBase<typename T1::elem_type,T1>& expr, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  arma_debug_check( (dim > 1), "index_min(): parameter 'dim' must be 0 or 1" );
  
  const unwrap_spmat<T1> U(expr.get_ref());
  const SpMat<eT>& X   = U.M;
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_index_min::apply(): dim = 0");
    
    out.set_size((X_n_rows > 0) ? 1 : 0, X_n_cols);
    
    if(X_n_rows == 0)  { return; }
    
    uword* out_mem = out.memptr();
    
    for(uword col=0; col < X_n_cols; ++col)
      {
      out_mem[col] = X.col(col).index_min();
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_index_min::apply(): dim = 1");
    
    out.set_size(X_n_rows, (X_n_cols > 0) ? 1 : 0);
    
    if(X_n_cols == 0)  { return; }
    
    uword* out_mem = out.memptr();
    
    const SpMat<eT> Xt = X.st();
    
    for(uword row=0; row < X_n_rows; ++row)
      {
      out_mem[row] = Xt.col(row).index_min();
      }
    }
  }



//! @}
