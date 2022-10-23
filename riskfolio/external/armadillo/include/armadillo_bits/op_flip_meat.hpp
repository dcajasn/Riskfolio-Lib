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


//! \addtogroup op_flip
//! @{



template<typename T1>
inline
void
op_flipud::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_flipud>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_Mat<T1>::value)
    {
    // allow detection of in-place operation
    
    const unwrap<T1> U(in.m);
    
    op_flipud::apply_direct(out, U.M);
    }
  else
    {
    const Proxy<T1> P(in.m);
    
    if(P.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_flipud::apply_proxy_noalias(tmp, P);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_flipud::apply_proxy_noalias(out, P);
      }
    }
  }



template<typename eT>
inline
void
op_flipud::apply_direct(Mat<eT>& out, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword X_n_rows_m1 = X_n_rows - 1;
  
  if(&out != &X)
    {
    out.set_size(X_n_rows, X_n_cols);
    
    if(X_n_cols == 1)
      {
      const eT*   X_mem =   X.memptr();
            eT* out_mem = out.memptr();
      
      for(uword row=0; row < X_n_rows; ++row)
        {
        out_mem[X_n_rows_m1 - row] = X_mem[row];
        }
      }
    else
      {
      for(uword col=0; col < X_n_cols; ++col)
        {
        const eT*   X_colmem =   X.colptr(col);
              eT* out_colmem = out.colptr(col);
        
        for(uword row=0; row < X_n_rows; ++row)
          {
          out_colmem[X_n_rows_m1 - row] = X_colmem[row];
          }
        }
      }
    }
  else  // in-place operation
    {
    const uword N = X_n_rows / 2;
    
    if(X_n_cols == 1)
      {
      eT* out_mem = out.memptr();
      
      for(uword row=0; row < N; ++row)
        {
        std::swap(out_mem[X_n_rows_m1 - row], out_mem[row]);
        }
      }
    else
      {
      for(uword col=0; col < X_n_cols; ++col)
        {
        eT* out_colmem = out.colptr(col);
        
        for(uword row=0; row < N; ++row)
          {
          std::swap(out_colmem[X_n_rows_m1 - row], out_colmem[row]);
          }
        }
      }
    }
  }



template<typename T1>
inline
void
op_flipud::apply_proxy_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  typedef typename Proxy<T1>::stored_type P_stored_type;
  
  if(is_Mat<P_stored_type>::value)
    {
    const unwrap<P_stored_type> U(P.Q);
    
    op_flipud::apply_direct(out, U.M);
    
    return;
    }
  
  const uword P_n_rows = P.get_n_rows();
  const uword P_n_cols = P.get_n_cols();
  
  const uword P_n_rows_m1 = P_n_rows - 1;
  
  out.set_size(P_n_rows, P_n_cols);
  
  if( ((T1::is_col) || (P_n_cols == 1)) && (Proxy<T1>::use_at == false) )
    {
    eT* out_mem = out.memptr();
    
    const typename Proxy<T1>::ea_type P_ea = P.get_ea();
    
    for(uword row=0; row < P_n_rows; ++row)
      {
      out_mem[P_n_rows_m1 - row] = P_ea[row];
      }
    }
  else
    {
    for(uword col=0; col < P_n_cols; ++col)
      {
      eT* out_colmem = out.colptr(col);
      
      for(uword row=0; row < P_n_rows; ++row)
        {
        out_colmem[P_n_rows_m1 - row] = P.at(row, col);
        }
      }
    }
  }



//



template<typename T1>
inline
void
op_fliplr::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_fliplr>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_Mat<T1>::value)
    {
    // allow detection of in-place operation
    
    const unwrap<T1> U(in.m);
    
    op_fliplr::apply_direct(out, U.M);
    }
  else
    {
    const Proxy<T1> P(in.m);
    
    if(P.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_fliplr::apply_proxy_noalias(tmp, P);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_fliplr::apply_proxy_noalias(out, P);
      }
    }
  }



template<typename eT>
inline
void
op_fliplr::apply_direct(Mat<eT>& out, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword X_n_cols_m1 = X_n_cols - 1;
  
  if(&out != &X)
    {
    out.set_size(X_n_rows, X_n_cols);
    
    if(X_n_rows == 1)
      {
      const eT*   X_mem =   X.memptr();
            eT* out_mem = out.memptr();
      
      for(uword col=0; col < X_n_cols; ++col)
        {
        out_mem[X_n_cols_m1 - col] = X_mem[col];
        }
      }
    else
      {
      for(uword col=0; col < X_n_cols; ++col)
        {
        out.col(X_n_cols_m1 - col) = X.col(col);
        }
      }
    }
  else  // in-place operation
    {
    const uword N = X_n_cols / 2;
    
    if(X_n_rows == 1)
      {
      eT* out_mem = out.memptr();
      
      for(uword col=0; col < N; ++col)
        {
        std::swap(out_mem[X_n_cols_m1 - col], out_mem[col]);
        }
      }
    else
      {
      for(uword col=0; col < N; ++col)
        {
        out.swap_cols(X_n_cols_m1 - col, col);
        }
      }
    }
  }



template<typename T1>
inline
void
op_fliplr::apply_proxy_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  typedef typename Proxy<T1>::stored_type P_stored_type;
  
  if(is_Mat<P_stored_type>::value)
    {
    const unwrap<P_stored_type> U(P.Q);
    
    op_fliplr::apply_direct(out, U.M);
    
    return;
    }
  
  const uword P_n_rows = P.get_n_rows();
  const uword P_n_cols = P.get_n_cols();
  
  const uword P_n_cols_m1 = P_n_cols - 1;
  
  out.set_size(P_n_rows, P_n_cols);
  
  if( ((T1::is_row) || (P_n_rows == 1)) && (Proxy<T1>::use_at == false) )
    {
    eT* out_mem = out.memptr();
    
    const typename Proxy<T1>::ea_type P_ea = P.get_ea();
    
    for(uword col=0; col < P_n_cols; ++col)
      {
      out_mem[P_n_cols_m1 - col] = P_ea[col];
      }
    }
  else
    {
    for(uword col=0; col < P_n_cols; ++col)
      {
      eT* out_colmem = out.colptr(P_n_cols_m1 - col);
      
      for(uword row=0; row < P_n_rows; ++row)
        {
        out_colmem[row] = P.at(row,col);
        }
      }
    }
  }



//! @}
