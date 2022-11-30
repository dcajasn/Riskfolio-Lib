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


//! \addtogroup op_htrans
//! @{



template<typename eT>
arma_hot
inline
void
op_htrans::apply_mat_noalias(Mat<eT>& out, const Mat<eT>& A, const typename arma_not_cx<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_strans::apply_mat_noalias(out, A);
  }



template<typename eT>
arma_hot
inline
void
op_htrans::apply_mat_noalias(Mat<eT>& out, const Mat<eT>& A, const typename arma_cx_only<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  out.set_size(A_n_cols, A_n_rows);
  
  if( (A_n_cols == 1) || (A_n_rows == 1) )
    {
    const uword n_elem = A.n_elem;
    
    const eT* A_mem   = A.memptr();
          eT* out_mem = out.memptr();
    
    for(uword i=0; i < n_elem; ++i)
      {
      out_mem[i] = std::conj(A_mem[i]);
      }
    }
  else
  if( (A_n_rows >= 512) && (A_n_cols >= 512) )
    {
    op_htrans::apply_mat_noalias_large(out, A);
    }
  else
    {
    eT* outptr = out.memptr();
    
    for(uword k=0; k < A_n_rows; ++k)
      {
      const eT* Aptr = &(A.at(k,0));
      
      for(uword j=0; j < A_n_cols; ++j)
        {
        (*outptr) = std::conj(*Aptr);
        
        Aptr += A_n_rows;
        outptr++;
        }
      }
    }
  }



template<typename T>
arma_hot
inline
void
op_htrans::block_worker(std::complex<T>* Y, const std::complex<T>* X, const uword X_n_rows, const uword Y_n_rows, const uword n_rows, const uword n_cols)
  {
  for(uword row = 0; row < n_rows; ++row)
    {
    const uword Y_offset = row * Y_n_rows;
    
    for(uword col = 0; col < n_cols; ++col)
      {
      const uword X_offset = col * X_n_rows;
      
      Y[col + Y_offset] = std::conj(X[row + X_offset]);
      }
    }
  }



template<typename T>
arma_hot
inline
void
op_htrans::apply_mat_noalias_large(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = A.n_rows;
  const uword n_cols = A.n_cols;
  
  const uword block_size = 64;
  
  const uword n_rows_base = block_size * (n_rows / block_size);
  const uword n_cols_base = block_size * (n_cols / block_size);
  
  const uword n_rows_extra = n_rows - n_rows_base;
  const uword n_cols_extra = n_cols - n_cols_base;
  
  const std::complex<T>* X =   A.memptr();
        std::complex<T>* Y = out.memptr();
  
  for(uword row = 0; row < n_rows_base; row += block_size)
    {
    const uword Y_offset = row * n_cols;
    
    for(uword col = 0; col < n_cols_base; col += block_size)
      {
      const uword X_offset = col * n_rows;
      
      op_htrans::block_worker(&Y[col + Y_offset], &X[row + X_offset], n_rows, n_cols, block_size, block_size);
      }
    
    const uword X_offset = n_cols_base * n_rows;
    
    op_htrans::block_worker(&Y[n_cols_base + Y_offset], &X[row + X_offset], n_rows, n_cols, block_size, n_cols_extra);
    }

  if(n_rows_extra == 0)  { return; }
  
  const uword Y_offset = n_rows_base * n_cols;
  
  for(uword col = 0; col < n_cols_base; col += block_size)
    {
    const uword X_offset = col * n_rows;
    
    op_htrans::block_worker(&Y[col + Y_offset], &X[n_rows_base + X_offset], n_rows, n_cols, n_rows_extra, block_size);
    }
  
  const uword X_offset = n_cols_base * n_rows;
  
  op_htrans::block_worker(&Y[n_cols_base + Y_offset], &X[n_rows_base + X_offset], n_rows, n_cols, n_rows_extra, n_cols_extra);
  }



template<typename eT>
arma_hot
inline
void
op_htrans::apply_mat_inplace(Mat<eT>& out, const typename arma_not_cx<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_strans::apply_mat_inplace(out);
  }



template<typename eT>
arma_hot
inline
void
op_htrans::apply_mat_inplace(Mat<eT>& out, const typename arma_cx_only<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword n_rows = out.n_rows;
  const uword n_cols = out.n_cols;
    
  if(n_rows == n_cols)
    {
    arma_extra_debug_print("doing in-place hermitian transpose of a square matrix");
    
    for(uword col=0; col < n_cols; ++col)
      {
      eT* coldata = out.colptr(col);
      
      out.at(col,col) = std::conj( out.at(col,col) );
      
      for(uword row=(col+1); row < n_rows; ++row)
        {
        const eT val1 = std::conj(coldata[row]);
        const eT val2 = std::conj(out.at(col,row));
        
        out.at(col,row) = val1;
        coldata[row]    = val2;
        }
      }
    }
  else
    {
    Mat<eT> tmp;
    
    op_htrans::apply_mat_noalias(tmp, out);
    
    out.steal_mem(tmp);
    }
  }



template<typename eT>
inline
void
op_htrans::apply_mat(Mat<eT>& out, const Mat<eT>& A, const typename arma_not_cx<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_strans::apply_mat(out, A);
  }



template<typename eT>
inline
void
op_htrans::apply_mat(Mat<eT>& out, const Mat<eT>& A, const typename arma_cx_only<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  if(&out != &A)
    {
    op_htrans::apply_mat_noalias(out, A);
    }
  else
    {
    op_htrans::apply_mat_inplace(out);
    }
  }



template<typename T1>
inline
void
op_htrans::apply_proxy(Mat<typename T1::elem_type>& out, const Proxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  if( (resolves_to_vector<T1>::yes) && (Proxy<T1>::use_at == false) )
    {
    out.set_size(n_cols, n_rows);
    
    eT* out_mem = out.memptr();
    
    const uword n_elem = P.get_n_elem();
    
    typename Proxy<T1>::ea_type Pea = P.get_ea();
    
    for(uword i=0; i < n_elem; ++i)
      {
      out_mem[i] = std::conj(Pea[i]);
      }
    }
  else
    {
    out.set_size(n_cols, n_rows);
    
    eT* outptr = out.memptr();
    
    for(uword k=0; k < n_rows; ++k)
      {
      for(uword j=0; j < n_cols; ++j)
        {
        (*outptr) = std::conj(P.at(k,j));
        
        outptr++;
        }
      }
    }
  }



template<typename T1>
inline
void
op_htrans::apply_direct(Mat<typename T1::elem_type>& out, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // allow detection of in-place transpose
  if(is_Mat<T1>::value || (arma_config::openmp && Proxy<T1>::use_mp))
    {
    const unwrap<T1> U(X);
    
    op_htrans::apply_mat(out, U.M);
    }
  else
    {
    const Proxy<T1> P(X);
    
    const bool is_alias = P.is_alias(out);
    
    if(is_Mat<typename Proxy<T1>::stored_type>::value)
      {
      const quasi_unwrap<typename Proxy<T1>::stored_type> U(P.Q);
      
      if(is_alias)
        {
        Mat<eT> tmp;
        
        op_htrans::apply_mat_noalias(tmp, U.M);
        
        out.steal_mem(tmp);
        }
      else
        {
        op_htrans::apply_mat_noalias(out, U.M);
        }
      }
    else
      {
      if(is_alias)
        {
        Mat<eT> tmp;
        
        op_htrans::apply_proxy(tmp, P);
        
        out.steal_mem(tmp);
        }
      else
        {
        op_htrans::apply_proxy(out, P);
        }
      }
    }
  }



template<typename T1>
inline
void
op_htrans::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans>& in, const typename arma_not_cx<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_strans::apply_direct(out, in.m);
  }



template<typename T1>
inline
void
op_htrans::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans>& in, const typename arma_cx_only<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_htrans::apply_direct(out, in.m);
  }



//
// op_htrans2



template<typename T1>
inline
void
op_htrans2::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans2>& in, const typename arma_not_cx<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_strans::apply_direct(out, in.m);
  
  arrayops::inplace_mul(out.memptr(), in.aux, out.n_elem);
  }



template<typename T1>
inline
void
op_htrans2::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans2>& in, const typename arma_cx_only<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  op_htrans::apply_direct(out, in.m);
  
  arrayops::inplace_mul(out.memptr(), in.aux, out.n_elem);
  }



//! @}
