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


//! \addtogroup op_trimat
//! @{



template<typename eT>
inline
void
op_trimat::fill_zeros(Mat<eT>& out, const bool upper)
  {
  arma_extra_debug_sigprint();
  
  const uword N = out.n_rows;
  
  if(upper)
    {
    // upper triangular: set all elements below the diagonal to zero
    
    for(uword i=0; i<N; ++i)
      {
      eT* data = out.colptr(i);
      
      arrayops::fill_zeros( &data[i+1], (N-(i+1)) );
      }
    }
  else
    {
    // lower triangular: set all elements above the diagonal to zero
    
    for(uword i=1; i<N; ++i)
      {
      eT* data = out.colptr(i);
      
      arrayops::fill_zeros( data, i );
      }
    }
  }



template<typename T1>
inline
void
op_trimat::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_trimat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const bool upper = (in.aux_uword_a == 0);
  
  // allow detection of in-place operation
  if(is_Mat<T1>::value || (arma_config::openmp && Proxy<T1>::use_mp))
    {
    const unwrap<T1> U(in.m);
    
    op_trimat::apply_unwrap(out, U.M, upper);
    }
  else
    {
    const Proxy<T1> P(in.m);
    
    const bool is_alias = P.is_alias(out);
    
    if(is_Mat<typename Proxy<T1>::stored_type>::value)
      {
      const quasi_unwrap<typename Proxy<T1>::stored_type> U(P.Q);
      
      if(is_alias)
        {
        Mat<eT> tmp;
        
        op_trimat::apply_unwrap(tmp, U.M, upper);
        
        out.steal_mem(tmp);
        }
      else
        {
        op_trimat::apply_unwrap(out, U.M, upper);
        }
      }
    else
      {
      if(is_alias)
        {
        Mat<eT> tmp;
        
        op_trimat::apply_proxy(tmp, P, upper);
        
        out.steal_mem(tmp);
        }
      else
        {
        op_trimat::apply_proxy(out, P, upper);
        }
      }
    }
  }



template<typename eT>
inline
void
op_trimat::apply_unwrap(Mat<eT>& out, const Mat<eT>& A, const bool upper)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.is_square() == false), "trimatu()/trimatl(): given matrix must be square sized" );
  
  if(&out != &A)
    {
    out.copy_size(A);
    
    const uword N = A.n_rows;
    
    if(upper)
      {
      // upper triangular: copy the diagonal and the elements above the diagonal
      for(uword i=0; i<N; ++i)
        {
        const eT* A_data   = A.colptr(i);
              eT* out_data = out.colptr(i);
        
        arrayops::copy( out_data, A_data, i+1 );
        }
      }
    else
      {
      // lower triangular: copy the diagonal and the elements below the diagonal
      for(uword i=0; i<N; ++i)
        {
        const eT* A_data   = A.colptr(i);
              eT* out_data = out.colptr(i);
        
        arrayops::copy( &out_data[i], &A_data[i], N-i );
        }
      }
    }
  
  op_trimat::fill_zeros(out, upper);
  }



template<typename T1>
inline
void
op_trimat::apply_proxy(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const bool upper)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (P.get_n_rows() != P.get_n_cols()), "trimatu()/trimatl(): given matrix must be square sized" );
  
  const uword N = P.get_n_rows();
  
  out.set_size(N,N);
  
  if(upper)
    {
    for(uword j=0; j < N;     ++j)
    for(uword i=0; i < (j+1); ++i)
      {
      out.at(i,j) = P.at(i,j);
      }
    }
  else
    {
    for(uword j=0; j<N; ++j)
    for(uword i=j; i<N; ++i)
      {
      out.at(i,j) = P.at(i,j);
      }
    }
  
  op_trimat::fill_zeros(out, upper);
  }



//



template<typename T1>
inline
void
op_trimatu_ext::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_trimatu_ext>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( (A.is_square() == false), "trimatu(): given matrix must be square sized" );
  
  const uword row_offset = in.aux_uword_a;
  const uword col_offset = in.aux_uword_b;
  
  const uword n_rows = A.n_rows;
  const uword n_cols = A.n_cols;
  
  arma_debug_check_bounds( ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)), "trimatu(): requested diagonal is out of bounds" );
  
  if(&out != &A)
    {
    out.copy_size(A);
    
    const uword N = (std::min)(n_rows - row_offset, n_cols - col_offset);
    
    for(uword i=0; i < n_cols; ++i)
      {
      const uword col = i + col_offset;
        
      if(i < N)
        {
        const uword end_row = i + row_offset;
        
        for(uword row=0; row <= end_row; ++row)
          {
          out.at(row,col) = A.at(row,col);
          }
        }
      else
        {
        if(col < n_cols)
          {
          arrayops::copy(out.colptr(col), A.colptr(col), n_rows);
          }
        }
      }
    }
  
  op_trimatu_ext::fill_zeros(out, row_offset, col_offset);
  }



template<typename eT>
inline
void
op_trimatu_ext::fill_zeros(Mat<eT>& out, const uword row_offset, const uword col_offset)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = out.n_rows;
  const uword n_cols = out.n_cols;
  
  const uword N = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  for(uword col=0; col < col_offset; ++col)
    {
    arrayops::fill_zeros(out.colptr(col), n_rows);
    }
  
  for(uword i=0; i < N; ++i)
    {
    const uword start_row = i + row_offset + 1;
    const uword col       = i + col_offset;
    
    for(uword row=start_row; row < n_rows; ++row)
      {
      out.at(row,col) = eT(0);
      }
    }
  }



//



template<typename T1>
inline
void
op_trimatl_ext::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_trimatl_ext>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( (A.is_square() == false), "trimatl(): given matrix must be square sized" );
  
  const uword row_offset = in.aux_uword_a;
  const uword col_offset = in.aux_uword_b;
  
  const uword n_rows = A.n_rows;
  const uword n_cols = A.n_cols;
  
  arma_debug_check_bounds( ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)), "trimatl(): requested diagonal is out of bounds" );
  
  if(&out != &A)
    {
    out.copy_size(A);
    
    const uword N = (std::min)(n_rows - row_offset, n_cols - col_offset);
    
    for(uword col=0; col < col_offset; ++col)
      {
      arrayops::copy( out.colptr(col), A.colptr(col), n_rows );
      }
    
    for(uword i=0; i<N; ++i)
      {
      const uword start_row = i + row_offset;
      const uword       col = i + col_offset;
      
      for(uword row=start_row; row < n_rows; ++row)
        {
        out.at(row,col) = A.at(row,col);
        }
      }
    }
  
  op_trimatl_ext::fill_zeros(out, row_offset, col_offset);
  }



template<typename eT>
inline
void
op_trimatl_ext::fill_zeros(Mat<eT>& out, const uword row_offset, const uword col_offset)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = out.n_rows;
  const uword n_cols = out.n_cols;
  
  const uword N = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  for(uword i=0; i < n_cols; ++i)
    {
    const uword col = i + col_offset;
      
    if(i < N)
      {
      const uword end_row = i + row_offset;
      
      for(uword row=0; row < end_row; ++row)
        {
        out.at(row,col) = eT(0);
        }
      }
    else
      {
      if(col < n_cols)
        {
        arrayops::fill_zeros(out.colptr(col), n_rows);
        }
      }
    }
  }



//! @}
