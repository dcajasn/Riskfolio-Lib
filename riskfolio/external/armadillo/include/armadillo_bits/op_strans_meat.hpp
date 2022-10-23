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


//! \addtogroup op_strans
//! @{



//! for tiny square matrices (size <= 4x4)
template<typename eT, typename TA>
arma_cold
inline
void
op_strans::apply_mat_noalias_tinysq(Mat<eT>& out, const TA& A)
  {
  const eT*   Am =   A.memptr();
        eT* outm = out.memptr();
  
  switch(A.n_rows)
    {
    case 1:
      {
      outm[0] = Am[0];
      }
      break;
      
    case 2:
      {
      outm[pos<false,0,0>::n2] = Am[pos<true,0,0>::n2];
      outm[pos<false,1,0>::n2] = Am[pos<true,1,0>::n2];
      
      outm[pos<false,0,1>::n2] = Am[pos<true,0,1>::n2];
      outm[pos<false,1,1>::n2] = Am[pos<true,1,1>::n2];
      }
      break;
    
    case 3:
      {
      outm[pos<false,0,0>::n3] = Am[pos<true,0,0>::n3];
      outm[pos<false,1,0>::n3] = Am[pos<true,1,0>::n3];
      outm[pos<false,2,0>::n3] = Am[pos<true,2,0>::n3];
      
      outm[pos<false,0,1>::n3] = Am[pos<true,0,1>::n3];
      outm[pos<false,1,1>::n3] = Am[pos<true,1,1>::n3];
      outm[pos<false,2,1>::n3] = Am[pos<true,2,1>::n3];
      
      outm[pos<false,0,2>::n3] = Am[pos<true,0,2>::n3];
      outm[pos<false,1,2>::n3] = Am[pos<true,1,2>::n3];
      outm[pos<false,2,2>::n3] = Am[pos<true,2,2>::n3];
      }
      break;
    
    case 4:
      {
      outm[pos<false,0,0>::n4] = Am[pos<true,0,0>::n4];
      outm[pos<false,1,0>::n4] = Am[pos<true,1,0>::n4];
      outm[pos<false,2,0>::n4] = Am[pos<true,2,0>::n4];
      outm[pos<false,3,0>::n4] = Am[pos<true,3,0>::n4];
      
      outm[pos<false,0,1>::n4] = Am[pos<true,0,1>::n4];
      outm[pos<false,1,1>::n4] = Am[pos<true,1,1>::n4];
      outm[pos<false,2,1>::n4] = Am[pos<true,2,1>::n4];
      outm[pos<false,3,1>::n4] = Am[pos<true,3,1>::n4];
      
      outm[pos<false,0,2>::n4] = Am[pos<true,0,2>::n4];
      outm[pos<false,1,2>::n4] = Am[pos<true,1,2>::n4];
      outm[pos<false,2,2>::n4] = Am[pos<true,2,2>::n4];
      outm[pos<false,3,2>::n4] = Am[pos<true,3,2>::n4];
      
      outm[pos<false,0,3>::n4] = Am[pos<true,0,3>::n4];
      outm[pos<false,1,3>::n4] = Am[pos<true,1,3>::n4];
      outm[pos<false,2,3>::n4] = Am[pos<true,2,3>::n4];
      outm[pos<false,3,3>::n4] = Am[pos<true,3,3>::n4];
      }
      break;
    
    default:
      ;
    }
  
  }



template<typename eT>
arma_hot
inline
void
op_strans::block_worker(eT* Y, const eT* X, const uword X_n_rows, const uword Y_n_rows, const uword n_rows, const uword n_cols)
  {
  for(uword row = 0; row < n_rows; ++row)
    {
    const uword Y_offset = row * Y_n_rows;
    
    for(uword col = 0; col < n_cols; ++col)
      {
      const uword X_offset = col * X_n_rows;
      
      Y[col + Y_offset] = X[row + X_offset];
      }
    }
  }



template<typename eT>
arma_hot
inline
void
op_strans::apply_mat_noalias_large(Mat<eT>& out, const Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = A.n_rows;
  const uword n_cols = A.n_cols;
  
  const uword block_size = 64;
  
  const uword n_rows_base = block_size * (n_rows / block_size);
  const uword n_cols_base = block_size * (n_cols / block_size);
  
  const uword n_rows_extra = n_rows - n_rows_base;
  const uword n_cols_extra = n_cols - n_cols_base;
  
  const eT* X =   A.memptr();
        eT* Y = out.memptr();
  
  for(uword row = 0; row < n_rows_base; row += block_size)
    {
    const uword Y_offset = row * n_cols;
    
    for(uword col = 0; col < n_cols_base; col += block_size)
      {
      const uword X_offset = col * n_rows;
      
      op_strans::block_worker(&Y[col + Y_offset], &X[row + X_offset], n_rows, n_cols, block_size, block_size);
      }
    
    const uword X_offset = n_cols_base * n_rows;
    
    op_strans::block_worker(&Y[n_cols_base + Y_offset], &X[row + X_offset], n_rows, n_cols, block_size, n_cols_extra);
    }

  if(n_rows_extra == 0)  { return; }
  
  const uword Y_offset = n_rows_base * n_cols;
  
  for(uword col = 0; col < n_cols_base; col += block_size)
    {
    const uword X_offset = col * n_rows;
    
    op_strans::block_worker(&Y[col + Y_offset], &X[n_rows_base + X_offset], n_rows, n_cols, n_rows_extra, block_size);
    }
  
  const uword X_offset = n_cols_base * n_rows;
  
  op_strans::block_worker(&Y[n_cols_base + Y_offset], &X[n_rows_base + X_offset], n_rows, n_cols, n_rows_extra, n_cols_extra);
  }



//! Immediate transpose of a dense matrix
template<typename eT, typename TA>
arma_hot
inline
void
op_strans::apply_mat_noalias(Mat<eT>& out, const TA& A)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_cols = A.n_cols;
  const uword A_n_rows = A.n_rows;
  
  out.set_size(A_n_cols, A_n_rows);
  
  if( (TA::is_row) || (TA::is_col) || (A_n_cols == 1) || (A_n_rows == 1) )
    {
    arrayops::copy( out.memptr(), A.memptr(), A.n_elem );
    }
  else
    {
    if( (A_n_rows <= 4) && (A_n_rows == A_n_cols) )
      {
      op_strans::apply_mat_noalias_tinysq(out, A);
      }
    else
    if( (A_n_rows >= 512) && (A_n_cols >= 512) )
      {
      op_strans::apply_mat_noalias_large(out, A);
      }
    else
      {
      eT* outptr = out.memptr();
      
      for(uword k=0; k < A_n_rows; ++k)
        {
        const eT* Aptr = &(A.at(k,0));
        
        uword j;
        for(j=1; j < A_n_cols; j+=2)
          {
          const eT tmp_i = (*Aptr);  Aptr += A_n_rows;
          const eT tmp_j = (*Aptr);  Aptr += A_n_rows;
          
          (*outptr) = tmp_i;  outptr++;
          (*outptr) = tmp_j;  outptr++;
          }
        
        if((j-1) < A_n_cols)
          {
          (*outptr) = (*Aptr);  outptr++;;
          }
        }
      }
    }
  }



template<typename eT>
arma_hot
inline
void
op_strans::apply_mat_inplace(Mat<eT>& out)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = out.n_rows;
  const uword n_cols = out.n_cols;
  
  if(n_rows == n_cols)
    {
    arma_extra_debug_print("op_strans::apply(): doing in-place transpose of a square matrix");
    
    const uword N = n_rows;
    
    for(uword k=0; k < N; ++k)
      {
      eT* colptr = &(out.at(k,k));
      eT* rowptr = colptr;
      
      colptr++;
      rowptr += N;
      
      uword j;
      
      for(j=(k+2); j < N; j+=2)
        {
        std::swap( (*rowptr), (*colptr) );  rowptr += N;  colptr++;
        std::swap( (*rowptr), (*colptr) );  rowptr += N;  colptr++;
        }
      
      if((j-1) < N)
        {
        std::swap( (*rowptr), (*colptr) );
        }
      }
    }
  else
    {
    Mat<eT> tmp;
    
    op_strans::apply_mat_noalias(tmp, out);
    
    out.steal_mem(tmp);
    }
  }



template<typename eT, typename TA>
inline
void
op_strans::apply_mat(Mat<eT>& out, const TA& A)
  {
  arma_extra_debug_sigprint();
  
  if(&out != &A)
    {
    op_strans::apply_mat_noalias(out, A);
    }
  else
    {
    op_strans::apply_mat_inplace(out);
    }
  }



template<typename T1>
inline
void
op_strans::apply_proxy(Mat<typename T1::elem_type>& out, const Proxy<T1>& P)
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
    
    uword i,j;
    for(i=0, j=1; j < n_elem; i+=2, j+=2)
      {
      const eT tmp_i = Pea[i];
      const eT tmp_j = Pea[j];
      
      out_mem[i] = tmp_i;
      out_mem[j] = tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] = Pea[i];
      }
    }
  else   // general matrix transpose
    {
    out.set_size(n_cols, n_rows);
    
    eT* outptr = out.memptr();
    
    for(uword k=0; k < n_rows; ++k)
      {
      uword j;
      for(j=1; j < n_cols; j+=2)
        {
        const uword i = j-1;
        
        const eT tmp_i = P.at(k,i);
        const eT tmp_j = P.at(k,j);
        
        (*outptr) = tmp_i;  outptr++;
        (*outptr) = tmp_j;  outptr++;
        }
      
      const uword i = j-1;
      
      if(i < n_cols)
        {
        (*outptr) = P.at(k,i);  outptr++;
        }
      }
    }
  }



template<typename T1>
inline
void
op_strans::apply_direct(Mat<typename T1::elem_type>& out, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // allow detection of in-place transpose
  if(is_Mat<T1>::value || (arma_config::openmp && Proxy<T1>::use_mp))
    {
    const unwrap<T1> U(X);
    
    op_strans::apply_mat(out, U.M);
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
        
        op_strans::apply_mat_noalias(tmp, U.M);
        
        out.steal_mem(tmp);
        }
      else
        {
        op_strans::apply_mat_noalias(out, U.M);
        }
      }
    else
      {
      if(is_alias)
        {
        Mat<eT> tmp;
        
        op_strans::apply_proxy(tmp, P);
        
        out.steal_mem(tmp);
        }
      else
        {
        op_strans::apply_proxy(out, P);
        }
      }
    }
  }



template<typename T1>
inline
void
op_strans::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_strans>& in)
  {
  arma_extra_debug_sigprint();
  
  op_strans::apply_direct(out, in.m);
  }



//
//
//



template<typename eT>
inline
void
op_strans_cube::apply_noalias(Cube<eT>& out, const Cube<eT>& X)
  {
  out.set_size(X.n_cols, X.n_rows, X.n_slices);
  
  for(uword s=0; s < X.n_slices; ++s)
    {
    Mat<eT> out_slice( out.slice_memptr(s), X.n_cols, X.n_rows, false, true );
    
    const Mat<eT> X_slice( const_cast<eT*>(X.slice_memptr(s)), X.n_rows, X.n_cols, false, true );
    
    op_strans::apply_mat_noalias(out_slice, X_slice);
    }
  }



//! @}
