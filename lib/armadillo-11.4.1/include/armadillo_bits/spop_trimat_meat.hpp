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


//! \addtogroup spop_trimat
//! @{



template<typename T1>
inline
void
spop_trimat::apply_noalias(SpMat<typename T1::elem_type>& out, const SpProxy<T1>& P, const bool upper)
  {
  arma_extra_debug_sigprint();
  
  typename SpProxy<T1>::const_iterator_type it = P.begin();
  
  const uword old_n_nonzero = P.get_n_nonzero();
        uword new_n_nonzero = 0;
  
  if(upper)
    {
    // upper triangular: count elements on the diagonal and above the diagonal
    
    for(uword i=0; i < old_n_nonzero; ++i)
      {
      new_n_nonzero += (it.row() <= it.col()) ? uword(1) : uword(0);
      ++it;
      }
    }
  else
    {
    // lower triangular: count elements on the diagonal and below the diagonal
    
    for(uword i=0; i < old_n_nonzero; ++i)
      {
      new_n_nonzero += (it.row() >= it.col()) ? uword(1) : uword(0);
      ++it;
      }
    }
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();  
  
  out.reserve(n_rows, n_cols, new_n_nonzero);
  
  uword new_index = 0;
  
  it = P.begin();
  
  if(upper)
    {
    // upper triangular: copy elements on the diagonal and above the diagonal
    
    for(uword i=0; i < old_n_nonzero; ++i)
      {
      const uword row = it.row();
      const uword col = it.col();
      
      if(row <= col)
        {
        access::rw(out.values[new_index])      = (*it);
        access::rw(out.row_indices[new_index]) = row;
        access::rw(out.col_ptrs[col + 1])++;
        ++new_index;
        }
      
      ++it;
      }
    }
  else
    {
    // lower triangular: copy elements on the diagonal and below the diagonal
    
    for(uword i=0; i < old_n_nonzero; ++i)
      {
      const uword row = it.row();
      const uword col = it.col();
      
      if(row >= col)
        {
        access::rw(out.values[new_index])      = (*it);
        access::rw(out.row_indices[new_index]) = row;
        access::rw(out.col_ptrs[col + 1])++;
        ++new_index;
        }
      
      ++it;
      }
    }
  
  for(uword i=0; i < n_cols; ++i)
    {
    access::rw(out.col_ptrs[i + 1]) += out.col_ptrs[i];
    }
  }



template<typename T1>
inline
void
spop_trimat::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_trimat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> P(in.m);
  
  arma_debug_check( (P.get_n_rows() != P.get_n_cols()), "trimatu()/trimatl(): given matrix must be square sized" );
  
  const bool upper = (in.aux_uword_a == 0);
  
  if(P.is_alias(out))
    {
    SpMat<eT> tmp;
    spop_trimat::apply_noalias(tmp, P, upper);
    out.steal_mem(tmp);
    }
  else
    {
    spop_trimat::apply_noalias(out, P, upper);
    }
  }



//



template<typename eT>
inline
void
spop_trimatu_ext::apply_noalias(SpMat<eT>& out, const SpMat<eT>& A, const uword row_offset, const uword col_offset)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = A.n_rows;
  const uword n_cols = A.n_cols;
  
  arma_debug_check_bounds( ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)), "trimatu(): requested diagonal is out of bounds" );
  
  if(A.n_nonzero == 0)  { out.zeros(n_rows, n_cols); return; }
  
  out.reserve(n_rows, n_cols, A.n_nonzero);  // upper bound on n_nonzero
  
  uword count = 0;
  
  const uword N = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  for(uword i=0; i < n_cols; ++i)
    {
    const uword col = i + col_offset;
    
    if(i < N)
      {
      typename SpMat<eT>::const_col_iterator it     = A.begin_col_no_sync(col);
      typename SpMat<eT>::const_col_iterator it_end = A.end_col_no_sync(col);
      
      const uword end_row = i + row_offset;
      
      for(; it != it_end; ++it)
        {
        const uword it_row = it.row();
        
        if(it_row <= end_row)
          {
          const uword it_col = it.col();
          
          access::rw(out.values[count])      = (*it);
          access::rw(out.row_indices[count]) = it_row;
          access::rw(out.col_ptrs[it_col + 1])++;
          ++count;
          }
        else
          {
          break;
          }
        }
      }
    else
      {
      if(col < n_cols)
        {
        typename SpMat<eT>::const_col_iterator it     = A.begin_col_no_sync(col);
        typename SpMat<eT>::const_col_iterator it_end = A.end_col_no_sync(col);
        
        for(; it != it_end; ++it)
          {
          const uword it_row = it.row();
          const uword it_col = it.col();
          
          access::rw(out.values[count])      = (*it);
          access::rw(out.row_indices[count]) = it_row;
          access::rw(out.col_ptrs[it_col + 1])++;
          ++count;
          }
        }
      }
    }
  
  for(uword i=0; i < n_cols; ++i)
    {
    access::rw(out.col_ptrs[i + 1]) += out.col_ptrs[i];
    }
  
  if(count < A.n_nonzero)  { out.mem_resize(count); }
  }



template<typename T1>
inline
void
spop_trimatu_ext::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_trimatu_ext>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> U(in.m);
  const SpMat<eT>&   A = U.M;
  
  arma_debug_check( (A.is_square() == false), "trimatu(): given matrix must be square sized" );
  
  const uword row_offset = in.aux_uword_a;
  const uword col_offset = in.aux_uword_b;
  
  if(U.is_alias(out))
    {
    SpMat<eT> tmp;
    spop_trimatu_ext::apply_noalias(tmp, A, row_offset, col_offset);
    out.steal_mem(tmp);
    }
  else
    {
    spop_trimatu_ext::apply_noalias(out, A, row_offset, col_offset);
    }
  }



//



template<typename eT>
inline
void
spop_trimatl_ext::apply_noalias(SpMat<eT>& out, const SpMat<eT>& A, const uword row_offset, const uword col_offset)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = A.n_rows;
  const uword n_cols = A.n_cols;
  
  arma_debug_check_bounds( ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)), "trimatl(): requested diagonal is out of bounds" );
  
  if(A.n_nonzero == 0)  { out.zeros(n_rows, n_cols); return; }
  
  out.reserve(n_rows, n_cols, A.n_nonzero);  // upper bound on n_nonzero
  
  uword count = 0;
  
  if(col_offset > 0)
    {
    typename SpMat<eT>::const_col_iterator it     = A.begin_col_no_sync(0);
    typename SpMat<eT>::const_col_iterator it_end = A.end_col_no_sync(col_offset-1);
    
    for(; it != it_end; ++it)
      {
      const uword it_row = it.row();
      const uword it_col = it.col();
      
      access::rw(out.values[count])      = (*it);
      access::rw(out.row_indices[count]) = it_row;
      access::rw(out.col_ptrs[it_col + 1])++;
      ++count;
      }
    }
  
  const uword N = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  for(uword i=0; i < N; ++i)
    {
    const uword start_row = i + row_offset;
    const uword       col = i + col_offset;
    
    typename SpMat<eT>::const_col_iterator it     = A.begin_col_no_sync(col);
    typename SpMat<eT>::const_col_iterator it_end = A.end_col_no_sync(col);
    
    for(; it != it_end; ++it)
      {
      const uword it_row = it.row();
      
      if(it_row >= start_row)
        {
        const uword it_col = it.col();
        
        access::rw(out.values[count])      = (*it);
        access::rw(out.row_indices[count]) = it_row;
        access::rw(out.col_ptrs[it_col + 1])++;
        ++count;
        }
      }
    }
  
  for(uword i=0; i < n_cols; ++i)
    {
    access::rw(out.col_ptrs[i + 1]) += out.col_ptrs[i];
    }
  
  if(count < A.n_nonzero)  { out.mem_resize(count); }
  }



template<typename T1>
inline
void
spop_trimatl_ext::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_trimatl_ext>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> U(in.m);
  const SpMat<eT>&   A = U.M;
  
  arma_debug_check( (A.is_square() == false), "trimatl(): given matrix must be square sized" );
  
  const uword row_offset = in.aux_uword_a;
  const uword col_offset = in.aux_uword_b;
  
  if(U.is_alias(out))
    {
    SpMat<eT> tmp;
    spop_trimatl_ext::apply_noalias(tmp, A, row_offset, col_offset);
    out.steal_mem(tmp);
    }
  else
    {
    spop_trimatl_ext::apply_noalias(out, A, row_offset, col_offset);
    }
  }



//! @}
