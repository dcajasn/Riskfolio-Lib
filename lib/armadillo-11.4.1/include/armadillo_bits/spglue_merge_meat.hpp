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


//! \addtogroup spglue_merge
//! @{



template<typename eT>
arma_hot
inline
void
spglue_merge::subview_merge(SpSubview<eT>& sv, const SpMat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  if(sv.n_elem == 0)  { return; }
  
  if(B.n_nonzero == 0)  { sv.zeros(); return; }
  
  SpMat<eT>& A = access::rw(sv.m);
  
  const uword merge_n_nonzero = A.n_nonzero - sv.n_nonzero + B.n_nonzero;
  
  const uword sv_row_start = sv.aux_row1;
  const uword sv_col_start = sv.aux_col1;
  
  const uword sv_row_end   = sv.aux_row1 + sv.n_rows - 1;
  const uword sv_col_end   = sv.aux_col1 + sv.n_cols - 1;
  
  
  if(A.n_nonzero == sv.n_nonzero)
    {
    // A is either all zeros or has all of its elements in the subview
    // so the merge is equivalent to overwrite of A
    
    SpMat<eT> tmp(arma_reserve_indicator(), A.n_rows, A.n_cols, B.n_nonzero);
    
    typename SpMat<eT>::const_iterator B_it     = B.begin();
    typename SpMat<eT>::const_iterator B_it_end = B.end();
    
    uword tmp_count = 0;
    
    for(; B_it != B_it_end; ++B_it)
      {
      access::rw(tmp.values[tmp_count])      = (*B_it);
      access::rw(tmp.row_indices[tmp_count]) = B_it.row() + sv_row_start;
      access::rw(tmp.col_ptrs[B_it.col() + sv_col_start + 1])++;
      ++tmp_count;
      }
    
    for(uword i=0; i < tmp.n_cols; ++i)
      {
      access::rw(tmp.col_ptrs[i + 1]) += tmp.col_ptrs[i];
      }
    
    A.steal_mem(tmp);
    
    access::rw(sv.n_nonzero) = B.n_nonzero;
    
    return;
    }
  
  
  if(sv.n_nonzero > (A.n_nonzero/2))
    {
    // A has most of its elements in the subview,
    // so regenerate A with zeros in the subview region
    // in order to increase merging efficiency
    
    sv.zeros();
    }
  
  
  SpMat<eT> out(arma_reserve_indicator(), A.n_rows, A.n_cols, merge_n_nonzero);
  
  typename SpMat<eT>::const_iterator x_it  = A.begin();
  typename SpMat<eT>::const_iterator x_end = A.end();
  
  typename SpMat<eT>::const_iterator y_it  = B.begin();
  typename SpMat<eT>::const_iterator y_end = B.end();
  
  uword count = 0;
  
  bool x_it_valid = (x_it != x_end);
  bool y_it_valid = (y_it != y_end);
  
  while(x_it_valid || y_it_valid)
    {
    eT out_val = eT(0);
    
    const uword x_it_row = (x_it_valid) ? uword(x_it.row()) : uword(0);
    const uword x_it_col = (x_it_valid) ? uword(x_it.col()) : uword(0);
    
    const uword y_it_row = (y_it_valid) ? uword(sv_row_start + y_it.row()) : uword(0);
    const uword y_it_col = (y_it_valid) ? uword(sv_col_start + y_it.col()) : uword(0);
    
    bool use_y_loc = false;
    
    if(x_it_valid && y_it_valid)
      {
      if( (x_it_row == y_it_row) && (x_it_col == y_it_col) )
        {
        out_val = (*y_it);
        
        ++x_it;
        ++y_it;
        }
      else
        {
        if((x_it_col < y_it_col) || ((x_it_col == y_it_col) && (x_it_row < y_it_row))) // if y is closer to the end
          {
          const bool x_inside_box = ((x_it_row >= sv_row_start) && (x_it_row <= sv_row_end)) && ((x_it_col >= sv_col_start) && (x_it_col <= sv_col_end));
          
          out_val = (x_inside_box) ? eT(0) : (*x_it);
          
          ++x_it;
          }
        else
          {
          out_val = (*y_it);
          
          ++y_it;
          
          use_y_loc = true;
          }
        }
      }
    else
    if(x_it_valid)
      {
      const bool x_inside_box = ((x_it_row >= sv_row_start) && (x_it_row <= sv_row_end)) && ((x_it_col >= sv_col_start) && (x_it_col <= sv_col_end));
      
      out_val = (x_inside_box) ? eT(0) : (*x_it);
      
      ++x_it;
      }
    else
    if(y_it_valid)
      {
      out_val = (*y_it);
      
      ++y_it;
      
      use_y_loc = true;
      }
    
    if(out_val != eT(0))
      {
      access::rw(out.values[count]) = out_val;
      
      const uword out_row = (use_y_loc == false) ? x_it_row : y_it_row;
      const uword out_col = (use_y_loc == false) ? x_it_col : y_it_col;
      
      access::rw(out.row_indices[count]) = out_row;
      access::rw(out.col_ptrs[out_col + 1])++;
      ++count;
      }
    
    x_it_valid = (x_it != x_end);
    y_it_valid = (y_it != y_end);
    }
  
  arma_check( (count != merge_n_nonzero), "internal error: spglue_merge::subview_merge(): count != merge_n_nonzero" );
  
  const uword out_n_cols = out.n_cols;
  
  uword* col_ptrs = access::rwp(out.col_ptrs);
  
  for(uword c = 1; c <= out_n_cols; ++c)
    {
    col_ptrs[c] += col_ptrs[c - 1];
    }
  
  A.steal_mem(out);
  
  access::rw(sv.n_nonzero) = B.n_nonzero;
  }



template<typename eT>
arma_hot
inline
void
spglue_merge::subview_merge(SpSubview<eT>& sv, const Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  if(sv.n_elem == 0)  { return; }
  
  const eT*   B_memptr = B.memptr();
  const uword B_n_elem = B.n_elem;
  
  uword B_n_nonzero = 0;
  
  for(uword i=0; i < B_n_elem; ++i)
    {
    B_n_nonzero += (B_memptr[i] != eT(0)) ? uword(1) : uword(0);
    }
  
  if(B_n_nonzero == 0)  { sv.zeros(); return; }
  
  SpMat<eT>& A = access::rw(sv.m);
  
  const uword merge_n_nonzero = A.n_nonzero - sv.n_nonzero + B_n_nonzero;
  
  const uword sv_row_start = sv.aux_row1;
  const uword sv_col_start = sv.aux_col1;
  
  const uword sv_row_end   = sv.aux_row1 + sv.n_rows - 1;
  const uword sv_col_end   = sv.aux_col1 + sv.n_cols - 1;
  
  
  if(A.n_nonzero == sv.n_nonzero)
    {
    // A is either all zeros or has all of its elements in the subview
    // so the merge is equivalent to overwrite of A
    
    SpMat<eT> tmp(arma_reserve_indicator(), A.n_rows, A.n_cols, B_n_nonzero);
    
    typename Mat<eT>::const_row_col_iterator B_it     = B.begin_row_col();
    typename Mat<eT>::const_row_col_iterator B_it_end = B.end_row_col();
    
    uword tmp_count = 0;
    
    for(; B_it != B_it_end; ++B_it)
      {
      const eT val = (*B_it);
      
      if(val != eT(0))
        {
        access::rw(tmp.values[tmp_count])      = val;
        access::rw(tmp.row_indices[tmp_count]) = B_it.row() + sv_row_start;
        access::rw(tmp.col_ptrs[B_it.col() + sv_col_start + 1])++;
        ++tmp_count;
        }
      }
    
    for(uword i=0; i < tmp.n_cols; ++i)
      {
      access::rw(tmp.col_ptrs[i + 1]) += tmp.col_ptrs[i];
      }
    
    A.steal_mem(tmp);
    
    access::rw(sv.n_nonzero) = B_n_nonzero;
    
    return;
    }
  
  
  if(sv.n_nonzero > (A.n_nonzero/2))
    {
    // A has most of its elements in the subview,
    // so regenerate A with zeros in the subview region
    // in order to increase merging efficiency
    
    sv.zeros();
    }
  
  
  SpMat<eT> out(arma_reserve_indicator(), A.n_rows, A.n_cols, merge_n_nonzero);
  
  typename SpMat<eT>::const_iterator x_it  = A.begin();
  typename SpMat<eT>::const_iterator x_end = A.end();
  
  typename Mat<eT>::const_row_col_iterator y_it  = B.begin_row_col();
  typename Mat<eT>::const_row_col_iterator y_end = B.end_row_col();
  
  uword count = 0;
  
  bool x_it_valid = (x_it != x_end);
  bool y_it_valid = (y_it != y_end);
    
  while(x_it_valid || y_it_valid)
    {
    eT out_val = eT(0);
    
    const uword x_it_row = (x_it_valid) ? uword(x_it.row()) : uword(0);
    const uword x_it_col = (x_it_valid) ? uword(x_it.col()) : uword(0);
    
    const uword y_it_row = (y_it_valid) ? uword(sv_row_start + y_it.row()) : uword(0);
    const uword y_it_col = (y_it_valid) ? uword(sv_col_start + y_it.col()) : uword(0);
    
    bool use_y_loc = false;
    
    if(x_it_valid && y_it_valid)
      {
      if( (x_it_row == y_it_row) && (x_it_col == y_it_col) )
        {
        out_val = (*y_it);
        
        ++x_it;
        ++y_it;
        }
      else
        {
        if((x_it_col < y_it_col) || ((x_it_col == y_it_col) && (x_it_row < y_it_row))) // if y is closer to the end
          {
          const bool x_inside_box = ((x_it_row >= sv_row_start) && (x_it_row <= sv_row_end)) && ((x_it_col >= sv_col_start) && (x_it_col <= sv_col_end));
          
          out_val = (x_inside_box) ? eT(0) : (*x_it);
          
          ++x_it;
          }
        else
          {
          out_val = (*y_it);
          
          ++y_it;
          
          use_y_loc = true;
          }
        }
      }
    else
    if(x_it_valid)
      {
      const bool x_inside_box = ((x_it_row >= sv_row_start) && (x_it_row <= sv_row_end)) && ((x_it_col >= sv_col_start) && (x_it_col <= sv_col_end));
      
      out_val = (x_inside_box) ? eT(0) : (*x_it);
      
      ++x_it;
      }
    else
    if(y_it_valid)
      {
      out_val = (*y_it);
      
      ++y_it;
      
      use_y_loc = true;
      }
    
    if(out_val != eT(0))
      {
      access::rw(out.values[count]) = out_val;
      
      const uword out_row = (use_y_loc == false) ? x_it_row : y_it_row;
      const uword out_col = (use_y_loc == false) ? x_it_col : y_it_col;
      
      access::rw(out.row_indices[count]) = out_row;
      access::rw(out.col_ptrs[out_col + 1])++;
      ++count;
      }
    
    x_it_valid = (x_it != x_end);
    y_it_valid = (y_it != y_end);
    }
  
  arma_check( (count != merge_n_nonzero), "internal error: spglue_merge::subview_merge(): count != merge_n_nonzero" );
  
  const uword out_n_cols = out.n_cols;
  
  uword* col_ptrs = access::rwp(out.col_ptrs);
  
  for(uword c = 1; c <= out_n_cols; ++c)
    {
    col_ptrs[c] += col_ptrs[c - 1];
    }
  
  A.steal_mem(out);
  
  access::rw(sv.n_nonzero) = B_n_nonzero;
  }



template<typename eT>
arma_hot
inline
void
spglue_merge::symmat_merge(SpMat<eT>& out, const SpMat<eT>& A, const SpMat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  out.reserve(A.n_rows, A.n_cols, 2*A.n_nonzero); // worst case scenario
  
  typename SpMat<eT>::const_iterator x_it  = A.begin();
  typename SpMat<eT>::const_iterator x_end = A.end();
  
  typename SpMat<eT>::const_iterator y_it  = B.begin();
  typename SpMat<eT>::const_iterator y_end = B.end();
  
  uword count = 0;
  
  while( (x_it != x_end) || (y_it != y_end) )
    {
    eT out_val;
    
    const uword x_it_col = x_it.col();
    const uword x_it_row = x_it.row();
    
    const uword y_it_col = y_it.col();
    const uword y_it_row = y_it.row();
    
    bool use_y_loc = false;
    
    if(x_it == y_it)
      {
      // this can only happen on the diagonal
      
      out_val = (*x_it);
      
      ++x_it;
      ++y_it;
      }
    else
      {
      if((x_it_col < y_it_col) || ((x_it_col == y_it_col) && (x_it_row < y_it_row))) // if y is closer to the end
        {
        out_val = (*x_it);
        
        ++x_it;
        }
      else
        {
        out_val = (*y_it);
        
        ++y_it;
        
        use_y_loc = true;
        }
      }
    
    access::rw(out.values[count]) = out_val;
    
    const uword out_row = (use_y_loc == false) ? x_it_row : y_it_row;
    const uword out_col = (use_y_loc == false) ? x_it_col : y_it_col;
    
    access::rw(out.row_indices[count]) = out_row;
    access::rw(out.col_ptrs[out_col + 1])++;
    ++count;
    }
  
  const uword out_n_cols = out.n_cols;
  
  uword* col_ptrs = access::rwp(out.col_ptrs);
  
  // Fix column pointers to be cumulative.
  for(uword c = 1; c <= out_n_cols; ++c)
    {
    col_ptrs[c] += col_ptrs[c - 1];
    }
  
  // quick resize without reallocating memory and copying data
  access::rw(         out.n_nonzero) = count;
  access::rw(     out.values[count]) = eT(0);
  access::rw(out.row_indices[count]) = uword(0);
  }



template<typename eT>
arma_hot
inline
void
spglue_merge::diagview_merge(SpMat<eT>& out, const SpMat<eT>& A, const SpMat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: assuming that B has non-zero elements only on the main diagonal
  
  out.reserve(A.n_rows, A.n_cols, A.n_nonzero + B.n_nonzero); // worst case scenario
  
  typename SpMat<eT>::const_iterator x_it  = A.begin();
  typename SpMat<eT>::const_iterator x_end = A.end();
  
  typename SpMat<eT>::const_iterator y_it  = B.begin();
  typename SpMat<eT>::const_iterator y_end = B.end();
  
  uword count = 0;
  
  while( (x_it != x_end) || (y_it != y_end) )
    {
    eT out_val = eT(0);
    
    const uword x_it_col = x_it.col();
    const uword x_it_row = x_it.row();
    
    const uword y_it_col = y_it.col();
    const uword y_it_row = y_it.row();
    
    bool use_y_loc = false;
    
    if(x_it == y_it)
      {
      // this can only happen on the diagonal
      
      out_val = (*y_it);
      
      ++x_it;
      ++y_it;
      }
    else
      {
      if((x_it_col < y_it_col) || ((x_it_col == y_it_col) && (x_it_row < y_it_row))) // if y is closer to the end
        {
        if(x_it_col != x_it_row)  { out_val = (*x_it); }  // don't take values from the main diagonal of A
        
        ++x_it;
        }
      else
        {
        if(y_it_col == y_it_row)  { out_val = (*y_it); use_y_loc = true; }  // take values only from the main diagonal of B
        
        ++y_it;
        }
      }
    
    if(out_val != eT(0))
      {
      access::rw(out.values[count]) = out_val;
      
      const uword out_row = (use_y_loc == false) ? x_it_row : y_it_row;
      const uword out_col = (use_y_loc == false) ? x_it_col : y_it_col;
      
      access::rw(out.row_indices[count]) = out_row;
      access::rw(out.col_ptrs[out_col + 1])++;
      ++count;
      }
    }
  
  const uword out_n_cols = out.n_cols;
  
  uword* col_ptrs = access::rwp(out.col_ptrs);
  
  // Fix column pointers to be cumulative.
  for(uword c = 1; c <= out_n_cols; ++c)
    {
    col_ptrs[c] += col_ptrs[c - 1];
    }
  
  // quick resize without reallocating memory and copying data
  access::rw(         out.n_nonzero) = count;
  access::rw(     out.values[count]) = eT(0);
  access::rw(out.row_indices[count]) = uword(0);
  }



//! @}
