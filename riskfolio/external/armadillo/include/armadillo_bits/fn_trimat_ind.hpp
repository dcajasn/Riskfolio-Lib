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


//! \addtogroup fn_trimat_ind
//! @{


arma_warn_unused
inline
uvec
trimatu_ind(const SizeMat& s, const sword k = 0)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = s.n_rows;
  const uword n_cols = s.n_cols;
  
  const uword row_offset = (k < 0) ? uword(-k) : uword(0);
  const uword col_offset = (k > 0) ? uword( k) : uword(0);
  
  arma_debug_check_bounds( ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)), "trimatu_ind(): requested diagonal is out of bounds" );
  
  const uword N = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  uvec   tmp(n_rows * n_cols, arma_nozeros_indicator());  // worst case scenario
  uword* tmp_mem = tmp.memptr();
  uword  count   = 0;
  
  for(uword i=0; i < n_cols; ++i)
    {
    const uword col = i + col_offset;
    
    if(i < N)
      {
      const uword end_row = i + row_offset;
      
      const uword index_offset = (n_rows * col);
      
      for(uword row=0; row <= end_row; ++row)
        {
        tmp_mem[count] = index_offset + row;
        ++count;
        }
      }
    else
      {
      if(col < n_cols)
        {
        const uword index_offset = (n_rows * col);
        
        for(uword row=0; row < n_rows; ++row)
          {
          tmp_mem[count] = index_offset + row;
          ++count;
          }
        }
      }
    }
  
  uvec out;
  
  out.steal_mem_col(tmp, count);
  
  return out;
  }



arma_warn_unused
inline
uvec
trimatl_ind(const SizeMat& s, const sword k = 0)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = s.n_rows;
  const uword n_cols = s.n_cols;
  
  const uword row_offset = (k < 0) ? uword(-k) : uword(0);
  const uword col_offset = (k > 0) ? uword( k) : uword(0);
  
  arma_debug_check_bounds( ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)), "trimatl_ind(): requested diagonal is out of bounds" );
  
  const uword N = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  uvec   tmp(n_rows * n_cols, arma_nozeros_indicator());  // worst case scenario
  uword* tmp_mem = tmp.memptr();
  uword  count   = 0;
  
  for(uword col=0; col < col_offset; ++col)
    {
    const uword index_offset = (n_rows * col);
    
    for(uword row=0; row < n_rows; ++row)
      {
      tmp_mem[count] = index_offset + row;
      ++count;
      }
    }
  
  for(uword i=0; i<N; ++i)
    {
    const uword start_row = i + row_offset;
    const uword       col = i + col_offset;
    
    const uword index_offset = (n_rows * col);
    
    for(uword row=start_row; row < n_rows; ++row)
      {
      tmp_mem[count] = index_offset + row;
      ++count;
      }
    }
  
  uvec out;
  
  out.steal_mem_col(tmp, count);
  
  return out;
  }



//! @}
