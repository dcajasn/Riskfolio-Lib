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



//! \addtogroup op_repelem
//! @{



template<typename obj>
inline
void
op_repelem::apply_noalias(Mat<typename obj::elem_type>& out, const obj& X, const uword copies_per_row, const uword copies_per_col)
  {
  arma_extra_debug_sigprint();
  
  typedef typename obj::elem_type eT;
  
  const uword X_n_rows = obj::is_row ? uword(1) : X.n_rows;
  const uword X_n_cols = obj::is_col ? uword(1) : X.n_cols;
  
  out.set_size(X_n_rows * copies_per_row, X_n_cols * copies_per_col);
  
  if(out.n_elem == 0)  { return; }
  
  for(uword col=0; col < X_n_cols; ++col)
    {
    const uword out_col_offset = col * copies_per_col;
    
    eT* out_colptr_first = out.colptr(out_col_offset);
    
    for(uword row=0; row < X_n_rows; ++row)
      {
      const uword out_row_offset = row * copies_per_row;
      
      const eT copy_value = X.at(row, col);
      
      for(uword row_copy=0; row_copy < copies_per_row; ++row_copy)
        {
        out_colptr_first[out_row_offset + row_copy] = copy_value;
        }
      
      if(copies_per_col != 1)
        {
        for(uword col_copy=1; col_copy < copies_per_col; ++col_copy)
          {
          eT* out_colptr = out.colptr(out_col_offset + col_copy);
          
          arrayops::copy(&out_colptr[out_row_offset], &out_colptr_first[out_row_offset], copies_per_row);
          }
        }
      }
    }
  }



template<typename T1>
inline
void
op_repelem::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_repelem>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword copies_per_row = in.aux_uword_a;
  const uword copies_per_col = in.aux_uword_b;
  
  const quasi_unwrap<T1> U(in.m);
  
  if(U.is_alias(out))
    {
    Mat<eT> tmp;
    
    op_repelem::apply_noalias(tmp, U.M, copies_per_row, copies_per_col);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_repelem::apply_noalias(out, U.M, copies_per_row, copies_per_col);
    }
  }



//! @}
