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


//! \addtogroup spop_repmat
//! @{



template<typename T1>
inline
void
spop_repmat::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_repmat>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> U(X.m);
  
  if(U.is_alias(out))
    {
    SpMat<eT> tmp;
    
    spop_repmat::apply_noalias(tmp, X.aux_uword_a, X.aux_uword_b, U.M);
    
    out.steal_mem(tmp);
    }
  else
    {
    spop_repmat::apply_noalias(out, X.aux_uword_a, X.aux_uword_b, U.M);
    }
  }



template<typename eT>
inline
void
spop_repmat::apply_noalias(SpMat<eT>& out, const uword A_n_rows, const uword A_n_cols, const SpMat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  const uword out_n_nonzero = A_n_rows * A_n_cols * B.n_nonzero;
  
  out.reserve(A_n_rows * B_n_rows, A_n_cols * B_n_cols, out_n_nonzero);
  
  if(out_n_nonzero == 0)  { return; }
  
  access::rw(out.col_ptrs[0]) = 0;
  
  uword count = 0;
  
  for(uword A_col=0; A_col < A_n_cols; ++A_col)
  for(uword B_col=0; B_col < B_n_cols; ++B_col)
    {
    for(uword A_row=0; A_row < A_n_rows; ++A_row)
      {
      const uword out_row = A_row * B_n_rows;
      
      for(uword B_i = B.col_ptrs[B_col]; B_i < B.col_ptrs[B_col+1]; ++B_i)
        {
        access::rw(out.values[count])      = B.values[B_i];
        access::rw(out.row_indices[count]) = out_row + B.row_indices[B_i];
        
        count++;
        }
      }
    
    access::rw(out.col_ptrs[A_col * B_n_cols + B_col + 1]) = count;
    }
  }



// template<typename T1>
// inline
// void
// spop_repmat::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_repmat>& in)
//   {
//   arma_extra_debug_sigprint();
//   
//   typedef typename T1::elem_type eT;
//   
//   const unwrap_spmat<T1> U(in.m);
//   const SpMat<eT>& X =   U.M;
//   
//   const uword X_n_rows = X.n_rows;
//   const uword X_n_cols = X.n_cols;
//   
//   const uword copies_per_row = in.aux_uword_a;
//   const uword copies_per_col = in.aux_uword_b;
//   
//   // out.set_size(X_n_rows * copies_per_row, X_n_cols * copies_per_col);
//   // 
//   // const uword out_n_rows = out.n_rows;
//   // const uword out_n_cols = out.n_cols;
//   // 
//   // if( (out_n_rows > 0) && (out_n_cols > 0) )
//   //   {
//   //   for(uword col = 0; col < out_n_cols; col += X_n_cols)
//   //   for(uword row = 0; row < out_n_rows; row += X_n_rows)
//   //     {
//   //     out.submat(row, col, row+X_n_rows-1, col+X_n_cols-1) = X;
//   //     }
//   //   }
//   
//   const uword out_n_rows = X_n_rows * copies_per_row;
//   const uword out_n_cols = X_n_cols * copies_per_col;
//   const uword out_nnz    = X.n_nonzero * copies_per_row * copies_per_col;
//   
//   if( (out_n_rows > 0) && (out_n_cols > 0) && (out_nnz > 0) )
//     {
//     umat    locs(2, out_nnz, arma_nozeros_indicator());
//     Col<eT> vals(   out_nnz, arma_nozeros_indicator());
//     
//     uword* locs_mem = locs.memptr();
//     eT*    vals_mem = vals.memptr();
//     
//     typename SpMat<eT>::const_iterator X_begin = X.begin();
//     typename SpMat<eT>::const_iterator X_end   = X.end();
//     typename SpMat<eT>::const_iterator X_it;
//     
//     for(uword col_offset = 0; col_offset < out_n_cols; col_offset += X_n_cols)
//     for(uword row_offset = 0; row_offset < out_n_rows; row_offset += X_n_rows)
//       {
//       for(X_it = X_begin; X_it != X_end; ++X_it)
//         {
//         const uword out_row = row_offset + X_it.row();
//         const uword out_col = col_offset + X_it.col();
//         
//         (*locs_mem) = out_row;  ++locs_mem;
//         (*locs_mem) = out_col;  ++locs_mem;
//         
//         (*vals_mem) = (*X_it);  ++vals_mem;
//         }
//       }
//     
//     out = SpMat<eT>(locs, vals, out_n_rows, out_n_cols);
//     }
//   else
//     {
//     out.zeros(out_n_rows, out_n_cols);
//     }
//   }



//! @}
