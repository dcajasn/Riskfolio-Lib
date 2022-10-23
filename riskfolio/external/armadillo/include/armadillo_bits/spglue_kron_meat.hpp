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


//! \addtogroup spglue_kron
//! @{



template<typename T1, typename T2>
inline
void
spglue_kron::apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_kron>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> UA(X.A);
  const unwrap_spmat<T2> UB(X.B);
  
  if(UA.is_alias(out) || UB.is_alias(out))
    {
    SpMat<eT> tmp;
    
    spglue_kron::apply_noalias(tmp, UA.M, UB.M);
    
    out.steal_mem(tmp);
    }
  else
    {
    spglue_kron::apply_noalias(out, UA.M, UB.M);
    }
  }



template<typename eT>
inline
void
spglue_kron::apply_noalias(SpMat<eT>& out, const SpMat<eT>& A, const SpMat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  const uword out_n_nonzero = A.n_nonzero * B.n_nonzero;
  
  out.reserve(A_n_rows * B_n_rows, A_n_cols * B_n_cols, out_n_nonzero);
  
  if(out_n_nonzero == 0)  { return; }
  
  access::rw(out.col_ptrs[0]) = 0;
  
  uword count = 0;
  
  for(uword A_col=0; A_col < A_n_cols; ++A_col)
  for(uword B_col=0; B_col < B_n_cols; ++B_col)
    {
    for(uword A_i = A.col_ptrs[A_col]; A_i < A.col_ptrs[A_col+1]; ++A_i)
      {
      const uword out_row = A.row_indices[A_i] * B_n_rows;
      
      const eT A_val = A.values[A_i];
      
      for(uword B_i = B.col_ptrs[B_col]; B_i < B.col_ptrs[B_col+1]; ++B_i)
        {
        access::rw(out.values[count])      = A_val * B.values[B_i];
        access::rw(out.row_indices[count]) = out_row + B.row_indices[B_i];
        
        count++;
        }
      }
    
    access::rw(out.col_ptrs[A_col * B_n_cols + B_col + 1]) = count;
    }
  }



// template<typename T1, typename T2>
// inline
// void
// spglue_kron::apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_kron>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   typedef typename T1::elem_type eT;
//   
//   const unwrap_spmat<T1> UA(X.A);
//   const unwrap_spmat<T2> UB(X.B);
//   
//   const SpMat<eT>& A = UA.M;
//   const SpMat<eT>& B = UB.M;
//   
//   umat    locs(2, A.n_nonzero * B.n_nonzero, arma_nozeros_indicator());
//   Col<eT> vals(   A.n_nonzero * B.n_nonzero, arma_nozeros_indicator());
//   
//   uword* locs_mem = locs.memptr();
//   eT*    vals_mem = vals.memptr();
//   
//   typename SpMat<eT>::const_iterator A_it     = A.begin();
//   typename SpMat<eT>::const_iterator A_it_end = A.end();
//   
//   typename SpMat<eT>::const_iterator B_it_start = B.begin();
//   typename SpMat<eT>::const_iterator B_it_end   = B.end();
//   
//   const uword B_n_rows = B.n_rows;
//   const uword B_n_cols = B.n_cols;
//   
//   uword i = 0;
//   
//   while(A_it != A_it_end)
//     {
//     typename SpMat<eT>::const_iterator B_it = B_it_start;
//     
//     const uword loc_row = A_it.row() * B_n_rows;
//     const uword loc_col = A_it.col() * B_n_cols;
//     
//     const eT A_val = (*A_it);
//     
//     while(B_it != B_it_end)
//       {
//       (*locs_mem) = loc_row + B_it.row();  locs_mem++;
//       (*locs_mem) = loc_col + B_it.col();  locs_mem++;
//       
//       vals_mem[i] = A_val * (*B_it);
//       
//       ++i;
//       ++B_it;
//       }
//     
//     ++A_it;
//     }
//   
//   out = SpMat<eT>(locs, vals, A.n_rows*B.n_rows, A.n_cols*B.n_cols);
//   }



//! @}
