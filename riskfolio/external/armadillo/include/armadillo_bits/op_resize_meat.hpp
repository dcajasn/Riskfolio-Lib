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



//! \addtogroup op_resize
//! @{



template<typename T1>
inline
void
op_resize::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_resize>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword new_n_rows = in.aux_uword_a;
  const uword new_n_cols = in.aux_uword_b;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  if(&out == &A)
    {
    op_resize::apply_mat_inplace(out, new_n_rows, new_n_cols);
    }
  else
    {
    op_resize::apply_mat_noalias(out, A, new_n_rows, new_n_cols);
    }
  }



template<typename eT>
inline
void
op_resize::apply_mat_inplace(Mat<eT>& A, const uword new_n_rows, const uword new_n_cols)
  {
  arma_extra_debug_sigprint();
  
  if( (A.n_rows == new_n_rows) && (A.n_cols == new_n_cols) )  { return; }
  
  if(A.is_empty())  { A.zeros(new_n_rows, new_n_cols); return; }
  
  Mat<eT> B;
  
  op_resize::apply_mat_noalias(B, A, new_n_rows, new_n_cols);
  
  A.steal_mem(B);
  }



template<typename eT>
inline
void
op_resize::apply_mat_noalias(Mat<eT>& out, const Mat<eT>& A, const uword new_n_rows, const uword new_n_cols)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(new_n_rows, new_n_cols);
  
  if( (new_n_rows > A.n_rows) || (new_n_cols > A.n_cols) )  { out.zeros(); }
  
  if( (out.n_elem > 0) && (A.n_elem > 0) )
    {
    const uword end_row = (std::min)(new_n_rows, A.n_rows) - 1;
    const uword end_col = (std::min)(new_n_cols, A.n_cols) - 1;
    
    out.submat(0, 0, end_row, end_col) = A.submat(0, 0, end_row, end_col);
    }
  }



//



template<typename T1>
inline
void
op_resize::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_resize>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword new_n_rows   = in.aux_uword_a;
  const uword new_n_cols   = in.aux_uword_b;
  const uword new_n_slices = in.aux_uword_c;
  
  const unwrap_cube<T1> tmp(in.m);
  const Cube<eT>& A   = tmp.M;
  
  if(&out == &A)
    {
    op_resize::apply_cube_inplace(out, new_n_rows, new_n_cols, new_n_slices);
    }
  else
    {
    op_resize::apply_cube_noalias(out, A, new_n_rows, new_n_cols, new_n_slices);
    }
  }



template<typename eT>
inline
void
op_resize::apply_cube_inplace(Cube<eT>& A, const uword new_n_rows, const uword new_n_cols, const uword new_n_slices)
  {
  arma_extra_debug_sigprint();
  
  if( (A.n_rows == new_n_rows) && (A.n_cols == new_n_cols) && (A.n_slices == new_n_slices) )  { return; }
  
  if(A.is_empty())  { A.zeros(new_n_rows, new_n_cols, new_n_slices); return; }
  
  Cube<eT> B;
  
  op_resize::apply_cube_noalias(B, A, new_n_rows, new_n_cols, new_n_slices);
  
  A.steal_mem(B);
  }



template<typename eT>
inline
void
op_resize::apply_cube_noalias(Cube<eT>& out, const Cube<eT>& A, const uword new_n_rows, const uword new_n_cols, const uword new_n_slices)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(new_n_rows, new_n_cols, new_n_slices);
  
  if( (new_n_rows > A.n_rows) || (new_n_cols > A.n_cols) || (new_n_slices > A.n_slices) )  { out.zeros(); }
  
  if( (out.n_elem > 0) && (A.n_elem > 0) )
    {
    const uword end_row   = (std::min)(new_n_rows,   A.n_rows)   - 1;
    const uword end_col   = (std::min)(new_n_cols,   A.n_cols)   - 1;
    const uword end_slice = (std::min)(new_n_slices, A.n_slices) - 1;
    
    out.subcube(0, 0, 0, end_row, end_col, end_slice) = A.subcube(0, 0, 0, end_row, end_col, end_slice);
    }
  }



//! @}
