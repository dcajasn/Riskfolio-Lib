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


//! \addtogroup spglue_join
//! @{



template<typename T1, typename T2>
inline
void
spglue_join_cols::apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_join_cols>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> UA(X.A);
  const unwrap_spmat<T2> UB(X.B);
  
  if(UA.is_alias(out) || UB.is_alias(out))
    {
    SpMat<eT> tmp;
    
    spglue_join_cols::apply_noalias(tmp, UA.M, UB.M);
    
    out.steal_mem(tmp);
    }
  else
    {
    spglue_join_cols::apply_noalias(out, UA.M, UB.M);
    }
  }



template<typename eT>
inline
void
spglue_join_cols::apply_noalias(SpMat<eT>& out, const SpMat<eT>& A, const SpMat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  arma_debug_check
    (
    ( (A_n_cols != B_n_cols) && ( (A_n_rows > 0) || (A_n_cols > 0) ) && ( (B_n_rows > 0) || (B_n_cols > 0) ) ),
    "join_cols() / join_vert(): number of columns must be the same"
    );
  
  out.set_size( A_n_rows + B_n_rows, (std::max)(A_n_cols, B_n_cols) );
  
  if( out.n_elem > 0 )
    {
    if(A.is_empty() == false)
      { 
      out.submat(0,        0,   A_n_rows-1, out.n_cols-1) = A;
      }
    
    if(B.is_empty() == false)
      {
      out.submat(A_n_rows, 0, out.n_rows-1, out.n_cols-1) = B;
      }
    }
  }



template<typename eT, typename T1, typename T2, typename T3>
inline
void
spglue_join_cols::apply(SpMat<eT>& out, const SpBase<eT,T1>& A_expr, const SpBase<eT,T2>& B_expr, const SpBase<eT,T3>& C_expr)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> UA(A_expr.get_ref());
  const unwrap_spmat<T2> UB(B_expr.get_ref());
  const unwrap_spmat<T3> UC(C_expr.get_ref());
  
  const SpMat<eT>& A = UA.M;
  const SpMat<eT>& B = UB.M;
  const SpMat<eT>& C = UC.M;
  
  const uword out_n_rows = A.n_rows + B.n_rows + C.n_rows;
  const uword out_n_cols = (std::max)((std::max)(A.n_cols, B.n_cols), C.n_cols);
  
  arma_debug_check( ((A.n_cols != out_n_cols) && ((A.n_rows > 0) || (A.n_cols > 0))), "join_cols() / join_vert(): number of columns must be the same" );
  arma_debug_check( ((B.n_cols != out_n_cols) && ((B.n_rows > 0) || (B.n_cols > 0))), "join_cols() / join_vert(): number of columns must be the same" );
  arma_debug_check( ((C.n_cols != out_n_cols) && ((C.n_rows > 0) || (C.n_cols > 0))), "join_cols() / join_vert(): number of columns must be the same" );
  
  out.set_size(out_n_rows, out_n_cols);
  
  if(out.n_elem == 0)  { return; }
  
  uword row_start  = 0;
  uword row_end_p1 = 0;
  
  if(A.n_elem > 0)  { row_end_p1 += A.n_rows; out.rows(row_start, row_end_p1 - 1) = A; }
  
  row_start = row_end_p1;
  
  if(B.n_elem > 0)  { row_end_p1 += B.n_rows; out.rows(row_start, row_end_p1 - 1) = B; }
  
  row_start = row_end_p1;
  
  if(C.n_elem > 0)  { row_end_p1 += C.n_rows; out.rows(row_start, row_end_p1 - 1) = C; }
  }



template<typename eT, typename T1, typename T2, typename T3, typename T4>
inline
void
spglue_join_cols::apply(SpMat<eT>& out, const SpBase<eT,T1>& A_expr, const SpBase<eT,T2>& B_expr, const SpBase<eT,T3>& C_expr, const SpBase<eT,T4>& D_expr)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> UA(A_expr.get_ref());
  const unwrap_spmat<T2> UB(B_expr.get_ref());
  const unwrap_spmat<T3> UC(C_expr.get_ref());
  const unwrap_spmat<T4> UD(D_expr.get_ref());
  
  const SpMat<eT>& A = UA.M;
  const SpMat<eT>& B = UB.M;
  const SpMat<eT>& C = UC.M;
  const SpMat<eT>& D = UD.M;
  
  const uword out_n_rows = A.n_rows + B.n_rows + C.n_rows + D.n_rows;
  const uword out_n_cols = (std::max)(((std::max)((std::max)(A.n_cols, B.n_cols), C.n_cols)), D.n_cols);
  
  arma_debug_check( ((A.n_cols != out_n_cols) && ((A.n_rows > 0) || (A.n_cols > 0))), "join_cols() / join_vert(): number of columns must be the same" );
  arma_debug_check( ((B.n_cols != out_n_cols) && ((B.n_rows > 0) || (B.n_cols > 0))), "join_cols() / join_vert(): number of columns must be the same" );
  arma_debug_check( ((C.n_cols != out_n_cols) && ((C.n_rows > 0) || (C.n_cols > 0))), "join_cols() / join_vert(): number of columns must be the same" );
  arma_debug_check( ((D.n_cols != out_n_cols) && ((D.n_rows > 0) || (D.n_cols > 0))), "join_cols() / join_vert(): number of columns must be the same" );
  
  out.set_size(out_n_rows, out_n_cols);
  
  if(out.n_elem == 0)  { return; }
  
  uword row_start  = 0;
  uword row_end_p1 = 0;
  
  if(A.n_elem > 0)  { row_end_p1 += A.n_rows; out.rows(row_start, row_end_p1 - 1) = A; }
  
  row_start = row_end_p1;
  
  if(B.n_elem > 0)  { row_end_p1 += B.n_rows; out.rows(row_start, row_end_p1 - 1) = B; }
  
  row_start = row_end_p1;
  
  if(C.n_elem > 0)  { row_end_p1 += C.n_rows; out.rows(row_start, row_end_p1 - 1) = C; }
  
  row_start = row_end_p1;
  
  if(D.n_elem > 0)  { row_end_p1 += D.n_rows; out.rows(row_start, row_end_p1 - 1) = D; }
  }



template<typename T1, typename T2>
inline
void
spglue_join_rows::apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_join_rows>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> UA(X.A);
  const unwrap_spmat<T2> UB(X.B);
  
  if(UA.is_alias(out) || UB.is_alias(out))
    {
    SpMat<eT> tmp;
    
    spglue_join_rows::apply_noalias(tmp, UA.M, UB.M);
    
    out.steal_mem(tmp);
    }
  else
    {
    spglue_join_rows::apply_noalias(out, UA.M, UB.M);
    }
  }



template<typename eT>
inline
void
spglue_join_rows::apply_noalias(SpMat<eT>& out, const SpMat<eT>& A, const SpMat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  const uword A_n_nz   = A.n_nonzero;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  const uword B_n_nz   = B.n_nonzero;
  
  arma_debug_check
    (
    ( (A_n_rows != B.n_rows) && ( (A_n_rows > 0) || (A_n_cols > 0) ) && ( (B_n_rows > 0) || (B_n_cols > 0) ) ),
    "join_rows() / join_horiz(): number of rows must be the same"
    );
  
  const uword C_n_rows = (std::max)(A_n_rows, B_n_rows);
  const uword C_n_cols = A_n_cols + B_n_cols;
  const uword C_n_nz   = A_n_nz + B_n_nz;
  
  if( ((C_n_rows * C_n_cols) == 0) || (C_n_nz == 0) )
    {
    out.zeros(C_n_rows, C_n_cols);
    return;
    }
  
  out.reserve(C_n_rows, C_n_cols, C_n_nz);
  
  arrayops::copy( access::rwp(out.values),          A.values, A_n_nz   );
  arrayops::copy( access::rwp(out.values) + A_n_nz, B.values, B_n_nz+1 );
  
  arrayops::copy( access::rwp(out.row_indices),          A.row_indices, A_n_nz   );
  arrayops::copy( access::rwp(out.row_indices) + A_n_nz, B.row_indices, B_n_nz+1 );
  
  arrayops::copy( access::rwp(out.col_ptrs),            A.col_ptrs, A_n_cols   );
  arrayops::copy( access::rwp(out.col_ptrs) + A_n_cols, B.col_ptrs, B_n_cols+2 );
  
  arrayops::inplace_plus( access::rwp(out.col_ptrs) + A_n_cols, A_n_nz, B_n_cols+1 );
  
  
  // // OLD METHOD
  // 
  // umat    locs(2, C_n_nz, arma_nozeros_indicator());
  // Col<eT> vals(   C_n_nz, arma_nozeros_indicator());
  // 
  // uword* locs_mem = locs.memptr();
  // eT*    vals_mem = vals.memptr();
  // 
  // typename SpMat<eT>::const_iterator A_it = A.begin();
  // 
  // for(uword i=0; i < A_n_nz; ++i)
  //   {
  //   const uword row = A_it.row();
  //   const uword col = A_it.col();
  //   
  //   (*locs_mem) = row;  locs_mem++;
  //   (*locs_mem) = col;  locs_mem++;
  //   
  //   (*vals_mem) = (*A_it); vals_mem++;
  //   
  //   ++A_it;
  //   }
  // 
  // typename SpMat<eT>::const_iterator B_it = B.begin();
  // 
  // for(uword i=0; i < B_n_nz; ++i)
  //   {
  //   const uword row =            B_it.row();
  //   const uword col = A_n_cols + B_it.col();
  //   
  //   (*locs_mem) = row;  locs_mem++;
  //   (*locs_mem) = col;  locs_mem++;
  //   
  //   (*vals_mem) = (*B_it); vals_mem++;
  //   
  //   ++B_it;
  //   }
  // 
  // // TODO: the first element of B within C will always have a larger index than the last element of A in C;
  // // TODO: so, is sorting really necessary here?
  // SpMat<eT> tmp(locs, vals, C_n_rows, C_n_cols, true, false);
  // 
  // out.steal_mem(tmp);
  }



template<typename eT, typename T1, typename T2, typename T3>
inline
void
spglue_join_rows::apply(SpMat<eT>& out, const SpBase<eT,T1>& A_expr, const SpBase<eT,T2>& B_expr, const SpBase<eT,T3>& C_expr)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> UA(A_expr.get_ref());
  const unwrap_spmat<T2> UB(B_expr.get_ref());
  const unwrap_spmat<T3> UC(C_expr.get_ref());
  
  const SpMat<eT>& A = UA.M;
  const SpMat<eT>& B = UB.M;
  const SpMat<eT>& C = UC.M;
  
  SpMat<eT> tmp;
  
  spglue_join_rows::apply_noalias(tmp, A,   B);
  spglue_join_rows::apply_noalias(out, tmp, C);
  }



template<typename eT, typename T1, typename T2, typename T3, typename T4>
inline
void
spglue_join_rows::apply(SpMat<eT>& out, const SpBase<eT,T1>& A_expr, const SpBase<eT,T2>& B_expr, const SpBase<eT,T3>& C_expr, const SpBase<eT,T4>& D_expr)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> UA(A_expr.get_ref());
  const unwrap_spmat<T2> UB(B_expr.get_ref());
  const unwrap_spmat<T3> UC(C_expr.get_ref());
  const unwrap_spmat<T4> UD(D_expr.get_ref());
  
  const SpMat<eT>& A = UA.M;
  const SpMat<eT>& B = UB.M;
  const SpMat<eT>& C = UC.M;
  const SpMat<eT>& D = UD.M;
  
  SpMat<eT> AB;
  SpMat<eT> ABC;
  
  spglue_join_rows::apply_noalias(AB,  A,   B);
  spglue_join_rows::apply_noalias(ABC, AB,  C);
  spglue_join_rows::apply_noalias(out, ABC, D);
  }



//! @}
