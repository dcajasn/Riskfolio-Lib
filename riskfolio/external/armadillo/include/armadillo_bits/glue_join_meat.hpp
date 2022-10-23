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


//! \addtogroup glue_join
//! @{



template<typename T1, typename T2>
inline
void
glue_join_cols::apply_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& A, const Proxy<T2>& B)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_rows = A.get_n_rows();
  const uword A_n_cols = A.get_n_cols();
  
  const uword B_n_rows = B.get_n_rows();
  const uword B_n_cols = B.get_n_cols();
  
  arma_debug_check
    (
    ( (A_n_cols != B_n_cols) && ( (A_n_rows > 0) || (A_n_cols > 0) ) && ( (B_n_rows > 0) || (B_n_cols > 0) ) ),
    "join_cols() / join_vert(): number of columns must be the same"
    );
  
  out.set_size( A_n_rows + B_n_rows, (std::max)(A_n_cols, B_n_cols) );
  
  if( out.n_elem > 0 )
    {
    if(A.get_n_elem() > 0)
      { 
      out.submat(0,        0,   A_n_rows-1, out.n_cols-1) = A.Q;
      }
    
    if(B.get_n_elem() > 0)
      {
      out.submat(A_n_rows, 0, out.n_rows-1, out.n_cols-1) = B.Q;
      }
    }
  }
  
  
  
  
template<typename T1, typename T2>
inline
void
glue_join_cols::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_join_cols>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  if( (A.is_alias(out) == false) && (B.is_alias(out) == false) )
    {
    glue_join_cols::apply_noalias(out, A, B);
    }
  else
    {
    Mat<eT> tmp;
    
    glue_join_cols::apply_noalias(tmp, A, B);
    
    out.steal_mem(tmp);
    }
  }



template<typename eT, typename T1, typename T2, typename T3>
inline
void
glue_join_cols::apply(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const Base<eT,T3>& C_expr)
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<T1> UA(A_expr.get_ref());
  const quasi_unwrap<T2> UB(B_expr.get_ref());
  const quasi_unwrap<T3> UC(C_expr.get_ref());
  
  const Mat<eT>& A = UA.M;
  const Mat<eT>& B = UB.M;
  const Mat<eT>& C = UC.M;
  
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
glue_join_cols::apply(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const Base<eT,T3>& C_expr, const Base<eT,T4>& D_expr)
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<T1> UA(A_expr.get_ref());
  const quasi_unwrap<T2> UB(B_expr.get_ref());
  const quasi_unwrap<T3> UC(C_expr.get_ref());
  const quasi_unwrap<T4> UD(D_expr.get_ref());
  
  const Mat<eT>& A = UA.M;
  const Mat<eT>& B = UB.M;
  const Mat<eT>& C = UC.M;
  const Mat<eT>& D = UD.M;
  
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
glue_join_rows::apply_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& A, const Proxy<T2>& B)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_rows = A.get_n_rows();
  const uword A_n_cols = A.get_n_cols();
  
  const uword B_n_rows = B.get_n_rows();
  const uword B_n_cols = B.get_n_cols();
  
  arma_debug_check
    (
    ( (A_n_rows != B_n_rows) && ( (A_n_rows > 0) || (A_n_cols > 0) ) && ( (B_n_rows > 0) || (B_n_cols > 0) ) ),
    "join_rows() / join_horiz(): number of rows must be the same"
    );
  
  out.set_size( (std::max)(A_n_rows, B_n_rows), A_n_cols + B_n_cols );
  
  if( out.n_elem > 0 )
    {
    if(A.get_n_elem() > 0)
      {
      out.submat(0, 0,        out.n_rows-1,   A_n_cols-1) = A.Q;
      }
    
    if(B.get_n_elem() > 0)
      {
      out.submat(0, A_n_cols, out.n_rows-1, out.n_cols-1) = B.Q;
      }
    }
  }
  
  
  
  
template<typename T1, typename T2>
inline
void
glue_join_rows::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_join_rows>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  if( (A.is_alias(out) == false) && (B.is_alias(out) == false) )
    {
    glue_join_rows::apply_noalias(out, A, B);
    }
  else
    {
    Mat<eT> tmp;
    
    glue_join_rows::apply_noalias(tmp, A, B);
    
    out.steal_mem(tmp);
    }
  }



template<typename eT, typename T1, typename T2, typename T3>
inline
void
glue_join_rows::apply(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const Base<eT,T3>& C_expr)
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<T1> UA(A_expr.get_ref());
  const quasi_unwrap<T2> UB(B_expr.get_ref());
  const quasi_unwrap<T3> UC(C_expr.get_ref());
  
  const Mat<eT>& A = UA.M;
  const Mat<eT>& B = UB.M;
  const Mat<eT>& C = UC.M;
  
  const uword out_n_rows = (std::max)((std::max)(A.n_rows, B.n_rows), C.n_rows);
  const uword out_n_cols = A.n_cols + B.n_cols + C.n_cols;
  
  arma_debug_check( ((A.n_rows != out_n_rows) && ((A.n_rows > 0) || (A.n_cols > 0))), "join_rows() / join_horiz(): number of rows must be the same" );
  arma_debug_check( ((B.n_rows != out_n_rows) && ((B.n_rows > 0) || (B.n_cols > 0))), "join_rows() / join_horiz(): number of rows must be the same" );
  arma_debug_check( ((C.n_rows != out_n_rows) && ((C.n_rows > 0) || (C.n_cols > 0))), "join_rows() / join_horiz(): number of rows must be the same" );
  
  out.set_size(out_n_rows, out_n_cols);
  
  if(out.n_elem == 0)  { return; }
  
  uword col_start  = 0;
  uword col_end_p1 = 0;
  
  if(A.n_elem > 0)  { col_end_p1 += A.n_cols; out.cols(col_start, col_end_p1 - 1) = A; }
  
  col_start = col_end_p1;
  
  if(B.n_elem > 0)  { col_end_p1 += B.n_cols; out.cols(col_start, col_end_p1 - 1) = B; }
  
  col_start = col_end_p1;
  
  if(C.n_elem > 0)  { col_end_p1 += C.n_cols; out.cols(col_start, col_end_p1 - 1) = C; }
  }



template<typename eT, typename T1, typename T2, typename T3, typename T4>
inline
void
glue_join_rows::apply(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const Base<eT,T3>& C_expr, const Base<eT,T4>& D_expr)
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<T1> UA(A_expr.get_ref());
  const quasi_unwrap<T2> UB(B_expr.get_ref());
  const quasi_unwrap<T3> UC(C_expr.get_ref());
  const quasi_unwrap<T4> UD(D_expr.get_ref());
  
  const Mat<eT>& A = UA.M;
  const Mat<eT>& B = UB.M;
  const Mat<eT>& C = UC.M;
  const Mat<eT>& D = UD.M;
  
  const uword out_n_rows = (std::max)(((std::max)((std::max)(A.n_rows, B.n_rows), C.n_rows)), D.n_rows);
  const uword out_n_cols = A.n_cols + B.n_cols + C.n_cols + D.n_cols;
  
  arma_debug_check( ((A.n_rows != out_n_rows) && ((A.n_rows > 0) || (A.n_cols > 0))), "join_rows() / join_horiz(): number of rows must be the same" );
  arma_debug_check( ((B.n_rows != out_n_rows) && ((B.n_rows > 0) || (B.n_cols > 0))), "join_rows() / join_horiz(): number of rows must be the same" );
  arma_debug_check( ((C.n_rows != out_n_rows) && ((C.n_rows > 0) || (C.n_cols > 0))), "join_rows() / join_horiz(): number of rows must be the same" );
  arma_debug_check( ((D.n_rows != out_n_rows) && ((D.n_rows > 0) || (D.n_cols > 0))), "join_rows() / join_horiz(): number of rows must be the same" );
  
  out.set_size(out_n_rows, out_n_cols);
  
  if(out.n_elem == 0)  { return; }
  
  uword col_start  = 0;
  uword col_end_p1 = 0;
  
  if(A.n_elem > 0)  { col_end_p1 += A.n_cols; out.cols(col_start, col_end_p1 - 1) = A; }
  
  col_start = col_end_p1;
  
  if(B.n_elem > 0)  { col_end_p1 += B.n_cols; out.cols(col_start, col_end_p1 - 1) = B; }
  
  col_start = col_end_p1;
  
  if(C.n_elem > 0)  { col_end_p1 += C.n_cols; out.cols(col_start, col_end_p1 - 1) = C; }
  
  col_start = col_end_p1;
  
  if(D.n_elem > 0)  { col_end_p1 += D.n_cols; out.cols(col_start, col_end_p1 - 1) = D; }
  }



template<typename T1, typename T2>
inline
void
glue_join_slices::apply(Cube<typename T1::elem_type>& out, const GlueCube<T1,T2,glue_join_slices>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;

  const unwrap_cube<T1> A_tmp(X.A);
  const unwrap_cube<T2> B_tmp(X.B);
  
  const Cube<eT>& A = A_tmp.M;
  const Cube<eT>& B = B_tmp.M;
  
  if(A.n_elem == 0)  { out = B; return; }
  if(B.n_elem == 0)  { out = A; return; }
  
  arma_debug_check( ( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) ), "join_slices(): size of slices must be the same" );
  
  if( (&out != &A) && (&out != &B) )
    {
    out.set_size(A.n_rows, A.n_cols, A.n_slices + B.n_slices);
    
    out.slices(0,          A.n_slices-1  ) = A;
    out.slices(A.n_slices, out.n_slices-1) = B;
    }
  else  // we have aliasing
    {
    Cube<eT> C(A.n_rows, A.n_cols, A.n_slices + B.n_slices, arma_nozeros_indicator());
    
    C.slices(0,          A.n_slices-1) = A;
    C.slices(A.n_slices, C.n_slices-1) = B;
    
    out.steal_mem(C);
    }
  
  }



//! @}
