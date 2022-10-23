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


//! \addtogroup op_diagmat
//! @{



template<typename T1>
inline
void
op_diagmat::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_diagmat>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_Mat<T1>::value)
    {
    // allow detection of in-place operation
    
    const unwrap<T1>   U(X.m);
    const Mat<eT>& A = U.M;
    
    if(&out != &A)  // no aliasing
      {
      const Proxy< Mat<eT> > P(A);
      
      op_diagmat::apply(out, P);
      }
    else  // we have aliasing
      {
      const uword n_rows = out.n_rows;
      const uword n_cols = out.n_cols;
      
      if((n_rows == 1) || (n_cols == 1))  // create diagonal matrix from vector
        {
        const eT*   out_mem = out.memptr();
        const uword N       = out.n_elem;
        
        Mat<eT> tmp(N,N, arma_zeros_indicator());
        
        for(uword i=0; i<N; ++i)  { tmp.at(i,i) = out_mem[i]; }
        
        out.steal_mem(tmp);
        }
      else  // create diagonal matrix from matrix
        {
        const uword N = (std::min)(n_rows, n_cols);
        
        for(uword i=0; i < n_cols; ++i)
          {
          if(i < N)
            {
            eT& out_ii = out.at(i,i);
            
            const eT val = out_ii;
            
            arrayops::fill_zeros(out.colptr(i), n_rows);
            
            out_ii = val;
            }
          else
            {
            arrayops::fill_zeros(out.colptr(i), n_rows);
            }
          }
        }
      }
    }
  else
    {
    const Proxy<T1> P(X.m);
    
    if(P.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_diagmat::apply(tmp, P);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_diagmat::apply(out, P);
      }
    }
  }



template<typename T1>
inline
void
op_diagmat::apply(Mat<typename T1::elem_type>& out, const Proxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  const uword n_elem = P.get_n_elem();
  
  if(n_elem == 0)  { out.reset(); return; }
  
  const bool P_is_vec = (T1::is_row) || (T1::is_col) || (n_rows == 1) || (n_cols == 1);
  
  if(P_is_vec)
    {
    out.zeros(n_elem, n_elem);
    
    if(Proxy<T1>::use_at == false)
      {
      typename Proxy<T1>::ea_type Pea = P.get_ea();
      
      for(uword i=0; i < n_elem; ++i)  { out.at(i,i) = Pea[i]; }
      }
    else
      {
      if(n_rows == 1)
        {
        for(uword i=0; i < n_elem; ++i)  { out.at(i,i) = P.at(0,i); }
        }
      else
        {
        for(uword i=0; i < n_elem; ++i)  { out.at(i,i) = P.at(i,0); }
        }
      }
    }
  else  // P represents a matrix 
    {
    out.zeros(n_rows, n_cols);
    
    const uword N = (std::min)(n_rows, n_cols);
    
    for(uword i=0; i<N; ++i)  { out.at(i,i) = P.at(i,i); }
    }
  }



template<typename T1, typename T2>
inline
void
op_diagmat::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2,glue_times>, op_diagmat>& X)
  {
  arma_extra_debug_sigprint();
  
  op_diagmat::apply_times(out, X.m.A, X.m.B);
  }



template<typename T1, typename T2>
inline
void
op_diagmat::apply_times(Mat<typename T1::elem_type>& actual_out, const T1& X, const T2& Y, const typename arma_not_cx<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  const partial_unwrap<T1> UA(X);
  const partial_unwrap<T2> UB(Y);
  
  const typename partial_unwrap<T1>::stored_type& A = UA.M;
  const typename partial_unwrap<T2>::stored_type& B = UB.M;
  
  arma_debug_assert_trans_mul_size< partial_unwrap<T1>::do_trans, partial_unwrap<T2>::do_trans >(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  const bool use_alpha = partial_unwrap<T1>::do_times || partial_unwrap<T2>::do_times;
  const eT       alpha = use_alpha ? (UA.get_val() * UB.get_val()) : eT(0);
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;

  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  // check if the multiplication results in a vector
  
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == false) )
    {
    if((A_n_rows == 1) || (B_n_cols == 1))
      {
      arma_extra_debug_print("trans_A = false; trans_B = false; vector result");
      
      const Mat<eT> C     = A*B;
      const eT*     C_mem = C.memptr();
      const uword   N     = C.n_elem;
      
      actual_out.zeros(N,N);
      
      for(uword i=0; i<N; ++i)  { actual_out.at(i,i) = (use_alpha) ? eT(alpha * C_mem[i]) : eT(C_mem[i]); }
      
      return;
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true ) && (partial_unwrap<T2>::do_trans == false) )
    {
    if((A_n_cols == 1) || (B_n_cols == 1))
      {
      arma_extra_debug_print("trans_A = true; trans_B = false; vector result");
      
      const Mat<eT> C     = trans(A)*B;
      const eT*     C_mem = C.memptr();
      const uword   N     = C.n_elem;
      
      actual_out.zeros(N,N);
      
      for(uword i=0; i<N; ++i)  { actual_out.at(i,i) = (use_alpha) ? eT(alpha * C_mem[i]) : eT(C_mem[i]); }
      
      return;
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == true ) )
    {
    if((A_n_rows == 1) || (B_n_rows == 1))
      {
      arma_extra_debug_print("trans_A = false; trans_B = true; vector result");
      
      const Mat<eT> C     = A*trans(B);
      const eT*     C_mem = C.memptr();
      const uword   N     = C.n_elem;
      
      actual_out.zeros(N,N);
      
      for(uword i=0; i<N; ++i)  { actual_out.at(i,i) = (use_alpha) ? eT(alpha * C_mem[i]) : eT(C_mem[i]); }
      
      return;
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true ) && (partial_unwrap<T2>::do_trans == true ) )
    {
    if((A_n_cols == 1) || (B_n_rows == 1))
      {
      arma_extra_debug_print("trans_A = true; trans_B = true; vector result");
      
      const Mat<eT> C     = trans(A)*trans(B);
      const eT*     C_mem = C.memptr();
      const uword   N     = C.n_elem;
      
      actual_out.zeros(N,N);
      
      for(uword i=0; i<N; ++i)  { actual_out.at(i,i) = (use_alpha) ? eT(alpha * C_mem[i]) : eT(C_mem[i]); }
      
      return;
      }
    }
  
  // if we got to this point, the multiplication results in a matrix

  const bool is_alias = (UA.is_alias(actual_out) || UB.is_alias(actual_out));
  
  Mat<eT>  tmp;
  Mat<eT>& out = (is_alias) ? tmp : actual_out;
  
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == false) )
    {
    arma_extra_debug_print("trans_A = false; trans_B = false; matrix result");
    
    out.zeros(A_n_rows, B_n_cols);
    
    const uword N = (std::min)(A_n_rows, B_n_cols);
    
    for(uword k=0; k < N; ++k)
      {
      eT acc1 = eT(0);
      eT acc2 = eT(0);
      
      const eT* B_colptr = B.colptr(k);
      
      // condition: A_n_cols = B_n_rows
      
      uword j;
      
      for(j=1; j < A_n_cols; j+=2)
        {
        const uword i = (j-1);
        
        const eT tmp_i = B_colptr[i];
        const eT tmp_j = B_colptr[j];
        
        acc1 += A.at(k, i) * tmp_i;
        acc2 += A.at(k, j) * tmp_j;
        }
      
      const uword i = (j-1);
      
      if(i < A_n_cols)
        {
        acc1 += A.at(k, i) * B_colptr[i];
        }
      
      const eT acc = acc1 + acc2;
      
      out.at(k,k) = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true ) && (partial_unwrap<T2>::do_trans == false) )
    {
    arma_extra_debug_print("trans_A = true; trans_B = false; matrix result");
    
    out.zeros(A_n_cols, B_n_cols);
    
    const uword N = (std::min)(A_n_cols, B_n_cols);
    
    for(uword k=0; k < N; ++k)
      {
      const eT* A_colptr = A.colptr(k);
      const eT* B_colptr = B.colptr(k);
      
      // condition: A_n_rows = B_n_rows
      
      const eT acc = op_dot::direct_dot(A_n_rows, A_colptr, B_colptr);
      
      out.at(k,k) = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == true ) )
    {
    arma_extra_debug_print("trans_A = false; trans_B = true; matrix result");
    
    out.zeros(A_n_rows, B_n_rows);
    
    const uword N = (std::min)(A_n_rows, B_n_rows);
    
    for(uword k=0; k < N; ++k)
      {
      eT acc = eT(0);
      
      // condition: A_n_cols = B_n_cols
      
      for(uword i=0; i < A_n_cols; ++i)
        {
        acc += A.at(k,i) * B.at(k,i);
        }
      
      out.at(k,k) = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true ) && (partial_unwrap<T2>::do_trans == true ) )
    {
    arma_extra_debug_print("trans_A = true; trans_B = true; matrix result");
    
    out.zeros(A_n_cols, B_n_rows);
    
    const uword N = (std::min)(A_n_cols, B_n_rows);
    
    for(uword k=0; k < N; ++k)
      {
      eT acc = eT(0);
      
      const eT* A_colptr = A.colptr(k);
      
      // condition: A_n_rows = B_n_cols
      
      for(uword i=0; i < A_n_rows; ++i)
        {
        acc += A_colptr[i] * B.at(k,i);
        }
      
      out.at(k,k) = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  
  if(is_alias)  { actual_out.steal_mem(tmp); }
  }



template<typename T1, typename T2>
inline
void
op_diagmat::apply_times(Mat<typename T1::elem_type>& actual_out, const T1& X, const T2& Y, const typename arma_cx_only<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type   T;
  typedef typename T1::elem_type eT;
  
  const partial_unwrap<T1> UA(X);
  const partial_unwrap<T2> UB(Y);
  
  const typename partial_unwrap<T1>::stored_type& A = UA.M;
  const typename partial_unwrap<T2>::stored_type& B = UB.M;
  
  arma_debug_assert_trans_mul_size< partial_unwrap<T1>::do_trans, partial_unwrap<T2>::do_trans >(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  const bool use_alpha = partial_unwrap<T1>::do_times || partial_unwrap<T2>::do_times;
  const eT       alpha = use_alpha ? (UA.get_val() * UB.get_val()) : eT(0);
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  // check if the multiplication results in a vector
  
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == false) )
    {
    if((A_n_rows == 1) || (B_n_cols == 1))
      {
      arma_extra_debug_print("trans_A = false; trans_B = false; vector result");
      
      const Mat<eT> C     = A*B;
      const eT*     C_mem = C.memptr();
      const uword   N     = C.n_elem;
      
      actual_out.zeros(N,N);
      
      for(uword i=0; i<N; ++i)  { actual_out.at(i,i) = (use_alpha) ? eT(alpha * C_mem[i]) : eT(C_mem[i]); }
      
      return;
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true ) && (partial_unwrap<T2>::do_trans == false) )
    {
    if((A_n_cols == 1) || (B_n_cols == 1))
      {
      arma_extra_debug_print("trans_A = true; trans_B = false; vector result");
      
      const Mat<eT> C     = trans(A)*B;
      const eT*     C_mem = C.memptr();
      const uword   N     = C.n_elem;
      
      actual_out.zeros(N,N);
      
      for(uword i=0; i<N; ++i)  { actual_out.at(i,i) = (use_alpha) ? eT(alpha * C_mem[i]) : eT(C_mem[i]); }
      
      return;
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == true ) )
    {
    if((A_n_rows == 1) || (B_n_rows == 1))
      {
      arma_extra_debug_print("trans_A = false; trans_B = true; vector result");
      
      const Mat<eT> C     = A*trans(B);
      const eT*     C_mem = C.memptr();
      const uword   N     = C.n_elem;
      
      actual_out.zeros(N,N);
      
      for(uword i=0; i<N; ++i)  { actual_out.at(i,i) = (use_alpha) ? eT(alpha * C_mem[i]) : eT(C_mem[i]); }
      
      return;
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true ) && (partial_unwrap<T2>::do_trans == true ) )
    {
    if((A_n_cols == 1) || (B_n_rows == 1))
      {
      arma_extra_debug_print("trans_A = true; trans_B = true; vector result");
      
      const Mat<eT> C     = trans(A)*trans(B);
      const eT*     C_mem = C.memptr();
      const uword   N     = C.n_elem;
      
      actual_out.zeros(N,N);
      
      for(uword i=0; i<N; ++i)  { actual_out.at(i,i) = (use_alpha) ? eT(alpha * C_mem[i]) : eT(C_mem[i]); }
      
      return;
      }
    }
  
  // if we got to this point, the multiplication results in a matrix

  const bool is_alias = (UA.is_alias(actual_out) || UB.is_alias(actual_out));
  
  Mat<eT>  tmp;
  Mat<eT>& out = (is_alias) ? tmp : actual_out;
  
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == false) )
    {
    arma_extra_debug_print("trans_A = false; trans_B = false; matrix result");
    
    out.zeros(A_n_rows, B_n_cols);
    
    const uword N = (std::min)(A_n_rows, B_n_cols);
    
    for(uword k=0; k < N; ++k)
      {
      T acc_real = T(0);
      T acc_imag = T(0);
      
      const eT* B_colptr = B.colptr(k);
      
      // condition: A_n_cols = B_n_rows
      
      for(uword i=0; i < A_n_cols; ++i)
        {
        // acc += A.at(k, i) * B_colptr[i];
        
        const std::complex<T>& xx = A.at(k, i);
        const std::complex<T>& yy = B_colptr[i];
        
        const T a = xx.real();
        const T b = xx.imag();
        
        const T c = yy.real();
        const T d = yy.imag();
        
        acc_real += (a*c) - (b*d);
        acc_imag += (a*d) + (b*c);
        }
      
      const eT acc = std::complex<T>(acc_real, acc_imag);
      
      out.at(k,k) = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true) && (partial_unwrap<T2>::do_trans == false) )
    {
    arma_extra_debug_print("trans_A = true; trans_B = false; matrix result");
    
    out.zeros(A_n_cols, B_n_cols);
    
    const uword N = (std::min)(A_n_cols, B_n_cols);
    
    for(uword k=0; k < N; ++k)
      {
      T acc_real = T(0);
      T acc_imag = T(0);
      
      const eT* A_colptr = A.colptr(k);
      const eT* B_colptr = B.colptr(k);
      
      // condition: A_n_rows = B_n_rows
      
      for(uword i=0; i < A_n_rows; ++i)
        {
        // acc += std::conj(A_colptr[i]) * B_colptr[i];
        
        const std::complex<T>& xx = A_colptr[i];
        const std::complex<T>& yy = B_colptr[i];
        
        const T a = xx.real();
        const T b = xx.imag();
        
        const T c = yy.real();
        const T d = yy.imag();
        
        // take into account the complex conjugate of xx
        
        acc_real += (a*c) + (b*d);
        acc_imag += (a*d) - (b*c);
        }
      
      const eT acc = std::complex<T>(acc_real, acc_imag);
      
      out.at(k,k) = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == true) )
    {
    arma_extra_debug_print("trans_A = false; trans_B = true; matrix result");
    
    out.zeros(A_n_rows, B_n_rows);
    
    const uword N = (std::min)(A_n_rows, B_n_rows);
    
    for(uword k=0; k < N; ++k)
      {
      T acc_real = T(0);
      T acc_imag = T(0);
      
      // condition: A_n_cols = B_n_cols
      
      for(uword i=0; i < A_n_cols; ++i)
        {
        // acc += A.at(k,i) * std::conj(B.at(k,i));
        
        const std::complex<T>& xx = A.at(k, i);
        const std::complex<T>& yy = B.at(k, i);
        
        const T a = xx.real();
        const T b = xx.imag();
        
        const T c =  yy.real();
        const T d = -yy.imag();  // take the conjugate
        
        acc_real += (a*c) - (b*d);
        acc_imag += (a*d) + (b*c);
        }
      
      const eT acc = std::complex<T>(acc_real, acc_imag);
      
      out.at(k,k) = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true) && (partial_unwrap<T2>::do_trans == true) )
    {
    arma_extra_debug_print("trans_A = true; trans_B = true; matrix result");
    
    out.zeros(A_n_cols, B_n_rows);
    
    const uword N = (std::min)(A_n_cols, B_n_rows);
    
    for(uword k=0; k < N; ++k)
      {
      T acc_real = T(0);
      T acc_imag = T(0);
      
      const eT* A_colptr = A.colptr(k);
      
      // condition: A_n_rows = B_n_cols
      
      for(uword i=0; i < A_n_rows; ++i)
        {
        // acc += std::conj(A_colptr[i]) * std::conj(B.at(k,i));
        
        const std::complex<T>& xx = A_colptr[i];
        const std::complex<T>& yy = B.at(k, i);
        
        const T a =  xx.real();
        const T b = -xx.imag();  // take the conjugate
        
        const T c =  yy.real();
        const T d = -yy.imag();  // take the conjugate
        
        acc_real += (a*c) - (b*d);
        acc_imag += (a*d) + (b*c);
        }
      
      const eT acc = std::complex<T>(acc_real, acc_imag);
      
      out.at(k,k) = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  
  if(is_alias)  { actual_out.steal_mem(tmp); }
  }



//
//
//



template<typename T1>
inline
void
op_diagmat2::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_diagmat2>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword row_offset = X.aux_uword_a;
  const uword col_offset = X.aux_uword_b;
  
  const Proxy<T1> P(X.m);
  
  if(P.is_alias(out))
    {
    Mat<eT> tmp;
    
    op_diagmat2::apply(tmp, P, row_offset, col_offset);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_diagmat2::apply(out, P, row_offset, col_offset);
    }
  }



template<typename T1>
inline
void
op_diagmat2::apply(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword row_offset, const uword col_offset)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  const uword n_elem = P.get_n_elem();
  
  if(n_elem == 0)  { out.reset(); return; }
  
  const bool P_is_vec = (T1::is_row) || (T1::is_col) || (n_rows == 1) || (n_cols == 1);
  
  if(P_is_vec)
    {
    const uword n_pad = (std::max)(row_offset, col_offset);
    
    out.zeros(n_elem + n_pad, n_elem + n_pad);
    
    if(Proxy<T1>::use_at == false)
      {
      typename Proxy<T1>::ea_type Pea = P.get_ea();
      
      for(uword i=0; i < n_elem; ++i)  { out.at(row_offset + i, col_offset + i) = Pea[i]; }
      }
    else
      {
      if(n_rows == 1)
        {
        for(uword i=0; i < n_elem; ++i)  { out.at(row_offset + i, col_offset + i) = P.at(0,i); }
        }
      else
        {
        for(uword i=0; i < n_elem; ++i)  { out.at(row_offset + i, col_offset + i) = P.at(i,0); }
        }
      }
    }
  else  // P represents a matrix 
    {
    arma_debug_check_bounds
      (
      ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)),
      "diagmat(): requested diagonal out of bounds"
      );
    
    out.zeros(n_rows, n_cols);
    
    const uword N = (std::min)(n_rows - row_offset, n_cols - col_offset);
    
    for(uword i=0; i<N; ++i)
      {
      const uword row = i + row_offset;
      const uword col = i + col_offset;
      
      out.at(row,col) = P.at(row,col);
      }
    }
  }



//! @}
