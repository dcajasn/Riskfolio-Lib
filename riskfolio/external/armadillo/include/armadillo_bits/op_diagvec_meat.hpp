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


//! \addtogroup op_diagvec
//! @{



template<typename T1>
inline
void
op_diagvec::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_diagvec>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(X.m);
  
  if(P.is_alias(out) == false)
    {
    op_diagvec::apply_proxy(out, P);
    }
  else
    {
    Mat<eT> tmp;
    
    op_diagvec::apply_proxy(tmp, P);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1>
inline
void
op_diagvec::apply_proxy(Mat<typename T1::elem_type>& out, const Proxy<T1>& P)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  const uword len = (std::min)(n_rows, n_cols);
  
  out.set_size(len, 1);
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < len; i+=2, j+=2)
    {
    const eT tmp_i = P.at(i, i);
    const eT tmp_j = P.at(j, j);
    
    out_mem[i] = tmp_i;
    out_mem[j] = tmp_j;
    }
  
  if(i < len)
    {
    out_mem[i] = P.at(i, i);
    }
  }



template<typename T1, typename T2>
inline
void
op_diagvec::apply(Mat<typename T1::elem_type>& actual_out, const Op< Glue<T1,T2,glue_times>, op_diagvec>& X, const typename arma_not_cx<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  const partial_unwrap<T1> UA(X.m.A);
  const partial_unwrap<T2> UB(X.m.B);
  
  const typename partial_unwrap<T1>::stored_type& A = UA.M;
  const typename partial_unwrap<T2>::stored_type& B = UB.M;
  
  arma_debug_assert_trans_mul_size< partial_unwrap<T1>::do_trans, partial_unwrap<T2>::do_trans >(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  if( (A.n_elem == 0) || (B.n_elem == 0) )  { actual_out.reset(); return; }
  
  const bool use_alpha = partial_unwrap<T1>::do_times || partial_unwrap<T2>::do_times;
  const eT       alpha = use_alpha ? (UA.get_val() * UB.get_val()) : eT(0);
  
  const bool is_alias  = (UA.is_alias(actual_out) || UB.is_alias(actual_out));
  
  Mat<eT>  tmp;
  Mat<eT>& out = (is_alias) ? tmp : actual_out;
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;

  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == false) )
    {
    arma_extra_debug_print("trans_A = false; trans_B = false;");
    
    const uword N = (std::min)(A_n_rows, B_n_cols);
    
    out.set_size(N,1);
    
    eT* out_mem = out.memptr();
    
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
      
      out_mem[k] = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true ) && (partial_unwrap<T2>::do_trans == false) )
    {
    arma_extra_debug_print("trans_A = true; trans_B = false;");
    
    const uword N = (std::min)(A_n_cols, B_n_cols);
    
    out.set_size(N,1);
    
    eT* out_mem = out.memptr();
    
    for(uword k=0; k < N; ++k)
      {
      const eT* A_colptr = A.colptr(k);
      const eT* B_colptr = B.colptr(k);
      
      // condition: A_n_rows = B_n_rows
      
      const eT acc = op_dot::direct_dot(A_n_rows, A_colptr, B_colptr);
      
      out_mem[k] = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == true ) )
    {
    arma_extra_debug_print("trans_A = false; trans_B = true;");
    
    const uword N = (std::min)(A_n_rows, B_n_rows);
    
    out.set_size(N,1);
    
    eT* out_mem = out.memptr();
    
    for(uword k=0; k < N; ++k)
      {
      eT acc = eT(0);
      
      // condition: A_n_cols = B_n_cols
      
      for(uword i=0; i < A_n_cols; ++i)
        {
        acc += A.at(k,i) * B.at(k,i);
        }
      
      out_mem[k] = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true ) && (partial_unwrap<T2>::do_trans == true ) )
    {
    arma_extra_debug_print("trans_A = true; trans_B = true;");
    
    const uword N = (std::min)(A_n_cols, B_n_rows);
    
    out.set_size(N,1);
    
    eT* out_mem = out.memptr();
    
    for(uword k=0; k < N; ++k)
      {
      eT acc = eT(0);
      
      const eT* A_colptr = A.colptr(k);
      
      // condition: A_n_rows = B_n_cols
      
      for(uword i=0; i < A_n_rows; ++i)
        {
        acc += A_colptr[i] * B.at(k,i);
        }
      
      out_mem[k] = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  
  if(is_alias)  { actual_out.steal_mem(tmp); }
  }



template<typename T1, typename T2>
inline
void
op_diagvec::apply(Mat<typename T1::elem_type>& actual_out, const Op< Glue<T1,T2,glue_times>, op_diagvec>& X, const typename arma_cx_only<typename T1::elem_type>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type   T;
  typedef typename T1::elem_type eT;
  
  const partial_unwrap<T1> UA(X.m.A);
  const partial_unwrap<T2> UB(X.m.B);
  
  const typename partial_unwrap<T1>::stored_type& A = UA.M;
  const typename partial_unwrap<T2>::stored_type& B = UB.M;
  
  arma_debug_assert_trans_mul_size< partial_unwrap<T1>::do_trans, partial_unwrap<T2>::do_trans >(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  if( (A.n_elem == 0) || (B.n_elem == 0) )  { actual_out.reset(); return; }
  
  const bool use_alpha = partial_unwrap<T1>::do_times || partial_unwrap<T2>::do_times;
  const eT       alpha = use_alpha ? (UA.get_val() * UB.get_val()) : eT(0);
  
  const bool is_alias  = (UA.is_alias(actual_out) || UB.is_alias(actual_out));
  
  Mat<eT>  tmp;
  Mat<eT>& out = (is_alias) ? tmp : actual_out;
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == false) )
    {
    arma_extra_debug_print("trans_A = false; trans_B = false;");
    
    const uword N = (std::min)(A_n_rows, B_n_cols);
    
    out.set_size(N,1);
    
    eT* out_mem = out.memptr();
    
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
      
      out_mem[k] = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true) && (partial_unwrap<T2>::do_trans == false) )
    {
    arma_extra_debug_print("trans_A = true; trans_B = false;");
    
    const uword N = (std::min)(A_n_cols, B_n_cols);
    
    out.set_size(N,1);
    
    eT* out_mem = out.memptr();
    
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
      
      out_mem[k] = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == true) )
    {
    arma_extra_debug_print("trans_A = false; trans_B = true;");
    
    const uword N = (std::min)(A_n_rows, B_n_rows);
    
    out.set_size(N,1);
    
    eT* out_mem = out.memptr();
    
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
      
      out_mem[k] = (use_alpha) ? eT(alpha * acc) : eT(acc);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true) && (partial_unwrap<T2>::do_trans == true) )
    {
    arma_extra_debug_print("trans_A = true; trans_B = true;");
    
    const uword N = (std::min)(A_n_cols, B_n_rows);
    
    out.set_size(N,1);
    
    eT* out_mem = out.memptr();
    
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
      
      out_mem[k] = (use_alpha) ? eT(alpha * acc) : eT(acc);
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
op_diagvec2::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_diagvec2>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword a = X.aux_uword_a;
  const uword b = X.aux_uword_b;
  
  const uword row_offset = (b >  0) ? a : 0;
  const uword col_offset = (b == 0) ? a : 0;
  
  const Proxy<T1> P(X.m);
  
  if(P.is_alias(out) == false)
    {
    op_diagvec2::apply_proxy(out, P, row_offset, col_offset);
    }
  else
    {
    Mat<eT> tmp;
    
    op_diagvec2::apply_proxy(tmp, P, row_offset, col_offset);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1>
inline
void
op_diagvec2::apply_proxy(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword row_offset, const uword col_offset)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  arma_debug_check_bounds
    (
    ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)),
    "diagvec(): requested diagonal is out of bounds"
    );
  
  const uword len = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  out.set_size(len, 1);
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < len; i+=2, j+=2)
    {
    const eT tmp_i = P.at( i + row_offset, i + col_offset );
    const eT tmp_j = P.at( j + row_offset, j + col_offset );
    
    out_mem[i] = tmp_i;
    out_mem[j] = tmp_j;
    }
  
  if(i < len)
    {
    out_mem[i] = P.at( i + row_offset, i + col_offset );
    }
  }



//! @}
