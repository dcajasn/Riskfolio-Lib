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


//! \addtogroup fn_trace
//! @{


template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
trace(const Base<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(X.get_ref());
  
  const uword N = (std::min)(P.get_n_rows(), P.get_n_cols());
  
  eT val1 = eT(0);
  eT val2 = eT(0);
  
  uword i,j;
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    val1 += P.at(i,i);
    val2 += P.at(j,j);
    }
  
  if(i < N)
    {
    val1 += P.at(i,i);
    }
  
  return val1 + val2;
  }



template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
trace(const Op<T1, op_diagmat>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const diagmat_proxy<T1> A(X.m);
  
  const uword N = (std::min)(A.n_rows, A.n_cols);
  
  eT val = eT(0);
  
  for(uword i=0; i<N; ++i)
    {
    val += A[i];
    }
  
  return val;
  }



//! speedup for trace(A*B); non-complex elements
template<typename T1, typename T2>
arma_warn_unused
inline
typename enable_if2< is_cx<typename T1::elem_type>::no, typename T1::elem_type>::result
trace(const Glue<T1, T2, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const partial_unwrap<T1> tmp1(X.A);
  const partial_unwrap<T2> tmp2(X.B);
  
  const typename partial_unwrap<T1>::stored_type& A = tmp1.M;
  const typename partial_unwrap<T2>::stored_type& B = tmp2.M;
  
  const bool use_alpha = partial_unwrap<T1>::do_times || partial_unwrap<T2>::do_times;
  const eT       alpha = use_alpha ? (tmp1.get_val() * tmp2.get_val()) : eT(0);
  
  arma_debug_assert_trans_mul_size< partial_unwrap<T1>::do_trans, partial_unwrap<T2>::do_trans >(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  if( (A.n_elem == 0) || (B.n_elem == 0) )  { return eT(0); }
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;

  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  eT acc = eT(0);
  
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == false) )
    {
    const uword N = (std::min)(A_n_rows, B_n_cols);
    
    eT acc1 = eT(0);
    eT acc2 = eT(0);
    
    for(uword k=0; k < N; ++k)
      {
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
      }
      
    acc = (acc1 + acc2);
    }
  else
  if( (partial_unwrap<T1>::do_trans == true ) && (partial_unwrap<T2>::do_trans == false) )
    {
    const uword N = (std::min)(A_n_cols, B_n_cols);
    
    for(uword k=0; k < N; ++k)
      {
      const eT* A_colptr = A.colptr(k);
      const eT* B_colptr = B.colptr(k);
      
      // condition: A_n_rows = B_n_rows
      acc += op_dot::direct_dot(A_n_rows, A_colptr, B_colptr);
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == true ) )
    {
    const uword N = (std::min)(A_n_rows, B_n_rows);
    
    for(uword k=0; k < N; ++k)
      {
      // condition: A_n_cols = B_n_cols
      for(uword i=0; i < A_n_cols; ++i)
        {
        acc += A.at(k,i) * B.at(k,i);
        }
      }
    }
  else
  if( (partial_unwrap<T1>::do_trans == true ) && (partial_unwrap<T2>::do_trans == true ) )
    {
    const uword N = (std::min)(A_n_cols, B_n_rows);
    
    for(uword k=0; k < N; ++k)
      {
      const eT* A_colptr = A.colptr(k);
      
      // condition: A_n_rows = B_n_cols
      for(uword i=0; i < A_n_rows; ++i)
        {
        acc += A_colptr[i] * B.at(k,i);
        }
      }
    }
  
  return (use_alpha) ? (alpha * acc) : acc;
  }



//! speedup for trace(A*B); complex elements
template<typename T1, typename T2>
arma_warn_unused
inline
typename enable_if2< is_cx<typename T1::elem_type>::yes, typename T1::elem_type>::result
trace(const Glue<T1, T2, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type   T;
  typedef typename T1::elem_type eT;
  
  const partial_unwrap<T1> tmp1(X.A);
  const partial_unwrap<T2> tmp2(X.B);
  
  const typename partial_unwrap<T1>::stored_type& A = tmp1.M;
  const typename partial_unwrap<T2>::stored_type& B = tmp2.M;
  
  const bool use_alpha = partial_unwrap<T1>::do_times || partial_unwrap<T2>::do_times;
  const eT       alpha = use_alpha ? (tmp1.get_val() * tmp2.get_val()) : eT(0);
  
  arma_debug_assert_trans_mul_size< partial_unwrap<T1>::do_trans, partial_unwrap<T2>::do_trans >(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  if( (A.n_elem == 0) || (B.n_elem == 0) )  { return eT(0); }
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  eT acc = eT(0);
  
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == false) )
    {
    const uword N = (std::min)(A_n_rows, B_n_cols);
    
    T acc_real = T(0);
    T acc_imag = T(0);
    
    for(uword k=0; k < N; ++k)
      {
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
      }
    
    acc = std::complex<T>(acc_real, acc_imag);
    }
  else
  if( (partial_unwrap<T1>::do_trans == true) && (partial_unwrap<T2>::do_trans == false) )
    {
    const uword N = (std::min)(A_n_cols, B_n_cols);
    
    T acc_real = T(0);
    T acc_imag = T(0);
    
    for(uword k=0; k < N; ++k)
      {
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
      }
    
    acc = std::complex<T>(acc_real, acc_imag);
    }
  else
  if( (partial_unwrap<T1>::do_trans == false) && (partial_unwrap<T2>::do_trans == true) )
    {
    const uword N = (std::min)(A_n_rows, B_n_rows);
    
    T acc_real = T(0);
    T acc_imag = T(0);
    
    for(uword k=0; k < N; ++k)
      {
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
      }
    
    acc = std::complex<T>(acc_real, acc_imag);
    }
  else
  if( (partial_unwrap<T1>::do_trans == true) && (partial_unwrap<T2>::do_trans == true) )
    {
    const uword N = (std::min)(A_n_cols, B_n_rows);
    
    T acc_real = T(0);
    T acc_imag = T(0);
    
    for(uword k=0; k < N; ++k)
      {
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
      }
    
    acc = std::complex<T>(acc_real, acc_imag);
    }
  
  return (use_alpha) ? eT(alpha * acc) : eT(acc);
  }



//! trace of sparse object; generic version
template<typename T1>
arma_warn_unused
inline
typename T1::elem_type
trace(const SpBase<typename T1::elem_type,T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> P(expr.get_ref());
  
  const uword N = (std::min)(P.get_n_rows(), P.get_n_cols());
  
  eT acc = eT(0);
  
  if( (is_SpMat<typename SpProxy<T1>::stored_type>::value) && (P.get_n_nonzero() >= 5*N) )
    {
    const unwrap_spmat<typename SpProxy<T1>::stored_type> U(P.Q);
    
    const SpMat<eT>& X = U.M;
    
    for(uword i=0; i < N; ++i)
      {
      acc += X.at(i,i);  // use binary search
      }
    }
  else
    {
    typename SpProxy<T1>::const_iterator_type it = P.begin();
    
    const uword P_n_nz = P.get_n_nonzero();
    
    for(uword i=0; i < P_n_nz; ++i)
      {
      if(it.row() == it.col())  { acc += (*it); }
      
      ++it;
      }
    }
  
  return acc;
  }



//! trace of sparse object; speedup for trace(A + B)
template<typename T1, typename T2>
arma_warn_unused
inline
typename T1::elem_type
trace(const SpGlue<T1, T2, spglue_plus>& expr)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> UA(expr.A);
  const unwrap_spmat<T2> UB(expr.B);
  
  arma_debug_assert_same_size(UA.M.n_rows, UA.M.n_cols, UB.M.n_rows, UB.M.n_cols, "addition");
  
  return (trace(UA.M) + trace(UB.M));
  }



//! trace of sparse object; speedup for trace(A - B)
template<typename T1, typename T2>
arma_warn_unused
inline
typename T1::elem_type
trace(const SpGlue<T1, T2, spglue_minus>& expr)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> UA(expr.A);
  const unwrap_spmat<T2> UB(expr.B);
  
  arma_debug_assert_same_size(UA.M.n_rows, UA.M.n_cols, UB.M.n_rows, UB.M.n_cols, "subtraction");
  
  return (trace(UA.M) - trace(UB.M));
  }



//! trace of sparse object; speedup for trace(A % B)
template<typename T1, typename T2>
arma_warn_unused
inline
typename T1::elem_type
trace(const SpGlue<T1, T2, spglue_schur>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> UA(expr.A);
  const unwrap_spmat<T2> UB(expr.B);
  
  const SpMat<eT>& A = UA.M;
  const SpMat<eT>& B = UB.M;
  
  arma_debug_assert_same_size(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "element-wise multiplication");
  
  const uword N = (std::min)(A.n_rows, A.n_cols);
  
  eT acc = eT(0);
  
  for(uword i=0; i<N; ++i)
    {
    acc += A.at(i,i) * B.at(i,i);
    }
  
  return acc;
  }



//! trace of sparse object; speedup for trace(A*B)
template<typename T1, typename T2>
arma_warn_unused
inline
typename T1::elem_type
trace(const SpGlue<T1, T2, spglue_times>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // better-than-nothing implementation
  
  const unwrap_spmat<T1> UA(expr.A);
  const unwrap_spmat<T2> UB(expr.B);
  
  const SpMat<eT>& A = UA.M;
  const SpMat<eT>& B = UB.M;
  
  arma_debug_assert_mul_size(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  if( (A.n_nonzero == 0) || (B.n_nonzero == 0) )  { return eT(0); }
  
  const uword N = (std::min)(A.n_rows, B.n_cols);
  
  eT acc = eT(0);
  
  // TODO: the threshold may need tuning for complex matrices
  if( (A.n_nonzero >= 5*N) || (B.n_nonzero >= 5*N) )
    {
    for(uword k=0; k < N; ++k)
      {
      typename SpMat<eT>::const_col_iterator B_it     = B.begin_col_no_sync(k);
      typename SpMat<eT>::const_col_iterator B_it_end = B.end_col_no_sync(k);
      
      while(B_it != B_it_end)
        {
        const eT    B_val = (*B_it);
        const uword i     = B_it.row();
        
        acc += A.at(k,i) * B_val;
        
        ++B_it;
        }
      }
    }
  else
    {
    const SpMat<eT> AB = A * B;
    
    acc = trace(AB);
    }
  
  return acc;
  }



//! trace of sparse object; speedup for trace(A.t()*B); non-complex elements
template<typename T1, typename T2>
arma_warn_unused
inline
typename enable_if2< is_cx<typename T1::elem_type>::no, typename T1::elem_type>::result
trace(const SpGlue<SpOp<T1, spop_htrans>, T2, spglue_times>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> UA(expr.A.m);
  const unwrap_spmat<T2> UB(expr.B);
  
  const SpMat<eT>& A = UA.M;
  const SpMat<eT>& B = UB.M;
  
  // NOTE: deliberately swapped A.n_rows and A.n_cols to take into account the requested transpose operation
  arma_debug_assert_mul_size(A.n_cols, A.n_rows, B.n_rows, B.n_cols, "matrix multiplication");
  
  if( (A.n_nonzero == 0) || (B.n_nonzero == 0) )  { return eT(0); }
  
  const uword N = (std::min)(A.n_cols, B.n_cols);
  
  eT acc = eT(0);
  
  if( (A.n_nonzero >= 5*N) || (B.n_nonzero >= 5*N) )
    {
    for(uword k=0; k < N; ++k)
      {
      typename SpMat<eT>::const_col_iterator B_it     = B.begin_col_no_sync(k);
      typename SpMat<eT>::const_col_iterator B_it_end = B.end_col_no_sync(k);
      
      while(B_it != B_it_end)
        {
        const eT    B_val = (*B_it);
        const uword i     = B_it.row();
        
        acc += A.at(i,k) * B_val;
        
        ++B_it;
        }
      }
    }
  else
    {
    const SpMat<eT> AtB = A.t() * B;
    
    acc = trace(AtB);
    }
  
  return acc;
  }



//! trace of sparse object; speedup for trace(A.t()*B); complex elements
template<typename T1, typename T2>
arma_warn_unused
inline
typename enable_if2< is_cx<typename T1::elem_type>::yes, typename T1::elem_type>::result
trace(const SpGlue<SpOp<T1, spop_htrans>, T2, spglue_times>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> UA(expr.A.m);
  const unwrap_spmat<T2> UB(expr.B);
  
  const SpMat<eT>& A = UA.M;
  const SpMat<eT>& B = UB.M;
  
  // NOTE: deliberately swapped A.n_rows and A.n_cols to take into account the requested transpose operation
  arma_debug_assert_mul_size(A.n_cols, A.n_rows, B.n_rows, B.n_cols, "matrix multiplication");
  
  if( (A.n_nonzero == 0) || (B.n_nonzero == 0) )  { return eT(0); }
  
  const uword N = (std::min)(A.n_cols, B.n_cols);
  
  eT acc = eT(0);
  
  // TODO: the threshold may need tuning for complex matrices
  if( (A.n_nonzero >= 5*N) || (B.n_nonzero >= 5*N) )
    {
    for(uword k=0; k < N; ++k)
      {
      typename SpMat<eT>::const_col_iterator B_it     = B.begin_col_no_sync(k);
      typename SpMat<eT>::const_col_iterator B_it_end = B.end_col_no_sync(k);
      
      while(B_it != B_it_end)
        {
        const eT    B_val = (*B_it);
        const uword i     = B_it.row();
        
        acc += std::conj(A.at(i,k)) * B_val;
        
        ++B_it;
        }
      }
    }
  else
    {
    const SpMat<eT> AtB = A.t() * B;
    
    acc = trace(AtB);
    }
  
  return acc;
  }



//! @}
