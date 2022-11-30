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


//! \addtogroup fn_misc
//! @{



template<typename out_type>
arma_warn_unused
inline
typename
enable_if2
  <
  is_Mat<out_type>::value,
  out_type
  >::result
linspace
  (
  const typename out_type::pod_type start,
  const typename out_type::pod_type end,
  const uword                       num = 100u
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename out_type::elem_type eT;
  typedef typename out_type::pod_type   T;
  
  out_type x;
  
  if(num == 1)
    {
    x.set_size(1);
    
    x[0] = eT(end);
    }
  else
  if(num >= 2)
    {
    x.set_size(num);
    
    eT* x_mem = x.memptr();
    
    const uword num_m1 = num - 1;
    
    if(is_non_integral<T>::value)
      {
      const T delta = (end-start)/T(num_m1);
      
      for(uword i=0; i<num_m1; ++i)
        {
        x_mem[i] = eT(start + i*delta);
        }
      
      x_mem[num_m1] = eT(end);
      }
    else
      {
      const double delta = (end >= start) ? double(end-start)/double(num_m1) : -double(start-end)/double(num_m1);
      
      for(uword i=0; i<num_m1; ++i)
        {
        x_mem[i] = eT(double(start) + i*delta);
        }
      
      x_mem[num_m1] = eT(end);
      }
    }
  
  return x;
  }



arma_warn_unused
inline
vec
linspace(const double start, const double end, const uword num = 100u)
  {
  arma_extra_debug_sigprint();
  return linspace<vec>(start, end, num);
  }



template<typename out_type>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_Mat<out_type>::value && is_real<typename out_type::pod_type>::value),
  out_type
  >::result
logspace
  (
  const typename out_type::pod_type A,
  const typename out_type::pod_type B,
  const uword                       N = 50u
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename out_type::elem_type eT;
  typedef typename out_type::pod_type   T;
  
  out_type x = linspace<out_type>(A,B,N);
  
  const uword n_elem = x.n_elem;
  
  eT* x_mem = x.memptr();
  
  for(uword i=0; i < n_elem; ++i)
    {
    x_mem[i] = std::pow(T(10), x_mem[i]);
    }
  
  return x;
  }



arma_warn_unused
inline
vec
logspace(const double A, const double B, const uword N = 50u)
  {
  arma_extra_debug_sigprint();
  return logspace<vec>(A, B, N);
  }



//
// log_exp_add

template<typename eT>
arma_warn_unused
inline
typename arma_real_only<eT>::result
log_add_exp(eT log_a, eT log_b)
  {
  if(log_a < log_b)
    {
    std::swap(log_a, log_b);
    }
  
  const eT negdelta = log_b - log_a;
  
  if( (negdelta < Datum<eT>::log_min) || (arma_isfinite(negdelta) == false) )
    {
    return log_a;
    }
  else
    {
    return (log_a + std::log1p(std::exp(negdelta)));
    }
  }



// for compatibility with earlier versions
template<typename eT>
arma_warn_unused
inline
typename arma_real_only<eT>::result
log_add(eT log_a, eT log_b)
  {
  return log_add_exp(log_a, log_b);
  }
  


//! kept for compatibility with old user code
template<typename eT>
arma_warn_unused
arma_inline
bool
is_finite(const eT x, const typename arma_scalar_only<eT>::result* junk = nullptr)
  {
  arma_ignore(junk);
  
  return arma_isfinite(x);
  }



//! kept for compatibility with old user code
template<typename T1>
arma_warn_unused
inline
bool
is_finite(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return X.is_finite();
  }



//! kept for compatibility with old user code
template<typename T1>
arma_warn_unused
inline
bool
is_finite(const SpBase<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return X.is_finite();
  }



//! kept for compatibility with old user code
template<typename T1>
arma_warn_unused
inline
bool
is_finite(const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return X.is_finite();
  }



template<typename eT>
inline
void
swap(Mat<eT>& A, Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  A.swap(B);
  }



template<typename eT>
inline
void
swap(Cube<eT>& A, Cube<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  A.swap(B);
  }



arma_warn_unused
inline
uvec
ind2sub(const SizeMat& s, const uword i)
  {
  arma_extra_debug_sigprint();
  
  const uword s_n_rows = s.n_rows;
  
  arma_debug_check( (i >= (s_n_rows * s.n_cols) ), "ind2sub(): index out of range" );
  
  const uword row = i % s_n_rows;
  const uword col = i / s_n_rows;
  
  uvec out(2, arma_nozeros_indicator());
  
  uword* out_mem = out.memptr();
  
  out_mem[0] = row;
  out_mem[1] = col;
  
  return out;
  }



template<typename T1>
arma_warn_unused
inline
typename enable_if2< (is_arma_type<T1>::value && is_same_type<uword,typename T1::elem_type>::yes), umat >::result
ind2sub(const SizeMat& s, const T1& indices)
  {
  arma_extra_debug_sigprint();
  
  const uword s_n_rows = s.n_rows;
  const uword s_n_elem = s_n_rows * s.n_cols;
  
  const Proxy<T1> P(indices);
  
  const uword P_n_rows = P.get_n_rows();
  const uword P_n_cols = P.get_n_cols();
  const uword P_n_elem = P.get_n_elem();
  
  const bool P_is_empty = (P_n_elem == 0);
  const bool P_is_vec   = ((P_n_rows == 1) || (P_n_cols == 1));
  
  arma_debug_check( ((P_is_empty == false) && (P_is_vec == false)), "ind2sub(): parameter 'indices' must be a vector" );
  
  umat out(2, P_n_elem, arma_nozeros_indicator());
  
  if(Proxy<T1>::use_at == false)
    {
    typename Proxy<T1>::ea_type Pea = P.get_ea();
    
    for(uword count=0; count < P_n_elem; ++count)
      {
      const uword i = Pea[count];
      
      arma_debug_check( (i >= s_n_elem), "ind2sub(): index out of range" );
      
      const uword row = i % s_n_rows;
      const uword col = i / s_n_rows;
      
      uword* out_colptr = out.colptr(count);
      
      out_colptr[0] = row;
      out_colptr[1] = col;
      }
    }
  else
    {
    if(P_n_rows == 1)
      {
      for(uword count=0; count < P_n_cols; ++count)
        {
        const uword i = P.at(0,count);
        
        arma_debug_check( (i >= s_n_elem), "ind2sub(): index out of range" );
        
        const uword row = i % s_n_rows;
        const uword col = i / s_n_rows;
        
        uword* out_colptr = out.colptr(count);
        
        out_colptr[0] = row;
        out_colptr[1] = col;
        }
      }
    else
    if(P_n_cols == 1)
      {
      for(uword count=0; count < P_n_rows; ++count)
        {
        const uword i = P.at(count,0);
        
        arma_debug_check( (i >= s_n_elem), "ind2sub(): index out of range" );
        
        const uword row = i % s_n_rows;
        const uword col = i / s_n_rows;
        
        uword* out_colptr = out.colptr(count);
        
        out_colptr[0] = row;
        out_colptr[1] = col;
        }
      }
    }
  
  return out;
  }



arma_warn_unused
inline
uvec
ind2sub(const SizeCube& s, const uword i)
  {
  arma_extra_debug_sigprint();
  
  const uword s_n_rows       = s.n_rows;
  const uword s_n_elem_slice = s_n_rows * s.n_cols;
  
  arma_debug_check( (i >= (s_n_elem_slice * s.n_slices) ), "ind2sub(): index out of range" );
  
  const uword slice  = i / s_n_elem_slice;
  const uword j      = i - (slice * s_n_elem_slice);
  const uword row    = j % s_n_rows;
  const uword col    = j / s_n_rows;
  
  uvec out(3, arma_nozeros_indicator());
  
  uword* out_mem = out.memptr();
  
  out_mem[0] = row;
  out_mem[1] = col;
  out_mem[2] = slice;
  
  return out;
  }



template<typename T1>
arma_warn_unused
inline
typename enable_if2< (is_arma_type<T1>::value && is_same_type<uword,typename T1::elem_type>::yes), umat >::result
ind2sub(const SizeCube& s, const T1& indices)
  {
  arma_extra_debug_sigprint();
  
  const uword s_n_rows       = s.n_rows;
  const uword s_n_elem_slice = s_n_rows * s.n_cols;
  const uword s_n_elem       = s.n_slices * s_n_elem_slice;
    
  const quasi_unwrap<T1> U(indices);
  
  arma_debug_check( ((U.M.is_empty() == false) && (U.M.is_vec() == false)), "ind2sub(): parameter 'indices' must be a vector" );
  
  const uword  U_n_elem = U.M.n_elem;
  const uword* U_mem    = U.M.memptr();
  
  umat out(3, U_n_elem, arma_nozeros_indicator());
  
  for(uword count=0; count < U_n_elem; ++count)
    {
    const uword i = U_mem[count];
    
    arma_debug_check( (i >= s_n_elem), "ind2sub(): index out of range" );
    
    const uword slice  = i / s_n_elem_slice;
    const uword j      = i - (slice * s_n_elem_slice);
    const uword row    = j % s_n_rows;
    const uword col    = j / s_n_rows;
    
    uword* out_colptr = out.colptr(count);
    
    out_colptr[0] = row;
    out_colptr[1] = col;
    out_colptr[2] = slice;
    }
  
  return out;
  }



arma_warn_unused
arma_inline
uword
sub2ind(const SizeMat& s, const uword row, const uword col)
  {
  arma_extra_debug_sigprint();
  
  const uword s_n_rows = s.n_rows;
  
  arma_debug_check( ((row >= s_n_rows) || (col >= s.n_cols)), "sub2ind(): subscript out of range" );
  
  return uword(row + col*s_n_rows);
  }



template<typename T1>
arma_warn_unused
inline
uvec
sub2ind(const SizeMat& s, const Base<uword,T1>& subscripts)
  {
  arma_extra_debug_sigprint();
  
  const uword s_n_rows = s.n_rows;
  const uword s_n_cols = s.n_cols;
  
  const quasi_unwrap<T1> U(subscripts.get_ref());
  
  arma_debug_check( (U.M.n_rows != 2), "sub2ind(): matrix of subscripts must have 2 rows" );
  
  const uword U_M_n_cols = U.M.n_cols;
  
  uvec out(U_M_n_cols, arma_nozeros_indicator());
  
        uword* out_mem = out.memptr();
  const uword* U_M_mem = U.M.memptr();
  
  for(uword count=0; count < U_M_n_cols; ++count)
    {
    const uword row = U_M_mem[0];
    const uword col = U_M_mem[1];
    
    U_M_mem += 2; // next column
    
    arma_debug_check( ((row >= s_n_rows) || (col >= s_n_cols)), "sub2ind(): subscript out of range" );
    
    out_mem[count] = uword(row + col*s_n_rows);
    }
  
  return out;
  }



arma_warn_unused
arma_inline
uword
sub2ind(const SizeCube& s, const uword row, const uword col, const uword slice)
  {
  arma_extra_debug_sigprint();
  
  const uword s_n_rows = s.n_rows;
  const uword s_n_cols = s.n_cols;
  
  arma_debug_check( ((row >= s_n_rows) || (col >= s_n_cols) || (slice >= s.n_slices)), "sub2ind(): subscript out of range" );
  
  return uword( (slice * s_n_rows * s_n_cols) + (col * s_n_rows) + row );
  }



template<typename T1>
arma_warn_unused
inline
uvec
sub2ind(const SizeCube& s, const Base<uword,T1>& subscripts)
  {
  arma_extra_debug_sigprint();
  
  const uword s_n_rows   = s.n_rows;
  const uword s_n_cols   = s.n_cols;
  const uword s_n_slices = s.n_slices;
  
  const quasi_unwrap<T1> U(subscripts.get_ref());
  
  arma_debug_check( (U.M.n_rows != 3), "sub2ind(): matrix of subscripts must have 3 rows" );
  
  const uword U_M_n_cols = U.M.n_cols;
  
  uvec out(U_M_n_cols, arma_nozeros_indicator());
  
        uword* out_mem = out.memptr();
  const uword* U_M_mem = U.M.memptr();
  
  for(uword count=0; count < U_M_n_cols; ++count)
    {
    const uword row   = U_M_mem[0];
    const uword col   = U_M_mem[1];
    const uword slice = U_M_mem[2];
    
    U_M_mem += 3; // next column
    
    arma_debug_check( ((row >= s_n_rows) || (col >= s_n_cols) || (slice >= s_n_slices)), "sub2ind(): subscript out of range" );
    
    out_mem[count] = uword( (slice * s_n_rows * s_n_cols) + (col * s_n_rows) + row );
    }
  
  return out;
  }



template<typename T1, typename T2>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const Glue<T1, T2, glue_affmul>
  >::result
affmul(const T1& A, const T2& B)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_affmul>(A,B);
  }



//! @}
