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



//! \addtogroup op_clamp
//! @{



template<typename T1>
inline
void
op_clamp::apply(Mat<typename T1::elem_type>& out, const mtOp<typename T1::elem_type, T1, op_clamp>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const eT min_val = in.aux;
  const eT max_val = in.aux_out_eT;
  
  arma_debug_check( (min_val > max_val), "clamp(): min_val must be less than max_val" );
  
  if(is_Mat<T1>::value)
    {
    const unwrap<T1> U(in.m);
    
    op_clamp::apply_direct(out, U.M, min_val, max_val);
    }
  else
    {
    const Proxy<T1> P(in.m);
    
    if(P.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_clamp::apply_proxy_noalias(tmp, P, min_val, max_val);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_clamp::apply_proxy_noalias(out, P, min_val, max_val);
      }
    }
  }



template<typename eT>
inline
void
op_clamp::apply_direct(Mat<eT>& out, const Mat<eT>& X, const eT min_val, const eT max_val)
  {
  arma_extra_debug_sigprint();
  
  if(&out != &X)
    {
    out.set_size(X.n_rows, X.n_cols);
    
    const uword N = out.n_elem;
    
    const eT*   X_mem =   X.memptr();
          eT* out_mem = out.memptr();
    
    for(uword i=0; i<N; ++i)
      {
      const eT val = X_mem[i];
      
      out_mem[i] = (val < min_val) ? min_val : ((val > max_val) ? max_val : val);
      }
    }
  else
    {
    arma_extra_debug_print("op_clamp::apply_direct(): inplace operation");
    
    arrayops::clamp(out.memptr(), out.n_elem, min_val, max_val);
    }
  }



template<typename T1>
inline
void
op_clamp::apply_proxy_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const typename T1::elem_type min_val, const typename T1::elem_type max_val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
  eT* out_mem = out.memptr();
  
  if(Proxy<T1>::use_at == false)
    {
    const uword N = P.get_n_elem();
    
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    for(uword i=0; i<N; ++i)
      {
      const eT val = A[i];
      
      out_mem[i] = (val < min_val) ? min_val : ((val > max_val) ? max_val : val);
      }
    }
  else
    {
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      const eT val = P.at(row,col);
      
      (*out_mem) = (val < min_val) ? min_val : ((val > max_val) ? max_val : val);
      
      out_mem++;
      }
    }
  }



//



template<typename T1>
inline
void
op_clamp::apply(Cube<typename T1::elem_type>& out, const mtOpCube<typename T1::elem_type, T1, op_clamp>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const eT min_val = in.aux;
  const eT max_val = in.aux_out_eT;
  
  arma_debug_check( (min_val > max_val), "clamp(): min_val must be less than max_val" );
  
  if(is_Cube<T1>::value)
    {
    const unwrap_cube<T1> U(in.m);
    
    op_clamp::apply_direct(out, U.M, min_val, max_val);
    }
  else
    {
    const ProxyCube<T1> P(in.m);
    
    if(P.is_alias(out))
      {
      Cube<eT> tmp;
      
      op_clamp::apply_proxy_noalias(tmp, P, min_val, max_val);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_clamp::apply_proxy_noalias(out, P, min_val, max_val);
      }
    }
  }



template<typename eT>
inline
void
op_clamp::apply_direct(Cube<eT>& out, const Cube<eT>& X, const eT min_val, const eT max_val)
  {
  arma_extra_debug_sigprint();
  
  if(&out != &X)
    {
    out.set_size(X.n_rows, X.n_cols, X.n_slices);
    
    const uword N = out.n_elem;
    
    const eT*   X_mem =   X.memptr();
          eT* out_mem = out.memptr();
    
    for(uword i=0; i<N; ++i)
      {
      const eT val = X_mem[i];
      
      out_mem[i] = (val < min_val) ? min_val : ((val > max_val) ? max_val : val);
      }
    }
  else
    {
    arma_extra_debug_print("op_clamp::apply_direct(): inplace operation");
    
    arrayops::clamp(out.memptr(), out.n_elem, min_val, max_val);
    }
  }



template<typename T1>
inline
void
op_clamp::apply_proxy_noalias(Cube<typename T1::elem_type>& out, const ProxyCube<T1>& P, const typename T1::elem_type min_val, const typename T1::elem_type max_val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword n_rows   = P.get_n_rows();
  const uword n_cols   = P.get_n_cols();
  const uword n_slices = P.get_n_slices();
  
  out.set_size(n_rows, n_cols, n_slices);
  
  eT* out_mem = out.memptr();
  
  if(ProxyCube<T1>::use_at == false)
    {
    const uword N = P.get_n_elem();
    
    typename ProxyCube<T1>::ea_type A = P.get_ea();
    
    for(uword i=0; i<N; ++i)
      {
      const eT val = A[i];
      
      out_mem[i] = (val < min_val) ? min_val : ((val > max_val) ? max_val : val);
      }
    }
  else
    {
    for(uword s=0; s < n_slices; ++s)
    for(uword c=0; c < n_cols;   ++c)
    for(uword r=0; r < n_rows;   ++r)
      {
      const eT val = P.at(r,c,s);
      
      (*out_mem) = (val < min_val) ? min_val : ((val > max_val) ? max_val : val);
      
      out_mem++;
      }
    }
  }



//



template<typename T1>
inline
void
op_clamp_cx::apply(Mat<typename T1::elem_type>& out, const mtOp<typename T1::elem_type, T1, op_clamp_cx>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_Mat<T1>::value)
    {
    const unwrap<T1> U(in.m);
    
    op_clamp_cx::apply_direct(out, U.M, in.aux, in.aux_out_eT);
    }
  else
    {
    const Proxy<T1> P(in.m);
    
    if(P.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_clamp_cx::apply_proxy_noalias(tmp, P, in.aux, in.aux_out_eT);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_clamp_cx::apply_proxy_noalias(out, P, in.aux, in.aux_out_eT);
      }
    }
  }



template<typename eT>
inline
void
op_clamp_cx::apply_direct(Mat<eT>& out, const Mat<eT>& X, const eT min_val, const eT max_val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const T min_val_real = std::real(min_val);
  const T min_val_imag = std::imag(min_val);
  
  const T max_val_real = std::real(max_val);
  const T max_val_imag = std::imag(max_val);
  
  arma_debug_check( (min_val_real > max_val_real), "clamp(): real(min_val) must be less than real(max_val)" );
  arma_debug_check( (min_val_imag > max_val_imag), "clamp(): imag(min_val) must be less than imag(max_val)" );
  
  if(&out != &X)
    {
    out.set_size(X.n_rows, X.n_cols);
    
    const uword N = out.n_elem;
    
    const eT*   X_mem =   X.memptr();
          eT* out_mem = out.memptr();
    
    for(uword i=0; i<N; ++i)
      {
      const eT& X_val = X_mem[i];
      
      T val_real = std::real(X_val);
      T val_imag = std::imag(X_val);
      
      val_real = (val_real < min_val_real) ? min_val_real : ((val_real > max_val_real) ? max_val_real : val_real);
      val_imag = (val_imag < min_val_imag) ? min_val_imag : ((val_imag > max_val_imag) ? max_val_imag : val_imag);
      
      out_mem[i] = std::complex<T>(val_real,val_imag);
      }
    }
  else
    {
    arma_extra_debug_print("op_clamp_cx::apply_direct(): inplace operation");
    
    arrayops::clamp(out.memptr(), out.n_elem, min_val, max_val);
    }
  }



template<typename T1>
inline
void
op_clamp_cx::apply_proxy_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const typename T1::elem_type min_val, const typename T1::elem_type max_val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const T min_val_real = std::real(min_val);
  const T min_val_imag = std::imag(min_val);
  
  const T max_val_real = std::real(max_val);
  const T max_val_imag = std::imag(max_val);
  
  arma_debug_check( (min_val_real > max_val_real), "clamp(): real(min_val) must be less than real(max_val)" );
  arma_debug_check( (min_val_imag > max_val_imag), "clamp(): imag(min_val) must be less than imag(max_val)" );
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
  eT* out_mem = out.memptr();
  
  if(Proxy<T1>::use_at == false)
    {
    const uword N = P.get_n_elem();
    
    typename Proxy<T1>::ea_type A = P.get_ea();
    
    for(uword i=0; i<N; ++i)
      {
      T val_real = std::real(A[i]);
      T val_imag = std::imag(A[i]);
      
      val_real = (val_real < min_val_real) ? min_val_real : ((val_real > max_val_real) ? max_val_real : val_real);
      val_imag = (val_imag < min_val_imag) ? min_val_imag : ((val_imag > max_val_imag) ? max_val_imag : val_imag);
      
      out_mem[i] = std::complex<T>(val_real,val_imag);
      }
    }
  else
    {
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      const eT val = P.at(row,col);
      
      T val_real = std::real(val);
      T val_imag = std::imag(val);
      
      val_real = (val_real < min_val_real) ? min_val_real : ((val_real > max_val_real) ? max_val_real : val_real);
      val_imag = (val_imag < min_val_imag) ? min_val_imag : ((val_imag > max_val_imag) ? max_val_imag : val_imag);
      
      (*out_mem) = std::complex<T>(val_real,val_imag);  out_mem++;
      }
    }
  }



//



template<typename T1>
inline
void
op_clamp_cx::apply(Cube<typename T1::elem_type>& out, const mtOpCube<typename T1::elem_type, T1, op_clamp_cx>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_Cube<T1>::value)
    {
    const unwrap_cube<T1> U(in.m);
    
    op_clamp_cx::apply_direct(out, U.M, in.aux, in.aux_out_eT);
    }
  else
    {
    const ProxyCube<T1> P(in.m);
    
    if(P.is_alias(out))
      {
      Cube<eT> tmp;
      
      op_clamp_cx::apply_proxy_noalias(tmp, P, in.aux, in.aux_out_eT);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_clamp_cx::apply_proxy_noalias(out, P, in.aux, in.aux_out_eT);
      }
    }
  }



template<typename eT>
inline
void
op_clamp_cx::apply_direct(Cube<eT>& out, const Cube<eT>& X, const eT min_val, const eT max_val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const T min_val_real = std::real(min_val);
  const T min_val_imag = std::imag(min_val);
  
  const T max_val_real = std::real(max_val);
  const T max_val_imag = std::imag(max_val);
  
  arma_debug_check( (min_val_real > max_val_real), "clamp(): real(min_val) must be less than real(max_val)" );
  arma_debug_check( (min_val_imag > max_val_imag), "clamp(): imag(min_val) must be less than imag(max_val)" );
  
  if(&out != &X)
    {
    out.set_size(X.n_rows, X.n_cols, X.n_slices);
    
    const uword N = out.n_elem;
    
    const eT*   X_mem =   X.memptr();
          eT* out_mem = out.memptr();
    
    for(uword i=0; i<N; ++i)
      {
      const eT& X_val = X_mem[i];
      
      T val_real = std::real(X_val);
      T val_imag = std::imag(X_val);
      
      val_real = (val_real < min_val_real) ? min_val_real : ((val_real > max_val_real) ? max_val_real : val_real);
      val_imag = (val_imag < min_val_imag) ? min_val_imag : ((val_imag > max_val_imag) ? max_val_imag : val_imag);
      
      out_mem[i] = std::complex<T>(val_real,val_imag);
      }
    }
  else
    {
    arma_extra_debug_print("op_clamp_cx::apply_direct(): inplace operation");
    
    arrayops::clamp(out.memptr(), out.n_elem, min_val, max_val);
    }
  }



template<typename T1>
inline
void
op_clamp_cx::apply_proxy_noalias(Cube<typename T1::elem_type>& out, const ProxyCube<T1>& P, const typename T1::elem_type min_val, const typename T1::elem_type max_val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const T min_val_real = std::real(min_val);
  const T min_val_imag = std::imag(min_val);
  
  const T max_val_real = std::real(max_val);
  const T max_val_imag = std::imag(max_val);
  
  arma_debug_check( (min_val_real > max_val_real), "clamp(): real(min_val) must be less than real(max_val)" );
  arma_debug_check( (min_val_imag > max_val_imag), "clamp(): imag(min_val) must be less than imag(max_val)" );
  
  const uword n_rows   = P.get_n_rows();
  const uword n_cols   = P.get_n_cols();
  const uword n_slices = P.get_n_slices();
  
  out.set_size(n_rows, n_cols, n_slices);
  
  eT* out_mem = out.memptr();
  
  if(ProxyCube<T1>::use_at == false)
    {
    const uword N = P.get_n_elem();
    
    typename ProxyCube<T1>::ea_type A = P.get_ea();
    
    for(uword i=0; i<N; ++i)
      {
      T val_real = std::real(A[i]);
      T val_imag = std::imag(A[i]);
      
      val_real = (val_real < min_val_real) ? min_val_real : ((val_real > max_val_real) ? max_val_real : val_real);
      val_imag = (val_imag < min_val_imag) ? min_val_imag : ((val_imag > max_val_imag) ? max_val_imag : val_imag);
      
      out_mem[i] = std::complex<T>(val_real,val_imag);
      }
    }
  else
    {
    for(uword s=0; s < n_slices; ++s)
    for(uword c=0; c < n_cols;   ++c)
    for(uword r=0; r < n_rows;   ++r)
      {
      const eT val = P.at(r,c,s);
      
      T val_real = std::real(val);
      T val_imag = std::imag(val);
      
      val_real = (val_real < min_val_real) ? min_val_real : ((val_real > max_val_real) ? max_val_real : val_real);
      val_imag = (val_imag < min_val_imag) ? min_val_imag : ((val_imag > max_val_imag) ? max_val_imag : val_imag);
      
      (*out_mem) = std::complex<T>(val_real,val_imag);  out_mem++;
      }
    }
  }



//! @}
