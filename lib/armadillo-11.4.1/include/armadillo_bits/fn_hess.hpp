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


//! \addtogroup fn_hess
//! @{


template<typename T1>
inline
bool
hess
  (
         Mat<typename T1::elem_type>&    H,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  Col<eT> tao;
  
  const bool status = auxlib::hess(H, X.get_ref(), tao);
  
  if(H.n_rows > 2)
    {
    for(uword i=0; i < H.n_rows-2; ++i)
      {
      H(span(i+2, H.n_rows-1), i).zeros();
      }
    }
  
  if(status == false)
    {
    H.soft_reset();
    arma_debug_warn_level(3, "hess(): decomposition failed");
    }
  
  return status;
  }



template<typename T1>
arma_warn_unused
inline
Mat<typename T1::elem_type>
hess
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> H;
  Col<eT> tao;
  
  const bool status = auxlib::hess(H, X.get_ref(), tao);
  
  if(H.n_rows > 2)
    {
    for(uword i=0; i < H.n_rows-2; ++i)
      {
      H(span(i+2, H.n_rows-1), i).zeros();
      }
    }
  
  if(status == false)
    {
    H.soft_reset();
    arma_stop_runtime_error("hess(): decomposition failed");
    }
  
  return H;
  }



template<typename T1> 
inline
bool
hess
  (
         Mat<typename T1::elem_type>&    U,
         Mat<typename T1::elem_type>&    H,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_debug_check( void_ptr(&U) == void_ptr(&H), "hess(): 'U' is an alias of 'H'" );
  
  typedef typename T1::elem_type eT;
  
  Col<eT> tao;
  
  const bool status = auxlib::hess(H, X.get_ref(), tao);
  
  if(H.n_rows == 0)
    {
    U.reset();
    }
  else
  if(H.n_rows == 1)
    {
    U.ones(1, 1);
    }
  else
  if(H.n_rows == 2)
    {
    U.eye(2, 2);
    }
  else
    {
    U.eye(size(H));
    
    Col<eT> v;
    
    for(uword i=0; i < H.n_rows-2; ++i)
      {
      // TODO: generate v in a more efficient manner; 
      // TODO: the .ones() operation is an overkill, as most of v is overwritten afterwards
      
      v.ones(H.n_rows-i-1);
      
      v(span(1, H.n_rows-i-2)) = H(span(i+2, H.n_rows-1), i);
      
      U(span::all, span(i+1, H.n_rows-1)) -= tao(i) * (U(span::all, span(i+1, H.n_rows-1)) * v * v.t());
      }
    
    U(span::all, H.n_rows-1) = U(span::all, H.n_rows-1) * (eT(1) - tao(H.n_rows-2));
    
    for(uword i=0; i < H.n_rows-2; ++i)
      {
      H(span(i+2, H.n_rows-1), i).zeros();
      }
    }
  
  if(status == false)
    {
    U.soft_reset();
    H.soft_reset();
    arma_debug_warn_level(3, "hess(): decomposition failed");
    }
  
  return status;
  }



//! @}
