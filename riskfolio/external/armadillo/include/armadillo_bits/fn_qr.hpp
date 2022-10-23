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


//! \addtogroup fn_qr
//! @{



//! QR decomposition
template<typename T1>
inline
bool
qr
  (
         Mat<typename T1::elem_type>&    Q,
         Mat<typename T1::elem_type>&    R,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_debug_check( (&Q == &R), "qr(): Q and R are the same object" );
  
  const bool status = auxlib::qr(Q, R, X);
  
  if(status == false)
    {
    Q.soft_reset();
    R.soft_reset();
    arma_debug_warn_level(3, "qr(): decomposition failed");
    }
  
  return status;
  }



//! economical QR decomposition
template<typename T1>
inline
bool
qr_econ
  (
         Mat<typename T1::elem_type>&    Q,
         Mat<typename T1::elem_type>&    R,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_debug_check( (&Q == &R), "qr_econ(): Q and R are the same object" );
  
  const bool status = auxlib::qr_econ(Q, R, X);
  
  if(status == false)
    {
    Q.soft_reset();
    R.soft_reset();
    arma_debug_warn_level(3, "qr_econ(): decomposition failed");
    }
  
  return status;
  }



//! QR decomposition with pivoting
template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
qr
  (
         Mat<typename T1::elem_type>&    Q,
         Mat<typename T1::elem_type>&    R,
         Mat<uword>&                     P,
  const Base<typename T1::elem_type,T1>& X,
  const char*                            P_mode = "matrix"
  )
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (&Q == &R), "qr(): Q and R are the same object" );
  
  const char sig = (P_mode != nullptr) ? P_mode[0] : char(0);
  
  arma_debug_check( ((sig != 'm') && (sig != 'v')), "qr(): argument 'P_mode' must be \"vector\" or \"matrix\"" );
  
  bool status = false;
  
  if(sig == 'v')
    {
    status = auxlib::qr_pivot(Q, R, P, X);
    }
  else
  if(sig == 'm')
    {
    Mat<uword> P_vec;
    
    status = auxlib::qr_pivot(Q, R, P_vec, X);
    
    if(status)
      {
      // construct P
      
      const uword N = P_vec.n_rows;
      
      P.zeros(N,N);
      
      for(uword row=0; row < N; ++row)  { P.at(P_vec[row], row) = uword(1); }
      }
    }
  
  if(status == false)
    {
    Q.soft_reset();
    R.soft_reset();
    P.soft_reset();
    arma_debug_warn_level(3, "qr(): decomposition failed");
    }
  
  return status;
  }



//! @}
