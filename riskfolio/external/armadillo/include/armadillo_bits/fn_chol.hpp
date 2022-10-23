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


//! \addtogroup fn_chol
//! @{



template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Op<T1, op_chol> >::result
chol
  (
  const Base<typename T1::elem_type,T1>& X,
  const char* layout = "upper"
  )
  {
  arma_extra_debug_sigprint();
  
  const char sig = (layout != nullptr) ? layout[0] : char(0);
  
  arma_debug_check( ((sig != 'u') && (sig != 'l')), "chol(): layout must be \"upper\" or \"lower\"" );
  
  return Op<T1, op_chol>(X.get_ref(), ((sig == 'u') ? 0 : 1), 0 );
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
chol
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& X,
  const char* layout = "upper"
  )
  {
  arma_extra_debug_sigprint();
  
  const char sig = (layout != nullptr) ? layout[0] : char(0);
  
  arma_debug_check( ((sig != 'u') && (sig != 'l')), "chol(): layout must be \"upper\" or \"lower\"" );
  
  const bool status = op_chol::apply_direct(out, X.get_ref(), ((sig == 'u') ? 0 : 1));
  
  if(status == false)
    {
    out.soft_reset();
    arma_debug_warn_level(3, "chol(): decomposition failed");
    }
  
  return status;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
chol
  (
         Mat<typename T1::elem_type>&    out,
         Mat<uword>&                     P,
  const Base<typename T1::elem_type,T1>& X,
  const char*                            layout = "upper",
  const char*                            P_mode = "matrix"
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const char sig_layout = (layout != nullptr) ? layout[0] : char(0);
  const char sig_P_mode = (P_mode != nullptr) ? P_mode[0] : char(0);
  
  arma_debug_check( ((sig_layout != 'u') && (sig_layout != 'l')), "chol(): argument 'layout' must be \"upper\" or \"lower\""   );
  arma_debug_check( ((sig_P_mode != 'm') && (sig_P_mode != 'v')), "chol(): argument 'P_mode' must be \"vector\" or \"matrix\"" );
  
  out = X.get_ref();
  
  arma_debug_check( (out.is_square() == false), "chol(): given matrix must be square sized", [&](){ out.soft_reset(); } );
  
  if(out.is_empty())
    {
    P.reset();
    return true;
    }
  
  if((arma_config::debug) && (auxlib::rudimentary_sym_check(out) == false))
    {
    if(is_cx<eT>::no )  { arma_debug_warn_level(1, "chol(): given matrix is not symmetric"); }
    if(is_cx<eT>::yes)  { arma_debug_warn_level(1, "chol(): given matrix is not hermitian"); }
    }
  
  bool status = false;
  
  if(sig_P_mode == 'v')
    {
    status = auxlib::chol_pivot(out, P, ((sig_layout == 'u') ? 0 : 1));
    }
  else
  if(sig_P_mode == 'm')
    {
    Mat<uword> P_vec;
    
    status = auxlib::chol_pivot(out, P_vec, ((sig_layout == 'u') ? 0 : 1));
    
    if(status)
      {
      // construct P
      
      const uword N = P_vec.n_rows;
      
      P.zeros(N,N);
      
      for(uword i=0; i < N; ++i)  { P.at(P_vec[i], i) = uword(1); }
      }
    }
  
  if(status == false)
    {
    out.soft_reset();
      P.soft_reset();
    arma_debug_warn_level(3, "chol(): decomposition failed");
    }
  
  return status;
  }



//! @}
