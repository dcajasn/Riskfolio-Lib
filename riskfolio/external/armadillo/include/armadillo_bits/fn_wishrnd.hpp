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


//! \addtogroup fn_wishrnd
//! @{



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_real<typename T1::elem_type>::value,
  const Op<T1, op_wishrnd>
  >::result
wishrnd(const Base<typename T1::elem_type, T1>& S, typename T1::elem_type df)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_wishrnd>(S.get_ref(), df, uword(1), uword(0));
  }



template<typename T1, typename T2>
arma_warn_unused
inline
typename
enable_if2
  <
  is_real<typename T1::elem_type>::value,
  const Op<T2, op_wishrnd>
  >::result
wishrnd(const Base<typename T1::elem_type, T1>& S, typename T1::elem_type df, const Base<typename T1::elem_type, T2>& D)
  {
  arma_extra_debug_sigprint();
  arma_ignore(S);
  
  return Op<T2, op_wishrnd>(D.get_ref(), df, uword(2), uword(0));
  }



template<typename T1>
inline
typename
enable_if2
  <
  is_real<typename T1::elem_type>::value,
  bool
  >::result
wishrnd(Mat<typename T1::elem_type>& W, const Base<typename T1::elem_type, T1>& S, typename T1::elem_type df)
  {
  arma_extra_debug_sigprint();
  
  const bool status = op_wishrnd::apply_direct(W, S.get_ref(), df, uword(1));
  
  if(status == false)
    {
    W.soft_reset();
    arma_debug_warn_level(3, "wishrnd(): given matrix is not symmetric positive definite");
    }
  
  return status;
  }



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  is_real<typename T1::elem_type>::value,
  bool
  >::result
wishrnd(Mat<typename T1::elem_type>& W, const Base<typename T1::elem_type, T1>& S, typename T1::elem_type df, const Base<typename T1::elem_type, T2>& D)
  {
  arma_extra_debug_sigprint();
  arma_ignore(S);
  
  const bool status = op_wishrnd::apply_direct(W, D.get_ref(), df, uword(2));
  
  if(status == false)
    {
    W.soft_reset();
    arma_debug_warn_level(3, "wishrnd(): problem with given 'D' matrix");
    }
  
  return status;
  }



//



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  is_real<typename T1::elem_type>::value,
  const Op<T1, op_iwishrnd>
  >::result
iwishrnd(const Base<typename T1::elem_type, T1>& T, typename T1::elem_type df)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_iwishrnd>(T.get_ref(), df, uword(1), uword(0));
  }



template<typename T1, typename T2>
arma_warn_unused
inline
typename
enable_if2
  <
  is_real<typename T1::elem_type>::value,
  const Op<T2, op_iwishrnd>
  >::result
iwishrnd(const Base<typename T1::elem_type, T1>& T, typename T1::elem_type df, const Base<typename T1::elem_type, T2>& Dinv)
  {
  arma_extra_debug_sigprint();
  arma_ignore(T);
  
  return Op<T2, op_iwishrnd>(Dinv.get_ref(), df, uword(2), uword(0));
  }



template<typename T1>
inline
typename
enable_if2
  <
  is_real<typename T1::elem_type>::value,
  bool
  >::result
iwishrnd(Mat<typename T1::elem_type>& W, const Base<typename T1::elem_type, T1>& T, typename T1::elem_type df)
  {
  arma_extra_debug_sigprint();
  
  const bool status = op_iwishrnd::apply_direct(W, T.get_ref(), df, uword(1));
  
  if(status == false)
    {
    W.soft_reset();
    arma_debug_warn_level(3, "iwishrnd(): given matrix is not symmetric positive definite and/or df is too low");
    }
  
  return status;
  }



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  is_real<typename T1::elem_type>::value,
  bool
  >::result
iwishrnd(Mat<typename T1::elem_type>& W, const Base<typename T1::elem_type, T1>& T, typename T1::elem_type df, const Base<typename T1::elem_type, T2>& Dinv)
  {
  arma_extra_debug_sigprint();
  arma_ignore(T);
  
  const bool status = op_iwishrnd::apply_direct(W, Dinv.get_ref(), df, uword(2));
  
  if(status == false)
    {
    W.soft_reset();
    arma_debug_warn_level(3, "wishrnd(): problem with given 'Dinv' matrix and/or df is too low");
    }
  
  return status;
  }



//! @}
