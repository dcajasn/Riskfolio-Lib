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


//! \addtogroup fn_quantile
//! @{


template<typename T1, typename T2>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no && is_real<typename T2::elem_type>::value,
  const mtGlue<typename T2::elem_type,T1,T2,glue_quantile_default>
  >::result
quantile(const T1& X, const Base<typename T2::elem_type,T2>& P)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<typename T2::elem_type,T1,T2,glue_quantile_default>(X, P.get_ref());
  }



template<typename T1, typename T2>
arma_warn_unused
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value && is_cx<typename T1::elem_type>::no && is_real<typename T2::elem_type>::value,
  const mtGlue<typename T2::elem_type,T1,T2,glue_quantile>
  >::result
quantile(const T1& X, const Base<typename T2::elem_type,T2>& P, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return mtGlue<typename T2::elem_type,T1,T2,glue_quantile>(X, P.get_ref(), dim);
  }


//! @}
