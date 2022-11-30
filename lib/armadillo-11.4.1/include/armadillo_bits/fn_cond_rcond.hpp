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


//! \addtogroup fn_cond
//! @{



template<typename T1>
arma_warn_unused
inline
typename enable_if2<is_supported_blas_type<typename T1::elem_type>::value, typename T1::pod_type>::result
cond(const Base<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return op_cond::apply(X.get_ref());
  }



template<typename T1>
arma_warn_unused
inline
typename enable_if2<is_supported_blas_type<typename T1::elem_type>::value, typename T1::pod_type>::result
rcond(const Base<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return op_rcond::apply(X.get_ref());
  }



// template<typename T1>
// arma_warn_unused
// inline
// typename enable_if2<is_supported_blas_type<typename T1::elem_type>::value, typename T1::pod_type>::result
// rcond(const SpBase<typename T1::elem_type, T1>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return sp_auxlib::rcond(X.get_ref());
//   }



//! @}
