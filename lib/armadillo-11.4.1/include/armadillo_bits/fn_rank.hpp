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


//! \addtogroup fn_rank
//! @{



template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, uword >::result
rank(const Base<typename T1::elem_type,T1>& expr, const typename T1::pod_type tol = 0)
  {
  arma_extra_debug_sigprint();
  
  uword out = uword(0);
  
  const bool status = op_rank::apply(out, expr.get_ref(), tol);
  
  if(status == false)  { arma_stop_runtime_error("rank(): failed"); return uword(0); }
  
  return out;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
rank(uword& out, const Base<typename T1::elem_type,T1>& expr, const typename T1::pod_type tol = 0)
  {
  arma_extra_debug_sigprint();
  
  out = uword(0);
  
  return op_rank::apply(out, expr.get_ref(), tol);
  }



//! @}
