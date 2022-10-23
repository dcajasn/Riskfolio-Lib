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


//! \addtogroup fn_inv
//! @{



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Op<T1, op_inv_gen_default> >::result
inv
  (
  const Base<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_inv_gen_default>(X.get_ref());
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
inv
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  const bool status = op_inv_gen_default::apply_direct(out, X.get_ref(), "inv()");
  
  if(status == false)
    {
    out.soft_reset();
    arma_debug_warn_level(3, "inv(): matrix is singular");
    }
  
  return status;
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Op<T1, op_inv_gen_full> >::result
inv
  (
  const Base<typename T1::elem_type,T1>& X,
  const inv_opts::opts&                  opts
  )
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_inv_gen_full>(X.get_ref(), opts.flags, uword(0));
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
inv
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& X,
  const inv_opts::opts&                  opts
  )
  {
  arma_extra_debug_sigprint();
  
  const bool status = op_inv_gen_full::apply_direct(out, X.get_ref(), "inv()", opts.flags);
  
  if(status == false)
    {
    out.soft_reset();
    arma_debug_warn_level(3, "inv(): matrix is singular");
    }
  
  return status;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
inv
  (
         Mat<typename T1::elem_type>&    out_inv,
             typename T1::pod_type&      out_rcond,
  const Base<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  op_inv_gen_state<T> inv_state;
  
  const bool status = op_inv_gen_rcond::apply_direct(out_inv, inv_state, X.get_ref());
  
  out_rcond = inv_state.rcond;
  
  if(status == false)
    {
    out_rcond = T(0);
    out_inv.soft_reset();
    arma_debug_warn_level(3, "inv(): matrix is singular");
    }
  
  return status;
  }



//! @}
