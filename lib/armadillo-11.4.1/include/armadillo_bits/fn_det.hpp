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


//! \addtogroup fn_det
//! @{



template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, typename T1::elem_type >::result
det(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  eT out_val = eT(0);
  
  const bool status = op_det::apply_direct(out_val, X.get_ref());
  
  if(status == false)
    {
    out_val = eT(0);
    arma_stop_runtime_error("det(): failed to find determinant");
    }
  
  return out_val;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
det(typename T1::elem_type& out_val, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const bool status = op_det::apply_direct(out_val, X.get_ref());
  
  if(status == false)
    {
    out_val = eT(0);
    arma_debug_warn_level(3, "det(): failed to find determinant");
    }
  
  return status;
  }



template<typename T>
arma_warn_unused
arma_inline
typename arma_scalar_only<T>::result
det(const T& x)
  {
  return x;
  }



//! @}
