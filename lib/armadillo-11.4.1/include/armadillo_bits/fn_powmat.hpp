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


//! \addtogroup fn_powmat
//! @{


template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Op<T1,op_powmat> >::result
powmat(const Base<typename T1::elem_type,T1>& X, const int y)
  {
  arma_extra_debug_sigprint();
  
  const uword aux_a = (y < int(0)) ? uword(-y) : uword(y);
  const uword aux_b = (y < int(0)) ? uword(1)  : uword(0);
  
  return Op<T1,op_powmat>(X.get_ref(), aux_a, aux_b);
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
powmat
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& X,
  const int                              y
  )
  {
  arma_extra_debug_sigprint();
  
  const uword y_val = (y < int(0)) ? uword(-y) : uword(y);
  const bool  y_neg = (y < int(0));
  
  const bool status = op_powmat::apply_direct(out, X.get_ref(), y_val, y_neg);
  
  if(status == false)
    {
    out.soft_reset();
    arma_debug_warn_level(3, "powmat(): transformation failed");
    }
  
  return status;
  }



template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const mtOp<std::complex<typename T1::pod_type>,T1,op_powmat_cx> >::result
powmat(const Base<typename T1::elem_type,T1>& X, const double y)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<typename T1::pod_type> out_eT;
  
  return mtOp<out_eT,T1,op_powmat_cx>('j', X.get_ref(), out_eT(y));
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
powmat
  (
         Mat< std::complex<typename T1::pod_type> >& out,
  const Base<typename T1::elem_type,T1>&             X,
  const double                                       y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const bool status = op_powmat_cx::apply_direct(out, X.get_ref(), T(y));
  
  if(status == false)
    {
    out.soft_reset();
    arma_debug_warn_level(3, "powmat(): transformation failed");
    }
  
  return status;
  }


//! @}
