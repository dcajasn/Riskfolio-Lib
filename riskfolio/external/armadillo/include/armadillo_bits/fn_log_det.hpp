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


//! \addtogroup fn_log_det
//! @{



//! log determinant of mat
template<typename T1>
inline
bool
log_det
  (
        typename T1::elem_type&          out_val,
        typename T1::pod_type&           out_sign,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const bool status = op_log_det::apply_direct(out_val, out_sign, X.get_ref());
  
  if(status == false)
    {
    out_val  = eT(Datum<T>::nan);
    out_sign = T(0);
    
    arma_debug_warn_level(3, "log_det(): failed to find determinant");
    }
  
  return status;
  }



template<typename T1>
inline
arma_warn_unused
std::complex<typename T1::pod_type>
log_det
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  eT out_val  = eT(0);
   T out_sign =  T(0);
  
  const bool status = op_log_det::apply_direct(out_val, out_sign, X.get_ref());
  
  if(status == false)
    {
    out_val  = eT(Datum<T>::nan);
    out_sign = T(0);
    
    arma_stop_runtime_error("log_det(): failed to find determinant");
    }
  
  return (out_sign >= T(1)) ? std::complex<T>(out_val) : (out_val + std::complex<T>(T(0),Datum<T>::pi));
  }



//



template<typename T1>
inline
bool
log_det_sympd
  (
        typename T1::pod_type&           out_val,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type T;
  
  out_val = T(0);
  
  const bool status = op_log_det_sympd::apply_direct(out_val, X.get_ref());
  
  if(status == false)
    {
    out_val = Datum<T>::nan;
    
    arma_debug_warn_level(3, "log_det_sympd(): given matrix is not symmetric positive definite");
    }
  
  return status;
  }



template<typename T1>
inline
arma_warn_unused
typename T1::pod_type
log_det_sympd
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = nullptr
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::pod_type T;
  
  T out_val = T(0);
  
  const bool status = op_log_det_sympd::apply_direct(out_val, X.get_ref());
  
  if(status == false)
    {
    out_val = Datum<T>::nan;
    
    arma_stop_runtime_error("log_det_sympd(): given matrix is not symmetric positive definite");
    }
  
  return out_val;
  }



//! @}
