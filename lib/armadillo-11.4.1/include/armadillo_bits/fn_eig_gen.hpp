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


//! \addtogroup fn_eig_gen
//! @{


template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_supported_blas_type<typename T1::pod_type>::value, Col< std::complex<typename T1::pod_type> > >::result
eig_gen
  (
  const Base<typename T1::elem_type, T1>& expr,
  const char* option = "nobalance"
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type     T;
  typedef typename std::complex<T> eT;
  
  const char sig = (option != nullptr) ? option[0] : char(0);
  
  arma_debug_check( ((sig != 'n') && (sig != 'b')), "eig_gen(): unknown option" );
  
  if( auxlib::crippled_lapack(expr) && (sig == 'b') )  { arma_debug_warn_level(1,  "eig_gen(): 'balance' option ignored due to linking with crippled lapack"); }
  
  Col<eT> eigvals;
  Mat<eT> eigvecs;
  
  const bool status = (sig == 'b') ? auxlib::eig_gen_balance(eigvals, eigvecs, false, expr.get_ref()) : auxlib::eig_gen(eigvals, eigvecs, false, expr.get_ref());
  
  if(status == false)
    {
    eigvals.soft_reset();
    arma_stop_runtime_error("eig_gen(): decomposition failed");
    }
  
  return eigvals;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::pod_type>::value, bool >::result
eig_gen
  (
         Col< std::complex<typename T1::pod_type> >& eigvals,
  const Base< typename T1::elem_type, T1>&           expr,
  const char* option = "nobalance"
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type     T;
  typedef typename std::complex<T> eT;
  
  const char sig = (option != nullptr) ? option[0] : char(0);
  
  arma_debug_check( ((sig != 'n') && (sig != 'b')), "eig_gen(): unknown option" );
  
  if( auxlib::crippled_lapack(expr) && (sig == 'b') )  { arma_debug_warn_level(1,  "eig_gen(): 'balance' option ignored due to linking with crippled lapack"); }
  
  Mat<eT> eigvecs;
  
  const bool status = (sig == 'b') ? auxlib::eig_gen_balance(eigvals, eigvecs, false, expr.get_ref()) : auxlib::eig_gen(eigvals, eigvecs, false, expr.get_ref());
  
  if(status == false)
    {
    eigvals.soft_reset();
    arma_debug_warn_level(3, "eig_gen(): decomposition failed");
    }
  
  return status;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::pod_type>::value, bool >::result
eig_gen
  (
        Col< std::complex<typename T1::pod_type> >& eigvals,
        Mat< std::complex<typename T1::pod_type> >& eigvecs,
  const Base<typename T1::elem_type, T1>&           expr,
  const char* option = "nobalance"
  )
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (void_ptr(&eigvals) == void_ptr(&eigvecs)), "eig_gen(): parameter 'eigval' is an alias of parameter 'eigvec'" );
  
  const char sig = (option != nullptr) ? option[0] : char(0);
  
  arma_debug_check( ((sig != 'n') && (sig != 'b')), "eig_gen(): unknown option" );
  
  if( auxlib::crippled_lapack(expr) && (sig == 'b') )  { arma_debug_warn_level(1,  "eig_gen(): 'balance' option ignored due to linking with crippled lapack"); }
  
  const bool status = (sig == 'b') ? auxlib::eig_gen_balance(eigvals, eigvecs, true, expr.get_ref()) : auxlib::eig_gen(eigvals, eigvecs, true, expr.get_ref());
  
  if(status == false)
    {
    eigvals.soft_reset();
    eigvecs.soft_reset();
    arma_debug_warn_level(3, "eig_gen(): decomposition failed");
    }
  
  return status;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::pod_type>::value, bool >::result
eig_gen
  (
        Col< std::complex<typename T1::pod_type> >&  eigvals,
        Mat< std::complex<typename T1::pod_type> >& leigvecs,
        Mat< std::complex<typename T1::pod_type> >& reigvecs,
  const Base<typename T1::elem_type, T1>&           expr,
  const char* option = "nobalance"
  )
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (void_ptr(&eigvals)  == void_ptr(&leigvecs)), "eig_gen(): parameter 'eigval' is an alias of parameter 'leigvec'" );
  arma_debug_check( (void_ptr(&eigvals)  == void_ptr(&reigvecs)), "eig_gen(): parameter 'eigval' is an alias of parameter 'reigvec'" );
  arma_debug_check( (void_ptr(&leigvecs) == void_ptr(&reigvecs)), "eig_gen(): parameter 'leigvec' is an alias of parameter 'reigvec'" );
  
  const char sig = (option != nullptr) ? option[0] : char(0);
  
  arma_debug_check( ((sig != 'n') && (sig != 'b')), "eig_gen(): unknown option" );
  
  if( auxlib::crippled_lapack(expr) && (sig == 'b') )  { arma_debug_warn_level(1,  "eig_gen(): 'balance' option ignored due to linking with crippled lapack"); }
  
  const bool status = (sig == 'b') ? auxlib::eig_gen_twosided_balance(eigvals, leigvecs, reigvecs, expr.get_ref()) : auxlib::eig_gen_twosided(eigvals, leigvecs, reigvecs, expr.get_ref());
  
  if(status == false)
    {
     eigvals.soft_reset();
    leigvecs.soft_reset();
    reigvecs.soft_reset();
    arma_debug_warn_level(3, "eig_gen(): decomposition failed");
    }
  
  return status;
  }



//! @}
