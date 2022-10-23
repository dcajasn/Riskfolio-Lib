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


//! \addtogroup fn_eigs_gen
//! @{


//! eigenvalues of general sparse matrix X
template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_real<typename T1::pod_type>::value, Col< std::complex<typename T1::pod_type> > >::result
eigs_gen
  (
  const SpBase<typename T1::elem_type, T1>& X,
  const uword                               n_eigvals,
  const char*                               form = "lm",
  const eigs_opts                           opts = eigs_opts()
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  Mat< std::complex<T> > eigvec;
  Col< std::complex<T> > eigval;
  
  sp_auxlib::form_type form_val = sp_auxlib::interpret_form_str(form);
  
  const bool status = sp_auxlib::eigs_gen(eigval, eigvec, X, n_eigvals, form_val, opts);
  
  if(status == false)
    {
    eigval.soft_reset();
    arma_stop_runtime_error("eigs_gen(): decomposition failed");
    }
  
  return eigval;
  }



//! this form is deprecated; use eigs_gen(X, n_eigvals, form, opts) instead
template<typename T1>
arma_deprecated
inline
typename enable_if2< is_real<typename T1::pod_type>::value, Col< std::complex<typename T1::pod_type> > >::result
eigs_gen
  (
  const SpBase<typename T1::elem_type, T1>& X,
  const uword                               n_eigvals,
  const char*                               form,
  const typename T1::pod_type               tol
  )
  {
  arma_extra_debug_sigprint();
  
  eigs_opts opts;
  opts.tol = tol;
  
  return eigs_gen(X, n_eigvals, form, opts);
  }



template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_real<typename T1::pod_type>::value, Col< std::complex<typename T1::pod_type> > >::result
eigs_gen
  (
  const SpBase<typename T1::elem_type, T1>& X,
  const uword                               n_eigvals,
  const std::complex<typename T1::pod_type> sigma,
  const eigs_opts                           opts = eigs_opts()
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  Mat< std::complex<T> > eigvec;
  Col< std::complex<T> > eigval;
  
  bool status = false;
  
  // If X is real and sigma is truly complex, treat X as complex.
  // The reason is that we are still not able to apply truly complex shifts to real matrices
  if( (is_real<typename T1::elem_type>::yes) && (std::imag(sigma) != T(0)) )
    {
    status = sp_auxlib::eigs_gen(eigval, eigvec, conv_to< SpMat< std::complex<T> > >::from(X), n_eigvals, sigma, opts);
    }
  else
    {
    status = sp_auxlib::eigs_gen(eigval, eigvec, X, n_eigvals, sigma, opts);
    }
  
  if(status == false)
    {
    eigval.soft_reset();
    arma_stop_runtime_error("eigs_gen(): decomposition failed");
    }
  
  return eigval;
  }



template<typename T1>
arma_warn_unused
inline
typename enable_if2< is_real<typename T1::pod_type>::value, Col< std::complex<typename T1::pod_type> > >::result
eigs_gen
  (
  const SpBase<typename T1::elem_type, T1>& X,
  const uword                               n_eigvals,
  const double                              sigma,
  const eigs_opts                           opts = eigs_opts()
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  Mat< std::complex<T> > eigvec;
  Col< std::complex<T> > eigval;
  
  const bool status = sp_auxlib::eigs_gen(eigval, eigvec, X, n_eigvals, std::complex<T>(T(sigma)), opts);
  
  if(status == false)
    {
    eigval.soft_reset();
    arma_stop_runtime_error("eigs_gen(): decomposition failed");
    }
  
  return eigval;
  }



//! eigenvalues of general sparse matrix X
template<typename T1>
inline
typename enable_if2< is_real<typename T1::pod_type>::value, bool >::result
eigs_gen
  (
           Col< std::complex<typename T1::pod_type> >& eigval,
  const SpBase<typename T1::elem_type, T1>&            X,
  const uword                                          n_eigvals,
  const char*                                          form = "lm",
  const eigs_opts                                      opts = eigs_opts()
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  Mat< std::complex<T> > eigvec;
  
  sp_auxlib::form_type form_val = sp_auxlib::interpret_form_str(form);
  
  const bool status = sp_auxlib::eigs_gen(eigval, eigvec, X, n_eigvals, form_val, opts);
  
  if(status == false)
    {
    eigval.soft_reset();
    arma_debug_warn_level(3, "eigs_gen(): decomposition failed");
    }
  
  return status;
  }



//! this form is deprecated; use eigs_gen(eigval, X, n_eigvals, form, opts) instead
template<typename T1>
arma_deprecated
inline
typename enable_if2< is_real<typename T1::pod_type>::value, bool >::result
eigs_gen
  (
           Col< std::complex<typename T1::pod_type> >& eigval,
  const SpBase<typename T1::elem_type, T1>&            X,
  const uword                                          n_eigvals,
  const char*                                          form,
  const typename T1::pod_type                          tol
  )
  {
  arma_extra_debug_sigprint();
  
  eigs_opts opts;
  opts.tol = tol;
  
  return eigs_gen(eigval, X, n_eigvals, form, opts);
  }



template<typename T1>
inline
typename enable_if2< is_real<typename T1::pod_type>::value, bool >::result
eigs_gen
  (
           Col< std::complex<typename T1::pod_type> >& eigval,
  const SpBase<typename T1::elem_type, T1>&            X,
  const uword                                          n_eigvals,
  const std::complex<typename T1::pod_type>            sigma,
  const eigs_opts                                      opts = eigs_opts()
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  Mat< std::complex<T> > eigvec;
  
  bool status = false;
  
  // If X is real and sigma is truly complex, treat X as complex.
  // The reason is that we are still not able to apply truly complex shifts to real matrices
  if( (is_real<typename T1::elem_type>::yes) && (std::imag(sigma) != T(0)) )
    {
    status = sp_auxlib::eigs_gen(eigval, eigvec, conv_to< SpMat< std::complex<T> > >::from(X), n_eigvals, sigma, opts);
    }
  else
    {
    status = sp_auxlib::eigs_gen(eigval, eigvec, X, n_eigvals, sigma, opts);
    }
  
  if(status == false)
    {
    eigval.soft_reset();
    arma_debug_warn_level(3, "eigs_gen(): decomposition failed");
    }
  
  return status;
  }



template<typename T1>
inline
typename enable_if2< is_real<typename T1::pod_type>::value, bool >::result
eigs_gen
  (
           Col< std::complex<typename T1::pod_type> >& eigval,
  const SpBase<typename T1::elem_type, T1>&            X,
  const uword                                          n_eigvals,
  const double                                         sigma,
  const eigs_opts                                      opts = eigs_opts()
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  Mat< std::complex<T> > eigvec;
  
  const bool status = sp_auxlib::eigs_gen(eigval, eigvec, X, n_eigvals, std::complex<T>(T(sigma)), opts);
  
  if(status == false)
    {
    eigval.soft_reset();
    arma_debug_warn_level(3, "eigs_gen(): decomposition failed");
    }
  
  return status;
  }



//! eigenvalues and eigenvectors of general sparse matrix X
template<typename T1>
inline
typename enable_if2< is_real<typename T1::pod_type>::value, bool >::result
eigs_gen
  (
         Col< std::complex<typename T1::pod_type> >& eigval,
         Mat< std::complex<typename T1::pod_type> >& eigvec,
  const SpBase<typename T1::elem_type, T1>&          X,
  const uword                                        n_eigvals,
  const char*                                        form = "lm",
  const eigs_opts                                    opts = eigs_opts()
  )
  {
  arma_extra_debug_sigprint();
  
  // typedef typename T1::pod_type T;
  
  arma_debug_check( void_ptr(&eigval) == void_ptr(&eigvec), "eigs_gen(): parameter 'eigval' is an alias of parameter 'eigvec'" );
  
  sp_auxlib::form_type form_val = sp_auxlib::interpret_form_str(form);
  
  const bool status = sp_auxlib::eigs_gen(eigval, eigvec, X, n_eigvals, form_val, opts);
  
  if(status == false)
    {
    eigval.soft_reset();
    eigvec.soft_reset();
    arma_debug_warn_level(3, "eigs_gen(): decomposition failed");
    }
  
  return status;
  }



//! this form is deprecated; use eigs_gen(eigval, eigvec, X, n_eigvals, form, opts) instead
template<typename T1>
arma_deprecated
inline
typename enable_if2< is_real<typename T1::pod_type>::value, bool >::result
eigs_gen
  (
         Col< std::complex<typename T1::pod_type> >& eigval,
         Mat< std::complex<typename T1::pod_type> >& eigvec,
  const SpBase<typename T1::elem_type, T1>&          X,
  const uword                                        n_eigvals,
  const char*                                        form,
  const typename T1::pod_type                        tol
  )
  {
  arma_extra_debug_sigprint();
  
  eigs_opts opts;
  opts.tol = tol;
  
  return eigs_gen(eigval, eigvec, X, n_eigvals, form, opts);
  }



template<typename T1>
inline
typename enable_if2< is_real<typename T1::pod_type>::value, bool >::result
eigs_gen
  (
         Col< std::complex<typename T1::pod_type> >& eigval,
         Mat< std::complex<typename T1::pod_type> >& eigvec,
  const SpBase<typename T1::elem_type, T1>&          X,
  const uword                                        n_eigvals,
  const std::complex<typename T1::pod_type>          sigma,
  const eigs_opts                                    opts = eigs_opts()
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  arma_debug_check( void_ptr(&eigval) == void_ptr(&eigvec), "eigs_gen(): parameter 'eigval' is an alias of parameter 'eigvec'" );
  
  bool status = false;
  
  // If X is real and sigma is truly complex, treat X as complex.
  // The reason is that we are still not able to apply truly complex shifts to real matrices
  if( (is_real<typename T1::elem_type>::yes) && (std::imag(sigma) != T(0)) )
    {
    status = sp_auxlib::eigs_gen(eigval, eigvec, conv_to< SpMat< std::complex<T> > >::from(X), n_eigvals, sigma, opts);
    }
  else
    {
    status = sp_auxlib::eigs_gen(eigval, eigvec, X, n_eigvals, sigma, opts);
    }
  
  if(status == false)
    {
    eigval.soft_reset();
    eigvec.soft_reset();
    arma_debug_warn_level(3, "eigs_gen(): decomposition failed");
    }
  
  return status;
  }



template<typename T1>
inline
typename enable_if2< is_real<typename T1::pod_type>::value, bool >::result
eigs_gen
  (
         Col< std::complex<typename T1::pod_type> >& eigval,
         Mat< std::complex<typename T1::pod_type> >& eigvec,
  const SpBase<typename T1::elem_type, T1>&          X,
  const uword                                        n_eigvals,
  const double                                       sigma,
  const eigs_opts                                    opts = eigs_opts()
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  arma_debug_check( void_ptr(&eigval) == void_ptr(&eigvec), "eigs_gen(): parameter 'eigval' is an alias of parameter 'eigvec'" );
  
  const bool status = sp_auxlib::eigs_gen(eigval, eigvec, X, n_eigvals, std::complex<T>(T(sigma)), opts);
  
  if(status == false)
    {
    eigval.soft_reset();
    eigvec.soft_reset();
    arma_debug_warn_level(3, "eigs_gen(): decomposition failed");
    }
  
  return status;
  }



//! @}
