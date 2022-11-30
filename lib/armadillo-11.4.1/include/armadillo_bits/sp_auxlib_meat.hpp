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


//! \addtogroup sp_auxlib
//! @{


inline
sp_auxlib::form_type
sp_auxlib::interpret_form_str(const char* form_str)
  {
  arma_extra_debug_sigprint();
  
  // the order of the 3 if statements below is important
  if( form_str    == nullptr )  { return form_none; }
  if( form_str[0] == char(0) )  { return form_none; }
  if( form_str[1] == char(0) )  { return form_none; }
  
  const char c1 = form_str[0];
  const char c2 = form_str[1];
  
  if(c1 == 'l')
    {
    if(c2 == 'm')  { return form_lm; }
    if(c2 == 'r')  { return form_lr; }
    if(c2 == 'i')  { return form_li; }
    if(c2 == 'a')  { return form_la; }
    }
  else
  if(c1 == 's')
    {
    if(c2 == 'm')  { return form_sm; }
    if(c2 == 'r')  { return form_sr; }
    if(c2 == 'i')  { return form_si; }
    if(c2 == 'a')  { return form_sa; }
    }
  
  return form_none;
  }



//! immediate eigendecomposition of symmetric real sparse object
template<typename eT, typename T1>
inline
bool
sp_auxlib::eigs_sym(Col<eT>& eigval, Mat<eT>& eigvec, const SpBase<eT, T1>& X, const uword n_eigvals, const form_type form_val, const eigs_opts& opts)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> U(X.get_ref());
  
  arma_debug_check( (U.M.is_square() == false), "eigs_sym(): given matrix must be square sized" );
  
  if((arma_config::debug) && (sp_auxlib::rudimentary_sym_check(U.M) == false))
    {
    if(is_cx<eT>::no )  { arma_debug_warn_level(1, "eigs_sym(): given matrix is not symmetric"); }
    if(is_cx<eT>::yes)  { arma_debug_warn_level(1, "eigs_sym(): given matrix is not hermitian"); }
    }
  
  if(arma_config::check_nonfinite && U.M.has_nonfinite())
    {
    arma_debug_warn_level(3, "eigs_sym(): detected non-finite elements");
    return false;
    }
  
  // TODO: investigate optional redirection of "sm" to ARPACK as it's capable of shift-invert;
  // TODO: in shift-invert mode, "sm" maps to "lm" of the shift-inverted matrix (with sigma = 0)
  
  #if   defined(ARMA_USE_NEWARP)
    {
    return sp_auxlib::eigs_sym_newarp(eigval, eigvec, U.M, n_eigvals, form_val, opts);
    }
  #elif defined(ARMA_USE_ARPACK)
    {
    constexpr eT sigma = eT(0);
    
    return sp_auxlib::eigs_sym_arpack<eT,false>(eigval, eigvec, U.M, n_eigvals, form_val, sigma, opts);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(n_eigvals);
    arma_ignore(form_val);
    arma_ignore(opts);
    
    arma_stop_logic_error("eigs_sym(): use of NEWARP or ARPACK must be enabled");
    return false;
    }
  #endif
  }



//! immediate eigendecomposition of symmetric real sparse object
template<typename eT, typename T1>
inline
bool
sp_auxlib::eigs_sym(Col<eT>& eigval, Mat<eT>& eigvec, const SpBase<eT, T1>& X, const uword n_eigvals, const eT sigma, const eigs_opts& opts)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> U(X.get_ref());
  
  arma_debug_check( (U.M.is_square() == false), "eigs_sym(): given matrix must be square sized" );
  
  if((arma_config::debug) && (sp_auxlib::rudimentary_sym_check(U.M) == false))
    {
    if(is_cx<eT>::no )  { arma_debug_warn_level(1, "eigs_sym(): given matrix is not symmetric"); }
    if(is_cx<eT>::yes)  { arma_debug_warn_level(1, "eigs_sym(): given matrix is not hermitian"); }
    }
  
  if(arma_config::check_nonfinite && U.M.has_nonfinite())
    {
    arma_debug_warn_level(3, "eigs_sym(): detected non-finite elements");
    return false;
    }
  
  #if   (defined(ARMA_USE_NEWARP) && defined(ARMA_USE_SUPERLU))
    {
    return sp_auxlib::eigs_sym_newarp(eigval, eigvec, U.M, n_eigvals, sigma, opts);
    }
  #elif (defined(ARMA_USE_ARPACK) && defined(ARMA_USE_SUPERLU))
    {
    constexpr form_type form_val = form_sigma;
    
    return sp_auxlib::eigs_sym_arpack<eT,true>(eigval, eigvec, U.M, n_eigvals, form_val, sigma, opts);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(n_eigvals);
    arma_ignore(sigma);
    arma_ignore(opts);
    
    arma_stop_logic_error("eigs_sym(): use of NEWARP or ARPACK as well as SuperLU must be enabled to use 'sigma'");
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
sp_auxlib::eigs_sym_newarp(Col<eT>& eigval, Mat<eT>& eigvec, const SpMat<eT>& X, const uword n_eigvals, const form_type form_val, const eigs_opts& opts)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_NEWARP)
    {
    arma_debug_check( (form_val != form_lm) && (form_val != form_sm) && (form_val != form_la) && (form_val != form_sa), "eigs_sym(): unknown form specified" );
    
    if(X.is_square() == false)  { return false; }
    
    const newarp::SparseGenMatProd<eT> op(X);
    
    arma_debug_check( (n_eigvals >= op.n_rows), "eigs_sym(): n_eigvals must be less than the number of rows in the matrix" );
    
    // If the matrix is empty, the case is trivial.
    if( (op.n_cols == 0) || (n_eigvals == 0) ) // We already know n_cols == n_rows.
      {
      eigval.reset();
      eigvec.reset();
      return true;
      }
    
    uword n   = op.n_rows;
    
    // Use max(2*k+1, 20) as default subspace dimension for the sym case; MATLAB uses max(2*k, 20), but we need to be backward-compatible.
    uword ncv_default = uword( ((2*n_eigvals+1)>(20)) ? (2*n_eigvals+1) : (20) );
    
    // Use opts.subdim only if it's within the limits, otherwise cap it.
    uword ncv = ncv_default;
    
    if(opts.subdim != 0)
      {
      if(opts.subdim < (n_eigvals + 1))
        {
        arma_debug_warn_level(1, "eigs_sym(): opts.subdim must be greater than k; using k+1 instead of ", opts.subdim);
        ncv = uword(n_eigvals + 1);
        }
      else
      if(opts.subdim > n)
        {
        arma_debug_warn_level(1, "eigs_sym(): opts.subdim cannot be greater than n_rows; using n_rows instead of ", opts.subdim);
        ncv = n;
        }
      else
        {
        ncv = uword(opts.subdim);
        }
      }
    
    // Re-check that we are within the limits
    if(ncv < (n_eigvals + 1)) { ncv = (n_eigvals + 1); }
    if(ncv > n              ) { ncv = n;               }
    
    eT tol = (std::max)(eT(opts.tol), std::numeric_limits<eT>::epsilon());
    
    uword maxiter = uword(opts.maxiter);
    
    // eigval.set_size(n_eigvals);
    // eigvec.set_size(n, n_eigvals);
    
    bool status = true;
    
    uword nconv = 0;
    
    try
      {
      if(form_val == form_lm)
        {
        newarp::SymEigsSolver< eT, newarp::EigsSelect::LARGEST_MAGN, newarp::SparseGenMatProd<eT> > eigs(op, n_eigvals, ncv);
        eigs.init();
        nconv  = eigs.compute(maxiter, tol);
        eigval = eigs.eigenvalues();
        eigvec = eigs.eigenvectors();
        }
      else
      if(form_val == form_sm)
        {
        newarp::SymEigsSolver< eT, newarp::EigsSelect::SMALLEST_MAGN, newarp::SparseGenMatProd<eT> > eigs(op, n_eigvals, ncv);
        eigs.init();
        nconv  = eigs.compute(maxiter, tol);
        eigval = eigs.eigenvalues();
        eigvec = eigs.eigenvectors();
        }
      else
      if(form_val == form_la)
        {
        newarp::SymEigsSolver< eT, newarp::EigsSelect::LARGEST_ALGE, newarp::SparseGenMatProd<eT> > eigs(op, n_eigvals, ncv);
        eigs.init();
        nconv  = eigs.compute(maxiter, tol);
        eigval = eigs.eigenvalues();
        eigvec = eigs.eigenvectors();
        }
      else
      if(form_val == form_sa)
        {
        newarp::SymEigsSolver< eT, newarp::EigsSelect::SMALLEST_ALGE, newarp::SparseGenMatProd<eT> > eigs(op, n_eigvals, ncv);
        eigs.init();
        nconv  = eigs.compute(maxiter, tol);
        eigval = eigs.eigenvalues();
        eigvec = eigs.eigenvectors();
        }
      }
    catch(const std::runtime_error&)
      {
      status = false;
      }
    
    if(status == true)
      {
      if(nconv == 0)  { status = false; }
      }
    
    return status;
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_ignore(n_eigvals);
    arma_ignore(form_val);
    arma_ignore(opts);
    
    return false;
    }
  #endif
  }



template<typename eT>
inline
bool
sp_auxlib::eigs_sym_newarp(Col<eT>& eigval, Mat<eT>& eigvec, const SpMat<eT>& X, const uword n_eigvals, const eT sigma, const eigs_opts& opts)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_NEWARP)
    {
    if(X.is_square() == false)  { return false; }
    
    const newarp::SparseGenRealShiftSolve<eT> op(X, sigma);
    
    if(op.valid == false)  { return false; }
    
    arma_debug_check( (n_eigvals >= op.n_rows), "eigs_sym(): n_eigvals must be less than the number of rows in the matrix" );
    
    // If the matrix is empty, the case is trivial.
    if( (op.n_cols == 0) || (n_eigvals == 0) ) // We already know n_cols == n_rows.
      {
      eigval.reset();
      eigvec.reset();
      return true;
      }
    
    uword n   = op.n_rows;
    
    // Use max(2*k+1, 20) as default subspace dimension for the sym case; MATLAB uses max(2*k, 20), but we need to be backward-compatible.
    uword ncv_default = uword( ((2*n_eigvals+1)>(20)) ? (2*n_eigvals+1) : (20) );
    
    // Use opts.subdim only if it's within the limits, otherwise cap it.
    uword ncv = ncv_default;
    
    if(opts.subdim != 0)
      {
      if(opts.subdim < (n_eigvals + 1))
        {
        arma_debug_warn_level(1, "eigs_sym(): opts.subdim must be greater than k; using k+1 instead of ", opts.subdim);
        ncv = uword(n_eigvals + 1);
        }
      else
      if(opts.subdim > n)
        {
        arma_debug_warn_level(1, "eigs_sym(): opts.subdim cannot be greater than n_rows; using n_rows instead of ", opts.subdim);
        ncv = n;
        }
      else
        {
        ncv = uword(opts.subdim);
        }
      }
    
    // Re-check that we are within the limits
    if(ncv < (n_eigvals + 1)) { ncv = (n_eigvals + 1); }
    if(ncv > n              ) { ncv = n;               }
    
    eT tol = (std::max)(eT(opts.tol), std::numeric_limits<eT>::epsilon());
    
    uword maxiter = uword(opts.maxiter);
    
    // eigval.set_size(n_eigvals);
    // eigvec.set_size(n, n_eigvals);
    
    bool status = true;
    
    uword nconv = 0;
    
    try
      {
      newarp::SymEigsShiftSolver< eT, newarp::EigsSelect::LARGEST_MAGN, newarp::SparseGenRealShiftSolve<eT> > eigs(op, n_eigvals, ncv, sigma);
      eigs.init();
      nconv  = eigs.compute(maxiter, tol);
      eigval = eigs.eigenvalues();
      eigvec = eigs.eigenvectors();
      }
    catch(const std::runtime_error&)
      {
      status = false;
      }
    
    if(status == true)
      {
      if(nconv == 0)  { status = false; }
      }
    
    return status;
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_ignore(n_eigvals);
    arma_ignore(sigma);
    arma_ignore(opts);
    
    return false;
    }
  #endif
  }



template<typename eT, bool use_sigma>
inline
bool
sp_auxlib::eigs_sym_arpack(Col<eT>& eigval, Mat<eT>& eigvec, const SpMat<eT>& X, const uword n_eigvals, const form_type form_val, const eT sigma, const eigs_opts& opts)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_ARPACK)
    {
    arma_debug_check( (form_val != form_lm) && (form_val != form_sm) && (form_val != form_la) && (form_val != form_sa) && (form_val != form_sigma), "eigs_sym(): unknown form specified" );
    
    if(X.is_square() == false)  { return false; }
    
    char  which_sm[3] = "SM";
    char  which_lm[3] = "LM";
    char  which_sa[3] = "SA";
    char  which_la[3] = "LA";
    char* which;
    
    switch(form_val)
      {
      case form_sm:  which = which_sm;  break;
      case form_lm:  which = which_lm;  break;
      case form_sa:  which = which_sa;  break;
      case form_la:  which = which_la;  break;
      
      default:       which = which_lm;  break;
      }
    
    // Make sure we aren't asking for every eigenvalue.
    // The _saupd() functions allow asking for one more eigenvalue than the _naupd() functions.
    arma_debug_check( (n_eigvals >= X.n_rows), "eigs_sym(): n_eigvals must be less than the number of rows in the matrix" );
    
    // If the matrix is empty, the case is trivial.
    if( (X.n_cols == 0) || (n_eigvals == 0) ) // We already know n_cols == n_rows.
      {
      eigval.reset();
      eigvec.reset();
      return true;
      }
    
    // Set up variables that get used for neupd().
    blas_int n, ncv, ncv_default, ldv, lworkl, info, maxiter;
    
    eT tol  = eT(opts.tol);
    maxiter = blas_int(opts.maxiter);
    
    podarray<eT> resid, v, workd, workl;
    podarray<blas_int> iparam, ipntr;
    podarray<eT> rwork; // Not used in this case.
    
    n = blas_int(X.n_rows); // The size of the matrix.
    
    // Use max(2*k+1, 20) as default subspace dimension for the sym case; MATLAB uses max(2*k, 20), but we need to be backward-compatible.
    ncv_default = blas_int( ((2*n_eigvals+1)>(20)) ? (2*n_eigvals+1) : (20) );
    
    // Use opts.subdim only if it's within the limits
    ncv = ncv_default;
    
    if(opts.subdim != 0)
      {
      if(opts.subdim < (n_eigvals + 1))
        {
        arma_debug_warn_level(1, "eigs_sym(): opts.subdim must be greater than k; using k+1 instead of ", opts.subdim);
        ncv = blas_int(n_eigvals + 1);
        }
      else
      if(blas_int(opts.subdim) > n)
        {
        arma_debug_warn_level(1, "eigs_sym(): opts.subdim cannot be greater than n_rows; using n_rows instead of ", opts.subdim);
        ncv = n;
        }
      else
        {
        ncv = blas_int(opts.subdim);
        }
      }
    
    if(use_sigma)
    //if(form_val == form_sigma)
      {
      run_aupd_shiftinvert(n_eigvals, sigma, X, true /* sym, not gen */, n, tol, maxiter, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info);
      }
    else
      {
      run_aupd_plain(n_eigvals, which, X, true /* sym, not gen */, n, tol, maxiter, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info);
      }
    
    if(info != 0)  { return false; }
    
    // The process has converged, and now we need to recover the actual eigenvectors using seupd()
    blas_int rvec = 1; // .TRUE
    blas_int nev  = blas_int(n_eigvals);
    
    char howmny = 'A';
    char bmat   = 'I'; // We are considering the standard eigenvalue problem.
    
    podarray<blas_int> select(ncv); // Logical array of dimension NCV.
    blas_int ldz = n;
    
    select.zeros();
    
    // seupd() will output directly into the eigval and eigvec objects.
    eigval.zeros(n_eigvals);
    eigvec.zeros(n, n_eigvals);
    
    arpack::seupd(&rvec, &howmny, select.memptr(), eigval.memptr(), eigvec.memptr(), &ldz, (eT*) &sigma, &bmat, &n, which, &nev, &tol, resid.memptr(), &ncv, v.memptr(), &ldv, iparam.memptr(), ipntr.memptr(), workd.memptr(), workl.memptr(), &lworkl, &info);
    
    // Check for errors.
    if(info != 0)  { arma_debug_warn_level(1, "eigs_sym(): ARPACK error ", info, " in seupd()"); return false; }
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_ignore(n_eigvals);
    arma_ignore(form_val);
    arma_ignore(sigma);
    arma_ignore(opts);
    
    return false;
    }
  #endif
  }



//! immediate eigendecomposition of non-symmetric real sparse object
template<typename T, typename T1>
inline
bool
sp_auxlib::eigs_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpBase<T, T1>& X, const uword n_eigvals, const form_type form_val, const eigs_opts& opts)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> U(X.get_ref());
  
  arma_debug_check( (U.M.is_square() == false), "eigs_gen(): given matrix must be square sized" );
  
  if(arma_config::check_nonfinite && U.M.has_nonfinite())
    {
    arma_debug_warn_level(3, "eigs_gen(): detected non-finite elements");
    return false;
    }
  
  // TODO: investigate optional redirection of "sm" to ARPACK as it's capable of shift-invert;
  // TODO: in shift-invert mode, "sm" maps to "lm" of the shift-inverted matrix (with sigma = 0)
  
  #if defined(ARMA_USE_NEWARP)
    {
    return sp_auxlib::eigs_gen_newarp(eigval, eigvec, U.M, n_eigvals, form_val, opts);
    }
  #elif defined(ARMA_USE_ARPACK)
    {
    constexpr std::complex<T> sigma = T(0);
    
    return sp_auxlib::eigs_gen_arpack<T,false>(eigval, eigvec, U.M, n_eigvals, form_val, sigma, opts);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(n_eigvals);
    arma_ignore(form_val);
    arma_ignore(opts);
    
    arma_stop_logic_error("eigs_gen(): use of NEWARP or ARPACK must be enabled");
    return false;
    }
  #endif
  }



//! immediate eigendecomposition of non-symmetric real sparse object
template<typename T, typename T1>
inline
bool
sp_auxlib::eigs_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpBase<T, T1>& X, const uword n_eigvals, const std::complex<T> sigma, const eigs_opts& opts)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> U(X.get_ref());
  
  arma_debug_check( (U.M.is_square() == false), "eigs_gen(): given matrix must be square sized" );
  
  if(arma_config::check_nonfinite && U.M.has_nonfinite())
    {
    arma_debug_warn_level(3, "eigs_gen(): detected non-finite elements");
    return false;
    }
  
  #if (defined(ARMA_USE_ARPACK) && defined(ARMA_USE_SUPERLU))
    {
    constexpr form_type form_val = form_sigma;
    
    return sp_auxlib::eigs_gen_arpack<T,true>(eigval, eigvec, U.M, n_eigvals, form_val, sigma, opts);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(n_eigvals);
    arma_ignore(sigma);
    arma_ignore(opts);
    
    arma_stop_logic_error("eigs_gen(): use of ARPACK and SuperLU must be enabled to use 'sigma'");
    return false;
    }
  #endif
  }



template<typename T>
inline
bool
sp_auxlib::eigs_gen_newarp(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpMat<T>& X, const uword n_eigvals, const form_type form_val, const eigs_opts& opts)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_NEWARP)
    {
    arma_debug_check( (form_val != form_lm) && (form_val != form_sm) && (form_val != form_lr) && (form_val != form_sr) && (form_val != form_li) && (form_val != form_si), "eigs_gen(): unknown form specified" );
    
    if(X.is_square() == false)  { return false; }
    
    const newarp::SparseGenMatProd<T> op(X);
    
    arma_debug_check( (n_eigvals + 1 >= op.n_rows), "eigs_gen(): n_eigvals + 1 must be less than the number of rows in the matrix" );
    
    // If the matrix is empty, the case is trivial.
    if( (op.n_cols == 0) || (n_eigvals == 0) ) // We already know n_cols == n_rows.
      {
      eigval.reset();
      eigvec.reset();
      return true;
      }
    
    uword n   = op.n_rows;
    
    // Use max(2*k+1, 20) as default subspace dimension for the gen case; same as MATLAB.
    uword ncv_default = uword( ((2*n_eigvals+1)>(20)) ? (2*n_eigvals+1) : (20) );
    
    // Use opts.subdim only if it's within the limits
    uword ncv = ncv_default;
    
    if(opts.subdim != 0)
      {
      if(opts.subdim < (n_eigvals + 3))
        {
        arma_debug_warn_level(1, "eigs_gen(): opts.subdim must be greater than k+2; using k+3 instead of ", opts.subdim);
        ncv = uword(n_eigvals + 3);
        }
      else
      if(opts.subdim > n)
        {
        arma_debug_warn_level(1, "eigs_gen(): opts.subdim cannot be greater than n_rows; using n_rows instead of ", opts.subdim);
        ncv = n;
        }
      else
        {
        ncv = uword(opts.subdim);
        }
      }
    
    // Re-check that we are within the limits
    if(ncv < (n_eigvals + 3)) { ncv = (n_eigvals + 3); }
    if(ncv > n              ) { ncv = n;               }
    
    T tol = (std::max)(T(opts.tol), std::numeric_limits<T>::epsilon());
    
    uword maxiter = uword(opts.maxiter);
    
    // eigval.set_size(n_eigvals);
    // eigvec.set_size(n, n_eigvals);
    
    bool status = true;
    
    uword nconv = 0;
    
    try
      {
      if(form_val == form_lm)
        {
        newarp::GenEigsSolver< T, newarp::EigsSelect::LARGEST_MAGN, newarp::SparseGenMatProd<T> > eigs(op, n_eigvals, ncv);
        eigs.init();
        nconv  = eigs.compute(maxiter, tol);
        eigval = eigs.eigenvalues();
        eigvec = eigs.eigenvectors();
        }
      else
      if(form_val == form_sm)
        {
        newarp::GenEigsSolver< T, newarp::EigsSelect::SMALLEST_MAGN, newarp::SparseGenMatProd<T> > eigs(op, n_eigvals, ncv);
        eigs.init();
        nconv  = eigs.compute(maxiter, tol);
        eigval = eigs.eigenvalues();
        eigvec = eigs.eigenvectors();
        }
      else
      if(form_val == form_lr)
        {
        newarp::GenEigsSolver< T, newarp::EigsSelect::LARGEST_REAL, newarp::SparseGenMatProd<T> > eigs(op, n_eigvals, ncv);
        eigs.init();
        nconv  = eigs.compute(maxiter, tol);
        eigval = eigs.eigenvalues();
        eigvec = eigs.eigenvectors();
        }
      else
      if(form_val == form_sr)
        {
        newarp::GenEigsSolver< T, newarp::EigsSelect::SMALLEST_REAL, newarp::SparseGenMatProd<T> > eigs(op, n_eigvals, ncv);
        eigs.init();
        nconv  = eigs.compute(maxiter, tol);
        eigval = eigs.eigenvalues();
        eigvec = eigs.eigenvectors();
        }
      else
      if(form_val == form_li)
        {
        newarp::GenEigsSolver< T, newarp::EigsSelect::LARGEST_IMAG, newarp::SparseGenMatProd<T> > eigs(op, n_eigvals, ncv);
        eigs.init();
        nconv  = eigs.compute(maxiter, tol);
        eigval = eigs.eigenvalues();
        eigvec = eigs.eigenvectors();
        }
      else
      if(form_val == form_si)
        {
        newarp::GenEigsSolver< T, newarp::EigsSelect::SMALLEST_IMAG, newarp::SparseGenMatProd<T> > eigs(op, n_eigvals, ncv);
        eigs.init();
        nconv  = eigs.compute(maxiter, tol);
        eigval = eigs.eigenvalues();
        eigvec = eigs.eigenvectors();
        }
      }
    catch(const std::runtime_error&)
      {
      status = false;
      }
    
    if(status == true)
      {
      if(nconv == 0)  { status = false; }
      }
    
    return status;
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_ignore(n_eigvals);
    arma_ignore(form_val);
    arma_ignore(opts);
    
    return false;
    }
  #endif
  }




template<typename T, bool use_sigma>
inline
bool
sp_auxlib::eigs_gen_arpack(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpMat<T>& X, const uword n_eigvals, const form_type form_val, const std::complex<T> sigma, const eigs_opts& opts)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_ARPACK)
    {
    arma_debug_check( (form_val != form_lm) && (form_val != form_sm) && (form_val != form_lr) && (form_val != form_sr) && (form_val != form_li) && (form_val != form_si) && (form_val != form_sigma), "eigs_gen(): unknown form specified" );
    
    if(X.is_square() == false)  { return false; }
    
    char which_lm[3] = "LM";
    char which_sm[3] = "SM";
    char which_lr[3] = "LR";
    char which_sr[3] = "SR";
    char which_li[3] = "LI";
    char which_si[3] = "SI";
    
    char* which;
    
    switch(form_val)
      {
      case form_lm:  which = which_lm;  break;
      case form_sm:  which = which_sm;  break;
      case form_lr:  which = which_lr;  break;
      case form_sr:  which = which_sr;  break;
      case form_li:  which = which_li;  break;
      case form_si:  which = which_si;  break;
      
      default:       which = which_lm;
      }
    
    // Make sure we aren't asking for every eigenvalue.
    arma_debug_check( (n_eigvals + 1 >= X.n_rows), "eigs_gen(): n_eigvals + 1 must be less than the number of rows in the matrix" );
    
    // If the matrix is empty, the case is trivial.
    if( (X.n_cols == 0) || (n_eigvals == 0) ) // We already know n_cols == n_rows.
      {
      eigval.reset();
      eigvec.reset();
      return true;
      }
    
    // Set up variables that get used for neupd().
    blas_int n, ncv, ncv_default, ldv, lworkl, info, maxiter;
    
    T tol   = T(opts.tol);
    maxiter = blas_int(opts.maxiter);
    
    podarray<T> resid, v, workd, workl;
    podarray<blas_int> iparam, ipntr;
    podarray<T> rwork; // Not used in the real case.
    
    n = blas_int(X.n_rows); // The size of the matrix.
    
    // Use max(2*k+1, 20) as default subspace dimension for the gen case; same as MATLAB.
    ncv_default = blas_int( ((2*n_eigvals+1)>(20)) ? (2*n_eigvals+1) : (20) );
    
    // Use opts.subdim only if it's within the limits
    ncv = ncv_default;
    
    if(opts.subdim != 0)
      {
      if(opts.subdim < (n_eigvals + 3))
        {
        arma_debug_warn_level(1, "eigs_gen(): opts.subdim must be greater than k+2; using k+3 instead of ", opts.subdim);
        ncv = blas_int(n_eigvals + 3);
        }
      else
      if(blas_int(opts.subdim) > n)
        {
        arma_debug_warn_level(1, "eigs_gen(): opts.subdim cannot be greater than n_rows; using n_rows instead of ", opts.subdim);
        ncv = n;
        }
      else
        {
        ncv = blas_int(opts.subdim);
        }
      }
    
    // WARNING!!!
    // We are still not able to apply truly complex shifts to real matrices,
    // in which case the OP that ARPACK wants is different (see [s/d]naupd).
    // Also, if sigma contains a non-zero imaginary part, retrieving the eigenvalues
    // becomes utterly messy (see [s/d]eupd, remark #3).
    // We should never get to the point in which the imaginary part of sigma is non-zero;
    // the user-facing functions currently convert X from real to complex if a complex sigma is detected.
    // The check here is just for extra safety, and as a reminder of what's missing.
    T sigmar = real(sigma);
    T sigmai = imag(sigma);
    
    if(use_sigma)
    //if(form_val == form_sigma)
      {
      if(sigmai != T(0))  { arma_stop_logic_error("eigs_gen(): complex 'sigma' not applicable to real matrix"); return false; }
    
      run_aupd_shiftinvert(n_eigvals, sigmar, X, false /* gen, not sym */, n, tol, maxiter, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info);
      }
    else
      {
      run_aupd_plain(n_eigvals, which, X, false /* gen, not sym */, n, tol, maxiter, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info);
      }
    
    if(info != 0)  { return false; }
    
    // The process has converged, and now we need to recover the actual eigenvectors using neupd().
    blas_int rvec = 1; // .TRUE
    blas_int nev  = blas_int(n_eigvals);
    
    char howmny = 'A';
    char bmat   = 'I'; // We are considering the standard eigenvalue problem.
    
    podarray<blas_int> select(ncv);      // logical array of dimension NCV
    podarray<T>        dr(nev + 1);      // real array of dimension NEV + 1
    podarray<T>        di(nev + 1);      // real array of dimension NEV + 1
    podarray<T>        z(n * (nev + 1)); // real N by NEV array if HOWMNY = 'A'
    blas_int ldz = n;
    podarray<T>        workev(3 * ncv);
    
    select.zeros();
    dr.zeros();
    di.zeros();
    z.zeros();
    workev.zeros();
    
    arpack::neupd(&rvec, &howmny, select.memptr(), dr.memptr(), di.memptr(), z.memptr(), &ldz, (T*) &sigmar, (T*) &sigmai, workev.memptr(), &bmat, &n, which, &nev, &tol, resid.memptr(), &ncv, v.memptr(), &ldv, iparam.memptr(), ipntr.memptr(), workd.memptr(), workl.memptr(), &lworkl, rwork.memptr(), &info);
    
    // Check for errors.
    if(info != 0)  { arma_debug_warn_level(1, "eigs_gen(): ARPACK error ", info, " in neupd()"); return false; }
    
    // Put it into the outputs.
    eigval.set_size(n_eigvals);
    eigvec.zeros(n, n_eigvals);
    
    for(uword i = 0; i < n_eigvals; ++i)
      {
      eigval[i] = std::complex<T>(dr[i], di[i]);
      }
    
    // Now recover the eigenvectors.
    for(uword i = 0; i < n_eigvals; ++i)
      {
      // ARPACK ?neupd lays things out kinda odd in memory;
      // so does LAPACK ?geev -- see auxlib::eig_gen()
      if((i < n_eigvals - 1) && (eigval[i] == std::conj(eigval[i + 1])))
        {
        for(uword j = 0; j < uword(n); ++j)
          {
          eigvec.at(j, i)     = std::complex<T>(z[n * i + j],  z[n * (i + 1) + j]);
          eigvec.at(j, i + 1) = std::complex<T>(z[n * i + j], -z[n * (i + 1) + j]);
          }
        ++i; // Skip the next one.
        }
      else
      if((i == n_eigvals - 1) && (std::complex<T>(eigval[i]).imag() != 0.0))
        {
        // We don't have the matched conjugate eigenvalue.
        for(uword j = 0; j < uword(n); ++j)
          {
          eigvec.at(j, i) = std::complex<T>(z[n * i + j], z[n * (i + 1) + j]);
          }
        }
      else
        {
        // The eigenvector is entirely real.
        for(uword j = 0; j < uword(n); ++j)
          {
          eigvec.at(j, i) = std::complex<T>(z[n * i + j], T(0));
          }
        }
      }
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_ignore(n_eigvals);
    arma_ignore(form_val);
    arma_ignore(sigma);
    arma_ignore(opts);
    
    return false;
    }
  #endif
  }



//! immediate eigendecomposition of non-symmetric complex sparse object
template<typename T, typename T1>
inline
bool
sp_auxlib::eigs_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpBase< std::complex<T>, T1>& X_expr, const uword n_eigvals, const form_type form_val, const eigs_opts& opts)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> U(X_expr.get_ref());
  
  arma_debug_check( (U.M.is_square() == false), "eigs_gen(): given matrix must be square sized" );
  
  if(arma_config::check_nonfinite && U.M.has_nonfinite())
    {
    arma_debug_warn_level(3, "eigs_gen(): detected non-finite elements");
    return false;
    }
  
  constexpr std::complex<T> sigma = T(0);
  
  return sp_auxlib::eigs_gen<T, false>(eigval, eigvec, U.M, n_eigvals, form_val, sigma, opts);
  }



//! immediate eigendecomposition of non-symmetric complex sparse object
template<typename T, typename T1>
inline
bool
sp_auxlib::eigs_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpBase< std::complex<T>, T1>& X, const uword n_eigvals, const std::complex<T> sigma, const eigs_opts& opts)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> U(X.get_ref());
  
  arma_debug_check( (U.M.is_square() == false), "eigs_gen(): given matrix must be square sized" );
  
  if(arma_config::check_nonfinite && U.M.has_nonfinite())
    {
    arma_debug_warn_level(3, "eigs_gen(): detected non-finite elements");
    return false;
    }
  
  #if (defined(ARMA_USE_ARPACK) && defined(ARMA_USE_SUPERLU))
    {
    constexpr form_type form_val = form_sigma;
    
    return sp_auxlib::eigs_gen<T, true>(eigval, eigvec, U.M, n_eigvals, form_val, sigma, opts);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(n_eigvals);
    arma_ignore(sigma);
    arma_ignore(opts);
    
    arma_stop_logic_error("eigs_gen(): use of ARPACK and SuperLU must be enabled to use 'sigma'");
    return false;
    }
  #endif
  }



template<typename T, bool use_sigma>
inline
bool
sp_auxlib::eigs_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpMat< std::complex<T> >& X, const uword n_eigvals, const form_type form_val, const std::complex<T> sigma, const eigs_opts& opts)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_ARPACK)
    {
    // typedef typename std::complex<T> eT;
    
    arma_debug_check( (form_val != form_lm) && (form_val != form_sm) && (form_val != form_lr) && (form_val != form_sr) && (form_val != form_li) && (form_val != form_si) && (form_val != form_sigma), "eigs_gen(): unknown form specified" );
    
    if(X.is_square() == false)  { return false; }
    
    char which_lm[3] = "LM";
    char which_sm[3] = "SM";
    char which_lr[3] = "LR";
    char which_sr[3] = "SR";
    char which_li[3] = "LI";
    char which_si[3] = "SI";
    
    char* which;
    
    switch(form_val)
      {
      case form_lm:  which = which_lm;  break;
      case form_sm:  which = which_sm;  break;
      case form_lr:  which = which_lr;  break;
      case form_sr:  which = which_sr;  break;
      case form_li:  which = which_li;  break;
      case form_si:  which = which_si;  break;
      
      default:       which = which_lm;
      }
    
    // Make sure we aren't asking for every eigenvalue.
    arma_debug_check( (n_eigvals + 1 >= X.n_rows), "eigs_gen(): n_eigvals + 1 must be less than the number of rows in the matrix" );
    
    // If the matrix is empty, the case is trivial.
    if( (X.n_cols == 0) || (n_eigvals == 0) ) // We already know n_cols == n_rows.
      {
      eigval.reset();
      eigvec.reset();
      return true;
      }
    
    // Set up variables that get used for neupd().
    blas_int n, ncv, ncv_default, ldv, lworkl, info, maxiter;
    
    T tol   = T(opts.tol);
    maxiter = blas_int(opts.maxiter);
    
    podarray< std::complex<T> > resid, v, workd, workl;
    podarray<blas_int> iparam, ipntr;
    podarray<T> rwork;
    
    n = blas_int(X.n_rows); // The size of the matrix.
    
    // Use max(2*k+1, 20) as default subspace dimension for the gen case; same as MATLAB.
    ncv_default = blas_int( ((2*n_eigvals+1)>(20)) ? (2*n_eigvals+1) : (20) );
    
    // Use opts.subdim only if it's within the limits
    ncv = ncv_default;
    
    if(opts.subdim != 0)
      {
      if(opts.subdim < (n_eigvals + 3))
        {
        arma_debug_warn_level(1, "eigs_gen(): opts.subdim must be greater than k+2; using k+3 instead of ", opts.subdim);
        ncv = blas_int(n_eigvals + 3);
        }
      else
      if(blas_int(opts.subdim) > n)
        {
        arma_debug_warn_level(1, "eigs_gen(): opts.subdim cannot be greater than n_rows; using n_rows instead of ", opts.subdim);
        ncv = n;
        }
      else
        {
        ncv = blas_int(opts.subdim);
        }
      }
    
    if(use_sigma)
    //if(form_val == form_sigma)
      {
      run_aupd_shiftinvert(n_eigvals, sigma, X, false /* gen, not sym */, n, tol, maxiter, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info);
      }
    else
      {
      run_aupd_plain(n_eigvals, which, X, false /* gen, not sym */, n, tol, maxiter, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info);
      }
    
    if(info != 0)  { return false; }
    
    // The process has converged, and now we need to recover the actual eigenvectors using neupd().
    blas_int rvec = 1; // .TRUE
    blas_int nev  = blas_int(n_eigvals);
    
    char howmny = 'A';
    char bmat   = 'I'; // We are considering the standard eigenvalue problem.
    
    podarray<blas_int>        select(ncv); // logical array of dimension NCV
    podarray<std::complex<T>> d(nev + 1);  // complex array of dimension NEV + 1
    podarray<std::complex<T>> z(n * nev);  // complex N by NEV array if HOWMNY = 'A'
    blas_int ldz = n;
    podarray<std::complex<T>> workev(2 * ncv);
    
    select.zeros();
    d.zeros();
    z.zeros();
    workev.zeros();
    
    // Prepare the outputs; neupd() will write directly to them.
    eigval.zeros(n_eigvals);
    eigvec.zeros(n, n_eigvals);
    
    arpack::neupd(&rvec, &howmny, select.memptr(), eigval.memptr(),
(std::complex<T>*) NULL, eigvec.memptr(), &ldz, (std::complex<T>*) &sigma, (std::complex<T>*) NULL, workev.memptr(), &bmat, &n, which, &nev, &tol, resid.memptr(), &ncv, v.memptr(), &ldv, iparam.memptr(), ipntr.memptr(), workd.memptr(), workl.memptr(), &lworkl, rwork.memptr(), &info);
    
    // Check for errors.
    if(info != 0)  { arma_debug_warn_level(1, "eigs_gen(): ARPACK error ", info, " in neupd()"); return false; }
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(eigval);
    arma_ignore(eigvec);
    arma_ignore(X);
    arma_ignore(n_eigvals);
    arma_ignore(form_val);
    arma_ignore(sigma);
    arma_ignore(opts);

    arma_stop_logic_error("eigs_gen(): use of ARPACK must be enabled for decomposition of complex matrices");
    return false;
    }
  #endif
  }



template<typename T1, typename T2>
inline
bool
sp_auxlib::spsolve_simple(Mat<typename T1::elem_type>& X, const SpBase<typename T1::elem_type, T1>& A_expr, const Base<typename T1::elem_type, T2>& B_expr, const superlu_opts& user_opts)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_SUPERLU)
    {
    typedef typename T1::elem_type eT;
    
    superlu::superlu_options_t  options;
    sp_auxlib::set_superlu_opts(options, user_opts);
    
    const unwrap_spmat<T1> tmp1(A_expr.get_ref());
    const SpMat<eT>& A =   tmp1.M;
    
    X = B_expr.get_ref();   // superlu::gssv() uses X as input (the B matrix) and as output (the solution)
    
    if(A.is_square() == false)
      {
      X.soft_reset();
      arma_stop_logic_error("spsolve(): solving under-determined / over-determined systems is currently not supported");
      return false;
      }
    
    arma_debug_check( (A.n_rows != X.n_rows), "spsolve(): number of rows in the given objects must be the same", [&](){ X.soft_reset(); } );
    
    if(A.is_empty() || X.is_empty())
      {
      X.zeros(A.n_cols, X.n_cols);
      return true;
      }
    
    if(A.n_nonzero == uword(0))  { X.soft_reset(); return false; }
    
    if(arma_config::check_nonfinite && (A.has_nonfinite() || X.has_nonfinite()))
      {
      arma_debug_warn_level(3, "spsolve(): detected non-finite elements");
      return false;
      }
    
    if(arma_config::debug)
      {
      bool overflow = false;
      
      overflow = (A.n_nonzero > INT_MAX);
      overflow = (A.n_rows > INT_MAX) || overflow;
      overflow = (A.n_cols > INT_MAX) || overflow;
      overflow = (X.n_rows > INT_MAX) || overflow;
      overflow = (X.n_cols > INT_MAX) || overflow;
      
      if(overflow)
        {
        arma_stop_runtime_error("spsolve(): integer overflow: matrix dimensions are too large for integer type used by SuperLU");
        return false;
        }
      }
    
    superlu_supermatrix_wrangler x;
    superlu_supermatrix_wrangler a;
    
    const bool status_x = wrap_to_supermatrix(x.get_ref(), X);
    const bool status_a = copy_to_supermatrix(a.get_ref(), A);
    
    if( (status_x == false) || (status_a == false) )  { X.soft_reset(); return false; }
    
    superlu_supermatrix_wrangler l;
    superlu_supermatrix_wrangler u;
    
    // paranoia: use SuperLU's memory allocation, in case it reallocs
    
    superlu_array_wrangler<int> perm_c(A.n_cols+1);  // extra paranoia: increase array length by 1
    superlu_array_wrangler<int> perm_r(A.n_rows+1);
    
    superlu_stat_wrangler stat;
    
    int info = 0; // Return code.
    
    arma_extra_debug_print("superlu::gssv()");
    superlu::gssv<eT>(&options, a.get_ptr(), perm_c.get_ptr(), perm_r.get_ptr(), l.get_ptr(), u.get_ptr(), x.get_ptr(), stat.get_ptr(), &info);
    
    
    // Process the return code.
    if( (info > 0) && (info <= int(A.n_cols)) )
      {
      // std::ostringstream tmp;
      // tmp << "spsolve(): could not solve system; LU factorisation completed, but detected zero in U(" << (info-1) << ',' << (info-1) << ')';
      // arma_debug_warn_level(1, tmp.str());
      }
    else
    if(info > int(A.n_cols))
      {
      arma_debug_warn_level(1, "spsolve(): memory allocation failure: could not allocate ", (info - int(A.n_cols)), " bytes");
      }
    else
    if(info < 0)
      {
      arma_debug_warn_level(1, "spsolve(): unknown SuperLU error code from gssv(): ", info);
      }
    
    // No need to extract the data from x, since it's using the same memory as X
    
    return (info == 0);
    }
  #else
    {
    arma_ignore(X);
    arma_ignore(A_expr);
    arma_ignore(B_expr);
    arma_ignore(user_opts);
    arma_stop_logic_error("spsolve(): use of SuperLU must be enabled");
    return false;
    }
  #endif
  }



template<typename T1, typename T2>
inline
bool
sp_auxlib::spsolve_refine(Mat<typename T1::elem_type>& X, typename T1::pod_type& out_rcond, const SpBase<typename T1::elem_type, T1>& A_expr, const Base<typename T1::elem_type, T2>& B_expr, const superlu_opts& user_opts)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_SUPERLU)
    {
    typedef typename T1::pod_type   T;
    typedef typename T1::elem_type eT;
    
    superlu::superlu_options_t  options;
    sp_auxlib::set_superlu_opts(options, user_opts);
    
    const unwrap_spmat<T1> tmp1(A_expr.get_ref());
    const SpMat<eT>& A =   tmp1.M;
    
    const unwrap<T2>          tmp2(B_expr.get_ref());
    const Mat<eT>& B_unwrap = tmp2.M;
    
    const bool B_is_modified = ( (user_opts.equilibrate) || (&B_unwrap == &X) );
    
    Mat<eT> B_copy;  if(B_is_modified)  { B_copy = B_unwrap; }
    
    const Mat<eT>& B = (B_is_modified) ?  B_copy : B_unwrap;
    
    if(A.is_square() == false)
      {
      X.soft_reset();
      arma_stop_logic_error("spsolve(): solving under-determined / over-determined systems is currently not supported");
      return false;
      }
    
    arma_debug_check( (A.n_rows != B.n_rows), "spsolve(): number of rows in the given objects must be the same", [&](){ X.soft_reset(); } );
    
    X.zeros(A.n_cols, B.n_cols);  // set the elements to zero, as we don't trust the SuperLU spaghetti code
    
    if(A.is_empty() || B.is_empty())  { return true; }
    
    if(A.n_nonzero == uword(0))  { X.soft_reset(); return false; }
    
    if(arma_config::check_nonfinite && (A.has_nonfinite() || X.has_nonfinite()))
      {
      arma_debug_warn_level(3, "spsolve(): detected non-finite elements");
      return false;
      }
    
    if(arma_config::debug)
      {
      bool overflow;
      
      overflow = (A.n_nonzero > INT_MAX);
      overflow = (A.n_rows > INT_MAX) || overflow;
      overflow = (A.n_cols > INT_MAX) || overflow;
      overflow = (B.n_rows > INT_MAX) || overflow;
      overflow = (B.n_cols > INT_MAX) || overflow;
      overflow = (X.n_rows > INT_MAX) || overflow;
      overflow = (X.n_cols > INT_MAX) || overflow;
      
      if(overflow)
        {
        arma_stop_runtime_error("spsolve(): integer overflow: matrix dimensions are too large for integer type used by SuperLU");
        return false;
        }
      }
    
    superlu_supermatrix_wrangler x;
    superlu_supermatrix_wrangler a;
    superlu_supermatrix_wrangler b;
    
    const bool status_x = wrap_to_supermatrix(x.get_ref(), X);
    const bool status_a = copy_to_supermatrix(a.get_ref(), A);  // NOTE: superlu::gssvx() modifies 'a' if equilibration is enabled
    const bool status_b = wrap_to_supermatrix(b.get_ref(), B);  // NOTE: superlu::gssvx() modifies 'b' if equilibration is enabled
    
    if( (status_x == false) || (status_a == false) || (status_b == false) )  { X.soft_reset(); return false; }
    
    superlu_supermatrix_wrangler l;
    superlu_supermatrix_wrangler u;
    
    // paranoia: use SuperLU's memory allocation, in case it reallocs
    
    superlu_array_wrangler<int> perm_c(A.n_cols+1);  // extra paranoia: increase array length by 1
    superlu_array_wrangler<int> perm_r(A.n_rows+1);
    superlu_array_wrangler<int>  etree(A.n_cols+1);
    
    superlu_array_wrangler<T>    R(A.n_rows+1);
    superlu_array_wrangler<T>    C(A.n_cols+1);
    superlu_array_wrangler<T> ferr(B.n_cols+1);
    superlu_array_wrangler<T> berr(B.n_cols+1);
    
    superlu::GlobalLU_t glu;
    arrayops::fill_zeros(reinterpret_cast<char*>(&glu), sizeof(superlu::GlobalLU_t));
    
    superlu::mem_usage_t  mu;
    arrayops::fill_zeros(reinterpret_cast<char*>(&mu), sizeof(superlu::mem_usage_t));
    
    superlu_stat_wrangler stat;
    
    char equed[8] = {};     // extra characters for paranoia
    T    rpg      = T(0);
    T    rcond    = T(0);
    int  info     = int(0); // Return code.
    
    char  work[8] = {};
    int  lwork    = int(0);  // 0 means superlu will allocate memory
    
    arma_extra_debug_print("superlu::gssvx()");
    superlu::gssvx<eT>(&options, a.get_ptr(), perm_c.get_ptr(), perm_r.get_ptr(), etree.get_ptr(), equed, R.get_ptr(), C.get_ptr(), l.get_ptr(), u.get_ptr(), &work[0], lwork, b.get_ptr(), x.get_ptr(), &rpg, &rcond, ferr.get_ptr(), berr.get_ptr(), &glu, &mu, stat.get_ptr(), &info);
    
    bool status = false;
    
    // Process the return code.
    if(info == 0)
      {
      status = true;
      }
    if( (info > 0) && (info <= int(A.n_cols)) )
      {
      // std::ostringstream tmp;
      // tmp << "spsolve(): could not solve system; LU factorisation completed, but detected zero in U(" << (info-1) << ',' << (info-1) << ')';
      // arma_debug_warn_level(1, tmp.str());
      }
    else
    if( (info == int(A.n_cols+1)) && (user_opts.allow_ugly) )
      {
      arma_debug_warn_level(2, "spsolve(): system is singular to working precision (rcond: ", rcond, ")");
      status = true;
      }
    else
    if(info > int(A.n_cols+1))
      {
      arma_debug_warn_level(1, "spsolve(): memory allocation failure");
      }
    else
    if(info < 0)
      {
      arma_debug_warn_level(1, "spsolve(): unknown SuperLU error code from gssvx(): ", info);
      }
    
    // No need to extract the data from x, since it's using the same memory as X
    
    out_rcond = rcond;
    
    return status;
    }
  #else
    {
    arma_ignore(X);
    arma_ignore(out_rcond);
    arma_ignore(A_expr);
    arma_ignore(B_expr);
    arma_ignore(user_opts);
    arma_stop_logic_error("spsolve(): use of SuperLU must be enabled");
    return false;
    }
  #endif
  }



// template<typename T1>
// inline
// typename T1::pod_type
// sp_auxlib::rcond(const SpBase<typename T1::elem_type, T1>& A_expr)
//   {
//   arma_extra_debug_sigprint();
//   
//   #if defined(ARMA_USE_SUPERLU)
//     {
//     typedef typename T1::pod_type   T;
//     typedef typename T1::elem_type eT;
//     
//     const unwrap_spmat<T1> tmp1(A_expr.get_ref());
//     const SpMat<eT>& A =   tmp1.M;
//     
//     arma_debug_check( (A.is_square() == false), "rcond(): matrix must be square sized" );
//     
//     superlu_opts superlu_opts_default;
//     superlu::superlu_options_t options;
//     sp_auxlib::set_superlu_opts(options, superlu_opts_default);
//     int lwork = 0;
//     
//     superlu::GlobalLU_t Glu;
//     arrayops::fill_zeros(reinterpret_cast<char*>(&Glu), sizeof(superlu::GlobalLU_t));
//     
//     superlu_supermatrix_wrangler a;
//     superlu_supermatrix_wrangler aC;
//     
//     const bool status_a = sp_auxlib::copy_to_supermatrix(a.get_ref(), A);
//     
//     if(status_a == false)  { arma_stop_runtime_error("rcond(): could not construct SuperLU matrix"); return T(0); }
//     
//     superlu_supermatrix_wrangler l;
//     superlu_supermatrix_wrangler u;
//     
//     superlu_array_wrangler<int> perm_c(A.n_cols+1);  // paranoia: increase array length by 1
//     superlu_array_wrangler<int> perm_r(A.n_rows+1);
//     superlu_array_wrangler<int>  etree(A.n_cols+1);
//     
//     superlu_stat_wrangler stat;
//     
//     int panel_size = superlu::sp_ispec_environ(1);
//     int relax      = superlu::sp_ispec_environ(2);
//     int slu_info   = 0; // Return code.
//     
//     arma_extra_debug_print("superlu::gstrf()");
//     superlu::get_permutation_c(options.ColPerm, a.get_ptr(), perm_c.get_ptr());
//     superlu::sp_preorder_mat(&options, a.get_ptr(), perm_c.get_ptr(), etree.get_ptr(), aC.get_ptr());
//     superlu::gstrf<T>(&options, aC.get_ptr(), relax, panel_size, etree.get_ptr(), NULL, lwork, perm_c.get_ptr(), perm_r.get_ptr(), l.get_ptr(), u.get_ptr(), &Glu, stat.get_ptr(), &slu_info);
//     
//     if(slu_info != 0)  { return T(0); }
//     
//     T a_norm_val = sp_auxlib::norm1<eT>(a.get_ptr());
//     T a_rcond    = sp_auxlib::lu_rcond<eT>(l.get_ptr(), u.get_ptr(), a_norm_val);
//     
//     if(arma_isnan(a_rcond))  { return T(0); }
//     
//     return a_rcond;
//     }
//   #else
//     {
//     typename T1::pod_type T;
//     
//     arma_ignore(A_expr);
//     arma_stop_logic_error("rcond(): use of SuperLU must be enabled");
//     
//     return T(0);
//     }
//   #endif
//   }



#if defined(ARMA_USE_SUPERLU)
  
  template<typename eT>
  inline
  typename get_pod_type<eT>::result
  sp_auxlib::norm1(superlu::SuperMatrix* A)
    {
    arma_extra_debug_sigprint();
    
    char norm_id = '1';
    
    return superlu::langs<eT>(&norm_id, A);
    }
  
  
  
  template<typename eT>
  inline
  typename get_pod_type<eT>::result
  sp_auxlib::lu_rcond(superlu::SuperMatrix* L, superlu::SuperMatrix* U, typename get_pod_type<eT>::result norm_val)
    {
    arma_extra_debug_sigprint();
    
    typedef typename get_pod_type<eT>::result T;
    
    char norm_id   = '1';
    T    rcond_out = T(0);
    int  info      = int(0);
    
    superlu_stat_wrangler stat;
    
    superlu::gscon<eT>(&norm_id, L, U, norm_val, &rcond_out, stat.get_ptr(), &info);
    
    return (info == 0) ? T(rcond_out) : T(0);
    }
  
  
  
  inline
  void
  sp_auxlib::set_superlu_opts(superlu::superlu_options_t& options, const superlu_opts& user_opts)
    {
    arma_extra_debug_sigprint();
    
    // default options as the starting point
    superlu::set_default_opts(&options);
    
    // our settings
    options.Trans           = superlu::NOTRANS;
    options.ConditionNumber = superlu::YES;
   
    // process user_opts
    
    if(user_opts.equilibrate == true)   { options.Equil = superlu::YES; }
    if(user_opts.equilibrate == false)  { options.Equil = superlu::NO;  }
    
    if(user_opts.symmetric == true)   { options.SymmetricMode = superlu::YES; }
    if(user_opts.symmetric == false)  { options.SymmetricMode = superlu::NO;  }
    
    options.DiagPivotThresh = user_opts.pivot_thresh;
    
    if(user_opts.permutation == superlu_opts::NATURAL)        { options.ColPerm = superlu::NATURAL;       }
    if(user_opts.permutation == superlu_opts::MMD_ATA)        { options.ColPerm = superlu::MMD_ATA;       }
    if(user_opts.permutation == superlu_opts::MMD_AT_PLUS_A)  { options.ColPerm = superlu::MMD_AT_PLUS_A; }
    if(user_opts.permutation == superlu_opts::COLAMD)         { options.ColPerm = superlu::COLAMD;        }
    
    if(user_opts.refine == superlu_opts::REF_NONE)    { options.IterRefine = superlu::NOREFINE;   }
    if(user_opts.refine == superlu_opts::REF_SINGLE)  { options.IterRefine = superlu::SLU_SINGLE; }
    if(user_opts.refine == superlu_opts::REF_DOUBLE)  { options.IterRefine = superlu::SLU_DOUBLE; }
    if(user_opts.refine == superlu_opts::REF_EXTRA)   { options.IterRefine = superlu::SLU_EXTRA;  }
    }
  
  
  
  template<typename eT>
  inline
  bool
  sp_auxlib::copy_to_supermatrix(superlu::SuperMatrix& out, const SpMat<eT>& A)
    {
    arma_extra_debug_sigprint();
    
    // We store in column-major CSC.
    out.Stype = superlu::SLU_NC;
    
         if(    is_float<eT>::value)  { out.Dtype = superlu::SLU_S; }
    else if(   is_double<eT>::value)  { out.Dtype = superlu::SLU_D; }
    else if( is_cx_float<eT>::value)  { out.Dtype = superlu::SLU_C; }
    else if(is_cx_double<eT>::value)  { out.Dtype = superlu::SLU_Z; }
    
    out.Mtype = superlu::SLU_GE; // Just a general matrix.  We don't know more now.
    
    // We have to actually create the object which stores the data.
    // This gets cleaned by destroy_supermatrix().
    // We have to use SuperLU's problematic memory allocation routines since they are
    // not guaranteed to be new and delete.  See the comments in def_superlu.hpp
    superlu::NCformat* nc = (superlu::NCformat*)superlu::malloc(sizeof(superlu::NCformat));
    
    if(nc == nullptr)  { return false; }
    
    A.sync();
    
    nc->nnz    = A.n_nonzero;
    nc->nzval  = (void*)          superlu::malloc(sizeof(eT)             * A.n_nonzero   );
    nc->colptr = (superlu::int_t*)superlu::malloc(sizeof(superlu::int_t) * (A.n_cols + 1));
    nc->rowind = (superlu::int_t*)superlu::malloc(sizeof(superlu::int_t) * A.n_nonzero   );
    
    if( (nc->nzval == nullptr) || (nc->colptr == nullptr) || (nc->rowind == nullptr) )  { return false; }
    
    // Fill the matrix.
    arrayops::copy((eT*) nc->nzval, A.values, A.n_nonzero);
    
    // // These have to be copied by hand, because the types may differ.
    // for(uword i = 0; i <= A.n_cols; ++i)  { nc->colptr[i] = (int_t) A.col_ptrs[i]; }
    // for(uword i = 0; i < A.n_nonzero; ++i) { nc->rowind[i] = (int_t) A.row_indices[i]; }
    
    arrayops::convert(nc->colptr, A.col_ptrs,    A.n_cols+1 );
    arrayops::convert(nc->rowind, A.row_indices, A.n_nonzero);
    
    out.nrow  = superlu::int_t(A.n_rows);
    out.ncol  = superlu::int_t(A.n_cols);
    out.Store = (void*) nc;
    
    return true;
    }
  
  
  
  // memory efficient implementation of out = A - shift*I, where A is a square matrix
  template<typename eT>
  inline
  bool
  sp_auxlib::copy_to_supermatrix_with_shift(superlu::SuperMatrix& out, const SpMat<eT>& A, const eT shift)
    {
    arma_extra_debug_sigprint();
    
    arma_debug_check( (A.is_square() == false), "sp_auxlib::copy_to_supermatrix_with_shift(): given matrix must be square sized" );
    
    if(shift == eT(0))
      {
      arma_extra_debug_print("sp_auxlib::copy_to_supermatrix_with_shift(): shift is zero; redirecting to sp_auxlib::copy_to_supermatrix()");
      return sp_auxlib::copy_to_supermatrix(out, A);
      }
    
    // We store in column-major CSC.
    out.Stype = superlu::SLU_NC;
    
         if(    is_float<eT>::value)  { out.Dtype = superlu::SLU_S; }
    else if(   is_double<eT>::value)  { out.Dtype = superlu::SLU_D; }
    else if( is_cx_float<eT>::value)  { out.Dtype = superlu::SLU_C; }
    else if(is_cx_double<eT>::value)  { out.Dtype = superlu::SLU_Z; }
    
    out.Mtype = superlu::SLU_GE; // Just a general matrix.  We don't know more now.
    
    // We have to actually create the object which stores the data.
    // This gets cleaned by destroy_supermatrix().
    superlu::NCformat* nc = (superlu::NCformat*)superlu::malloc(sizeof(superlu::NCformat));
    
    if(nc == nullptr)  { return false; }
    
    A.sync();
    
    uword n_nonzero_diag_old = 0;
    uword n_nonzero_diag_new = 0;
    
    const uword n_search_cols = (std::min)(A.n_rows, A.n_cols);
    
    for(uword j=0; j < n_search_cols; ++j)
      {
      const uword      col_offset = A.col_ptrs[j    ];
      const uword next_col_offset = A.col_ptrs[j + 1];
      
      const uword* start_ptr = &(A.row_indices[     col_offset]);
      const uword*   end_ptr = &(A.row_indices[next_col_offset]);
      
      const uword wanted_row = j;
      
      const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, wanted_row);  // binary search
      
      if( (pos_ptr != end_ptr) && ((*pos_ptr) == wanted_row) )
        {
        // element on the main diagonal is non-zero
        ++n_nonzero_diag_old;
        
        const uword offset = uword(pos_ptr - start_ptr);
        const uword index  = offset + col_offset;
        
        const eT new_val = A.values[index] - shift;
        
        if(new_val != eT(0))  { ++n_nonzero_diag_new; }
        }
      else
        {
        // element on the main diagonal is zero, but sigma is non-zero,
        // so the number of new non-zero elments on the diagonal is increased
        ++n_nonzero_diag_new;
        }
      }
    
    const uword out_n_nonzero = A.n_nonzero - n_nonzero_diag_old + n_nonzero_diag_new;
    
    arma_extra_debug_print( arma_str::format("A.n_nonzero:        %d") % A.n_nonzero        );
    arma_extra_debug_print( arma_str::format("n_nonzero_diag_old: %d") % n_nonzero_diag_old );
    arma_extra_debug_print( arma_str::format("n_nonzero_diag_new: %d") % n_nonzero_diag_new );
    arma_extra_debug_print( arma_str::format("out_n_nonzero:      %d") % out_n_nonzero      );
    
    nc->nnz    = out_n_nonzero;
    nc->nzval  = (void*)          superlu::malloc(sizeof(eT)             * out_n_nonzero );
    nc->colptr = (superlu::int_t*)superlu::malloc(sizeof(superlu::int_t) * (A.n_cols + 1));
    nc->rowind = (superlu::int_t*)superlu::malloc(sizeof(superlu::int_t) * out_n_nonzero );
    
    if( (nc->nzval == nullptr) || (nc->colptr == nullptr) || (nc->rowind == nullptr) )  { return false; }
    
    // fill the matrix column by column, and insert diagonal elements when necessary
    
    nc->colptr[0] = 0;
    
    eT*             values_current = (eT*) nc->nzval;
    superlu::int_t* rowind_current = nc->rowind;
    
    uword count = 0;
    
    for(uword j=0; j < A.n_cols; ++j)
      {
      const uword idx_start = A.col_ptrs[j    ];
      const uword idx_end   = A.col_ptrs[j + 1];
      
      const eT* values_start = values_current;
      
      uword i = idx_start;
      
      // elements in the upper triangular part, excluding the main diagonal
      for(; (i < idx_end) && (A.row_indices[i] < j); ++i)
        {
        (*values_current) = A.values[i];
        (*rowind_current) = superlu::int_t(A.row_indices[i]);
        
        ++values_current;
        ++rowind_current;
        
        ++count;
        }
      
      // elements on the main diagonal
      if( (i < idx_end) && (A.row_indices[i] == j) )
        {
        // A(j,j) is non-zero
        
        const eT new_diag_val = A.values[i] - shift;
        
        if(new_diag_val != eT(0))
          {
          (*values_current) = new_diag_val;
          (*rowind_current) = superlu::int_t(j);
          
          ++values_current;
          ++rowind_current;
          
          ++count;
          }
        
        ++i;
        }
      else
        {
        // A(j,j) is zero, so insert a new element
        
        if(j < n_search_cols)   
          {
          (*values_current) = -shift;
          (*rowind_current) = superlu::int_t(j);
          
          ++values_current;
          ++rowind_current;
          
          ++count;
          }
        }
      
      // elements in the lower triangular part, excluding the main diagonal
      for(; i < idx_end; ++i)
        {
        (*values_current) = A.values[i];
        (*rowind_current) = superlu::int_t(A.row_indices[i]);
        
        ++values_current;
        ++rowind_current;
        
        ++count;
        }
      
      // number of non-zero elements in the j-th column of out
      const uword nnz_col = values_current - values_start;
      nc->colptr[j + 1] = superlu::int_t(nc->colptr[j] + nnz_col);
      }
    
    arma_extra_debug_print( arma_str::format("count: %d") % count );
    
    arma_check( (count != out_n_nonzero), "internal error: sp_auxlib::copy_to_supermatrix_with_shift(): count != out_n_nonzero" );
    
    out.nrow  = superlu::int_t(A.n_rows);
    out.ncol  = superlu::int_t(A.n_cols);
    out.Store = (void*) nc;
    
    return true;
    }
  
  
  
//   // for debugging only
//   template<typename eT>
//   inline
//   void
//   sp_auxlib::copy_to_spmat(SpMat<eT>& out, const superlu::SuperMatrix& A)
//     {
//     arma_extra_debug_sigprint();
//     
//     bool type_matched = false;
//     
//          if(    is_float<eT>::value)  { type_matched = (A.Dtype == superlu::SLU_S); }
//     else if(   is_double<eT>::value)  { type_matched = (A.Dtype == superlu::SLU_D); }
//     else if( is_cx_float<eT>::value)  { type_matched = (A.Dtype == superlu::SLU_C); }
//     else if(is_cx_double<eT>::value)  { type_matched = (A.Dtype == superlu::SLU_Z); }
//     
//     arma_debug_check( (type_matched == false),      "copy_to_spmat(): type mismatch"  );
//     arma_debug_check( (A.Mtype != superlu::SLU_GE), "copy_to_spmat(): unknown layout" );
//     
//     const superlu::NCformat* nc = (const superlu::NCformat*)(A.Store);
//     
//     if(nc == nullptr)  { out.reset(); return; }
//     
//     if( (nc->nzval == nullptr) || (nc->colptr == nullptr) || (nc->rowind == nullptr) )  { out.reset(); return; }
//     
//     const uword A_n_rows    = uword(A.nrow );
//     const uword A_n_cols    = uword(A.ncol );
//     const uword A_n_nonzero = uword(nc->nnz);
//     
//     if(A_n_nonzero == 0)  { out.zeros(A_n_rows, A_n_cols); return; }
//     
//     out.reserve(A_n_rows, A_n_cols, A_n_nonzero);
//     
//     arrayops::copy(access::rwp(out.values), (const eT*)(nc->nzval), A_n_nonzero);
//     
//     arrayops::convert(access::rwp(out.col_ptrs),    nc->colptr, A_n_cols+1 );
//     arrayops::convert(access::rwp(out.row_indices), nc->rowind, A_n_nonzero);
//     }
  
  
  
  template<typename eT>
  inline
  bool
  sp_auxlib::wrap_to_supermatrix(superlu::SuperMatrix& out, const Mat<eT>& A)
    {
    arma_extra_debug_sigprint();
    
    // NOTE: this function re-uses memory from matrix A
    
    // This is being stored as a dense matrix.
    out.Stype = superlu::SLU_DN;
    
         if(    is_float<eT>::value)  { out.Dtype = superlu::SLU_S; }
    else if(   is_double<eT>::value)  { out.Dtype = superlu::SLU_D; }
    else if( is_cx_float<eT>::value)  { out.Dtype = superlu::SLU_C; }
    else if(is_cx_double<eT>::value)  { out.Dtype = superlu::SLU_Z; }
    
    out.Mtype = superlu::SLU_GE;
    
    // We have to create the object that stores the data.
    superlu::DNformat* dn = (superlu::DNformat*)superlu::malloc(sizeof(superlu::DNformat));
    
    if(dn == nullptr)  { return false; }
    
    dn->lda   = A.n_rows;
    dn->nzval = (void*) A.memptr();  // re-use memory instead of copying
    
    out.nrow  = A.n_rows;
    out.ncol  = A.n_cols;
    out.Store = (void*) dn;
    
    return true;
    }
  
  
  
  inline
  void
  sp_auxlib::destroy_supermatrix(superlu::SuperMatrix& out)
    {
    arma_extra_debug_sigprint();
    
    // Clean up.
    if(out.Stype == superlu::SLU_NC)
      {
      superlu::destroy_compcol_mat(&out);
      }
    else
    if(out.Stype == superlu::SLU_NCP)
      {
      superlu::destroy_compcolperm_mat(&out);
      }
    else
    if(out.Stype == superlu::SLU_DN)
      {
      // superlu::destroy_dense_mat(&out);
      
      // since dn->nzval is set to re-use memory from a Mat object (which manages its own memory),
      // we cannot simply call superlu::destroy_dense_mat().
      // Only the out.Store structure can be freed.
      
      superlu::DNformat* dn = (superlu::DNformat*) out.Store;
      
      if(dn != nullptr)  { superlu::free(dn); }
      }
    else
    if(out.Stype == superlu::SLU_SC)
      {
      superlu::destroy_supernode_mat(&out);
      }
    else
      {
      // Uh, crap.
      
      std::ostringstream tmp;
      
      tmp << "sp_auxlib::destroy_supermatrix(): unhandled Stype" << std::endl;
      tmp << "Stype  val: " << out.Stype << std::endl;
      tmp << "Stype name: ";
      
      if(out.Stype == superlu::SLU_NC)      { tmp << "SLU_NC";     }
      if(out.Stype == superlu::SLU_NCP)     { tmp << "SLU_NCP";    }
      if(out.Stype == superlu::SLU_NR)      { tmp << "SLU_NR";     }
      if(out.Stype == superlu::SLU_SC)      { tmp << "SLU_SC";     }
      if(out.Stype == superlu::SLU_SCP)     { tmp << "SLU_SCP";    }
      if(out.Stype == superlu::SLU_SR)      { tmp << "SLU_SR";     }
      if(out.Stype == superlu::SLU_DN)      { tmp << "SLU_DN";     }
      if(out.Stype == superlu::SLU_NR_loc)  { tmp << "SLU_NR_loc"; }
      
      arma_debug_warn_level(1, tmp.str());
      arma_stop_runtime_error("internal error: sp_auxlib::destroy_supermatrix()");
      }
    }
  
#endif



template<typename eT, typename T>
inline
void
sp_auxlib::run_aupd_plain
  (
  const uword n_eigvals, char* which, 
  const SpMat<T>& X, const bool sym,
  blas_int& n, eT& tol, blas_int& maxiter,
  podarray<T>& resid, blas_int& ncv, podarray<T>& v, blas_int& ldv,
  podarray<blas_int>& iparam, podarray<blas_int>& ipntr,
  podarray<T>& workd, podarray<T>& workl, blas_int& lworkl, podarray<eT>& rwork,
  blas_int& info
  )
  {
  #if defined(ARMA_USE_ARPACK)
    {
    // ARPACK provides a "reverse communication interface" which is an
    // entertainingly archaic FORTRAN software engineering technique that
    // basically means that we call saupd()/naupd() and it tells us with some
    // return code what we need to do next (usually a matrix-vector product) and
    // then call it again.  So this results in some type of iterative process
    // where we call saupd()/naupd() many times.
    
    blas_int ido = 0; // This must be 0 for the first call.
    char bmat = 'I'; // We are considering the standard eigenvalue problem.
    n = X.n_rows; // The size of the matrix (should already be set outside).
    blas_int nev = n_eigvals;
    
    resid.zeros(n);
    
    // Two contraints on NCV: (NCV > NEV) for sym problems or
    // (NCV > NEV + 2) for gen problems and (NCV <= N)
    // 
    // We're calling either arpack::saupd() or arpack::naupd(),
    // which have slighly different minimum constraint and recommended value for NCV:
    // http://www.caam.rice.edu/software/ARPACK/UG/node136.html
    // http://www.caam.rice.edu/software/ARPACK/UG/node138.html
    
    if(ncv < (nev + (sym ? 1 : 3))) { ncv = (nev + (sym ? 1 : 3)); }
    if(ncv > n                    ) { ncv = n;                     }
    
    v.zeros(n * ncv); // Array N by NCV (output).
    rwork.zeros(ncv); // Work array of size NCV for complex calls.
    ldv = n; // "Leading dimension of V exactly as declared in the calling program."
    
    // IPARAM: integer array of length 11.
    iparam.zeros(11);
    iparam(0) = 1; // Exact shifts (not provided by us).
    iparam(2) = maxiter; // Maximum iterations; all the examples use 300, but they were written in the ancient times.
    iparam(6) = 1; // Mode 1: A * x = lambda * x.
    
    // IPNTR: integer array of length 14 (output).
    ipntr.zeros(14);
    
    // Real work array used in the basic Arnoldi iteration for reverse communication.
    workd.zeros(3 * n);
    
    // lworkl must be at least 3 * NCV^2 + 6 * NCV.
    lworkl = 3 * (ncv * ncv) + 6 * ncv;
    
    // Real work array of length lworkl.
    workl.zeros(lworkl);
    
    info = 0; // Set to 0 initially to use random initial vector.
    
    // All the parameters have been set or created.  Time to loop a lot.
    while(ido != 99)
      {
      // Call saupd() or naupd() with the current parameters.
      if(sym)
        {
        arma_extra_debug_print("arpack::saupd()");
        arpack::saupd(&ido, &bmat, &n, which, &nev, &tol, resid.memptr(), &ncv, v.memptr(), &ldv, iparam.memptr(), ipntr.memptr(), workd.memptr(), workl.memptr(), &lworkl, &info);
        }
      else
        {
        arma_extra_debug_print("arpack::naupd()");
        arpack::naupd(&ido, &bmat, &n, which, &nev, &tol, resid.memptr(), &ncv, v.memptr(), &ldv, iparam.memptr(), ipntr.memptr(), workd.memptr(), workl.memptr(), &lworkl, rwork.memptr(), &info);
        }
      
      // What do we do now?
      switch (ido)
        {
        case -1:
          // fallthrough
        case 1:
          {
          // We need to calculate the matrix-vector multiplication y = OP * x
          // where x is of length n and starts at workd(ipntr(0)), and y is of
          // length n and starts at workd(ipntr(1)).
          
          // operator*(sp_mat, vec) doesn't properly put the result into the
          // right place so we'll just reimplement it here for now...
          
          // Set the output to point at the right memory.  We have to subtract
          // one from FORTRAN pointers...
          Col<T> out(workd.memptr() + ipntr(1) - 1, n, false /* don't copy */);
          // Set the input to point at the right memory.
          Col<T> in(workd.memptr() + ipntr(0) - 1, n, false /* don't copy */);
          
          out.zeros();
          
          typename SpMat<T>::const_iterator X_it     = X.begin();
          typename SpMat<T>::const_iterator X_it_end = X.end();
          
          while(X_it != X_it_end)
            {
            const uword X_it_row = X_it.row();
            const uword X_it_col = X_it.col();
            
            out[X_it_row] += (*X_it) * in[X_it_col];
            ++X_it;
            }
          
          // No need to modify memory further since it was all done in-place.
          
          break;
          }
        case 99:
          // Nothing to do here, things have converged.
          break;
        default:
          {
          return; // Parent frame can look at the value of info.
          }
        }
      }
    
    // The process has ended; check the return code.
    if( (info != 0) && (info != 1) )
      {
      // Print warnings if there was a failure.
      
      if(sym)
        {
        arma_debug_warn_level(1, "eigs_sym(): ARPACK error ", info, " in saupd()");
        }
      else
        {
        arma_debug_warn_level(1, "eigs_gen(): ARPACK error ", info, " in naupd()");
        }
      
      return; // Parent frame can look at the value of info.
      }
    }
  #else
    {
    arma_ignore(n_eigvals);
    arma_ignore(which);
    arma_ignore(X);
    arma_ignore(sym);
    arma_ignore(n);
    arma_ignore(tol);
    arma_ignore(maxiter);
    arma_ignore(resid);
    arma_ignore(ncv);
    arma_ignore(v);
    arma_ignore(ldv);
    arma_ignore(iparam);
    arma_ignore(ipntr);
    arma_ignore(workd);
    arma_ignore(workl);
    arma_ignore(lworkl);
    arma_ignore(rwork);
    arma_ignore(info);
    }
  #endif
  }



// Here 'sigma' is 'T', but should be 'eT'.
// Applying complex shifts to real matrices is currently not directly implemented
template<typename eT, typename T>
inline
void
sp_auxlib::run_aupd_shiftinvert
  (
  const uword n_eigvals, const T sigma,
  const SpMat<T>& X, const bool sym,
  blas_int& n, eT& tol, blas_int& maxiter,
  podarray<T>& resid, blas_int& ncv, podarray<T>& v, blas_int& ldv,
  podarray<blas_int>& iparam, podarray<blas_int>& ipntr,
  podarray<T>& workd, podarray<T>& workl, blas_int& lworkl, podarray<eT>& rwork,
  blas_int& info
  )
  {
  // TODO: inconsistent use of type names: T can be complex while eT can be real
  
  #if (defined(ARMA_USE_ARPACK) && defined(ARMA_USE_SUPERLU))
    {
    char which_lm[3] = "LM";
    
    char* which = which_lm;  // NOTE: which_lm is the assumed operation when using shift-invert
    
    blas_int ido = 0; // This must be 0 for the first call.
    char bmat = 'I'; // We are considering the standard eigenvalue problem.
    n = X.n_rows; // The size of the matrix (should already be set outside).
    blas_int nev = n_eigvals;
    
    resid.zeros(n);
    
    // Two contraints on NCV: (NCV > NEV) for sym problems or
    // (NCV > NEV + 2) for gen problems and (NCV <= N)
    // 
    // We're calling either arpack::saupd() or arpack::naupd(),
    // which have slighly different minimum constraint and recommended value for NCV:
    // http://www.caam.rice.edu/software/ARPACK/UG/node136.html
    // http://www.caam.rice.edu/software/ARPACK/UG/node138.html
    
    if(ncv < (nev + (sym ? 1 : 3))) { ncv = (nev + (sym ? 1 : 3)); }
    if(ncv > n                    ) { ncv = n;                     }
    
    v.zeros(n * ncv); // Array N by NCV (output).
    rwork.zeros(ncv); // Work array of size NCV for complex calls.
    ldv = n; // "Leading dimension of V exactly as declared in the calling program."
    
    // IPARAM: integer array of length 11.
    iparam.zeros(11);
    iparam(0) = 1; // Exact shifts (not provided by us).
    iparam(2) = maxiter; // Maximum iterations; all the examples use 300, but they were written in the ancient times.
    // iparam(6) = 1; // Mode 1: A * x = lambda * x.
    
    // Change IPARAM for shift-invert
    iparam(6) = 3; // Mode 3:  A * x = lambda * M * x, M symmetric semi-definite. OP = inv[A - sigma*M]*M  (A complex)  or  Real_Part{ inv[A - sigma*M]*M }  (A real)  and  B = M.
    
    // IPNTR: integer array of length 14 (output).
    ipntr.zeros(14);
    
    // Real work array used in the basic Arnoldi iteration for reverse communication.
    workd.zeros(3 * n);
    
    // lworkl must be at least 3 * NCV^2 + 6 * NCV.
    lworkl = 3 * (ncv * ncv) + 6 * ncv;
    
    // Real work array of length lworkl.
    workl.zeros(lworkl);
    
    info = 0; // Set to 0 initially to use random initial vector.
    
    superlu_opts superlu_opts_default;
    superlu::superlu_options_t options;
    sp_auxlib::set_superlu_opts(options, superlu_opts_default);
    int lwork = 0;
    superlu::trans_t trans = superlu::NOTRANS;
    
    superlu::GlobalLU_t Glu; /* Not needed on return. */
    arrayops::fill_zeros(reinterpret_cast<char*>(&Glu), sizeof(superlu::GlobalLU_t));
    
    superlu_supermatrix_wrangler x;
    superlu_supermatrix_wrangler xC;
    
    const bool status_x = sp_auxlib::copy_to_supermatrix_with_shift(x.get_ref(), X, sigma);
    
    if(status_x == false)
      {
      arma_stop_runtime_error("run_aupd_shiftinvert(): could not construct SuperLU matrix");
      info = blas_int(-1);
      return;
      }
    
    // // for debugging only
    // if(true)
    //   {
    //   cout << "*** testing output of copy_to_supermatrix_with_shift()" << endl;
    //   cout << "*** sigma: " << sigma << endl;
    //   
    //   SpMat<T> Y(X);
    //   Y.diag() -= sigma;
    //   
    //   SpMat<T> Z;
    //   
    //   sp_auxlib::copy_to_spmat(Z, x.get_ref());
    //   
    //   cout << "*** size(Y): " << arma::size(Y) << endl;
    //   cout << "*** size(Z): " << arma::size(Z) << endl;
    //   cout << "*** accu(abs(Y)): " << accu(abs(Y)) << endl;
    //   cout << "*** accu(abs(Z)): " << accu(abs(Z)) << endl;
    //   
    //   if(arma::size(Y) == arma::size(Z))
    //     {
    //     cout << "*** error: " << accu(abs(Y-Z)) << endl;
    //     }
    //   }
    
    superlu_supermatrix_wrangler l;
    superlu_supermatrix_wrangler u;
    
    superlu_array_wrangler<int> perm_c(X.n_cols+1);  // paranoia: increase array length by 1
    superlu_array_wrangler<int> perm_r(X.n_rows+1);
    superlu_array_wrangler<int>  etree(X.n_cols+1);
    
    superlu_stat_wrangler stat;
    
    int panel_size = superlu::sp_ispec_environ(1);
    int relax      = superlu::sp_ispec_environ(2);
    int slu_info   = 0; // Return code.
    
    arma_extra_debug_print("superlu::gstrf()");
    superlu::get_permutation_c(options.ColPerm, x.get_ptr(), perm_c.get_ptr());
    superlu::sp_preorder_mat(&options, x.get_ptr(), perm_c.get_ptr(), etree.get_ptr(), xC.get_ptr());
    superlu::gstrf<T>(&options, xC.get_ptr(), relax, panel_size, etree.get_ptr(), NULL, lwork, perm_c.get_ptr(), perm_r.get_ptr(), l.get_ptr(), u.get_ptr(), &Glu, stat.get_ptr(), &slu_info);
    
    if(slu_info != 0)
      {
      arma_debug_warn_level(2, "matrix is singular to working precision");
      info = blas_int(-1);
      return;
      }
    
    eT x_norm_val = sp_auxlib::norm1<T>(x.get_ptr());
    eT x_rcond    = sp_auxlib::lu_rcond<T>(l.get_ptr(), u.get_ptr(), x_norm_val);
    
    if( (x_rcond < std::numeric_limits<eT>::epsilon()) || arma_isnan(x_rcond) )
      {
      arma_debug_warn_level(2, "matrix is singular to working precision (rcond: ", x_rcond, ")");
      info = blas_int(-1);
      return;
      }
    
    // All the parameters have been set or created.  Time to loop a lot.
    while(ido != 99)
      {
      // Call saupd() or naupd() with the current parameters.
      if(sym)
        {
        arma_extra_debug_print("arpack::saupd()");
        arpack::saupd(&ido, &bmat, &n, which, &nev, &tol, resid.memptr(), &ncv, v.memptr(), &ldv, iparam.memptr(), ipntr.memptr(), workd.memptr(), workl.memptr(), &lworkl, &info);
        }
      else
        {
        arma_extra_debug_print("arpack::naupd()");
        arpack::naupd(&ido, &bmat, &n, which, &nev, &tol, resid.memptr(), &ncv, v.memptr(), &ldv, iparam.memptr(), ipntr.memptr(), workd.memptr(), workl.memptr(), &lworkl, rwork.memptr(), &info);
        }
      
      // What do we do now?
      switch (ido)
        {
        case -1:
          // fallthrough
        case 1:
          {
          // We need to calculate the matrix-vector multiplication y = OP * x
          // where x is of length n and starts at workd(ipntr(0)), and y is of
          // length n and starts at workd(ipntr(1)).
          
          // operator*(sp_mat, vec) doesn't properly put the result into the
          // right place so we'll just reimplement it here for now...
          
          // Set the output to point at the right memory.  We have to subtract
          // one from FORTRAN pointers...
          Col<T> out(workd.memptr() + ipntr(1) - 1, n, false /* don't copy */);
          // Set the input to point at the right memory.
          Col<T> in(workd.memptr() + ipntr(0) - 1, n, false /* don't copy */);
          
          // Consider getting the LU factorization from ZGSTRF, and then
          // solve the system L*U*out = in (possibly with permutation matrix?)
          // Instead of "spsolve(out,X,in)" we call gstrf above and gstrs below
          
          out = in;
          superlu_supermatrix_wrangler out_slu;
          
          const bool status_out_slu = sp_auxlib::wrap_to_supermatrix(out_slu.get_ref(), out);
          
          if(status_out_slu == false)  { arma_stop_runtime_error("run_aupd_shiftinvert(): could not construct SuperLU matrix"); return; }
          
          arma_extra_debug_print("superlu::gstrs()");
          superlu::gstrs<T>(trans, l.get_ptr(), u.get_ptr(), perm_c.get_ptr(), perm_r.get_ptr(), out_slu.get_ptr(), stat.get_ptr(), &info);
          
          // No need to modify memory further since it was all done in-place.
          
          break;
          }
        case 99:
          // Nothing to do here, things have converged.
          break;
        default:
          {
          return; // Parent frame can look at the value of info.
          }
        }
      }
    
    // The process has ended; check the return code.
    if( (info != 0) && (info != 1) )
      {
      // Print warnings if there was a failure.
      
      if(sym)
        {
        arma_debug_warn_level(2, "eigs_sym(): ARPACK error ", info, " in saupd()");
        }
      else
        {
        arma_debug_warn_level(2, "eigs_gen(): ARPACK error ", info, " in naupd()");
        }
      
      return; // Parent frame can look at the value of info.
      }
    }
  #else
    {
    arma_ignore(n_eigvals);
    arma_ignore(sigma);
    arma_ignore(X);
    arma_ignore(sym);
    arma_ignore(n);
    arma_ignore(tol);
    arma_ignore(maxiter);
    arma_ignore(resid);
    arma_ignore(ncv);
    arma_ignore(v);
    arma_ignore(ldv);
    arma_ignore(iparam);
    arma_ignore(ipntr);
    arma_ignore(workd);
    arma_ignore(workl);
    arma_ignore(lworkl);
    arma_ignore(rwork);
    arma_ignore(info);
    }
  #endif
  }



template<typename eT>
inline
bool
sp_auxlib::rudimentary_sym_check(const SpMat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  if(X.n_rows != X.n_cols)  { return false; }
  
  const eT tol = eT(10000) * std::numeric_limits<eT>::epsilon();  // allow some leeway
  
  typename SpMat<eT>::const_iterator it     = X.begin();
  typename SpMat<eT>::const_iterator it_end = X.end();
  
  const uword n_check_limit = (std::max)( uword(2), uword(X.n_nonzero/100) );
  
  uword n_check = 1;
  
  while( (it != it_end) && (n_check <= n_check_limit) )
    {
    const uword it_row = it.row();
    const uword it_col = it.col();
    
    if(it_row != it_col)
      {
      const eT A = (*it);
      const eT B = X.at( it_col, it_row );  // deliberately swapped
      
      const eT C = (std::max)(std::abs(A), std::abs(B));
      
      const eT delta = std::abs(A - B);
     
      if(( (delta <= tol) || (delta <= (C * tol)) ) == false)  { return false; }
      
      ++n_check;
      }
    
    ++it;
    }
  
  return true;
  }



template<typename T>
inline
bool
sp_auxlib::rudimentary_sym_check(const SpMat< std::complex<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: the function name is a misnomer, as it checks for hermitian complex matrices;
  // NOTE: for simplicity of use, the function name is the same as for real matrices
  
  typedef typename std::complex<T> eT;
  
  if(X.n_rows != X.n_cols)  { return false; }
  
  const T tol = T(10000) * std::numeric_limits<T>::epsilon();  // allow some leeway
  
  typename SpMat<eT>::const_iterator it     = X.begin();
  typename SpMat<eT>::const_iterator it_end = X.end();
  
  const uword n_check_limit = (std::max)( uword(2), uword(X.n_nonzero/100) );
  
  uword n_check = 1;
  
  while( (it != it_end) && (n_check <= n_check_limit) )
    {
    const uword it_row = it.row();
    const uword it_col = it.col();
    
    if(it_row != it_col)
      {
      const eT A = (*it);
      const eT B = X.at( it_col, it_row );  // deliberately swapped
      
      const T C_real = (std::max)(std::abs(A.real()), std::abs(B.real()));
      const T C_imag = (std::max)(std::abs(A.imag()), std::abs(B.imag()));
      
      const T delta_real = std::abs(A.real() - B.real());
      const T delta_imag = std::abs(A.imag() + B.imag());  // take into account the conjugate
      
      const bool okay_real = ( (delta_real <= tol) || (delta_real <= (C_real * tol)) );
      const bool okay_imag = ( (delta_imag <= tol) || (delta_imag <= (C_imag * tol)) );
      
      if( (okay_real == false) || (okay_imag == false) )  { return false; }
      
      ++n_check;
      }
    else
      {
      const eT A = (*it);
      
      if(std::abs(A.imag()) > tol)  { return false; }
      }
    
    ++it;
    }
  
  return true;
  }



#if defined(ARMA_USE_SUPERLU)

inline
superlu_supermatrix_wrangler::~superlu_supermatrix_wrangler()
  {
  arma_extra_debug_sigprint_this(this);
  
  if(used == false)  { return; }
  
  char* m_char   = reinterpret_cast<char*>(&m);
  bool  all_zero = true;
  
  for(size_t i=0; i < sizeof(superlu::SuperMatrix); ++i)
    {
    if(m_char[i] != char(0))  { all_zero = false; break; }
    }
  
  if(all_zero == false)  { sp_auxlib::destroy_supermatrix(m); }
  }

inline
superlu_supermatrix_wrangler::superlu_supermatrix_wrangler()
  {
  arma_extra_debug_sigprint_this(this);
  
  arrayops::fill_zeros(reinterpret_cast<char*>(&m), sizeof(superlu::SuperMatrix));
  }

inline
superlu::SuperMatrix&
superlu_supermatrix_wrangler::get_ref()
  {
  used = true;
  
  return m;
  }

inline
superlu::SuperMatrix*
superlu_supermatrix_wrangler::get_ptr()
  {
  used = true;
  
  return &m;
  }


//


inline
superlu_stat_wrangler::~superlu_stat_wrangler()
  {
  arma_extra_debug_sigprint_this(this);
  
  superlu::free_stat(&stat);
  }

inline
superlu_stat_wrangler::superlu_stat_wrangler()
  {
  arma_extra_debug_sigprint_this(this);
  
  arrayops::fill_zeros(reinterpret_cast<char*>(&stat), sizeof(superlu::SuperLUStat_t));
  
  superlu::init_stat(&stat);
  }

inline
superlu::SuperLUStat_t*
superlu_stat_wrangler::get_ptr()
  {
  return &stat;
  }


//


template<typename eT>
inline
superlu_array_wrangler<eT>::~superlu_array_wrangler()
  {
  arma_extra_debug_sigprint_this(this);
  
  if(mem != nullptr)
    {
    superlu::free(mem);
    mem = nullptr;
    }
  }

template<typename eT>
inline
superlu_array_wrangler<eT>::superlu_array_wrangler(const uword n_elem)
  {
  arma_extra_debug_sigprint_this(this);
  
  mem = (eT*)(superlu::malloc(n_elem * sizeof(eT)));
  
  arma_check_bad_alloc( (mem == nullptr), "superlu::malloc(): out of memory" );
  
  arrayops::fill_zeros(mem, n_elem);
  }

template<typename eT>
inline
eT*
superlu_array_wrangler<eT>::get_ptr()
  {
  return mem;
  }

#endif


//! @}
