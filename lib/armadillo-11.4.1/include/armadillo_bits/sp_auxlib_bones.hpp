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


//! wrapper for accesing external functions in ARPACK and SuperLU
class sp_auxlib
  {
  public:
  
  enum form_type
    {
    form_none, form_lm, form_sm, form_lr, form_la, form_sr, form_li, form_si, form_sa, form_sigma
    };
  
  inline static form_type interpret_form_str(const char* form_str);
  
  //
  // eigs_sym() for real matrices
  
  template<typename eT, typename T1>
  inline static bool eigs_sym(Col<eT>& eigval, Mat<eT>& eigvec, const SpBase<eT, T1>& X, const uword n_eigvals, const form_type form_val, const eigs_opts& opts);
  
  template<typename eT, typename T1>
  inline static bool eigs_sym(Col<eT>& eigval, Mat<eT>& eigvec, const SpBase<eT, T1>& X, const uword n_eigvals, const eT sigma, const eigs_opts& opts);
  
  template<typename eT>
  inline static bool eigs_sym_newarp(Col<eT>& eigval, Mat<eT>& eigvec, const SpMat<eT>& X, const uword n_eigvals, const form_type form_val, const eigs_opts& opts);

  template<typename eT>
  inline static bool eigs_sym_newarp(Col<eT>& eigval, Mat<eT>& eigvec, const SpMat<eT>& X, const uword n_eigvals, const eT sigma, const eigs_opts& opts);

  template<typename eT, bool use_sigma>
  inline static bool eigs_sym_arpack(Col<eT>& eigval, Mat<eT>& eigvec, const SpMat<eT>& X, const uword n_eigvals, const form_type form_val, const eT sigma, const eigs_opts& opts);
  
  //
  // eigs_gen() for real matrices
  
  template<typename T, typename T1>
  inline static bool eigs_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpBase<T, T1>& X, const uword n_eigvals, const form_type form_val, const eigs_opts& opts);
  
  template<typename T, typename T1>
  inline static bool eigs_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpBase<T, T1>& X, const uword n_eigvals, const std::complex<T> sigma, const eigs_opts& opts);
  
  template<typename T>
  inline static bool eigs_gen_newarp(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpMat<T>& X, const uword n_eigvals, const form_type form_val, const eigs_opts& opts);
  
  template<typename T, bool use_sigma>
  inline static bool eigs_gen_arpack(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpMat<T>& X, const uword n_eigvals, const form_type form_val, const std::complex<T> sigma, const eigs_opts& opts);
  
  //
  // eigs_gen() for complex matrices
  
  template<typename T, typename T1>
  inline static bool eigs_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpBase< std::complex<T>, T1>& X, const uword n_eigvals, const form_type form_val, const eigs_opts& opts);
  
  template<typename T, typename T1>
  inline static bool eigs_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpBase< std::complex<T>, T1>& X, const uword n_eigvals, const std::complex<T> sigma, const eigs_opts& opts);
  
  template<typename T, bool use_sigma>
  inline static bool eigs_gen(Col< std::complex<T> >& eigval, Mat< std::complex<T> >& eigvec, const SpMat< std::complex<T> >& X, const uword n_eigvals, const form_type form_val, const std::complex<T> sigma, const eigs_opts& opts);
  
  //
  // spsolve() via SuperLU
  
  template<typename T1, typename T2>
  inline static bool spsolve_simple(Mat<typename T1::elem_type>& out, const SpBase<typename T1::elem_type, T1>& A, const Base<typename T1::elem_type, T2>& B, const superlu_opts& user_opts);
  
  template<typename T1, typename T2>
  inline static bool spsolve_refine(Mat<typename T1::elem_type>& out, typename T1::pod_type& out_rcond, const SpBase<typename T1::elem_type, T1>& A, const Base<typename T1::elem_type, T2>& B, const superlu_opts& user_opts);
  
  // //
  // // rcond() via SuperLU
  // 
  // template<typename T1>
  // sinline static typename T1::pod_type rcond(const SpBase<typename T1::elem_type, T1>& A);
  
  //
  // support functions
  
  #if defined(ARMA_USE_SUPERLU)
    
    template<typename eT>
    inline static typename get_pod_type<eT>::result norm1(superlu::SuperMatrix* A);
    
    template<typename eT>
    inline static typename get_pod_type<eT>::result lu_rcond(superlu::SuperMatrix* L, superlu::SuperMatrix* U, typename get_pod_type<eT>::result norm_val);
    
    inline static void set_superlu_opts(superlu::superlu_options_t& options, const superlu_opts& user_opts);
    
    template<typename eT>
    inline static bool copy_to_supermatrix(superlu::SuperMatrix& out, const SpMat<eT>& A);
    
    template<typename eT>
    inline static bool copy_to_supermatrix_with_shift(superlu::SuperMatrix& out, const SpMat<eT>& A, const eT shift);
    
    // // for debugging only
    // template<typename eT>
    // inline static void copy_to_spmat(SpMat<eT>& out, const superlu::SuperMatrix& A);
    
    template<typename eT>
    inline static bool wrap_to_supermatrix(superlu::SuperMatrix& out, const Mat<eT>& A);
    
    inline static void destroy_supermatrix(superlu::SuperMatrix& out);
  
  #endif
  
  
  
  private:
  
  // calls arpack saupd()/naupd() because the code is so similar for each
  // all of the extra variables are later used by seupd()/neupd(), but those
  // functions are very different and we can't combine their code
  
  template<typename eT, typename T>
  inline static void run_aupd_plain
    (
    const uword n_eigvals, char* which,
    const SpMat<T>& X, const bool sym,
    blas_int& n, eT& tol, blas_int& maxiter,
    podarray<T>& resid, blas_int& ncv, podarray<T>& v, blas_int& ldv,
    podarray<blas_int>& iparam, podarray<blas_int>& ipntr,
    podarray<T>& workd, podarray<T>& workl, blas_int& lworkl, podarray<eT>& rwork,
    blas_int& info
    );
  
  template<typename eT, typename T>
  inline static void run_aupd_shiftinvert
    (
    const uword n_eigvals, const T sigma,
    const SpMat<T>& X, const bool sym,
    blas_int& n, eT& tol, blas_int& maxiter,
    podarray<T>& resid, blas_int& ncv, podarray<T>& v, blas_int& ldv,
    podarray<blas_int>& iparam, podarray<blas_int>& ipntr,
    podarray<T>& workd, podarray<T>& workl, blas_int& lworkl, podarray<eT>& rwork,
    blas_int& info
    );
  
  
  template<typename eT>
  inline static bool rudimentary_sym_check(const SpMat<eT>& X);
  
  template<typename T>
  inline static bool rudimentary_sym_check(const SpMat< std::complex<T> >& X);
  };



#if defined(ARMA_USE_SUPERLU)

class superlu_supermatrix_wrangler
  {
  private:
  
  bool used = false;
  
  arma_aligned superlu::SuperMatrix m;
  
  public:
  
  inline ~superlu_supermatrix_wrangler();
  inline  superlu_supermatrix_wrangler();
  
  inline superlu_supermatrix_wrangler(const superlu_supermatrix_wrangler&) = delete;
  inline void operator=              (const superlu_supermatrix_wrangler&) = delete;
  
  inline superlu::SuperMatrix& get_ref();
  inline superlu::SuperMatrix* get_ptr();
  };


class superlu_stat_wrangler
  {
  private:
  
  arma_aligned superlu::SuperLUStat_t stat;
  
  public:
  
  inline ~superlu_stat_wrangler();
  inline  superlu_stat_wrangler();
  
  inline superlu_stat_wrangler(const superlu_stat_wrangler&) = delete;
  inline void operator=       (const superlu_stat_wrangler&) = delete;
  
  inline superlu::SuperLUStat_t* get_ptr();
  };


template<typename eT>
class superlu_array_wrangler
  {
  private:
  
  arma_aligned eT* mem = nullptr;
  
  public:
  
  inline ~superlu_array_wrangler();
  inline  superlu_array_wrangler(const uword n_elem);
  
  inline superlu_array_wrangler()                              = delete;
  inline superlu_array_wrangler(const superlu_array_wrangler&) = delete;
  inline void operator=        (const superlu_array_wrangler&) = delete;
  
  inline eT* get_ptr();
  };

#endif



//! @}

