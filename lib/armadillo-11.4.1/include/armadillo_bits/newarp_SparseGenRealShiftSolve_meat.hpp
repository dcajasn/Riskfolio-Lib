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


namespace newarp
{


template<typename eT>
inline
SparseGenRealShiftSolve<eT>::SparseGenRealShiftSolve(const SpMat<eT>& mat_obj, const eT shift)
  #if defined(ARMA_USE_SUPERLU)
    : perm_c(mat_obj.n_cols + 1)
    , perm_r(mat_obj.n_rows + 1)
    , n_rows(mat_obj.n_rows)
    , n_cols(mat_obj.n_cols)
  #else
    : n_rows(0)
    , n_cols(0)
  #endif
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_SUPERLU)
    {
    // Derived from sp_auxlib::run_aupd_shiftinvert()
    superlu_opts superlu_opts_default;
    superlu::superlu_options_t options;
    sp_auxlib::set_superlu_opts(options, superlu_opts_default);
    
    superlu::GlobalLU_t Glu;
    arrayops::fill_zeros(reinterpret_cast<char*>(&Glu), sizeof(superlu::GlobalLU_t));
    
    superlu_supermatrix_wrangler x;
    superlu_supermatrix_wrangler xC;
    superlu_array_wrangler<int> etree(mat_obj.n_cols+1);
    
    // Copy A-shift*I to x
    const bool status_x = sp_auxlib::copy_to_supermatrix_with_shift(x.get_ref(), mat_obj, shift);
    
    if(status_x == false)  { arma_stop_runtime_error("newarp::SparseGenRealShiftSolve::SparseGenRealShiftSolve(): could not construct SuperLU matrix"); return; }
    
    int panel_size = superlu::sp_ispec_environ(1);
    int relax      = superlu::sp_ispec_environ(2);
    int slu_info   = 0; // Return code
    int lwork      = 0; // lwork = 0: allocate space internally by system malloc
    
    superlu_stat_wrangler stat;
    
    arma_extra_debug_print("superlu::gstrf()");
    superlu::get_permutation_c(options.ColPerm, x.get_ptr(), perm_c.get_ptr());
    superlu::sp_preorder_mat(&options, x.get_ptr(), perm_c.get_ptr(), etree.get_ptr(), xC.get_ptr());
    superlu::gstrf<eT>(&options, xC.get_ptr(), relax, panel_size, etree.get_ptr(), NULL, lwork, perm_c.get_ptr(), perm_r.get_ptr(), l.get_ptr(), u.get_ptr(), &Glu, stat.get_ptr(), &slu_info);
    
    if(slu_info != 0)
      {
      arma_debug_warn_level(2, "matrix is singular to working precision");
      return;
      }
    
    eT x_norm_val = sp_auxlib::norm1<eT>(x.get_ptr());
    eT x_rcond    = sp_auxlib::lu_rcond<eT>(l.get_ptr(), u.get_ptr(), x_norm_val);
    
    if( (x_rcond < std::numeric_limits<eT>::epsilon()) || arma_isnan(x_rcond) )
      {
      if(x_rcond == eT(0))  { arma_debug_warn_level(2, "matrix is singular to working precision");                        }
      else                  { arma_debug_warn_level(2, "matrix is singular to working precision (rcond: ", x_rcond, ")"); }
      return;
      }
    
    valid = true;
    }
  #else
    {
    arma_ignore(mat_obj);
    arma_ignore(shift);
    }
  #endif
  }



// Perform the shift-solve operation \f$y=(A-\sigma I)^{-1}x\f$.
// y_out = inv(A - sigma * I) * x_in
template<typename eT>
inline
void
SparseGenRealShiftSolve<eT>::perform_op(eT* x_in, eT* y_out) const
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_SUPERLU)
    {
    const Col<eT> x(x_in , n_cols, false, true);
          Col<eT> y(y_out, n_rows, false, true);
    
    // Derived from sp_auxlib::run_aupd_shiftinvert()
    y = x;
    superlu_supermatrix_wrangler out_slu;
    
    const bool status_out_slu = sp_auxlib::wrap_to_supermatrix(out_slu.get_ref(), y);
    
    if(status_out_slu == false)  { arma_stop_runtime_error("newarp::SparseGenRealShiftSolve::perform_op(): could not construct SuperLU matrix"); return; }
    
    superlu_stat_wrangler stat;
    int info = 0;
    
    arma_extra_debug_print("superlu::gstrs()");
    superlu::gstrs<eT>(superlu::NOTRANS, l.get_ptr(), u.get_ptr(), perm_c.get_ptr(), perm_r.get_ptr(), out_slu.get_ptr(), stat.get_ptr(), &info);
    
    if(info != 0)  { arma_stop_runtime_error("newarp::SparseGenRealShiftSolve::perform_op(): could not solve linear equation"); return; }
    
    // No need to modify memory further since it was all done in-place.
    }
  #else
    {
    arma_ignore(x_in);
    arma_ignore(y_out);
    }
  #endif
  }


}  // namespace newarp
