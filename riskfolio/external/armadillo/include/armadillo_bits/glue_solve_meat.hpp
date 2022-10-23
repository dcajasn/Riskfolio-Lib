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


//! \addtogroup glue_solve
//! @{



//
// glue_solve_gen_default


template<typename T1, typename T2>
inline
void
glue_solve_gen_default::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve_gen_default>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_solve_gen_default::apply(out, X.A, X.B);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("solve(): solution not found");
    }
  }



template<typename eT, typename T1, typename T2>
inline
bool
glue_solve_gen_default::apply(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  return glue_solve_gen_full::apply<eT,T1,T2,false>( out, A_expr, B_expr, uword(0));
  }



//
// glue_solve_gen_full


template<typename T1, typename T2>
inline
void
glue_solve_gen_full::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve_gen_full>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_solve_gen_full::apply( out, X.A, X.B, X.aux_uword );
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("solve(): solution not found");
    }
  }



template<typename eT, typename T1, typename T2, const bool has_user_flags>
inline
bool
glue_solve_gen_full::apply(Mat<eT>& actual_out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const uword flags)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  if(has_user_flags == true )  { arma_extra_debug_print("glue_solve_gen_full::apply(): has_user_flags = true" ); }
  if(has_user_flags == false)  { arma_extra_debug_print("glue_solve_gen_full::apply(): has_user_flags = false"); }
  
  const bool fast         = has_user_flags && bool(flags & solve_opts::flag_fast        );
  const bool equilibrate  = has_user_flags && bool(flags & solve_opts::flag_equilibrate );
  const bool no_approx    = has_user_flags && bool(flags & solve_opts::flag_no_approx   );
  const bool no_band      = has_user_flags && bool(flags & solve_opts::flag_no_band     );
  const bool no_sympd     = has_user_flags && bool(flags & solve_opts::flag_no_sympd    );
  const bool allow_ugly   = has_user_flags && bool(flags & solve_opts::flag_allow_ugly  );
  const bool likely_sympd = has_user_flags && bool(flags & solve_opts::flag_likely_sympd);
  const bool refine       = has_user_flags && bool(flags & solve_opts::flag_refine      );
  const bool no_trimat    = has_user_flags && bool(flags & solve_opts::flag_no_trimat   );
  const bool force_approx = has_user_flags && bool(flags & solve_opts::flag_force_approx);
  
  if(has_user_flags)
    {
    arma_extra_debug_print("glue_solve_gen_full::apply(): enabled flags:");
    
    if(fast        )  { arma_extra_debug_print("fast");         }
    if(equilibrate )  { arma_extra_debug_print("equilibrate");  }
    if(no_approx   )  { arma_extra_debug_print("no_approx");    }
    if(no_band     )  { arma_extra_debug_print("no_band");      }
    if(no_sympd    )  { arma_extra_debug_print("no_sympd");     }
    if(allow_ugly  )  { arma_extra_debug_print("allow_ugly");   }
    if(likely_sympd)  { arma_extra_debug_print("likely_sympd"); }
    if(refine      )  { arma_extra_debug_print("refine");       }
    if(no_trimat   )  { arma_extra_debug_print("no_trimat");    }
    if(force_approx)  { arma_extra_debug_print("force_approx"); }
    
    arma_debug_check( (fast     && equilibrate ), "solve(): options 'fast' and 'equilibrate' are mutually exclusive"      );
    arma_debug_check( (fast     && refine      ), "solve(): options 'fast' and 'refine' are mutually exclusive"           );
    arma_debug_check( (no_sympd && likely_sympd), "solve(): options 'no_sympd' and 'likely_sympd' are mutually exclusive" );
    }
  
  Mat<eT> A = A_expr.get_ref();
  
  if(force_approx)
    {
    arma_extra_debug_print("glue_solve_gen_full::apply(): forced approximate solution");
    
    arma_debug_check( no_approx, "solve(): options 'no_approx' and 'force_approx' are mutually exclusive" );
    
    if(fast)          { arma_debug_warn_level(2,  "solve(): option 'fast' ignored for forced approximate solution"         ); }
    if(equilibrate)   { arma_debug_warn_level(2,  "solve(): option 'equilibrate' ignored for forced approximate solution"  ); }
    if(refine)        { arma_debug_warn_level(2,  "solve(): option 'refine' ignored for forced approximate solution"       ); }
    if(likely_sympd)  { arma_debug_warn_level(2,  "solve(): option 'likely_sympd' ignored for forced approximate solution" ); }
    
    return auxlib::solve_approx_svd(actual_out, A, B_expr.get_ref());  // A is overwritten
    }
  
  // A_expr and B_expr can be used more than once (sympd optimisation fails or approximate solution required),
  // so ensure they are not overwritten in case we have aliasing
  
  bool is_alias = true;  // assume we have aliasing until we can prove otherwise
  
  if(is_Mat<T1>::value && is_Mat<T2>::value)
    {
    const quasi_unwrap<T1> UA( A_expr.get_ref() );
    const quasi_unwrap<T2> UB( B_expr.get_ref() );
    
    is_alias = UA.is_alias(actual_out) || UB.is_alias(actual_out);
    }
  
  Mat<eT>  tmp;
  Mat<eT>& out = (is_alias) ? tmp : actual_out;
  
  T    rcond  = T(0);
  bool status = false;
  
  if(A.n_rows == A.n_cols)
    {
    arma_extra_debug_print("glue_solve_gen_full::apply(): detected square system");
    
    uword KL = 0;
    uword KU = 0;
    
    const bool is_band  = arma_config::optimise_band && ((no_band || auxlib::crippled_lapack(A)) ? false : band_helper::is_band(KL, KU, A, uword(32)));
    
    const bool is_triu = (no_trimat || refine || equilibrate || likely_sympd || is_band           ) ? false : trimat_helper::is_triu(A);
    const bool is_tril = (no_trimat || refine || equilibrate || likely_sympd || is_band || is_triu) ? false : trimat_helper::is_tril(A);
    
    const bool try_sympd = arma_config::optimise_sympd && ((no_sympd || auxlib::crippled_lapack(A) || is_band || is_triu || is_tril) ? false : (likely_sympd ? true : sympd_helper::guess_sympd(A, uword(16))));
    
    if(fast)
      {
      // fast mode: solvers without refinement and without rcond estimate
      
      arma_extra_debug_print("glue_solve_gen_full::apply(): fast mode");
      
      if(is_band)
        {
        if( (KL == 1) && (KU == 1) )
          {
          arma_extra_debug_print("glue_solve_gen_full::apply(): fast + tridiagonal");
          
          status = auxlib::solve_tridiag_fast(out, A, B_expr.get_ref());
          }
        else
          {
          arma_extra_debug_print("glue_solve_gen_full::apply(): fast + band");
          
          status = auxlib::solve_band_fast(out, A, KL, KU, B_expr.get_ref());
          }
        }
      else
      if(is_triu || is_tril)
        {
        if(is_triu)  { arma_extra_debug_print("glue_solve_gen_full::apply(): fast + upper triangular matrix"); }
        if(is_tril)  { arma_extra_debug_print("glue_solve_gen_full::apply(): fast + lower triangular matrix"); }
        
        const uword layout = (is_triu) ? uword(0) : uword(1);
        
        status = auxlib::solve_trimat_fast(out, A, B_expr.get_ref(), layout);
        }
      else
      if(try_sympd)
        {
        arma_extra_debug_print("glue_solve_gen_full::apply(): fast + try_sympd");
        
        status = auxlib::solve_sympd_fast(out, A, B_expr.get_ref());  // A is overwritten
        
        if(status == false)
          {
          // auxlib::solve_sympd_fast() may have failed because A isn't really sympd
          
          arma_extra_debug_print("glue_solve_gen_full::apply(): auxlib::solve_sympd_fast() failed; retrying");
          
          A = A_expr.get_ref();
          
          status = auxlib::solve_square_fast(out, A, B_expr.get_ref());  // A is overwritten
          }
        }
      else
        {
        arma_extra_debug_print("glue_solve_gen_full::apply(): fast + dense");
        
        status = auxlib::solve_square_fast(out, A, B_expr.get_ref());  // A is overwritten
        }
      }
    else
    if(refine || equilibrate)
      {
      // refine mode: solvers with refinement and with rcond estimate
      
      arma_extra_debug_print("glue_solve_gen_full::apply(): refine mode");
      
      if(is_band)
        {
        arma_extra_debug_print("glue_solve_gen_full::apply(): refine + band");
        
        status = auxlib::solve_band_refine(out, rcond, A, KL, KU, B_expr, equilibrate);
        }
      else
      if(try_sympd)
        {
        arma_extra_debug_print("glue_solve_gen_full::apply(): refine + try_sympd");
        
        status = auxlib::solve_sympd_refine(out, rcond, A, B_expr.get_ref(), equilibrate);  // A is overwritten
        
        if( (status == false) && (rcond == T(0)) )
          {
          // auxlib::solve_sympd_refine() may have failed because A isn't really sympd;
          // in that case rcond is set to zero
          
          arma_extra_debug_print("glue_solve_gen_full::apply(): auxlib::solve_sympd_refine() failed; retrying");
          
          A = A_expr.get_ref();
          
          status = auxlib::solve_square_refine(out, rcond, A, B_expr.get_ref(), equilibrate);  // A is overwritten
          }
        }
      else
        {
        arma_extra_debug_print("glue_solve_gen_full::apply(): refine + dense");
        
        status = auxlib::solve_square_refine(out, rcond, A, B_expr, equilibrate);  // A is overwritten
        }
      }
    else
      {
      // default mode: solvers without refinement but with rcond estimate
      
      arma_extra_debug_print("glue_solve_gen_full::apply(): default mode");
      
      if(is_band)
        {
        arma_extra_debug_print("glue_solve_gen_full::apply(): rcond + band");
        
        status = auxlib::solve_band_rcond(out, rcond, A, KL, KU, B_expr.get_ref());
        }
      else
      if(is_triu || is_tril)
        {
        if(is_triu)  { arma_extra_debug_print("glue_solve_gen_full::apply(): rcond + upper triangular matrix"); }
        if(is_tril)  { arma_extra_debug_print("glue_solve_gen_full::apply(): rcond + lower triangular matrix"); }
        
        const uword layout = (is_triu) ? uword(0) : uword(1);
        
        status = auxlib::solve_trimat_rcond(out, rcond, A, B_expr.get_ref(), layout);
        }
      else
      if(try_sympd)
        {
        bool sympd_state = false;
        
        status = auxlib::solve_sympd_rcond(out, sympd_state, rcond, A, B_expr.get_ref());  // A is overwritten
        
        if( (status == false) && (sympd_state == false) )
          {
          arma_extra_debug_print("glue_solve_gen_full::apply(): auxlib::solve_sympd_rcond() failed; retrying");
          
          A = A_expr.get_ref();
          
          status = auxlib::solve_square_rcond(out, rcond, A, B_expr.get_ref());  // A is overwritten
          }
        }
      else
        {
        status = auxlib::solve_square_rcond(out, rcond, A, B_expr.get_ref());  // A is overwritten
        }
      }
    }
  else
    {
    arma_extra_debug_print("glue_solve_gen_full::apply(): detected non-square system");
    
    if(equilibrate)   { arma_debug_warn_level(2,  "solve(): option 'equilibrate' ignored for non-square matrix"  ); }
    if(refine)        { arma_debug_warn_level(2,  "solve(): option 'refine' ignored for non-square matrix"       ); }
    if(likely_sympd)  { arma_debug_warn_level(2,  "solve(): option 'likely_sympd' ignored for non-square matrix" ); }
    
    if(fast)
      {
      status = auxlib::solve_rect_fast(out, A, B_expr.get_ref());  // A is overwritten
      }
    else
      {
      status = auxlib::solve_rect_rcond(out, rcond, A, B_expr.get_ref());  // A is overwritten
      }
    }
  
  
  if( (status == true) && (fast == false) && (allow_ugly == false) && ((rcond < std::numeric_limits<T>::epsilon()) || arma_isnan(rcond)) )
    {
    status = false;
    }
  
  
  if( (status == false) && (no_approx == false) )
    {
    arma_extra_debug_print("glue_solve_gen_full::apply(): solving rank deficient system");
    
    if(rcond == T(0))
      {
      arma_debug_warn_level(2, "solve(): system is singular; attempting approx solution");
      }
    else
      {
      arma_debug_warn_level(2, "solve(): system is singular (rcond: ", rcond, "); attempting approx solution");
      }
    
    // TODO: conditionally recreate A: have a separate state flag which indicates whether A was previously overwritten
    
    A = A_expr.get_ref();  // as A may have been overwritten
    
    status = auxlib::solve_approx_svd(out, A, B_expr.get_ref());  // A is overwritten
    }
  
  if(is_alias)  { actual_out.steal_mem(out); }
  
  return status;
  }



//
// glue_solve_tri_default


template<typename T1, typename T2>
inline
void
glue_solve_tri_default::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve_tri_default>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_solve_tri_default::apply( out, X.A, X.B, X.aux_uword );
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("solve(): solution not found");
    }
  }



template<typename eT, typename T1, typename T2>
inline
bool
glue_solve_tri_default::apply(Mat<eT>& actual_out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const uword flags)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const bool triu = bool(flags & solve_opts::flag_triu);
  const bool tril = bool(flags & solve_opts::flag_tril);
  
  arma_extra_debug_print("glue_solve_tri_default::apply(): enabled flags:");
  
  if(triu)  { arma_extra_debug_print("triu"); }
  if(tril)  { arma_extra_debug_print("tril"); }
  
  const quasi_unwrap<T1> UA(A_expr.get_ref());
  const Mat<eT>& A     = UA.M;
  
  arma_debug_check( (A.is_square() == false), "solve(): matrix marked as triangular must be square sized" );
  
  const uword layout = (triu) ? uword(0) : uword(1);
  
  bool is_alias = true;
  
  if(is_Mat<T2>::value)
    {
    const quasi_unwrap<T2> UB(B_expr.get_ref());
    
    is_alias = UA.is_alias(actual_out) || UB.is_alias(actual_out);
    }
  
  T    rcond  = T(0);
  bool status = false;
  
  Mat<eT>  tmp;
  Mat<eT>& out = (is_alias) ? tmp : actual_out;
  
  status = auxlib::solve_trimat_rcond(out, rcond, A, B_expr.get_ref(), layout);  // A is not modified
  
  
  if( (status == true) && ( (rcond < std::numeric_limits<T>::epsilon()) || arma_isnan(rcond) ) )
    {
    status = false;
    }
  
  
  if(status == false)
    {
    arma_extra_debug_print("glue_solve_tri_default::apply(): solving rank deficient system");
    
    if(rcond == T(0))
      {
      arma_debug_warn_level(2, "solve(): system is singular; attempting approx solution");
      }
    else
      {
      arma_debug_warn_level(2, "solve(): system is singular (rcond: ", rcond, "); attempting approx solution");
      }
    
    Mat<eT> triA = (triu) ? trimatu(A) : trimatl(A);  // trimatu() and trimatl() return the same type
    
    status = auxlib::solve_approx_svd(out, triA, B_expr.get_ref());  // triA is overwritten
    }
  
  
  if(is_alias)  { actual_out.steal_mem(out); }
  
  return status;
  }



//
// glue_solve_tri_full


template<typename T1, typename T2>
inline
void
glue_solve_tri_full::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve_tri_full>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_solve_tri_full::apply( out, X.A, X.B, X.aux_uword );
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("solve(): solution not found");
    }
  }



template<typename eT, typename T1, typename T2>
inline
bool
glue_solve_tri_full::apply(Mat<eT>& actual_out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const uword flags)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  const bool fast         = bool(flags & solve_opts::flag_fast        );
  const bool equilibrate  = bool(flags & solve_opts::flag_equilibrate );
  const bool no_approx    = bool(flags & solve_opts::flag_no_approx   );
  const bool triu         = bool(flags & solve_opts::flag_triu        );
  const bool tril         = bool(flags & solve_opts::flag_tril        );
  const bool allow_ugly   = bool(flags & solve_opts::flag_allow_ugly  );
  const bool likely_sympd = bool(flags & solve_opts::flag_likely_sympd);
  const bool refine       = bool(flags & solve_opts::flag_refine      );
  const bool no_trimat    = bool(flags & solve_opts::flag_no_trimat   );
  const bool force_approx = bool(flags & solve_opts::flag_force_approx);
  
  arma_extra_debug_print("glue_solve_tri_full::apply(): enabled flags:");
  
  if(fast        )  { arma_extra_debug_print("fast");         }
  if(equilibrate )  { arma_extra_debug_print("equilibrate");  }
  if(no_approx   )  { arma_extra_debug_print("no_approx");    }
  if(triu        )  { arma_extra_debug_print("triu");         }
  if(tril        )  { arma_extra_debug_print("tril");         }
  if(allow_ugly  )  { arma_extra_debug_print("allow_ugly");   }
  if(likely_sympd)  { arma_extra_debug_print("likely_sympd"); }
  if(refine      )  { arma_extra_debug_print("refine");       }
  if(no_trimat   )  { arma_extra_debug_print("no_trimat");    }
  if(force_approx)  { arma_extra_debug_print("force_approx"); }
  
  if(no_trimat || equilibrate || refine || force_approx)
    {
    const uword mask = ~(solve_opts::flag_triu | solve_opts::flag_tril);
    
    return glue_solve_gen_full::apply(actual_out, ((triu) ? trimatu(A_expr.get_ref()) : trimatl(A_expr.get_ref())), B_expr, (flags & mask));
    }
  
  if(likely_sympd)  { arma_debug_warn_level(2, "solve(): option 'likely_sympd' ignored for triangular matrix"); }
  
  const quasi_unwrap<T1> UA(A_expr.get_ref());
  const Mat<eT>& A     = UA.M;
  
  arma_debug_check( (A.is_square() == false), "solve(): matrix marked as triangular must be square sized" );
  
  const uword layout = (triu) ? uword(0) : uword(1);
  
  bool is_alias = true;
  
  if(is_Mat<T2>::value)
    {
    const quasi_unwrap<T2> UB(B_expr.get_ref());
    
    is_alias = UA.is_alias(actual_out) || UB.is_alias(actual_out);
    }
  
  T    rcond  = T(0);
  bool status = false;
  
  Mat<eT>  tmp;
  Mat<eT>& out = (is_alias) ? tmp : actual_out;
  
  if(fast)
    {
    status = auxlib::solve_trimat_fast(out, A, B_expr.get_ref(), layout);  // A is not modified
    }
  else
    {
    status = auxlib::solve_trimat_rcond(out, rcond, A, B_expr.get_ref(), layout);  // A is not modified
    }
  
  
  if( (status == true) && (fast == false) && (allow_ugly == false) && ((rcond < std::numeric_limits<T>::epsilon()) || arma_isnan(rcond)) )
    {
    status = false;
    }
  
  
  if( (status == false) && (no_approx == false) )
    {
    arma_extra_debug_print("glue_solve_tri_full::apply(): solving rank deficient system");
    
    if(rcond == T(0))
      {
      arma_debug_warn_level(2, "solve(): system is singular; attempting approx solution");
      }
    else
      {
      arma_debug_warn_level(2, "solve(): system is singular (rcond: ", rcond, "); attempting approx solution");
      }
    
    Mat<eT> triA = (triu) ? trimatu(A) : trimatl(A);  // trimatu() and trimatl() return the same type
    
    status = auxlib::solve_approx_svd(out, triA, B_expr.get_ref());  // triA is overwritten
    }
  
  
  if(is_alias)  { actual_out.steal_mem(out); }
  
  return status;
  }



//! @}
