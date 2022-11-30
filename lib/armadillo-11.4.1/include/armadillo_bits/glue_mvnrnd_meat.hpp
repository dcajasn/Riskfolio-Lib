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


//! \addtogroup glue_mvnrnd
//! @{


// implementation based on:
// James E. Gentle.
// Generation of Random Numbers.
// Computational Statistics, pp. 305-331, 2009.
// http://dx.doi.org/10.1007/978-0-387-98144-4_7


template<typename T1, typename T2>
inline
void
glue_mvnrnd_vec::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_mvnrnd_vec>& expr)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_mvnrnd::apply_direct(out, expr.A, expr.B, uword(1));
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("mvnrnd(): given covariance matrix is not symmetric positive semi-definite");
    }
  }



template<typename T1, typename T2>
inline
void
glue_mvnrnd::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_mvnrnd>& expr)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_mvnrnd::apply_direct(out, expr.A, expr.B, expr.aux_uword);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("mvnrnd(): given covariance matrix is not symmetric positive semi-definite");
    }
  }



template<typename T1, typename T2>
inline
bool
glue_mvnrnd::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& M, const Base<typename T1::elem_type,T2>& C, const uword N)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> UM(M.get_ref());
  const quasi_unwrap<T2> UC(C.get_ref());
  
  arma_debug_check( (UM.M.is_colvec() == false) && (UM.M.is_empty() == false),  "mvnrnd(): given mean must be a column vector"     );
  arma_debug_check( (UC.M.is_square() == false),  "mvnrnd(): given covariance matrix must be square sized"                         );
  arma_debug_check( (UM.M.n_rows != UC.M.n_rows), "mvnrnd(): number of rows in given mean vector and covariance matrix must match" );
  
  if( UM.M.is_empty() || UC.M.is_empty() )
    {
    out.set_size(0,N);
    return true;
    }
  
  if((arma_config::debug) && (auxlib::rudimentary_sym_check(UC.M) == false))
    {
    arma_debug_warn_level(1, "mvnrnd(): given matrix is not symmetric");
    }
  
  bool status = false;
  
  if(UM.is_alias(out) || UC.is_alias(out))
    {
    Mat<eT> tmp;
    
    status = glue_mvnrnd::apply_noalias(tmp, UM.M, UC.M, N);
    
    out.steal_mem(tmp);
    }
  else
    {
    status = glue_mvnrnd::apply_noalias(out, UM.M, UC.M, N);
    }
  
  return status;
  }



template<typename eT>
inline
bool
glue_mvnrnd::apply_noalias(Mat<eT>& out, const Mat<eT>& M, const Mat<eT>& C, const uword N)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> D;
  
  const bool chol_status = op_chol::apply_direct(D, C, 1);  // '1' means "lower triangular"
  
  if(chol_status == false)
    {
    // C is not symmetric positive definite, so find approximate square root of C
    
    Col<eT> eigval;  // NOTE: eT is constrained to be real (ie. float or double) in fn_mvnrnd.hpp
    Mat<eT> eigvec;
    
    const bool eig_status = eig_sym_helper(eigval, eigvec, C, 'd', "mvnrnd()");
    
    if(eig_status == false)  { return false; }
    
          eT*   eigval_mem    = eigval.memptr();
    const uword eigval_n_elem = eigval.n_elem;
    
    // since we're doing an approximation, tolerate tiny negative eigenvalues
    
    const eT tol = eT(-100) * Datum<eT>::eps * norm(C, "fro");
    
    if(arma_isfinite(tol) == false)  { return false; }
    
    for(uword i=0; i<eigval_n_elem; ++i)
      {
      const eT val = eigval_mem[i];
      
      if( (val < tol) || (arma_isfinite(val) == false) )  { return false; }
      }
    
    for(uword i=0; i<eigval_n_elem; ++i)  { if(eigval_mem[i] < eT(0))  { eigval_mem[i] = eT(0); } }
    
    Mat<eT> DD = eigvec * diagmat(sqrt(eigval));
    
    D.steal_mem(DD);
    }
  
  out = D * randn< Mat<eT> >(M.n_rows, N);
    
  if(N == 1)
    {
    out += M;
    }
  else
  if(N > 1)
    {
    out.each_col() += M;
    }
  
  return true;
  }



//! @}
