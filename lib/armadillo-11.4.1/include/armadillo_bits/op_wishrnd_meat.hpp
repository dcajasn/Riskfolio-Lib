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


//! \addtogroup op_wishrnd
//! @{


// implementation based on:
// Yu-Cheng Ku and Peter Bloomfield.
// Generating Random Wishart Matrices with Fractional Degrees of Freedom in OX.
// Oxmetrics User Conference, 2010.
  

template<typename T1>
inline
void
op_wishrnd::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_wishrnd>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const eT    df   = expr.aux;
  const uword mode = expr.aux_uword_a;
  
  const bool status = op_wishrnd::apply_direct(out, expr.m, df, mode);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("wishrnd(): given matrix is not symmetric positive definite");
    }
  }



template<typename T1>
inline
bool
op_wishrnd::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type df, const uword mode)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> U(X.get_ref());
  
  bool status = false;
  
  if(U.is_alias(out))
    {
    Mat<eT> tmp;
    
    if(mode == 1)  { status = op_wishrnd::apply_noalias_mode1(tmp, U.M, df); }
    if(mode == 2)  { status = op_wishrnd::apply_noalias_mode2(tmp, U.M, df); }
    
    out.steal_mem(tmp);
    }
  else
    {
    if(mode == 1)  { status = op_wishrnd::apply_noalias_mode1(out, U.M, df); }
    if(mode == 2)  { status = op_wishrnd::apply_noalias_mode2(out, U.M, df); }
    }
  
  return status;
  }



template<typename eT>
inline
bool
op_wishrnd::apply_noalias_mode1(Mat<eT>& out, const Mat<eT>& S, const eT df)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (S.is_square() == false), "wishrnd(): given matrix must be square sized" );
  
  if(S.is_empty())  { out.reset(); return true; }
  
  if(auxlib::rudimentary_sym_check(S) == false)  { return false; }
  
  Mat<eT> D;
  
  const bool status = op_chol::apply_direct(D, S, 0);
  
  if(status == false)  { return false; }
  
  return op_wishrnd::apply_noalias_mode2(out, D, df);
  }



template<typename eT>
inline
bool
op_wishrnd::apply_noalias_mode2(Mat<eT>& out, const Mat<eT>& D, const eT df)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (df <= eT(0)),            "df must be greater than zero"                 );
  arma_debug_check( (D.is_square() == false), "wishrnd(): given matrix must be square sized" );
  
  if(D.is_empty())  { out.reset(); return true; }
  
  const uword N = D.n_rows;
  
  if(df < eT(N))
    {
    arma_extra_debug_print("simple generator");
    
    const uword df_floor = uword(std::floor(df));
    
    const Mat<eT> tmp = (randn< Mat<eT> >(df_floor, N)) * D;
    
    out = tmp.t() * tmp;
    }
  else
    {
    arma_extra_debug_print("standard generator");
    
    op_chi2rnd_varying_df<eT> chi2rnd_generator;
    
    Mat<eT> A(N, N, arma_zeros_indicator());
    
    for(uword i=0; i<N; ++i)
      {
      A.at(i,i) = std::sqrt( chi2rnd_generator(df - eT(i)) );
      }
    
    for(uword i=1; i < N; ++i)
      {
      arma_rng::randn<eT>::fill( A.colptr(i), i );
      }
    
    const Mat<eT> tmp = A * D;
    
    A.reset();
    
    out = tmp.t() * tmp;
    }
  
  return true;
  }



//



template<typename T1>
inline
void
op_iwishrnd::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_iwishrnd>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const eT    df   = expr.aux;
  const uword mode = expr.aux_uword_a;
  
  const bool status = op_iwishrnd::apply_direct(out, expr.m, df, mode);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("iwishrnd(): given matrix is not symmetric positive definite and/or df is too low");
    }
  }



template<typename T1>
inline
bool
op_iwishrnd::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type df, const uword mode)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> U(X.get_ref());
  
  bool status = false;
  
  if(U.is_alias(out))
    {
    Mat<eT> tmp;
    
    if(mode == 1)  { status = op_iwishrnd::apply_noalias_mode1(tmp, U.M, df); }
    if(mode == 2)  { status = op_iwishrnd::apply_noalias_mode2(tmp, U.M, df); }
    
    out.steal_mem(tmp);
    }
  else
    {
    if(mode == 1)  { status = op_iwishrnd::apply_noalias_mode1(out, U.M, df); }
    if(mode == 2)  { status = op_iwishrnd::apply_noalias_mode2(out, U.M, df); }
    }
  
  return status;
  }



template<typename eT>
inline
bool
op_iwishrnd::apply_noalias_mode1(Mat<eT>& out, const Mat<eT>& T, const eT df)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (T.is_square() == false), "iwishrnd(): given matrix must be square sized" );
  
  if(T.is_empty())  { out.reset(); return true; }
  
  if(auxlib::rudimentary_sym_check(T) == false)  { return false; }
  
  Mat<eT> Tinv;
  Mat<eT> Dinv;
  
  const bool inv_status = auxlib::inv_sympd(Tinv, T);
  
  if(inv_status == false)  { return false; }
  
  const bool chol_status = op_chol::apply_direct(Dinv, Tinv, 0);
  
  if(chol_status == false)  { return false; }
  
  return op_iwishrnd::apply_noalias_mode2(out, Dinv, df);
  }



template<typename eT>
inline
bool
op_iwishrnd::apply_noalias_mode2(Mat<eT>& out, const Mat<eT>& Dinv, const eT df)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (df <= eT(0)),               "df must be greater than zero"                  );
  arma_debug_check( (Dinv.is_square() == false), "iwishrnd(): given matrix must be square sized" );
  
  if(Dinv.is_empty())  { out.reset(); return true; }
  
  Mat<eT> tmp;
  
  const bool wishrnd_status = op_wishrnd::apply_noalias_mode2(tmp, Dinv, df);
  
  if(wishrnd_status == false)  { return false; }
  
  const bool inv_status1 = auxlib::inv_sympd(out, tmp);
  
  const bool inv_status2 = (inv_status1) ? bool(true) : bool(auxlib::inv(out, tmp));
  
  if(inv_status2 == false)  { return false; }
  
  return true;
  }



//! @}
