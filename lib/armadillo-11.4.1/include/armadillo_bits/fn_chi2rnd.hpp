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


//! \addtogroup fn_chi2rnd
//! @{



arma_warn_unused
inline
double
chi2rnd(const double df)
  {
  arma_extra_debug_sigprint();
  
  op_chi2rnd_varying_df<double> generator;
  
  return generator(df);
  }



template<typename eT>
arma_warn_unused
inline
typename arma_real_only<eT>::result
chi2rnd(const eT df)
  {
  arma_extra_debug_sigprint();
  
  op_chi2rnd_varying_df<eT> generator;
  
  return generator(df);
  }



template<typename T1>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_real<typename T1::elem_type>::value),
  const Op<T1, op_chi2rnd>
  >::result
chi2rnd(const T1& expr)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_chi2rnd>(expr);
  }



template<typename obj_type>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_Mat<obj_type>::value && is_real<typename obj_type::elem_type>::value),
  obj_type
  >::result
chi2rnd(const typename obj_type::elem_type df, const uword n_rows, const uword n_cols)
  {
  arma_extra_debug_sigprint();
  
  if(is_Col<obj_type>::value)
    {
    arma_debug_check( (n_cols != 1), "chi2rnd(): incompatible size" );
    }
  else
  if(is_Row<obj_type>::value)
    {
    arma_debug_check( (n_rows != 1), "chi2rnd(): incompatible size" );
    }
  
  obj_type out(n_rows, n_cols, arma_nozeros_indicator());
  
  op_chi2rnd::fill_constant_df(out, df);
  
  return out;
  }



template<typename obj_type>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_Mat<obj_type>::value && is_real<typename obj_type::elem_type>::value),
  obj_type
  >::result
chi2rnd(const typename obj_type::elem_type df, const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  return chi2rnd<obj_type>(df, s.n_rows, s.n_cols);
  }



template<typename obj_type>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_Mat<obj_type>::value && is_real<typename obj_type::elem_type>::value),
  obj_type
  >::result
chi2rnd(const typename obj_type::elem_type df, const uword n_elem)
  {
  arma_extra_debug_sigprint();
  
  if(is_Row<obj_type>::value)
    {
    return chi2rnd<obj_type>(df, 1, n_elem);
    }
  else
    {
    return chi2rnd<obj_type>(df, n_elem, 1);
    }
  }



arma_warn_unused
inline
mat
chi2rnd(const double df, const uword n_rows, const uword n_cols)
  {
  arma_extra_debug_sigprint();
  
  return chi2rnd<mat>(df, n_rows, n_cols);
  }



arma_warn_unused
inline
mat
chi2rnd(const double df, const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  return chi2rnd<mat>(df, s.n_rows, s.n_cols);
  }



arma_warn_unused
inline
vec
chi2rnd(const double df, const uword n_elem)
  {
  arma_extra_debug_sigprint();
  
  return chi2rnd<vec>(df, n_elem, 1);
  }



//! @}
