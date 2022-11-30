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


//! \addtogroup fn_randg
//! @{



template<typename obj_type>
arma_warn_unused
inline
obj_type
randg(const uword n_rows, const uword n_cols, const distr_param& param = distr_param(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename obj_type::elem_type eT;
  
  if(is_Col<obj_type>::value)
    {
    arma_debug_check( (n_cols != 1), "randg(): incompatible size" );
    }
  else
  if(is_Row<obj_type>::value)
    {
    arma_debug_check( (n_rows != 1), "randg(): incompatible size" );
    }
  
  double a = double(1);
  double b = double(1);
  
  param.get_double_vals(a,b);
  
  arma_debug_check( ((a <= double(0)) || (b <= double(0))), "randg(): incorrect distribution parameters; a and b must be greater than zero" );
  
  obj_type out(n_rows, n_cols, arma_nozeros_indicator());
  
  arma_rng::randg<eT>::fill(out.memptr(), out.n_elem, a, b);
  
  return out;
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
randg(const SizeMat& s, const distr_param& param = distr_param(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return randg<obj_type>(s.n_rows, s.n_cols, param);
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
randg(const uword n_elem, const distr_param& param = distr_param(), const arma_empty_class junk1 = arma_empty_class(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk2 = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  const uword n_rows = (is_Row<obj_type>::value) ? uword(1) : n_elem;
  const uword n_cols = (is_Row<obj_type>::value) ? n_elem   : uword(1);
  
  return randg<obj_type>(n_rows, n_cols, param);
  }



arma_warn_unused
inline
mat
randg(const uword n_rows, const uword n_cols, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randg<mat>(n_rows, n_cols, param);
  }



arma_warn_unused
inline
mat
randg(const SizeMat& s, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randg<mat>(s.n_rows, s.n_cols, param);
  }



arma_warn_unused
inline
vec
randg(const uword n_elem, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randg<vec>(n_elem, uword(1), param);
  }



arma_warn_unused
inline
double
randg(const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return as_scalar( randg<vec>(uword(1), uword(1), param) );
  }



template<typename eT>
arma_warn_unused
inline
typename arma_real_or_cx_only<eT>::result
randg(const distr_param& param = distr_param())
  {
  return eT( as_scalar( randg< Col<eT> >(uword(1), uword(1), param) ) );
  }



template<typename cube_type>
arma_warn_unused
inline
cube_type
randg(const uword n_rows, const uword n_cols, const uword n_slices, const distr_param& param = distr_param(), const typename arma_Cube_only<cube_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename cube_type::elem_type eT;
  
  double a = double(1);
  double b = double(1);
  
  param.get_double_vals(a,b);
  
  arma_debug_check( ((a <= double(0)) || (b <= double(0))), "randg(): incorrect distribution parameters; a and b must be greater than zero" );
  
  cube_type out(n_rows, n_cols, n_slices, arma_nozeros_indicator());
  
  arma_rng::randg<eT>::fill(out.memptr(), out.n_elem, a, b);
  
  return out;
  }



template<typename cube_type>
arma_warn_unused
inline
cube_type
randg(const SizeCube& s, const distr_param& param = distr_param(), const typename arma_Cube_only<cube_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return randg<cube_type>(s.n_rows, s.n_cols, s.n_slices, param);
  }



arma_warn_unused
inline
cube
randg(const uword n_rows, const uword n_cols, const uword n_slices, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randg<cube>(n_rows, n_cols, n_slices, param);
  }



arma_warn_unused
inline
cube
randg(const SizeCube& s, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randg<cube>(s.n_rows, s.n_cols, s.n_slices, param);
  }



//! @}
