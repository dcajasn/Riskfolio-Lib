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


//! \addtogroup fn_randi
//! @{



template<typename obj_type>
arma_warn_unused
inline
obj_type
randi(const uword n_rows, const uword n_cols, const distr_param& param = distr_param(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename obj_type::elem_type eT;
  
  if(is_Col<obj_type>::value)
    {
    arma_debug_check( (n_cols != 1), "randi(): incompatible size" );
    }
  else
  if(is_Row<obj_type>::value)
    {
    arma_debug_check( (n_rows != 1), "randi(): incompatible size" );
    }
  
  int a = 0;
  int b = arma_rng::randi<eT>::max_val();
  
  param.get_int_vals(a,b);
  
  arma_debug_check( (a > b), "randi(): incorrect distribution parameters; a must be less than b" );
  
  obj_type out(n_rows, n_cols, arma_nozeros_indicator());
  
  arma_rng::randi<eT>::fill(out.memptr(), out.n_elem, a, b);
  
  return out;
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
randi(const SizeMat& s, const distr_param& param = distr_param(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return randi<obj_type>(s.n_rows, s.n_cols, param);
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
randi(const uword n_elem, const distr_param& param = distr_param(), const arma_empty_class junk1 = arma_empty_class(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk2 = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  if(is_Row<obj_type>::value)
    {
    return randi<obj_type>(1, n_elem, param);
    }
  else
    {
    return randi<obj_type>(n_elem, 1, param);
    }
  }



arma_warn_unused
inline
imat
randi(const uword n_rows, const uword n_cols, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randi<imat>(n_rows, n_cols, param);
  }



arma_warn_unused
inline
imat
randi(const SizeMat& s, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randi<imat>(s.n_rows, s.n_cols, param);
  }



arma_warn_unused
inline
ivec
randi(const uword n_elem, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randi<ivec>(n_elem, uword(1), param);
  }



arma_warn_unused
inline
sword
randi(const distr_param& param)
  {
  return as_scalar( randi<ivec>(uword(1), uword(1), param) );
  }



template<typename eT>
arma_warn_unused
inline
typename arma_scalar_only<eT>::result
randi(const distr_param& param)
  {
  return eT( as_scalar( randi< Col<eT> >(uword(1), uword(1), param) ) );
  }



arma_warn_unused
inline
sword
randi()
  {
  return sword( arma_rng::randi<sword>() );
  }



template<typename eT>
arma_warn_unused
inline
typename arma_scalar_only<eT>::result
randi()
  {
  return eT( arma_rng::randi<eT>() );
  }



template<typename cube_type>
arma_warn_unused
inline
cube_type
randi(const uword n_rows, const uword n_cols, const uword n_slices, const distr_param& param = distr_param(), const typename arma_Cube_only<cube_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename cube_type::elem_type eT;
  
  int a = 0;
  int b = arma_rng::randi<eT>::max_val();
  
  param.get_int_vals(a,b);
  
  arma_debug_check( (a > b), "randi(): incorrect distribution parameters; a must be less than b" );
  
  cube_type out(n_rows, n_cols, n_slices, arma_nozeros_indicator());
  
  arma_rng::randi<eT>::fill(out.memptr(), out.n_elem, a, b);
  
  return out;
  }



template<typename cube_type>
arma_warn_unused
inline
cube_type
randi(const SizeCube& s, const distr_param& param = distr_param(), const typename arma_Cube_only<cube_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return randi<cube_type>(s.n_rows, s.n_cols, s.n_slices, param);
  }



arma_warn_unused
inline
icube
randi(const uword n_rows, const uword n_cols, const uword n_slices, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randi<icube>(n_rows, n_cols, n_slices, param);
  }



arma_warn_unused
inline
icube
randi(const SizeCube& s, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randi<icube>(s.n_rows, s.n_cols, s.n_slices, param);
  }



//! @}
