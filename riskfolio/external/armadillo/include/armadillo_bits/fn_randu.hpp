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


//! \addtogroup fn_randu
//! @{



// scalars

arma_warn_unused
inline
double
randu()
  {
  return double(arma_rng::randu<double>());
  }



template<typename eT>
arma_warn_unused
inline
typename arma_real_or_cx_only<eT>::result
randu()
  {
  return eT(arma_rng::randu<eT>());
  }



arma_warn_unused
inline
double
randu(const distr_param& param)
  {
  arma_extra_debug_sigprint();
  
  if(param.state == 0)  { return double(arma_rng::randu<double>()); }
  
  double a = double(0);
  double b = double(1);
  
  param.get_double_vals(a,b);
  
  arma_debug_check( (a >= b), "randu(): incorrect distribution parameters; a must be less than b" );
  
  const double val = double(arma_rng::randu<double>());
  
  return ((val * (b - a)) + a);
  }



template<typename eT>
arma_warn_unused
inline
typename arma_real_or_cx_only<eT>::result
randu(const distr_param& param)
  {
  arma_extra_debug_sigprint();
  
  if(param.state == 0)  { return eT(arma_rng::randu<eT>()); }
  
  double a = double(0);
  double b = double(1);
  
  param.get_double_vals(a,b);
  
  arma_debug_check( (a >= b), "randu(): incorrect distribution parameters; a must be less than b" );
  
  eT val = eT(0);
  
  arma_rng::randu<eT>::fill(&val, 1, a, b);  // using fill() as eT can be complex
  
  return val;
  }



// vectors

arma_warn_unused
inline
vec
randu(const uword n_elem, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  vec out(n_elem, arma_nozeros_indicator());
  
  if(param.state == 0)
    {
    arma_rng::randu<double>::fill(out.memptr(), n_elem);
    }
  else
    {
    double a = double(0);
    double b = double(1);
    
    param.get_double_vals(a,b);
    
    arma_debug_check( (a >= b), "randu(): incorrect distribution parameters; a must be less than b" );
    
    arma_rng::randu<double>::fill(out.memptr(), n_elem, a, b);
    }
  
  return out;
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
randu(const uword n_elem, const distr_param& param = distr_param(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename obj_type::elem_type eT;
  
  const uword n_rows = (is_Row<obj_type>::value) ? uword(1) : n_elem;
  const uword n_cols = (is_Row<obj_type>::value) ? n_elem   : uword(1);
  
  obj_type out(n_rows, n_cols, arma_nozeros_indicator());
  
  if(param.state == 0)
    {
    arma_rng::randu<eT>::fill(out.memptr(), out.n_elem);
    }
  else
    {
    double a = double(0);
    double b = double(1);
    
    param.get_double_vals(a,b);
    
    arma_debug_check( (a >= b), "randu(): incorrect distribution parameters; a must be less than b" );
    
    arma_rng::randu<eT>::fill(out.memptr(), out.n_elem, a, b);
    }
  
  return out;
  }



// matrices

arma_warn_unused
inline
mat
randu(const uword n_rows, const uword n_cols, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  mat out(n_rows, n_cols, arma_nozeros_indicator());
  
  if(param.state == 0)
    {
    arma_rng::randu<double>::fill(out.memptr(), out.n_elem);
    }
  else
    {
    double a = double(0);
    double b = double(1);
    
    param.get_double_vals(a,b);
    
    arma_debug_check( (a >= b), "randu(): incorrect distribution parameters; a must be less than b" );
    
    arma_rng::randu<double>::fill(out.memptr(), out.n_elem, a, b);
    }
  
  return out;
  }



arma_warn_unused
inline
mat
randu(const SizeMat& s, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randu(s.n_rows, s.n_cols, param);
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
randu(const uword n_rows, const uword n_cols, const distr_param& param = distr_param(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename obj_type::elem_type eT;
  
  if(is_Col<obj_type>::value)  { arma_debug_check( (n_cols != 1), "randu(): incompatible size" ); }
  if(is_Row<obj_type>::value)  { arma_debug_check( (n_rows != 1), "randu(): incompatible size" ); }
  
  obj_type out(n_rows, n_cols, arma_nozeros_indicator());
  
  if(param.state == 0)
    {
    arma_rng::randu<eT>::fill(out.memptr(), out.n_elem);
    }
  else
    {
    double a = double(0);
    double b = double(1);
    
    param.get_double_vals(a,b);
    
    arma_debug_check( (a >= b), "randu(): incorrect distribution parameters; a must be less than b" );
    
    arma_rng::randu<eT>::fill(out.memptr(), out.n_elem, a, b);
    }
  
  return out;
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
randu(const SizeMat& s, const distr_param& param = distr_param(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return randu<obj_type>(s.n_rows, s.n_cols, param);
  }



// cubes


arma_warn_unused
inline
cube
randu(const uword n_rows, const uword n_cols, const uword n_slices, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  cube out(n_rows, n_cols, n_slices, arma_nozeros_indicator());
  
  if(param.state == 0)
    {
    arma_rng::randu<double>::fill(out.memptr(), out.n_elem);
    }
  else
    {
    double a = double(0);
    double b = double(1);
    
    param.get_double_vals(a,b);
    
    arma_debug_check( (a >= b), "randu(): incorrect distribution parameters; a must be less than b" );
    
    arma_rng::randu<double>::fill(out.memptr(), out.n_elem, a, b);
    }
  
  return out;
  }



arma_warn_unused
inline
cube
randu(const SizeCube& s, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randu(s.n_rows, s.n_cols, s.n_slices, param);
  }



template<typename cube_type>
arma_warn_unused
inline
cube_type
randu(const uword n_rows, const uword n_cols, const uword n_slices, const distr_param& param = distr_param(), const typename arma_Cube_only<cube_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename cube_type::elem_type eT;
  
  cube_type out(n_rows, n_cols, n_slices, arma_nozeros_indicator());
  
  if(param.state == 0)
    {
    arma_rng::randu<eT>::fill(out.memptr(), out.n_elem);
    }
  else
    {
    double a = double(0);
    double b = double(1);
    
    param.get_double_vals(a,b);
    
    arma_debug_check( (a >= b), "randu(): incorrect distribution parameters; a must be less than b" );
    
    arma_rng::randu<eT>::fill(out.memptr(), out.n_elem, a, b);
    }
  
  return out;
  }



template<typename cube_type>
arma_warn_unused
inline
cube_type
randu(const SizeCube& s, const distr_param& param = distr_param(), const typename arma_Cube_only<cube_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return randu<cube_type>(s.n_rows, s.n_cols, s.n_slices, param);
  }



//! @}
