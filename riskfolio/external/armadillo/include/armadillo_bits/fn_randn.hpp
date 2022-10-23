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


//! \addtogroup fn_randn
//! @{



// scalars

arma_warn_unused
inline
double
randn()
  {
  return double(arma_rng::randn<double>());
  }



template<typename eT>
arma_warn_unused
inline
typename arma_real_or_cx_only<eT>::result
randn()
  {
  return eT(arma_rng::randn<eT>());
  }



arma_warn_unused
inline
double
randn(const distr_param& param)
  {
  arma_extra_debug_sigprint();
  
  if(param.state == 0)  { return double(arma_rng::randn<double>()); }
  
  double mu = double(0);
  double sd = double(1);
  
  param.get_double_vals(mu,sd);
  
  arma_debug_check( (sd <= double(0)), "randn(): incorrect distribution parameters; standard deviation must be > 0" );
  
  const double val = double(arma_rng::randn<double>());
  
  return ((val * sd) + mu);
  }



template<typename eT>
arma_warn_unused
inline
typename arma_real_or_cx_only<eT>::result
randn(const distr_param& param)
  {
  arma_extra_debug_sigprint();
  
  if(param.state == 0)  { return eT(arma_rng::randn<eT>()); }
  
  double mu = double(0);
  double sd = double(1);
  
  param.get_double_vals(mu,sd);
  
  arma_debug_check( (sd <= double(0)), "randn(): incorrect distribution parameters; standard deviation must be > 0" );
  
  eT val = eT(0);
  
  arma_rng::randn<eT>::fill(&val, 1, mu, sd);  // using fill() as eT can be complex
  
  return val;
  }



// vectors

arma_warn_unused
inline
vec
randn(const uword n_elem, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  vec out(n_elem, arma_nozeros_indicator());
  
  if(param.state == 0)
    {
    arma_rng::randn<double>::fill(out.memptr(), n_elem);
    }
  else
    {
    double mu = double(0);
    double sd = double(1);
    
    param.get_double_vals(mu,sd);
    
    arma_debug_check( (sd <= double(0)), "randn(): incorrect distribution parameters; standard deviation must be > 0" );
    
    arma_rng::randn<double>::fill(out.memptr(), n_elem, mu, sd);
    }
  
  return out;
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
randn(const uword n_elem, const distr_param& param = distr_param(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename obj_type::elem_type eT;
  
  const uword n_rows = (is_Row<obj_type>::value) ? uword(1) : n_elem;
  const uword n_cols = (is_Row<obj_type>::value) ? n_elem   : uword(1);
  
  obj_type out(n_rows, n_cols, arma_nozeros_indicator());
  
  if(param.state == 0)
    {
    arma_rng::randn<eT>::fill(out.memptr(), out.n_elem);
    }
  else
    {
    double mu = double(0);
    double sd = double(1);
    
    param.get_double_vals(mu,sd);
    
    arma_debug_check( (sd <= double(0)), "randn(): incorrect distribution parameters; standard deviation must be > 0" );
    
    arma_rng::randn<eT>::fill(out.memptr(), out.n_elem, mu, sd);
    }
  
  return out;
  }



// matrices

arma_warn_unused
inline
mat
randn(const uword n_rows, const uword n_cols, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  mat out(n_rows, n_cols, arma_nozeros_indicator());
  
  if(param.state == 0)
    {
    arma_rng::randn<double>::fill(out.memptr(), out.n_elem);
    }
  else
    {
    double mu = double(0);
    double sd = double(1);
    
    param.get_double_vals(mu,sd);
    
    arma_debug_check( (sd <= double(0)), "randn(): incorrect distribution parameters; standard deviation must be > 0" );
    
    arma_rng::randn<double>::fill(out.memptr(), out.n_elem, mu, sd);
    }
  
  return out;
  }



arma_warn_unused
inline
mat
randn(const SizeMat& s, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randn(s.n_rows, s.n_cols, param);
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
randn(const uword n_rows, const uword n_cols, const distr_param& param = distr_param(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename obj_type::elem_type eT;
  
  if(is_Col<obj_type>::value)  { arma_debug_check( (n_cols != 1), "randn(): incompatible size" ); }
  if(is_Row<obj_type>::value)  { arma_debug_check( (n_rows != 1), "randn(): incompatible size" ); }
  
  obj_type out(n_rows, n_cols, arma_nozeros_indicator());
  
  if(param.state == 0)
    {
    arma_rng::randn<eT>::fill(out.memptr(), out.n_elem);
    }
  else
    {
    double mu = double(0);
    double sd = double(1);
    
    param.get_double_vals(mu,sd);
    
    arma_debug_check( (sd <= double(0)), "randn(): incorrect distribution parameters; standard deviation must be > 0" );
    
    arma_rng::randn<eT>::fill(out.memptr(), out.n_elem, mu, sd);
    }
  
  return out;
  }



template<typename obj_type>
arma_warn_unused
inline
obj_type
randn(const SizeMat& s, const distr_param& param = distr_param(), const typename arma_Mat_Col_Row_only<obj_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return randn<obj_type>(s.n_rows, s.n_cols, param);
  }



// cubes


arma_warn_unused
inline
cube
randn(const uword n_rows, const uword n_cols, const uword n_slices, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  cube out(n_rows, n_cols, n_slices, arma_nozeros_indicator());
  
  if(param.state == 0)
    {
    arma_rng::randn<double>::fill(out.memptr(), out.n_elem);
    }
  else
    {
    double mu = double(0);
    double sd = double(1);
    
    param.get_double_vals(mu,sd);
    
    arma_debug_check( (sd <= double(0)), "randn(): incorrect distribution parameters; standard deviation must be > 0" );
    
    arma_rng::randn<double>::fill(out.memptr(), out.n_elem, mu, sd);
    }
  
  return out;
  }



arma_warn_unused
inline
cube
randn(const SizeCube& s, const distr_param& param = distr_param())
  {
  arma_extra_debug_sigprint();
  
  return randn(s.n_rows, s.n_cols, s.n_slices, param);
  }



template<typename cube_type>
arma_warn_unused
inline
cube_type
randn(const uword n_rows, const uword n_cols, const uword n_slices, const distr_param& param = distr_param(), const typename arma_Cube_only<cube_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename cube_type::elem_type eT;
  
  cube_type out(n_rows, n_cols, n_slices, arma_nozeros_indicator());
  
  if(param.state == 0)
    {
    arma_rng::randn<eT>::fill(out.memptr(), out.n_elem);
    }
  else
    {
    double mu = double(0);
    double sd = double(1);
    
    param.get_double_vals(mu,sd);
    
    arma_debug_check( (sd <= double(0)), "randn(): incorrect distribution parameters; standard deviation must be > 0" );
    
    arma_rng::randn<eT>::fill(out.memptr(), out.n_elem, mu, sd);
    }
  
  return out;
  }



template<typename cube_type>
arma_warn_unused
inline
cube_type
randn(const SizeCube& s, const distr_param& param = distr_param(), const typename arma_Cube_only<cube_type>::result* junk = nullptr)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return randn<cube_type>(s.n_rows, s.n_cols, s.n_slices, param);
  }



//! @}
