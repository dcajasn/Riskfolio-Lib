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


//! \addtogroup BaseCube
//! @{



template<typename elem_type, typename derived>
arma_inline
const derived&
BaseCube<elem_type,derived>::get_ref() const
  {
  return static_cast<const derived&>(*this);
  }



template<typename elem_type, typename derived>
arma_cold
inline
void
BaseCube<elem_type,derived>::print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<derived> tmp( (*this).get_ref() );
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = get_cout_stream().width();
    
    get_cout_stream() << extra_text << '\n';
    
    get_cout_stream().width(orig_width);
    }
  
  arma_ostream::print(get_cout_stream(), tmp.M, true);
  }



template<typename elem_type, typename derived>
arma_cold
inline
void
BaseCube<elem_type,derived>::print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<derived> tmp( (*this).get_ref() );
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = user_stream.width();
    
    user_stream << extra_text << '\n';
    
    user_stream.width(orig_width);
    }
  
  arma_ostream::print(user_stream, tmp.M, true);
  }
  


template<typename elem_type, typename derived>
arma_cold
inline
void
BaseCube<elem_type,derived>::raw_print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<derived> tmp( (*this).get_ref() );
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = get_cout_stream().width();
    
    get_cout_stream() << extra_text << '\n';
    
    get_cout_stream().width(orig_width);
    }
  
  arma_ostream::print(get_cout_stream(), tmp.M, false);
  }



template<typename elem_type, typename derived>
arma_cold
inline
void
BaseCube<elem_type,derived>::raw_print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<derived> tmp( (*this).get_ref() );
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = user_stream.width();
    
    user_stream << extra_text << '\n';
    
    user_stream.width(orig_width);
    }
  
  arma_ostream::print(user_stream, tmp.M, false);
  }



template<typename elem_type, typename derived>
arma_cold
inline
void
BaseCube<elem_type,derived>::brief_print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<derived> tmp( (*this).get_ref() );
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = get_cout_stream().width();
    
    get_cout_stream() << extra_text << '\n';
    
    get_cout_stream().width(orig_width);
    }
  
  arma_ostream::brief_print(get_cout_stream(), tmp.M);
  }



template<typename elem_type, typename derived>
arma_cold
inline
void
BaseCube<elem_type,derived>::brief_print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_cube<derived> tmp( (*this).get_ref() );
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = user_stream.width();
    
    user_stream << extra_text << '\n';
    
    user_stream.width(orig_width);
    }
  
  arma_ostream::brief_print(user_stream, tmp.M);
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
elem_type
BaseCube<elem_type,derived>::min() const
  {
  return op_min::min( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
elem_type
BaseCube<elem_type,derived>::max() const
  {
  return op_max::max( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
uword
BaseCube<elem_type,derived>::index_min() const
  {
  const ProxyCube<derived> P( (*this).get_ref() );
  
  uword index = 0;
  
  if(P.get_n_elem() == 0)
    {
    arma_debug_check(true, "index_min(): object has no elements");
    }
  else
    {
    op_min::min_with_index(P, index);
    }
  
  return index;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
uword
BaseCube<elem_type,derived>::index_max() const
  {
  const ProxyCube<derived> P( (*this).get_ref() );
  
  uword index = 0;
  
  if(P.get_n_elem() == 0)
    {
    arma_debug_check(true, "index_max(): object has no elements");
    }
  else
    {
    op_max::max_with_index(P, index);
    }
  
  return index;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
BaseCube<elem_type,derived>::is_zero(const typename get_pod_type<elem_type>::result tol) const
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<elem_type>::result T;
  
  arma_debug_check( (tol < T(0)), "is_zero(): parameter 'tol' must be >= 0" );
  
  if(ProxyCube<derived>::use_at || is_Cube<typename ProxyCube<derived>::stored_type>::value)
    {
    const unwrap_cube<derived> U( (*this).get_ref() );
    
    return arrayops::is_zero( U.M.memptr(), U.M.n_elem, tol );
    }
  
  const ProxyCube<derived> P( (*this).get_ref() );
  
  const uword n_elem = P.get_n_elem();
  
  if(n_elem == 0)  { return false; }
  
  const typename ProxyCube<derived>::ea_type Pea = P.get_ea();
  
  if(is_cx<elem_type>::yes)
    {
    for(uword i=0; i<n_elem; ++i)
      {
      const elem_type val = Pea[i];
      
      const T val_real = access::tmp_real(val);
      const T val_imag = access::tmp_imag(val);
      
      if(eop_aux::arma_abs(val_real) > tol)  { return false; }
      if(eop_aux::arma_abs(val_imag) > tol)  { return false; }
      }
    }
  else  // not complex
    {
    for(uword i=0; i < n_elem; ++i)
      {
      if(eop_aux::arma_abs(Pea[i]) > tol)  { return false; }
      }
    }
  
  return true;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
BaseCube<elem_type,derived>::is_empty() const
  {
  arma_extra_debug_sigprint();
  
  const ProxyCube<derived> P( (*this).get_ref() );
  
  return (P.get_n_elem() == uword(0));
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
BaseCube<elem_type,derived>::is_finite() const
  {
  arma_extra_debug_sigprint();
  
  const ProxyCube<derived> P( (*this).get_ref() );
  
  if(is_Cube<typename ProxyCube<derived>::stored_type>::value)
    {
    const unwrap_cube<typename ProxyCube<derived>::stored_type> U(P.Q);
    
    return arrayops::is_finite( U.M.memptr(), U.M.n_elem );
    }
  
  const uword n_r = P.get_n_rows();
  const uword n_c = P.get_n_cols();
  const uword n_s = P.get_n_slices();
  
  for(uword s=0; s<n_s; ++s)
  for(uword c=0; c<n_c; ++c)
  for(uword r=0; r<n_r; ++r)
    {
    if( arma_isfinite(P.at(r,c,s)) == false )  { return false; }
    }
  
  return true;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
BaseCube<elem_type,derived>::has_inf() const
  {
  arma_extra_debug_sigprint();
  
  const ProxyCube<derived> P( (*this).get_ref() );
  
  if(is_Cube<typename ProxyCube<derived>::stored_type>::value)
    {
    const unwrap_cube<typename ProxyCube<derived>::stored_type> U(P.Q);
    
    return arrayops::has_inf( U.M.memptr(), U.M.n_elem );
    }
  
  const uword n_r = P.get_n_rows();
  const uword n_c = P.get_n_cols();
  const uword n_s = P.get_n_slices();
  
  for(uword s=0; s<n_s; ++s)
  for(uword c=0; c<n_c; ++c)
  for(uword r=0; r<n_r; ++r)
    {
    if(arma_isinf(P.at(r,c,s)))  { return true; }
    }
  
  return false;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
BaseCube<elem_type,derived>::has_nan() const
  {
  arma_extra_debug_sigprint();
  
  const ProxyCube<derived> P( (*this).get_ref() );
  
  if(is_Cube<typename ProxyCube<derived>::stored_type>::value)
    {
    const unwrap_cube<typename ProxyCube<derived>::stored_type> U(P.Q);
    
    return arrayops::has_nan( U.M.memptr(), U.M.n_elem );
    }
  
  const uword n_r = P.get_n_rows();
  const uword n_c = P.get_n_cols();
  const uword n_s = P.get_n_slices();
  
  for(uword s=0; s<n_s; ++s)
  for(uword c=0; c<n_c; ++c)
  for(uword r=0; r<n_r; ++r)
    {
    if(arma_isnan(P.at(r,c,s)))  { return true; }
    }
  
  return false;
  }



//
// extra functions defined in BaseCube_eval_Cube

template<typename elem_type, typename derived>
arma_inline
arma_warn_unused
const derived&
BaseCube_eval_Cube<elem_type, derived>::eval() const
  {
  arma_extra_debug_sigprint();
  
  return static_cast<const derived&>(*this);
  }



//
// extra functions defined in BaseCube_eval_expr

template<typename elem_type, typename derived>
inline
arma_warn_unused
Cube<elem_type>
BaseCube_eval_expr<elem_type, derived>::eval() const
  {
  arma_extra_debug_sigprint();
  
  return Cube<elem_type>( static_cast<const derived&>(*this) );
  }



//! @}
