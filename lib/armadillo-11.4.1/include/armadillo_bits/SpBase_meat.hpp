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


//! \addtogroup SpBase
//! @{



template<typename elem_type, typename derived>
arma_inline
const derived&
SpBase<elem_type,derived>::get_ref() const
  {
  return static_cast<const derived&>(*this);
  }



template<typename elem_type, typename derived>
arma_inline
bool
SpBase<elem_type,derived>::is_alias(const SpMat<elem_type>& X) const
  {
  return (*this).get_ref().is_alias(X);
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
const SpOp<derived, spop_htrans>
SpBase<elem_type,derived>::t() const
  {
  return SpOp<derived,spop_htrans>( (*this).get_ref() );
  }


template<typename elem_type, typename derived>
inline
arma_warn_unused
const SpOp<derived, spop_htrans>
SpBase<elem_type,derived>::ht() const
  {
  return SpOp<derived, spop_htrans>( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
const SpOp<derived, spop_strans>
SpBase<elem_type,derived>::st() const
  {
  return SpOp<derived, spop_strans>( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
arma_cold
inline
void
SpBase<elem_type,derived>::print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
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
SpBase<elem_type,derived>::print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
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
SpBase<elem_type,derived>::raw_print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
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
SpBase<elem_type,derived>::raw_print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
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
SpBase<elem_type, derived>::print_dense(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );

  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = get_cout_stream().width();
    
    get_cout_stream() << extra_text << '\n';
    
    get_cout_stream().width(orig_width);
    }
  
  arma_ostream::print_dense(get_cout_stream(), tmp.M, true);
  }



template<typename elem_type, typename derived>
arma_cold
inline
void
SpBase<elem_type, derived>::print_dense(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );

  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = user_stream.width();
    
    user_stream << extra_text << '\n';
    
    user_stream.width(orig_width);
    }
  
  arma_ostream::print_dense(user_stream, tmp.M, true);
  }



template<typename elem_type, typename derived>
arma_cold
inline
void
SpBase<elem_type, derived>::raw_print_dense(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );

  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = get_cout_stream().width();
    
    get_cout_stream() << extra_text << '\n';
    
    get_cout_stream().width(orig_width);
    }
  
  arma_ostream::print_dense(get_cout_stream(), tmp.M, false);
  }



template<typename elem_type, typename derived>
arma_cold
inline
void
SpBase<elem_type, derived>::raw_print_dense(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );

  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = user_stream.width();
    
    user_stream << extra_text << '\n';
    
    user_stream.width(orig_width);
    }
  
  arma_ostream::print_dense(user_stream, tmp.M, false);
  }



template<typename elem_type, typename derived>
arma_cold
inline
void
SpBase<elem_type,derived>::brief_print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
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
SpBase<elem_type,derived>::brief_print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
  if(extra_text.length() != 0)
    {
    const std::streamsize orig_width = user_stream.width();
    
    user_stream << extra_text << '\n';
    
    user_stream.width(orig_width);
    }
  
  arma_ostream::brief_print(user_stream, tmp.M);
  }
  


//
// extra functions defined in SpBase_eval_SpMat

template<typename elem_type, typename derived>
inline
arma_warn_unused
const derived&
SpBase_eval_SpMat<elem_type, derived>::eval() const
  {
  arma_extra_debug_sigprint();
  
  return static_cast<const derived&>(*this);
  }



//
// extra functions defined in SpBase_eval_expr

template<typename elem_type, typename derived>
inline
arma_warn_unused
SpMat<elem_type>
SpBase_eval_expr<elem_type, derived>::eval() const
  {
  arma_extra_debug_sigprint();
  
  return SpMat<elem_type>( static_cast<const derived&>(*this) );
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
elem_type
SpBase<elem_type, derived>::min() const
  {
  return spop_min::min( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
elem_type
SpBase<elem_type, derived>::max() const
  {
  return spop_max::max( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
inline
elem_type
SpBase<elem_type, derived>::min(uword& index_of_min_val) const
  {
  const SpProxy<derived> P( (*this).get_ref() );
  
  return spop_min::min_with_index(P, index_of_min_val);
  }



template<typename elem_type, typename derived>
inline
elem_type
SpBase<elem_type, derived>::max(uword& index_of_max_val) const
  {
  const SpProxy<derived> P( (*this).get_ref() );
  
  return spop_max::max_with_index(P, index_of_max_val);
  }



template<typename elem_type, typename derived>
inline
elem_type
SpBase<elem_type, derived>::min(uword& row_of_min_val, uword& col_of_min_val) const
  {
  const SpProxy<derived> P( (*this).get_ref() );
  
  uword index = 0;
  
  const elem_type val = spop_min::min_with_index(P, index);
  
  const uword local_n_rows = P.get_n_rows();
  
  row_of_min_val = index % local_n_rows;
  col_of_min_val = index / local_n_rows;
  
  return val;
  }



template<typename elem_type, typename derived>
inline
elem_type
SpBase<elem_type, derived>::max(uword& row_of_max_val, uword& col_of_max_val) const
  {
  const SpProxy<derived> P( (*this).get_ref() );
  
  uword index = 0;
  
  const elem_type val = spop_max::max_with_index(P, index);
  
  const uword local_n_rows = P.get_n_rows();
  
  row_of_max_val = index % local_n_rows;
  col_of_max_val = index / local_n_rows;
  
  return val;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
uword
SpBase<elem_type,derived>::index_min() const
  {
  const SpProxy<derived> P( (*this).get_ref() );
  
  uword index = 0;
  
  if(P.get_n_elem() == 0)
    {
    arma_debug_check(true, "index_min(): object has no elements");
    }
  else
    {
    spop_min::min_with_index(P, index);
    }
  
  return index;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
uword
SpBase<elem_type,derived>::index_max() const
  {
  const SpProxy<derived> P( (*this).get_ref() );
  
  uword index = 0;
  
  if(P.get_n_elem() == 0)
    {
    arma_debug_check(true, "index_max(): object has no elements");
    }
  else
    {
    spop_max::max_with_index(P, index);
    }
  
  return index;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_symmetric() const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
  return tmp.M.is_symmetric();
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_symmetric(const typename get_pod_type<elem_type>::result tol) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
  return tmp.M.is_symmetric(tol);
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_hermitian() const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
  return tmp.M.is_hermitian();
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_hermitian(const typename get_pod_type<elem_type>::result tol) const
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<derived> tmp( (*this).get_ref() );
  
  return tmp.M.is_hermitian(tol);
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_zero(const typename get_pod_type<elem_type>::result tol) const
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<elem_type>::result T;
  
  arma_debug_check( (tol < T(0)), "is_zero(): parameter 'tol' must be >= 0" );
  
  const SpProxy<derived> P( (*this).get_ref() );
  
  if(P.get_n_elem() == 0)  { return false; }
  
  if(P.get_n_nonzero() == 0)  { return true; }
  
  if(is_SpMat<typename SpProxy<derived>::stored_type>::value)
    {
    const unwrap_spmat<typename SpProxy<derived>::stored_type> U(P.Q);
    
    return arrayops::is_zero(U.M.values, U.M.n_nonzero, tol);
    }
  
  typename SpProxy<derived>::const_iterator_type it     = P.begin();
  typename SpProxy<derived>::const_iterator_type it_end = P.end();
  
  if(is_cx<elem_type>::yes)
    {
    while(it != it_end)
      {
      const elem_type val = (*it);
      
      const T val_real = access::tmp_real(val);
      const T val_imag = access::tmp_imag(val);
      
      if(eop_aux::arma_abs(val_real) > tol)  { return false; }
      if(eop_aux::arma_abs(val_imag) > tol)  { return false; }
      
      ++it;
      }
    }
  else  // not complex
    {
    while(it != it_end)
      {
      if(eop_aux::arma_abs(*it) > tol)  { return false; }
      
      ++it;
      }
    }
  
  return true;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_trimatu() const
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<derived> P( (*this).get_ref() );
  
  if(P.get_n_rows() != P.get_n_cols())  { return false; }
  
  typename SpProxy<derived>::const_iterator_type it     = P.begin();
  typename SpProxy<derived>::const_iterator_type it_end = P.end();
  
  while(it != it_end)
    {
    if(it.row() > it.col())  { return false; }
    ++it;
    }
  
  return true;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_trimatl() const
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<derived> P( (*this).get_ref() );
  
  if(P.get_n_rows() != P.get_n_cols())  { return false; }
  
  typename SpProxy<derived>::const_iterator_type it     = P.begin();
  typename SpProxy<derived>::const_iterator_type it_end = P.end();
  
  while(it != it_end)
    {
    if(it.row() < it.col())  { return false; }
    ++it;
    }
  
  return true;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_diagmat() const
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<derived> P( (*this).get_ref() );
  
  typename SpProxy<derived>::const_iterator_type it     = P.begin();
  typename SpProxy<derived>::const_iterator_type it_end = P.end();
  
  while(it != it_end)
    {
    if(it.row() != it.col())  { return false; }
    ++it;
    }
  
  return true;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_empty() const
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<derived> P( (*this).get_ref() );
  
  return (P.get_n_elem() == uword(0));
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_square() const
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<derived> P( (*this).get_ref() );
  
  return (P.get_n_rows() == P.get_n_cols());
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_vec() const
  {
  arma_extra_debug_sigprint();
  
  if( (SpProxy<derived>::is_row) || (SpProxy<derived>::is_col) || (SpProxy<derived>::is_xvec) )  { return true; }
  
  const SpProxy<derived> P( (*this).get_ref() );
  
  return ( (P.get_n_rows() == uword(1)) || (P.get_n_cols() == uword(1)) );
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_colvec() const
  {
  arma_extra_debug_sigprint();
  
  if(SpProxy<derived>::is_col)  { return true; }
  
  const SpProxy<derived> P( (*this).get_ref() );
  
  return (P.get_n_cols() == uword(1));
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_rowvec() const
  {
  arma_extra_debug_sigprint();
  
  if(SpProxy<derived>::is_row)  { return true; }
  
  const SpProxy<derived> P( (*this).get_ref() );
  
  return (P.get_n_rows() == uword(1));
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::is_finite() const
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<derived> P( (*this).get_ref() );
  
  if(is_SpMat<typename SpProxy<derived>::stored_type>::value)
    {
    const unwrap_spmat<typename SpProxy<derived>::stored_type> U(P.Q);
    
    return U.M.is_finite();
    }
  else
    {
    typename SpProxy<derived>::const_iterator_type it     = P.begin();
    typename SpProxy<derived>::const_iterator_type it_end = P.end();
    
    while(it != it_end)
      {
      if(arma_isfinite(*it) == false)  { return false; }
      ++it;
      }
    }
  
  return true;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::has_inf() const
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<derived> P( (*this).get_ref() );
  
  if(is_SpMat<typename SpProxy<derived>::stored_type>::value)
    {
    const unwrap_spmat<typename SpProxy<derived>::stored_type> U(P.Q);
    
    return U.M.has_inf();
    }
  else
    {
    typename SpProxy<derived>::const_iterator_type it     = P.begin();
    typename SpProxy<derived>::const_iterator_type it_end = P.end();
    
    while(it != it_end)
      {
      if(arma_isinf(*it))  { return true; }
      ++it;
      }
    }
  
  return false;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
SpBase<elem_type,derived>::has_nan() const
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<derived> P( (*this).get_ref() );
  
  if(is_SpMat<typename SpProxy<derived>::stored_type>::value)
    {
    const unwrap_spmat<typename SpProxy<derived>::stored_type> U(P.Q);
    
    return U.M.has_nan();
    }
  else
    {
    typename SpProxy<derived>::const_iterator_type it     = P.begin();
    typename SpProxy<derived>::const_iterator_type it_end = P.end();
    
    while(it != it_end)
      {
      if(arma_isnan(*it))  { return true; }
      ++it;
      }
    }
  
  return false;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
const SpOp<derived,spop_vectorise_col>
SpBase<elem_type, derived>::as_col() const
  {
  return SpOp<derived,spop_vectorise_col>( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
const SpOp<derived,spop_vectorise_row>
SpBase<elem_type, derived>::as_row() const
  {
  return SpOp<derived,spop_vectorise_row>( (*this).get_ref() );
  }



//! @}
