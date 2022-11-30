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


//! \addtogroup Base
//! @{



template<typename elem_type, typename derived>
arma_inline
const derived&
Base<elem_type,derived>::get_ref() const
  {
  return static_cast<const derived&>(*this);
  }



template<typename elem_type, typename derived>
arma_cold
inline
void
Base<elem_type,derived>::print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<derived> tmp( (*this).get_ref() );
  
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
Base<elem_type,derived>::print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<derived> tmp( (*this).get_ref() );
  
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
Base<elem_type,derived>::raw_print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<derived> tmp( (*this).get_ref() );
  
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
Base<elem_type,derived>::raw_print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<derived> tmp( (*this).get_ref() );
  
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
Base<elem_type,derived>::brief_print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<derived> tmp( (*this).get_ref() );
  
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
Base<elem_type,derived>::brief_print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<derived> tmp( (*this).get_ref() );
  
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
Base<elem_type,derived>::min() const
  {
  return op_min::min( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
elem_type
Base<elem_type,derived>::max() const
  {
  return op_max::max( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
inline
elem_type
Base<elem_type,derived>::min(uword& index_of_min_val) const
  {
  const Proxy<derived> P( (*this).get_ref() );
  
  return op_min::min_with_index(P, index_of_min_val);
  }



template<typename elem_type, typename derived>
inline
elem_type
Base<elem_type,derived>::max(uword& index_of_max_val) const
  {
  const Proxy<derived> P( (*this).get_ref() );
  
  return op_max::max_with_index(P, index_of_max_val);
  }



template<typename elem_type, typename derived>
inline
elem_type
Base<elem_type,derived>::min(uword& row_of_min_val, uword& col_of_min_val) const
  {
  const Proxy<derived> P( (*this).get_ref() );
  
  uword index = 0;
  
  const elem_type val = op_min::min_with_index(P, index);
  
  const uword local_n_rows = P.get_n_rows();
  
  row_of_min_val = index % local_n_rows;
  col_of_min_val = index / local_n_rows;
  
  return val;
  }



template<typename elem_type, typename derived>
inline
elem_type
Base<elem_type,derived>::max(uword& row_of_max_val, uword& col_of_max_val) const
  {
  const Proxy<derived> P( (*this).get_ref() );
  
  uword index = 0;
  
  const elem_type val = op_max::max_with_index(P, index);
  
  const uword local_n_rows = P.get_n_rows();
  
  row_of_max_val = index % local_n_rows;
  col_of_max_val = index / local_n_rows;
  
  return val;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
uword
Base<elem_type,derived>::index_min() const
  {
  const Proxy<derived> P( (*this).get_ref() );
  
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
Base<elem_type,derived>::index_max() const
  {
  const Proxy<derived> P( (*this).get_ref() );
  
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
Base<elem_type,derived>::is_symmetric() const
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<derived> U( (*this).get_ref() );
  
  const Mat<elem_type>& A = U.M;
  
  if(A.n_rows != A.n_cols)  { return false; }
  if(A.n_elem <= 1       )  { return true;  }
  
  const uword N   = A.n_rows;
  const uword Nm1 = N-1;
  
  const elem_type* A_col = A.memptr();
  
  for(uword j=0; j < Nm1; ++j)
    {
    const uword jp1 = j+1;
    
    const elem_type* A_row = &(A.at(j,jp1));
    
    for(uword i=jp1; i < N; ++i)
      {
      if(A_col[i] != (*A_row))  { return false; }
      
      A_row += N;
      }
    
    A_col += N;
    }
  
  return true;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::is_symmetric(const typename get_pod_type<elem_type>::result tol) const
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<elem_type>::result T;
  
  if(tol == T(0))  { return (*this).is_symmetric(); }
  
  arma_debug_check( (tol < T(0)), "is_symmetric(): parameter 'tol' must be >= 0" );
  
  const quasi_unwrap<derived> U( (*this).get_ref() );
  
  const Mat<elem_type>& A = U.M;
  
  if(A.n_rows != A.n_cols)  { return false; }
  if(A.n_elem <= 1       )  { return true;  }
  
  const T norm_A = as_scalar( arma::max(sum(abs(A), 1), 0) );
  
  if(norm_A == T(0))  { return true; }
  
  const T norm_A_Ast = as_scalar( arma::max(sum(abs(A - A.st()), 1), 0) );
  
  return ( (norm_A_Ast / norm_A) <= tol );
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::is_hermitian() const
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<elem_type>::result T;
  
  const quasi_unwrap<derived> U( (*this).get_ref() );
  
  const Mat<elem_type>& A = U.M;
  
  if(A.n_rows != A.n_cols)  { return false; }
  if(A.n_elem == 0       )  { return true;  }
  
  const uword N = A.n_rows;
  
  const elem_type* A_col = A.memptr();
  
  for(uword j=0; j < N; ++j)
    {
    if( access::tmp_imag(A_col[j]) != T(0) )  { return false; }
    
    A_col += N;
    }
  
  A_col = A.memptr();
  
  const uword Nm1 = N-1;
  
  for(uword j=0; j < Nm1; ++j)
    {
    const uword jp1 = j+1;
    
    const elem_type* A_row = &(A.at(j,jp1));
    
    for(uword i=jp1; i < N; ++i)
      {
      if(A_col[i] != access::alt_conj(*A_row))  { return false; }
      
      A_row += N;
      }
    
    A_col += N;
    }
  
  return true;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::is_hermitian(const typename get_pod_type<elem_type>::result tol) const
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<elem_type>::result T;
  
  if(tol == T(0))  { return (*this).is_hermitian(); }
  
  arma_debug_check( (tol < T(0)), "is_hermitian(): parameter 'tol' must be >= 0" );
  
  const quasi_unwrap<derived> U( (*this).get_ref() );
  
  const Mat<elem_type>& A = U.M;
  
  if(A.n_rows != A.n_cols)  { return false; }
  if(A.n_elem == 0       )  { return true;  }
  
  const T norm_A = as_scalar( arma::max(sum(abs(A), 1), 0) );
  
  if(norm_A == T(0))  { return true; }
  
  const T norm_A_At = as_scalar( arma::max(sum(abs(A - A.t()), 1), 0) );
  
  return ( (norm_A_At / norm_A) <= tol );
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::is_zero(const typename get_pod_type<elem_type>::result tol) const
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<elem_type>::result T;
  
  arma_debug_check( (tol < T(0)), "is_zero(): parameter 'tol' must be >= 0" );
  
  if(Proxy<derived>::use_at || is_Mat<typename Proxy<derived>::stored_type>::value)
    {
    const quasi_unwrap<derived> U( (*this).get_ref() );
    
    return arrayops::is_zero( U.M.memptr(), U.M.n_elem, tol );
    }
  
  const Proxy<derived> P( (*this).get_ref() );
  
  const uword n_elem = P.get_n_elem();
  
  if(n_elem == 0)  { return false; }
  
  const typename Proxy<derived>::ea_type Pea = P.get_ea();
  
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
    for(uword i=0; i<n_elem; ++i)
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
Base<elem_type,derived>::is_trimatu() const
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<derived> U( (*this).get_ref() );
  
  if(U.M.n_rows != U.M.n_cols)  { return false; }
  
  if(U.M.n_elem <= 1)  { return true; }
  
  return trimat_helper::is_triu(U.M);
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::is_trimatl() const
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<derived> U( (*this).get_ref() );
  
  if(U.M.n_rows != U.M.n_cols)  { return false; }
  
  if(U.M.n_elem <= 1)  { return true; }
  
  return trimat_helper::is_tril(U.M);
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::is_diagmat() const
  {
  arma_extra_debug_sigprint();
  
  const quasi_unwrap<derived> U( (*this).get_ref() );
  
  const Mat<elem_type>& A = U.M;
  
  if(A.n_elem <= 1)  { return true; }
  
  // NOTE: we're NOT assuming the matrix has a square size
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const elem_type* A_mem = A.memptr();
  
  if(A_mem[1] != elem_type(0))  { return false; }
  
  // if we got to this point, do a thorough check
  
  for(uword A_col=0; A_col < A_n_cols; ++A_col)
    {
    for(uword A_row=0; A_row < A_n_rows; ++A_row)
      {
      if( (A_mem[A_row] != elem_type(0)) && (A_row != A_col) )  { return false; }
      }
    
    A_mem += A_n_rows;
    }
  
  return true;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::is_empty() const
  {
  arma_extra_debug_sigprint();
  
  const Proxy<derived> P( (*this).get_ref() );
  
  return (P.get_n_elem() == uword(0));
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::is_square() const
  {
  arma_extra_debug_sigprint();
  
  const Proxy<derived> P( (*this).get_ref() );
  
  return (P.get_n_rows() == P.get_n_cols());
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::is_vec() const
  {
  arma_extra_debug_sigprint();
  
  if( (Proxy<derived>::is_row) || (Proxy<derived>::is_col) || (Proxy<derived>::is_xvec) )  { return true; }
  
  const Proxy<derived> P( (*this).get_ref() );
  
  return ( (P.get_n_rows() == uword(1)) || (P.get_n_cols() == uword(1)) );
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::is_colvec() const
  {
  arma_extra_debug_sigprint();
  
  if(Proxy<derived>::is_col)  { return true; }
  
  const Proxy<derived> P( (*this).get_ref() );
  
  return (P.get_n_cols() == uword(1));
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::is_rowvec() const
  {
  arma_extra_debug_sigprint();
  
  if(Proxy<derived>::is_row)  { return true; }
  
  const Proxy<derived> P( (*this).get_ref() );
  
  return (P.get_n_rows() == uword(1));
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::is_finite() const
  {
  arma_extra_debug_sigprint();
  
  const Proxy<derived> P( (*this).get_ref() );
  
  if(is_Mat<typename Proxy<derived>::stored_type>::value)
    {
    const quasi_unwrap<typename Proxy<derived>::stored_type> U(P.Q);
    
    return arrayops::is_finite( U.M.memptr(), U.M.n_elem );
    }
  
  if(Proxy<derived>::use_at == false)
    {
    const typename Proxy<derived>::ea_type Pea = P.get_ea();
    
    const uword n_elem = P.get_n_elem();
    
    for(uword i=0; i<n_elem; ++i)
      {
      if(arma_isfinite(Pea[i]) == false)  { return false; }
      }
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    for(uword col=0; col<n_cols; ++col)
    for(uword row=0; row<n_rows; ++row)
      {
      if(arma_isfinite(P.at(row,col)) == false)  { return false; }
      }
    }
  
  return true;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::has_inf() const
  {
  arma_extra_debug_sigprint();
  
  const Proxy<derived> P( (*this).get_ref() );
  
  if(is_Mat<typename Proxy<derived>::stored_type>::value)
    {
    const quasi_unwrap<typename Proxy<derived>::stored_type> U(P.Q);
    
    return arrayops::has_inf( U.M.memptr(), U.M.n_elem );
    }
  
  if(Proxy<derived>::use_at == false)
    {
    const typename Proxy<derived>::ea_type Pea = P.get_ea();
    
    const uword n_elem = P.get_n_elem();
    
    for(uword i=0; i<n_elem; ++i)
      {
      if(arma_isinf(Pea[i]))  { return true; }
      }
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    for(uword col=0; col<n_cols; ++col)
    for(uword row=0; row<n_rows; ++row)
      {
      if(arma_isinf(P.at(row,col)))  { return true; }
      }
    }
  
  return false;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base<elem_type,derived>::has_nan() const
  {
  arma_extra_debug_sigprint();
  
  const Proxy<derived> P( (*this).get_ref() );
  
  if(is_Mat<typename Proxy<derived>::stored_type>::value)
    {
    const quasi_unwrap<typename Proxy<derived>::stored_type> U(P.Q);
    
    return arrayops::has_nan( U.M.memptr(), U.M.n_elem );
    }
  
  if(Proxy<derived>::use_at == false)
    {
    const typename Proxy<derived>::ea_type Pea = P.get_ea();
    
    const uword n_elem = P.get_n_elem();
    
    for(uword i=0; i<n_elem; ++i)
      {
      if(arma_isnan(Pea[i]))  { return true; }
      }
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    for(uword col=0; col<n_cols; ++col)
    for(uword row=0; row<n_rows; ++row)
      {
      if(arma_isnan(P.at(row,col)))  { return true; }
      }
    }
  
  return false;
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
const Op<derived,op_vectorise_col>
Base<elem_type, derived>::as_col() const
  {
  return Op<derived,op_vectorise_col>( (*this).get_ref() );
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
const Op<derived,op_vectorise_row>
Base<elem_type, derived>::as_row() const
  {
  return Op<derived,op_vectorise_row>( (*this).get_ref() );
  }



//
// extra functions defined in Base_extra_yes

template<typename elem_type, typename derived>
inline
arma_warn_unused
const Op<derived,op_inv_gen_default>
Base_extra_yes<elem_type, derived>::i() const
  {
  return Op<derived,op_inv_gen_default>(static_cast<const derived&>(*this));
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base_extra_yes<elem_type,derived>::is_sympd() const
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<elem_type>::result T;
  
  Mat<elem_type> X = static_cast<const derived&>(*this);
  
  // default value for tol
  const T tol = T(100) * std::numeric_limits<T>::epsilon() * norm(X, "fro");
  
  if(X.is_hermitian(tol) == false)  { return false; }
  
  if(X.is_empty())  { return false; }
  
  X.diag() -= elem_type(tol);
  
  return auxlib::chol_simple(X);
  }



template<typename elem_type, typename derived>
inline
arma_warn_unused
bool
Base_extra_yes<elem_type,derived>::is_sympd(typename get_pod_type<elem_type>::result tol) const
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<elem_type>::result T;
  
  arma_debug_check( (tol < T(0)), "is_sympd(): parameter 'tol' must be >= 0" );
  
  Mat<elem_type> X = static_cast<const derived&>(*this);
  
  if(X.is_hermitian(tol) == false)  { return false; }
  
  if(X.is_empty())  { return false; }
  
  X.diag() -= elem_type(tol);
  
  return auxlib::chol_simple(X);
  }



//
// extra functions defined in Base_eval_Mat

template<typename elem_type, typename derived>
arma_inline
arma_warn_unused
const derived&
Base_eval_Mat<elem_type, derived>::eval() const
  {
  arma_extra_debug_sigprint();
  
  return static_cast<const derived&>(*this);
  }



//
// extra functions defined in Base_eval_expr

template<typename elem_type, typename derived>
inline
arma_warn_unused
Mat<elem_type>
Base_eval_expr<elem_type, derived>::eval() const
  {
  arma_extra_debug_sigprint();
  
  return Mat<elem_type>( static_cast<const derived&>(*this) );
  }



//
// extra functions defined in Base_trans_cx

template<typename derived>
arma_inline
arma_warn_unused
const Op<derived,op_htrans>
Base_trans_cx<derived>::t() const
  {
  return Op<derived,op_htrans>( static_cast<const derived&>(*this) );
  }



template<typename derived>
arma_inline
arma_warn_unused
const Op<derived,op_htrans>
Base_trans_cx<derived>::ht() const
  {
  return Op<derived,op_htrans>( static_cast<const derived&>(*this) );
  }



template<typename derived>
arma_inline
arma_warn_unused
const Op<derived,op_strans>
Base_trans_cx<derived>::st() const
  {
  return Op<derived,op_strans>( static_cast<const derived&>(*this) );
  }



//
// extra functions defined in Base_trans_default

template<typename derived>
arma_inline
arma_warn_unused
const Op<derived,op_htrans>
Base_trans_default<derived>::t() const
  {
  return Op<derived,op_htrans>( static_cast<const derived&>(*this) );
  }



template<typename derived>
arma_inline
arma_warn_unused
const Op<derived,op_htrans>
Base_trans_default<derived>::ht() const
  {
  return Op<derived,op_htrans>( static_cast<const derived&>(*this) );
  }



template<typename derived>
arma_inline
arma_warn_unused
const Op<derived,op_htrans>
Base_trans_default<derived>::st() const
  {
  return Op<derived,op_htrans>( static_cast<const derived&>(*this) );
  }



//! @}
