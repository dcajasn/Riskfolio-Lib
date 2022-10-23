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
struct Base_extra_yes
  {
  inline arma_warn_unused const Op<derived,op_inv_gen_default> i() const;   //!< matrix inverse
  
  inline arma_warn_unused bool is_sympd() const;
  inline arma_warn_unused bool is_sympd(typename get_pod_type<elem_type>::result tol) const;
  };


template<typename elem_type, typename derived>
struct Base_extra_no
  {
  };


template<typename elem_type, typename derived, bool condition>
struct Base_extra {};

template<typename elem_type, typename derived>
struct Base_extra<elem_type, derived, true>  { typedef Base_extra_yes<elem_type, derived> result; };

template<typename elem_type, typename derived>
struct Base_extra<elem_type, derived, false> { typedef Base_extra_no<elem_type, derived>  result; };



template<typename elem_type, typename derived>
struct Base_eval_Mat
  {
  arma_inline arma_warn_unused const derived& eval() const;
  };


template<typename elem_type, typename derived>
struct Base_eval_expr
  {
  inline arma_warn_unused Mat<elem_type> eval() const;   //!< force the immediate evaluation of a delayed expression
  };


template<typename elem_type, typename derived, bool condition>
struct Base_eval {};

template<typename elem_type, typename derived>
struct Base_eval<elem_type, derived, true>  { typedef Base_eval_Mat<elem_type, derived>  result; };

template<typename elem_type, typename derived>
struct Base_eval<elem_type, derived, false> { typedef Base_eval_expr<elem_type, derived> result; };



template<typename derived>
struct Base_trans_cx
  {
  arma_inline arma_warn_unused const Op<derived,op_htrans>  t() const;
  arma_inline arma_warn_unused const Op<derived,op_htrans> ht() const;
  arma_inline arma_warn_unused const Op<derived,op_strans> st() const;  // simple transpose: no complex conjugates
  };


template<typename derived>
struct Base_trans_default
  {
  arma_inline arma_warn_unused const Op<derived,op_htrans>  t() const;
  arma_inline arma_warn_unused const Op<derived,op_htrans> ht() const;
  arma_inline arma_warn_unused const Op<derived,op_htrans> st() const;  // return op_htrans instead of op_strans, as it's handled better by matrix multiplication code
  };


template<typename derived, bool condition>
struct Base_trans {};

template<typename derived>
struct Base_trans<derived, true>  { typedef Base_trans_cx<derived>      result; };

template<typename derived>
struct Base_trans<derived, false> { typedef Base_trans_default<derived> result; };



//! Class for static polymorphism, modelled after the "Curiously Recurring Template Pattern" (CRTP).
//! Used for type-safe downcasting in functions that restrict their input(s) to be classes that are
//! derived from Base (eg. Mat, Op, Glue, diagview, subview).
//! A Base object can be converted to a Mat object by the unwrap class.

template<typename elem_type, typename derived>
struct Base
  : public Base_extra<elem_type, derived, is_supported_blas_type<elem_type>::value>::result
  , public Base_eval<elem_type, derived, is_Mat<derived>::value>::result
  , public Base_trans<derived, is_cx<elem_type>::value>::result
  {
  arma_inline const derived& get_ref() const;
  
  arma_cold inline void print(                           const std::string extra_text = "") const;
  arma_cold inline void print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  arma_cold inline void raw_print(                           const std::string extra_text = "") const;
  arma_cold inline void raw_print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  arma_cold inline void brief_print(                           const std::string extra_text = "") const;
  arma_cold inline void brief_print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  inline arma_warn_unused elem_type min() const;
  inline arma_warn_unused elem_type max() const;
  
  inline elem_type min(uword& index_of_min_val) const;
  inline elem_type max(uword& index_of_max_val) const;
  
  inline elem_type min(uword& row_of_min_val, uword& col_of_min_val) const;
  inline elem_type max(uword& row_of_max_val, uword& col_of_max_val) const;
  
  inline arma_warn_unused uword index_min() const;
  inline arma_warn_unused uword index_max() const;
  
  inline arma_warn_unused bool is_symmetric() const;
  inline arma_warn_unused bool is_symmetric(const typename get_pod_type<elem_type>::result tol) const;
  
  inline arma_warn_unused bool is_hermitian() const;
  inline arma_warn_unused bool is_hermitian(const typename get_pod_type<elem_type>::result tol) const;
  
  inline arma_warn_unused bool is_zero(const typename get_pod_type<elem_type>::result tol = 0) const;
  
  inline arma_warn_unused bool is_trimatu() const;
  inline arma_warn_unused bool is_trimatl() const;
  inline arma_warn_unused bool is_diagmat() const;
  inline arma_warn_unused bool is_empty()   const;
  inline arma_warn_unused bool is_square()  const;
  inline arma_warn_unused bool is_vec()     const;
  inline arma_warn_unused bool is_colvec()  const;
  inline arma_warn_unused bool is_rowvec()  const;
  inline arma_warn_unused bool is_finite()  const;
  inline arma_warn_unused bool has_inf()    const;
  inline arma_warn_unused bool has_nan()    const;
  
  inline arma_warn_unused const Op<derived,op_vectorise_col> as_col() const;
  inline arma_warn_unused const Op<derived,op_vectorise_row> as_row() const;
  };



//! @}
