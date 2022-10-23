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


//! \addtogroup SpToDOp
//! @{



//! Class for storing data required for delayed unary operations on a sparse
//! matrix that produce a dense matrix; the data for storage may include
//! the operand (eg. the matrix to which the operation is to be applied) and the unary operator (eg. inverse).
//! The operand is stored as a reference (which can be optimised away),
//! while the operator is "stored" through the template definition (op_type).
//! The operands can be 'SpMat', 'SpRow', 'SpCol', 'SpOp', and 'SpGlue'.
//! Note that as 'SpGlue' can be one of the operands, more than one matrix can be stored.
//!
//! For example, we could have:
//! SpToDOp< SpGlue< SpMat, SpMat, sp_glue_times >, op_sp_plus >

template<typename T1, typename op_type>
class SpToDOp : public Base< typename T1::elem_type, SpToDOp<T1, op_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  inline explicit SpToDOp(const T1& in_m);
  inline          SpToDOp(const T1& in_m, const elem_type in_aux);
  inline         ~SpToDOp();
  
  arma_aligned const T1&       m;            //!< the operand; must be derived from SpBase
  arma_aligned       elem_type aux;          //!< auxiliary data, using the element type as used by T1
  
  static constexpr bool is_row  = op_type::template traits<T1>::is_row;
  static constexpr bool is_col  = op_type::template traits<T1>::is_col;
  static constexpr bool is_xvec = op_type::template traits<T1>::is_xvec;
  };



//! @}
