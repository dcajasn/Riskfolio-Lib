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


//! \addtogroup SpOp
//! @{



template<typename T1, typename op_type>
class SpOp : public SpBase< typename T1::elem_type, SpOp<T1, op_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static constexpr bool is_row  = op_type::template traits<T1>::is_row;
  static constexpr bool is_col  = op_type::template traits<T1>::is_col;
  static constexpr bool is_xvec = op_type::template traits<T1>::is_xvec;
  
  inline explicit SpOp(const T1& in_m);
  inline          SpOp(const T1& in_m, const elem_type in_aux);
  inline          SpOp(const T1& in_m, const uword     in_aux_uword_a, const uword in_aux_uword_b);
  inline         ~SpOp();
  
  arma_inline bool is_alias(const SpMat<elem_type>& X) const;
  
  arma_aligned const T1&       m;            //!< the operand; must be derived from SpBase
  arma_aligned       elem_type aux;          //!< auxiliary data, using the element type as used by T1
  arma_aligned       uword     aux_uword_a;  //!< auxiliary data, uword format
  arma_aligned       uword     aux_uword_b;  //!< auxiliary data, uword format
  };



//! @}
