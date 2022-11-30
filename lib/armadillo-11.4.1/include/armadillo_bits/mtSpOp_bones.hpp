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


//! \addtogroup mtSpOp
//! @{

// Class for delayed multi-type sparse operations.  These are operations where
// the resulting type is different than the stored type.



template<typename out_eT, typename T1, typename spop_type>
class mtSpOp : public SpBase< out_eT, mtSpOp<out_eT, T1, spop_type> >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  
  typedef typename T1::elem_type                in_eT;
  
  static constexpr bool is_row  = spop_type::template traits<T1>::is_row;
  static constexpr bool is_col  = spop_type::template traits<T1>::is_col;
  static constexpr bool is_xvec = spop_type::template traits<T1>::is_xvec;
  
  inline explicit  mtSpOp(const T1& in_m);
  inline           mtSpOp(const T1& in_m, const uword aux_uword_a, const uword aux_uword_b);
  inline           mtSpOp(const char junk, const T1& in_m, const out_eT in_aux);
  inline          ~mtSpOp();
  
  template<typename eT2>
  arma_inline bool is_alias(const SpMat<eT2>& X) const;
  
  arma_aligned const T1&    m;            //!< the operand; must be derived from SpBase
  arma_aligned       out_eT aux_out_eT;   //!< auxiliary data, using the element type as specified by the out_eT template parameter
  arma_aligned       uword  aux_uword_a;
  arma_aligned       uword  aux_uword_b;
  };



//! @}
