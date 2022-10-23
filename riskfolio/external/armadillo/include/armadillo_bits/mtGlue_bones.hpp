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


//! \addtogroup mtGlue
//! @{



template<typename out_eT, typename T1, typename T2, typename glue_type>
class mtGlue : public Base< out_eT, mtGlue<out_eT, T1, T2, glue_type> >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  
  static constexpr bool is_row  = glue_type::template traits<T1,T2>::is_row;
  static constexpr bool is_col  = glue_type::template traits<T1,T2>::is_col;
  static constexpr bool is_xvec = glue_type::template traits<T1,T2>::is_xvec;
  
  arma_inline  mtGlue(const T1& in_A, const T2& in_B);
  arma_inline  mtGlue(const T1& in_A, const T2& in_B, const uword in_aux_uword);
  arma_inline ~mtGlue();
  
  arma_aligned const T1&   A;         //!< first operand;  must be derived from Base
  arma_aligned const T2&   B;         //!< second operand; must be derived from Base
  arma_aligned       uword aux_uword; //!< storage of auxiliary data, uword format
  };



//! @}
