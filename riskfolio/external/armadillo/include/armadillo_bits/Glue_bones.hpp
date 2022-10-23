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


//! \addtogroup Glue
//! @{



template<typename T1, typename T2, typename glue_type, bool condition>
struct Glue_traits {};
  

template<typename T1, typename T2, typename glue_type>
struct Glue_traits<T1, T2, glue_type, true>
  {
  static constexpr bool is_row  = glue_type::template traits<T1,T2>::is_row;
  static constexpr bool is_col  = glue_type::template traits<T1,T2>::is_col;
  static constexpr bool is_xvec = glue_type::template traits<T1,T2>::is_xvec;
  };

template<typename T1, typename T2, typename glue_type>
struct Glue_traits<T1, T2, glue_type, false>
  {
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  };


template<typename T1, typename T2, typename glue_type>
class Glue
  : public Base< typename T1::elem_type, Glue<T1, T2, glue_type> >
  , public Glue_traits<T1, T2, glue_type, has_nested_glue_traits<glue_type>::value>
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  inline  Glue(const T1& in_A, const T2& in_B);
  inline  Glue(const T1& in_A, const T2& in_B, const uword in_aux_uword);
  inline ~Glue();
  
  const T1&   A;          //!< first operand;  must be derived from Base
  const T2&   B;          //!< second operand; must be derived from Base
        uword aux_uword;  //!< storage of auxiliary data, uword format
  };



//! @}
