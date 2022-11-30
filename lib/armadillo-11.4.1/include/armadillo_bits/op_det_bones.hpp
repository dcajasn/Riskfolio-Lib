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


//! \addtogroup op_det
//! @{



class op_det
  : public traits_op_default
  {
  public:
  
  template<const uword row, const uword col>
  struct pos
    {
    static constexpr uword n2 = row + col*2;
    static constexpr uword n3 = row + col*3;
    static constexpr uword n4 = row + col*4;
    };
  
  template<typename T1>
  inline static bool apply_direct(typename T1::elem_type& out_val, const Base<typename T1::elem_type,T1>& expr);
  
  template<typename T1>
  inline static typename T1::elem_type apply_diagmat(const Base<typename T1::elem_type,T1>& expr);
  
  template<typename T1>
  inline static typename T1::elem_type apply_trimat(const Base<typename T1::elem_type,T1>& expr);
  
  template<typename eT>
  arma_cold inline static eT apply_tiny_2x2(const Mat<eT>& X);
  
  template<typename eT>
  arma_cold inline static eT apply_tiny_3x3(const Mat<eT>& X);
  
  template<typename eT>
  arma_cold inline static eT apply_tiny_4x4(const Mat<eT>& X);
  };



//! @}
