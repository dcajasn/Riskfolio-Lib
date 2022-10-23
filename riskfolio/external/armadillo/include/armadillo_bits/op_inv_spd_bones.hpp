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


//! \addtogroup op_inv_spd
//! @{



class op_inv_spd_default
  : public traits_op_default
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_spd_default>& in);
  
  template<typename T1>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr);
  };



class op_inv_spd_full
  : public traits_op_default
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_spd_full>& in);
  
  template<typename T1, const bool has_user_flags = true>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, const uword flags);
  
  template<typename eT>
  arma_cold inline static bool apply_tiny_2x2(Mat<eT>& X);
  
  template<typename eT>
  arma_cold inline static bool apply_tiny_3x3(Mat<eT>& X);
  
  template<typename eT>
  arma_cold inline static bool apply_tiny_4x4(Mat<eT>& X);
  };



template<typename T>
struct op_inv_spd_state
  {
  T    rcond   = T(0);
  bool is_diag = false;
  };



class op_inv_spd_rcond
  : public traits_op_default
  {
  public:
  
  template<typename T1>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out_inv, op_inv_spd_state<typename T1::pod_type>& out_state, const Base<typename T1::elem_type,T1>& expr);
  };



//! @}
