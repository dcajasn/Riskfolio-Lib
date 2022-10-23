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



//! \addtogroup spop_vectorise
//! @{


class spop_vectorise_col
  : public traits_op_col
  {
  public:
  
  template<typename T1> inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_vectorise_col>& in);
  
  template<typename T1> inline static void apply_direct(SpMat<typename T1::elem_type>& out, const T1& expr);
  };



class spop_vectorise_row
  : public traits_op_row
  {
  public:
  
  template<typename T1> inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_vectorise_row>& in);
  
  template<typename T1> inline static void apply_direct(SpMat<typename T1::elem_type>& out, const T1& expr);
  };



class spop_vectorise_all
  : public traits_op_xvec
  {
  public:
  
  template<typename T1> inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_vectorise_all>& in);
  };



//! @}
