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


//! \addtogroup glue_quantile
//! @{


class glue_quantile
  : public traits_glue_default
  {
  public:
  
  template<typename eTa, typename eTb>
  inline static void worker(eTb* out_mem, Col<eTa>& Y, const Mat<eTb>& P);
  
  
  template<typename eTa, typename eTb>
  inline static void apply_noalias(Mat<eTb>& out, const Mat<eTa>& X, const Mat<eTb>& P, const uword dim);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T2::elem_type>& out, const mtGlue<typename T2::elem_type,T1,T2,glue_quantile>& expr);
  };



class glue_quantile_default
  {
  public:
  
  template<typename T1, typename T2>
  struct traits
    {
    static constexpr bool is_row  = T1::is_row;
    static constexpr bool is_col  = T1::is_col;
    static constexpr bool is_xvec = T1::is_xvec;
    };
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T2::elem_type>& out, const mtGlue<typename T2::elem_type,T1,T2,glue_quantile_default>& expr);
  };


//! @}
