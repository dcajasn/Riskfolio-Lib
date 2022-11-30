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


//! \addtogroup spglue_relational
//! @{



class spglue_rel_lt
  : public traits_glue_or
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(SpMat<uword>& out, const mtSpGlue<uword, T1, T2, spglue_rel_lt>& X);
  
  template<typename T1, typename T2>
  inline static void apply_noalias(SpMat<uword>& out, const SpProxy<T1>& PA, const SpProxy<T2>& PB);
  };



class spglue_rel_gt
  : public traits_glue_or
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(SpMat<uword>& out, const mtSpGlue<uword, T1, T2, spglue_rel_gt>& X);
  
  template<typename T1, typename T2>
  inline static void apply_noalias(SpMat<uword>& out, const SpProxy<T1>& PA, const SpProxy<T2>& PB);
  };



class spglue_rel_and
  : public traits_glue_or
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(SpMat<uword>& out, const mtSpGlue<uword, T1, T2, spglue_rel_and>& X);
  
  template<typename T1, typename T2>
  inline static void apply_noalias(SpMat<uword>& out, const SpProxy<T1>& PA, const SpProxy<T2>& PB);
  };



class spglue_rel_or
  : public traits_glue_or
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(SpMat<uword>& out, const mtSpGlue<uword, T1, T2, spglue_rel_or>& X);
  
  template<typename T1, typename T2>
  inline static void apply_noalias(SpMat<uword>& out, const SpProxy<T1>& PA, const SpProxy<T2>& PB);
  };



//! @}
