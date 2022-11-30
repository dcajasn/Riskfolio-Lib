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


//! \addtogroup spglue_max
//! @{



class spglue_max
  : public traits_glue_or
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_max>& X);
  
  template<typename eT, typename T1, typename T2>
  inline static void apply_noalias(SpMat<eT>& out, const SpProxy<T1>& pa, const SpProxy<T2>& pb);
  
  template<typename eT>
  inline static void apply_noalias(SpMat<eT>& out, const SpMat<eT>& A, const SpMat<eT>& B);
  
  template<typename eT, typename T1, typename T2>
  inline static void dense_sparse_max(Mat<eT>& out, const Base<eT,T1>& X, const SpBase<eT,T2>& Y);
  
  template<typename eT>
  inline
  static
  typename enable_if2<is_cx<eT>::no, eT>::result
  elem_max(const eT& a, const eT& b);
  
  template<typename eT>
  inline
  static
  typename enable_if2<is_cx<eT>::yes, eT>::result
  elem_max(const eT& a, const eT& b);
  };



//! @}
