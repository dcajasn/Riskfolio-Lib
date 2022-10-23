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


//! \addtogroup spglue_schur
//! @{



class spglue_schur
  : public traits_glue_or
  {
  public:
  
  template<typename T1, typename T2>
  arma_hot inline static void apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_schur>& X);
  
  template<typename eT, typename T1, typename T2>
  arma_hot inline static void apply_noalias(SpMat<eT>& out, const SpProxy<T1>& pa, const SpProxy<T2>& pb);
  
  template<typename eT>
  arma_hot inline static void apply_noalias(SpMat<eT>& out, const SpMat<eT>& A, const SpMat<eT>& B);
  };



class spglue_schur_misc
  : public traits_glue_or
  {
  public:
  
  template<typename T1, typename T2>
  inline static void dense_schur_sparse(SpMat<typename T1::elem_type>& out, const T1& x, const T2& y);
  };



class spglue_schur_mixed
  : public traits_glue_or
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(SpMat<typename eT_promoter<T1,T2>::eT>& out, const mtSpGlue<typename eT_promoter<T1,T2>::eT, T1, T2, spglue_schur_mixed>& expr);
  
  template<typename T1, typename T2>
  inline static void dense_schur_sparse(SpMat< typename promote_type<typename T1::elem_type, typename T2::elem_type >::result>& out, const T1& X, const T2& Y);
  };



//! @}
