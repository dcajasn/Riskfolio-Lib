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


//! \addtogroup spop_norm
//! @{


class spop_norm
  : public traits_op_default
  {
  public:
  
  template<typename eT> inline static typename get_pod_type<eT>::result mat_norm_1(const SpMat<eT>& X);

  template<typename eT> inline static typename get_pod_type<eT>::result mat_norm_2(const SpMat<eT>& X, const typename arma_real_only<eT>::result* junk = nullptr);
  template<typename eT> inline static typename get_pod_type<eT>::result mat_norm_2(const SpMat<eT>& X, const typename   arma_cx_only<eT>::result* junk = nullptr);

  template<typename eT> inline static typename get_pod_type<eT>::result mat_norm_inf(const SpMat<eT>& X);
  
  template<typename eT> inline static typename get_pod_type<eT>::result vec_norm_k(const eT* mem, const uword N, const uword k);
  };


//! @}
