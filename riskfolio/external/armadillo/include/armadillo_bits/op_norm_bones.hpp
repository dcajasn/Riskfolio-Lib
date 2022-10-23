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


//! \addtogroup op_norm
//! @{


class op_norm
  : public traits_op_default
  {
  public:
  
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_1(const Proxy<T1>& P, const typename  arma_not_cx<typename T1::elem_type>::result* junk = nullptr);
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_1(const Proxy<T1>& P, const typename arma_cx_only<typename T1::elem_type>::result* junk = nullptr);
  template<typename eT> arma_hot inline static eT                    vec_norm_1_direct_std(const Mat<eT>& X);
  template<typename eT> arma_hot inline static eT                    vec_norm_1_direct_mem(const uword N, const eT* A);
  
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_2(const Proxy<T1>& P, const typename  arma_not_cx<typename T1::elem_type>::result* junk = nullptr);
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_2(const Proxy<T1>& P, const typename arma_cx_only<typename T1::elem_type>::result* junk = nullptr);
  template<typename eT> arma_hot inline static eT                    vec_norm_2_direct_std(const Mat<eT>& X);
  template<typename eT> arma_hot inline static eT                    vec_norm_2_direct_mem(const uword N, const eT* A);
  template<typename eT> arma_hot inline static eT                    vec_norm_2_direct_robust(const Mat<eT>& X);
  
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_k(const Proxy<T1>& P, const int k);
  
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_max(const Proxy<T1>& P);
  template<typename T1> arma_hot inline static typename T1::pod_type vec_norm_min(const Proxy<T1>& P);
  
  template<typename eT> inline static typename get_pod_type<eT>::result mat_norm_1(const Mat<eT>& X);
  template<typename eT> inline static typename get_pod_type<eT>::result mat_norm_2(const Mat<eT>& X);
  
  template<typename eT> inline static typename get_pod_type<eT>::result mat_norm_inf(const Mat<eT>& X);
  };



//! @}
