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


//! \addtogroup op_cond
//! @{


class op_cond
  : public traits_op_default
  {
  public:
  
  template<typename T1> static inline typename T1::pod_type apply(const Base<typename T1::elem_type, T1>& X);
  
  template<typename eT> static inline typename get_pod_type<eT>::result apply_diag(const Mat<eT>& A);
  template<typename eT> static inline typename get_pod_type<eT>::result apply_sym (      Mat<eT>& A);
  template<typename eT> static inline typename get_pod_type<eT>::result apply_gen (      Mat<eT>& A);
  };


//! @}
