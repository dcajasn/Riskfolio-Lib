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


//! \addtogroup op_chi2rnd
//! @{


class op_chi2rnd
  : public traits_op_passthru
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_chi2rnd>& in);
  
  template<typename T1>
  inline static void apply_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& P);
  
  template<typename eT>
  inline static void fill_constant_df(Mat<eT>& out, const eT df);
  };



template<typename eT>
class op_chi2rnd_varying_df
  {
  public:
  
  arma_aligned std::mt19937_64 motor;
  
  inline ~op_chi2rnd_varying_df();
  inline  op_chi2rnd_varying_df();
  
  inline eT operator()(const eT df);
  };


//! @}
