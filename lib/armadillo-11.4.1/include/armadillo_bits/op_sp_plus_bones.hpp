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


//! \addtogroup op_sp_plus
//! @{



// Add a scalar to a sparse matrix; this will return a dense matrix.
class op_sp_plus
  : public traits_op_passthru
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const SpToDOp<T1,op_sp_plus>& in);

  // force apply into an SpMat<>
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpToDOp<T1,op_sp_plus>& in);

  // used for the optimization of sparse % (sparse + scalar)
  template<typename eT, typename T2, typename T3>
  inline static void apply_inside_schur(SpMat<eT>& out, const T2& x, const SpToDOp<T3, op_sp_plus>& y);

  // used for the optimization of sparse / (sparse + scalar)
  template<typename eT, typename T2, typename T3>
  inline static void apply_inside_div(SpMat<eT>& out, const T2& x, const SpToDOp<T3, op_sp_plus>& y);
  };



//! @}
