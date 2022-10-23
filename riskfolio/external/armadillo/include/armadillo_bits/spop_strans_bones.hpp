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


//! \addtogroup spop_strans
//! @{


//! simple transpose operation (no complex conjugates) for sparse matrices 

class spop_strans
  {
  public:
  
  template<typename T1>
  struct traits
    {
    static constexpr bool is_row  = T1::is_col;  // deliberately swapped
    static constexpr bool is_col  = T1::is_row;
    static constexpr bool is_xvec = T1::is_xvec;
    };
  
  template<typename eT>
  inline static void apply_noalias(SpMat<eT>& B, const SpMat<eT>& A);
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_strans>& in);
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_htrans>& in);
  };



//! @}
