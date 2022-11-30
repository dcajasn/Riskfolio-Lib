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


namespace newarp
{


//! Define matrix operations on existing matrix objects
template<typename eT>
class SparseGenRealShiftSolve
  {
  private:
  
  #if defined(ARMA_USE_SUPERLU)
    // The following objects are read-only in perform_op()
    mutable superlu_supermatrix_wrangler l;
    mutable superlu_supermatrix_wrangler u;
    mutable superlu_array_wrangler<int>  perm_c;
    mutable superlu_array_wrangler<int>  perm_r;
  #endif
  
  
  public:
  
  bool valid = false;
  
  const uword n_rows;  // number of rows of the underlying matrix
  const uword n_cols;  // number of columns of the underlying matrix
  
  inline SparseGenRealShiftSolve(const SpMat<eT>& mat_obj, const eT shift);
  
  inline void perform_op(eT* x_in, eT* y_out) const;
  };


}  // namespace newarp
