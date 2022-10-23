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


template<typename eT, int SelectionRule, typename OpType>
inline
void
SymEigsShiftSolver<eT, SelectionRule, OpType>::sort_ritzpair()
  {
  arma_extra_debug_sigprint();

  // First transform back the Ritz values, and then sort
  for(uword i = 0; i < this->nev; i++)
    {
    this->ritz_val(i) = eT(1.0) / this->ritz_val(i) + sigma;
    }
  SymEigsSolver<eT, SelectionRule, OpType>::sort_ritzpair();
  }



template<typename eT, int SelectionRule, typename OpType>
inline
SymEigsShiftSolver<eT, SelectionRule, OpType>::SymEigsShiftSolver(const OpType& op_, uword nev_, uword ncv_, const eT sigma_)
  : SymEigsSolver<eT, SelectionRule, OpType>::SymEigsSolver(op_, nev_, ncv_)
  , sigma(sigma_)
  {
  arma_extra_debug_sigprint();
  }


}  // namespace newarp
