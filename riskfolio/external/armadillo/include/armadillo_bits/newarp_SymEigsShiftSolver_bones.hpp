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


//! This class implements the eigen solver for real symmetric matrices in the shift-and-invert mode.
template<typename eT, int SelectionRule, typename OpType>
class SymEigsShiftSolver : public SymEigsSolver<eT, SelectionRule, OpType>
  {
  private:

  const eT sigma;

  // Sort the first nev Ritz pairs in ascending algebraic order
  // This is used to return the final results
  void sort_ritzpair();


  public:

  //! Constructor to create a solver object.
  inline SymEigsShiftSolver(const OpType& op_, uword nev_, uword ncv_, const eT sigma_);
  };


}  // namespace newarp
