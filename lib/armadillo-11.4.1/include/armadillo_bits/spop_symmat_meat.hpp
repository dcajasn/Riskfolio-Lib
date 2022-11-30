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


//! \addtogroup spop_symmat
//! @{



template<typename T1>
inline
void
spop_symmat::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_symmat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> U(in.m);
  const SpMat<eT>& X   = U.M;
  
  arma_debug_check( (X.n_rows != X.n_cols), "symmatu()/symmatl(): given matrix must be square sized" );
  
  if(X.n_nonzero == uword(0))  { out.zeros(X.n_rows, X.n_cols); return; }
  
  const bool upper = (in.aux_uword_a == 0);
  
  const SpMat<eT> A = (upper) ? trimatu(X) : trimatl(X);  // in this case trimatu() and trimatl() return the same type
  const SpMat<eT> B = A.st();
  
  spglue_merge::symmat_merge(out, A, B);
  }



template<typename T1>
inline
void
spop_symmat_cx::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_symmat_cx>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> U(in.m);
  const SpMat<eT>& X   = U.M;
  
  arma_debug_check( (X.n_rows != X.n_cols), "symmatu()/symmatl(): given matrix must be square sized" );
  
  if(X.n_nonzero == uword(0))  { out.zeros(X.n_rows, X.n_cols); return; }
  
  const bool upper   = (in.aux_uword_a == 0);
  const bool do_conj = (in.aux_uword_b == 1);
  
  const SpMat<eT> A = (upper) ? trimatu(X) : trimatl(X);  // in this case trimatu() and trimatl() return the same type
  
  if(do_conj)
    {
    const SpMat<eT> B = A.t();
    
    spglue_merge::symmat_merge(out, A, B);
    }
  else
    {
    const SpMat<eT> B = A.st();
    
    spglue_merge::symmat_merge(out, A, B);
    }
  }



//! @}
