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



template<typename eT>
inline
void
spop_strans::apply_noalias(SpMat<eT>& B, const SpMat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  B.reserve(A.n_cols, A.n_rows, A.n_nonzero);  // deliberately swapped
  
  if(A.n_nonzero == 0)  { return; }
  
  // This follows the TRANSP algorithm described in
  // 'Sparse Matrix Multiplication Package (SMMP)'
  // (R.E. Bank and C.C. Douglas, 2001)
  
  const uword m = A.n_rows;
  const uword n = A.n_cols;

  const eT* a = A.values;
        eT* b = access::rwp(B.values);

  const uword* ia = A.col_ptrs;     
  const uword* ja = A.row_indices;
  
        uword* ib = access::rwp(B.col_ptrs);
        uword* jb = access::rwp(B.row_indices);
  
  // // ib is already zeroed, as B is freshly constructed
  // 
  // for(uword i=0; i < (m+1); ++i)
  //   {
  //   ib[i] = 0;
  //   }
  
  for(uword i=0; i < n; ++i)
    {
    for(uword j = ia[i]; j < ia[i+1]; ++j)
      {
      ib[ ja[j] + 1 ]++;
      }
    }
  
  for(uword i=0; i < m; ++i)
    {
    ib[i+1] += ib[i];
    }
 
  for(uword i=0; i < n; ++i)
    {
    for(uword j = ia[i]; j < ia[i+1]; ++j)
      {
      const uword jj = ja[j];
      
      const uword ib_jj = ib[jj];
      
      jb[ib_jj] = i;
      
      b[ib_jj] = a[j];
      
      ib[jj]++;
      }
    }
  
  for(uword i = m-1; i >= 1; --i)
    {
    ib[i] = ib[i-1];
    }
  
  ib[0] = 0;
  }



template<typename T1>
inline
void
spop_strans::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_strans>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> U(in.m);
  
  if(U.is_alias(out))
    {
    SpMat<eT> tmp;
    
    spop_strans::apply_noalias(tmp, U.M);
    
    out.steal_mem(tmp);
    }
  else
    {
    spop_strans::apply_noalias(out, U.M);
    }
  }



//! for transpose of non-complex matrices, redirected from spop_htrans::apply()
template<typename T1>
inline
void
spop_strans::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_htrans>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> U(in.m);
  
  if(U.is_alias(out))
    {
    SpMat<eT> tmp;
    
    spop_strans::apply_noalias(tmp, U.M);
    
    out.steal_mem(tmp);
    }
  else
    {
    spop_strans::apply_noalias(out, U.M);
    }
  }



//! @}
