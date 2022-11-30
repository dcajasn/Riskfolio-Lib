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


//! \addtogroup spop_reverse
//! @{



template<typename eT>
inline
void
spop_reverse::apply_spmat(SpMat<eT>& out, const SpMat<eT>& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword X_n_rows_m1 = X_n_rows - 1;
  const uword X_n_cols_m1 = X_n_cols - 1;
  
  const uword N = X.n_nonzero;
  
  if(N == uword(0))
    {
    out.zeros(X_n_rows, X_n_cols);
    return;
    }
  
  umat locs(2, N, arma_nozeros_indicator());
  
  uword* locs_mem = locs.memptr();
  
  typename SpMat<eT>::const_iterator it = X.begin();
  
  if(dim == 0)
    {
    for(uword i=0; i < N; ++i)
      {
      const uword row = it.row();
      const uword col = it.col();
      
      (*locs_mem) = X_n_rows_m1 - row;  locs_mem++;
      (*locs_mem) =               col;  locs_mem++;
      
      ++it;
      }
    }
  else
  if(dim == 1)
    {
    for(uword i=0; i < N; ++i)
      {
      const uword row = it.row();
      const uword col = it.col();
      
      (*locs_mem) =               row;  locs_mem++;
      (*locs_mem) = X_n_cols_m1 - col;  locs_mem++;
      
      ++it;
      }
    }
  
  const Col<eT> vals(const_cast<eT*>(X.values), N, false);
  
  SpMat<eT> tmp(locs, vals, X_n_rows, X_n_cols, true, false);
  
  out.steal_mem(tmp);
  }



template<typename T1>
inline
void
spop_reverse::apply_proxy(SpMat<typename T1::elem_type>& out, const T1& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> P(X);
  
  const uword P_n_rows = P.get_n_rows();
  const uword P_n_cols = P.get_n_cols();
  
  const uword P_n_rows_m1 = P_n_rows - 1;
  const uword P_n_cols_m1 = P_n_cols - 1;
  
  const uword N = P.get_n_nonzero();
  
  if(N == uword(0))
    {
    out.zeros(P_n_rows, P_n_cols);
    return;
    }
  
  umat    locs(2, N, arma_nozeros_indicator());
  Col<eT> vals(   N, arma_nozeros_indicator());
  
  uword* locs_mem = locs.memptr();
  eT*    vals_mem = vals.memptr();
  
  typename SpProxy<T1>::const_iterator_type it = P.begin();
  
  if(dim == 0)
    {
    for(uword i=0; i < N; ++i)
      {
      const uword row = it.row();
      const uword col = it.col();
      
      (*locs_mem) = P_n_rows_m1 - row;  locs_mem++;
      (*locs_mem) =               col;  locs_mem++;
      
      (*vals_mem) = (*it);  vals_mem++;
      
      ++it;
      }
    }
  else
  if(dim == 1)
    {
    for(uword i=0; i < N; ++i)
      {
      const uword row = it.row();
      const uword col = it.col();
      
      (*locs_mem) =               row;  locs_mem++;
      (*locs_mem) = P_n_cols_m1 - col;  locs_mem++;
      
      (*vals_mem) = (*it);  vals_mem++;
      
      ++it;
      }
    }
  
  SpMat<eT> tmp(locs, vals, P_n_rows, P_n_cols, true, false);
  
  out.steal_mem(tmp);
  }



template<typename T1>
inline
void
spop_reverse::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_reverse>& in)
  {
  arma_extra_debug_sigprint();
  
  const uword dim = in.aux_uword_a;
  
  arma_debug_check( (dim > 1), "reverse(): parameter 'dim' must be 0 or 1" );
  
  if(is_SpMat<T1>::value)
    {
    const unwrap_spmat<T1> tmp(in.m);
    
    spop_reverse::apply_spmat(out, tmp.M, dim);
    }
  else
    {
    spop_reverse::apply_proxy(out, in.m, dim);
    }
  }



//! @}
