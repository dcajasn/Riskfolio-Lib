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


//! \addtogroup spop_normalise
//! @{



template<typename T1>
inline
void
spop_normalise::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_normalise>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword p   = expr.aux_uword_a;
  const uword dim = expr.aux_uword_b;
  
  arma_debug_check( (p   == 0), "normalise(): parameter 'p' must be greater than zero" );
  arma_debug_check( (dim >  1), "normalise(): parameter 'dim' must be 0 or 1"          );
  
  const unwrap_spmat<T1> U(expr.m);
  
  const SpMat<eT>& X = U.M;
  
  X.sync();
  
  if( X.is_empty() || (X.n_nonzero == 0) )  { out.zeros(X.n_rows, X.n_cols); return; }
  
  if(dim == 0)
    {
    spop_normalise::apply_direct(out, X, p);
    }
  else
  if(dim == 1)
    {
    SpMat<eT> tmp1;
    SpMat<eT> tmp2;
    
    spop_strans::apply_noalias(tmp1, X);
    
    spop_normalise::apply_direct(tmp2, tmp1, p);
    
    spop_strans::apply_noalias(out, tmp2);
    }
  }



template<typename eT>
inline
void
spop_normalise::apply_direct(SpMat<eT>& out, const SpMat<eT>& X, const uword p)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  SpMat<eT> tmp(arma_reserve_indicator(), X.n_rows, X.n_cols, X.n_nonzero);
  
  bool has_zero = false;
  
  podarray<T> norm_vals(X.n_cols);
  
  T* norm_vals_mem = norm_vals.memptr();
  
  for(uword col=0; col < X.n_cols; ++col)
    {
    const uword      col_offset = X.col_ptrs[col    ];
    const uword next_col_offset = X.col_ptrs[col + 1];
    
    const eT* start_ptr = &X.values[     col_offset];
    const eT*   end_ptr = &X.values[next_col_offset];
    
    const uword n_elem = end_ptr - start_ptr;
    
    const Col<eT> fake_vec(const_cast<eT*>(start_ptr), n_elem, false, false);
    
    const T norm_val = norm(fake_vec, p);
    
    norm_vals_mem[col] = (norm_val != T(0)) ? norm_val : T(1);
    }
  
  const uword N = X.n_nonzero;
  
  typename SpMat<eT>::const_iterator it = X.begin();
  
  for(uword i=0; i < N; ++i)
    {
    const uword row = it.row();
    const uword col = it.col();
    
    const eT val = (*it) / norm_vals_mem[col];
    
    if(val == eT(0))  { has_zero = true; }
    
    access::rw(tmp.values[i])      = val;
    access::rw(tmp.row_indices[i]) = row;
    access::rw(tmp.col_ptrs[col + 1])++;
    
    ++it;
    }
  
  for(uword c=0; c < tmp.n_cols; ++c)
    {
    access::rw(tmp.col_ptrs[c + 1]) += tmp.col_ptrs[c];
    }
  
  if(has_zero)  { tmp.remove_zeros(); }
  
  out.steal_mem(tmp);
  }



//! @}
