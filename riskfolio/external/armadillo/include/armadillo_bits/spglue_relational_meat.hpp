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


//! \addtogroup spglue_relational
//! @{



template<typename T1, typename T2>
inline
void
spglue_rel_lt::apply(SpMat<uword>& out, const mtSpGlue<uword, T1, T2, spglue_rel_lt>& X)
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<T1> PA(X.A);
  const SpProxy<T2> PB(X.B);
  
  const bool is_alias = PA.is_alias(out) || PB.is_alias(out);
  
  if(is_alias == false)
    {
    spglue_rel_lt::apply_noalias(out, PA, PB);
    }
  else
    {
    SpMat<uword> tmp;
    
    spglue_rel_lt::apply_noalias(tmp, PA, PB);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1, typename T2>
inline
void
spglue_rel_lt::apply_noalias(SpMat<uword>& out, const SpProxy<T1>& PA, const SpProxy<T2>& PB)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  arma_debug_assert_same_size(PA.get_n_rows(), PA.get_n_cols(), PB.get_n_rows(), PB.get_n_cols(), "operator<");
  
  const uword max_n_nonzero = PA.get_n_nonzero() + PB.get_n_nonzero();
  
  // Resize memory to upper bound
  out.reserve(PA.get_n_rows(), PA.get_n_cols(), max_n_nonzero);
  
  // Now iterate across both matrices.
  typename SpProxy<T1>::const_iterator_type x_it  = PA.begin();
  typename SpProxy<T1>::const_iterator_type x_end = PA.end();
  
  typename SpProxy<T2>::const_iterator_type y_it  = PB.begin();
  typename SpProxy<T2>::const_iterator_type y_end = PB.end();
  
  uword count = 0;
  
  while( (x_it != x_end) || (y_it != y_end) )
    {
    uword out_val;
    
    const uword x_it_col = x_it.col();
    const uword x_it_row = x_it.row();
    
    const uword y_it_col = y_it.col();
    const uword y_it_row = y_it.row();
    
    bool use_y_loc = false;
    
    if(x_it == y_it)
      {
      out_val = ((*x_it) < (*y_it)) ? uword(1) : uword(0);
      
      ++x_it;
      ++y_it;
      }
    else
      {
      if((x_it_col < y_it_col) || ((x_it_col == y_it_col) && (x_it_row < y_it_row))) // if y is closer to the end
        {
        out_val = ((*x_it) < eT(0)) ? uword(1) : uword(0);
        
        ++x_it;
        }
      else
        {
        out_val = (eT(0) < (*y_it)) ? uword(1) : uword(0);
        
        ++y_it;
        
        use_y_loc = true;
        }
      }
    
    if(out_val != uword(0))
      {
      access::rw(out.values[count]) = out_val;
      
      const uword out_row = (use_y_loc == false) ? x_it_row : y_it_row;
      const uword out_col = (use_y_loc == false) ? x_it_col : y_it_col;
      
      access::rw(out.row_indices[count]) = out_row;
      access::rw(out.col_ptrs[out_col + 1])++;
      ++count;
      }
    
    arma_check( (count > max_n_nonzero), "internal error: spglue_rel_lt::apply_noalias(): count > max_n_nonzero" );
    }
  
  const uword out_n_cols = out.n_cols;
  
  uword* col_ptrs = access::rwp(out.col_ptrs);
  
  // Fix column pointers to be cumulative.
  for(uword c = 1; c <= out_n_cols; ++c)
    {
    col_ptrs[c] += col_ptrs[c - 1];
    }
  
  if(count < max_n_nonzero)
    {
    if(count <= (max_n_nonzero/2))
      {
      out.mem_resize(count);
      }
    else
      {
      // quick resize without reallocating memory and copying data
      access::rw(         out.n_nonzero) = count;
      access::rw(     out.values[count]) = eT(0);
      access::rw(out.row_indices[count]) = uword(0);
      }
    }
  }



//



template<typename T1, typename T2>
inline
void
spglue_rel_gt::apply(SpMat<uword>& out, const mtSpGlue<uword, T1, T2, spglue_rel_gt>& X)
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<T1> PA(X.A);
  const SpProxy<T2> PB(X.B);
  
  const bool is_alias = PA.is_alias(out) || PB.is_alias(out);
  
  if(is_alias == false)
    {
    spglue_rel_gt::apply_noalias(out, PA, PB);
    }
  else
    {
    SpMat<uword> tmp;
    
    spglue_rel_gt::apply_noalias(tmp, PA, PB);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1, typename T2>
inline
void
spglue_rel_gt::apply_noalias(SpMat<uword>& out, const SpProxy<T1>& PA, const SpProxy<T2>& PB)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  arma_debug_assert_same_size(PA.get_n_rows(), PA.get_n_cols(), PB.get_n_rows(), PB.get_n_cols(), "operator>");
  
  const uword max_n_nonzero = PA.get_n_nonzero() + PB.get_n_nonzero();
  
  // Resize memory to upper bound
  out.reserve(PA.get_n_rows(), PA.get_n_cols(), max_n_nonzero);
  
  // Now iterate across both matrices.
  typename SpProxy<T1>::const_iterator_type x_it  = PA.begin();
  typename SpProxy<T1>::const_iterator_type x_end = PA.end();
  
  typename SpProxy<T2>::const_iterator_type y_it  = PB.begin();
  typename SpProxy<T2>::const_iterator_type y_end = PB.end();
  
  uword count = 0;
  
  while( (x_it != x_end) || (y_it != y_end) )
    {
    uword out_val;
    
    const uword x_it_col = x_it.col();
    const uword x_it_row = x_it.row();
    
    const uword y_it_col = y_it.col();
    const uword y_it_row = y_it.row();
    
    bool use_y_loc = false;
    
    if(x_it == y_it)
      {
      out_val = ((*x_it) > (*y_it)) ? uword(1) : uword(0);
      
      ++x_it;
      ++y_it;
      }
    else
      {
      if((x_it_col < y_it_col) || ((x_it_col == y_it_col) && (x_it_row < y_it_row))) // if y is closer to the end
        {
        out_val = ((*x_it) > eT(0)) ? uword(1) : uword(0);
        
        ++x_it;
        }
      else
        {
        out_val = (eT(0) > (*y_it)) ? uword(1) : uword(0);
        
        ++y_it;
        
        use_y_loc = true;
        }
      }
    
    if(out_val != uword(0))
      {
      access::rw(out.values[count]) = out_val;
      
      const uword out_row = (use_y_loc == false) ? x_it_row : y_it_row;
      const uword out_col = (use_y_loc == false) ? x_it_col : y_it_col;
      
      access::rw(out.row_indices[count]) = out_row;
      access::rw(out.col_ptrs[out_col + 1])++;
      ++count;
      }
    
    arma_check( (count > max_n_nonzero), "internal error: spglue_rel_gt::apply_noalias(): count > max_n_nonzero" );
    }
  
  const uword out_n_cols = out.n_cols;
  
  uword* col_ptrs = access::rwp(out.col_ptrs);
  
  // Fix column pointers to be cumulative.
  for(uword c = 1; c <= out_n_cols; ++c)
    {
    col_ptrs[c] += col_ptrs[c - 1];
    }
  
  if(count < max_n_nonzero)
    {
    if(count <= (max_n_nonzero/2))
      {
      out.mem_resize(count);
      }
    else
      {
      // quick resize without reallocating memory and copying data
      access::rw(         out.n_nonzero) = count;
      access::rw(     out.values[count]) = eT(0);
      access::rw(out.row_indices[count]) = uword(0);
      }
    }
  }



//



template<typename T1, typename T2>
inline
void
spglue_rel_and::apply(SpMat<uword>& out, const mtSpGlue<uword, T1, T2, spglue_rel_and>& X)
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<T1> PA(X.A);
  const SpProxy<T2> PB(X.B);
  
  const bool is_alias = PA.is_alias(out) || PB.is_alias(out);
  
  if(is_alias == false)
    {
    spglue_rel_and::apply_noalias(out, PA, PB);
    }
  else
    {
    SpMat<uword> tmp;
    
    spglue_rel_and::apply_noalias(tmp, PA, PB);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1, typename T2>
inline
void
spglue_rel_and::apply_noalias(SpMat<uword>& out, const SpProxy<T1>& PA, const SpProxy<T2>& PB)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  arma_debug_assert_same_size(PA.get_n_rows(), PA.get_n_cols(), PB.get_n_rows(), PB.get_n_cols(), "operator&&");
  
  if( (PA.get_n_nonzero() == 0) || (PB.get_n_nonzero() == 0) )
    {
    out.zeros(PA.get_n_rows(), PA.get_n_cols());
    return;
    }
  
  const uword max_n_nonzero = (std::min)(PA.get_n_nonzero(), PB.get_n_nonzero());
  
  // Resize memory to upper bound
  out.reserve(PA.get_n_rows(), PA.get_n_cols(), max_n_nonzero);
  
  // Now iterate across both matrices.
  typename SpProxy<T1>::const_iterator_type x_it  = PA.begin();
  typename SpProxy<T1>::const_iterator_type x_end = PA.end();
  
  typename SpProxy<T2>::const_iterator_type y_it  = PB.begin();
  typename SpProxy<T2>::const_iterator_type y_end = PB.end();
  
  uword count = 0;
  
  while( (x_it != x_end) || (y_it != y_end) )
    {
    const uword x_it_row = x_it.row();
    const uword x_it_col = x_it.col();
    
    const uword y_it_row = y_it.row();
    const uword y_it_col = y_it.col();
    
    if(x_it == y_it)
      {
      access::rw(out.values[count]) = uword(1);
      
      access::rw(out.row_indices[count]) = x_it_row;
      access::rw(out.col_ptrs[x_it_col + 1])++;
      ++count;
      
      ++x_it;
      ++y_it;
      }
    else
      {
      if((x_it_col < y_it_col) || ((x_it_col == y_it_col) && (x_it_row < y_it_row))) // if y is closer to the end
        {
        ++x_it;
        }
      else
        {
        ++y_it;
        }
      }
    
    arma_check( (count > max_n_nonzero), "internal error: spglue_rel_and::apply_noalias(): count > max_n_nonzero" );
    }
  
  const uword out_n_cols = out.n_cols;
  
  uword* col_ptrs = access::rwp(out.col_ptrs);
  
  // Fix column pointers to be cumulative.
  for(uword c = 1; c <= out_n_cols; ++c)
    {
    col_ptrs[c] += col_ptrs[c - 1];
    }
  
  if(count < max_n_nonzero)
    {
    if(count <= (max_n_nonzero/2))
      {
      out.mem_resize(count);
      }
    else
      {
      // quick resize without reallocating memory and copying data
      access::rw(         out.n_nonzero) = count;
      access::rw(     out.values[count]) = eT(0);
      access::rw(out.row_indices[count]) = uword(0);
      }
    }
  }



//



template<typename T1, typename T2>
inline
void
spglue_rel_or::apply(SpMat<uword>& out, const mtSpGlue<uword, T1, T2, spglue_rel_or>& X)
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<T1> PA(X.A);
  const SpProxy<T2> PB(X.B);
  
  const bool is_alias = PA.is_alias(out) || PB.is_alias(out);
  
  if(is_alias == false)
    {
    spglue_rel_or::apply_noalias(out, PA, PB);
    }
  else
    {
    SpMat<uword> tmp;
    
    spglue_rel_or::apply_noalias(tmp, PA, PB);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1, typename T2>
inline
void
spglue_rel_or::apply_noalias(SpMat<uword>& out, const SpProxy<T1>& PA, const SpProxy<T2>& PB)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  arma_debug_assert_same_size(PA.get_n_rows(), PA.get_n_cols(), PB.get_n_rows(), PB.get_n_cols(), "operator||");
  
  const uword max_n_nonzero = PA.get_n_nonzero() + PB.get_n_nonzero();
  
  // Resize memory to upper bound
  out.reserve(PA.get_n_rows(), PA.get_n_cols(), max_n_nonzero);
  
  // Now iterate across both matrices.
  typename SpProxy<T1>::const_iterator_type x_it  = PA.begin();
  typename SpProxy<T1>::const_iterator_type x_end = PA.end();
  
  typename SpProxy<T2>::const_iterator_type y_it  = PB.begin();
  typename SpProxy<T2>::const_iterator_type y_end = PB.end();
  
  uword count = 0;
  
  while( (x_it != x_end) || (y_it != y_end) )
    {
    const uword x_it_col = x_it.col();
    const uword x_it_row = x_it.row();
    
    const uword y_it_col = y_it.col();
    const uword y_it_row = y_it.row();
    
    bool use_y_loc = false;
    
    if(x_it == y_it)
      {
      ++x_it;
      ++y_it;
      }
    else
      {
      if((x_it_col < y_it_col) || ((x_it_col == y_it_col) && (x_it_row < y_it_row))) // if y is closer to the end
        {
        ++x_it;
        }
      else
        {
        ++y_it;
        
        use_y_loc = true;
        }
      }
    
    access::rw(out.values[count]) = uword(1);
    
    const uword out_row = (use_y_loc == false) ? x_it_row : y_it_row;
    const uword out_col = (use_y_loc == false) ? x_it_col : y_it_col;
    
    access::rw(out.row_indices[count]) = out_row;
    access::rw(out.col_ptrs[out_col + 1])++;
    ++count;
    
    arma_check( (count > max_n_nonzero), "internal error: spglue_rel_or::apply_noalias(): count > max_n_nonzero" );
    }
  
  const uword out_n_cols = out.n_cols;
  
  uword* col_ptrs = access::rwp(out.col_ptrs);
  
  // Fix column pointers to be cumulative.
  for(uword c = 1; c <= out_n_cols; ++c)
    {
    col_ptrs[c] += col_ptrs[c - 1];
    }
  
  if(count < max_n_nonzero)
    {
    if(count <= (max_n_nonzero/2))
      {
      out.mem_resize(count);
      }
    else
      {
      // quick resize without reallocating memory and copying data
      access::rw(         out.n_nonzero) = count;
      access::rw(     out.values[count]) = eT(0);
      access::rw(out.row_indices[count]) = uword(0);
      }
    }
  }



//! @}
