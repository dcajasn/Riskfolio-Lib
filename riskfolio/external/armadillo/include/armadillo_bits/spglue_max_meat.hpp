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


//! \addtogroup spglue_max
//! @{



template<typename T1, typename T2>
inline
void
spglue_max::apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_max>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> pa(X.A);
  const SpProxy<T2> pb(X.B);
  
  const bool is_alias = pa.is_alias(out) || pb.is_alias(out);
  
  if(is_alias == false)
    {
    spglue_max::apply_noalias(out, pa, pb);
    }
  else
    {
    SpMat<eT> tmp;
    
    spglue_max::apply_noalias(tmp, pa, pb);
    
    out.steal_mem(tmp);
    }
  }



template<typename eT, typename T1, typename T2>
inline
void
spglue_max::apply_noalias(SpMat<eT>& out, const SpProxy<T1>& pa, const SpProxy<T2>& pb)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(pa.get_n_rows(), pa.get_n_cols(), pb.get_n_rows(), pb.get_n_cols(), "element-wise max()");
  
  const uword max_n_nonzero = pa.get_n_nonzero() + pb.get_n_nonzero();
  
  // Resize memory to upper bound
  out.reserve(pa.get_n_rows(), pa.get_n_cols(), max_n_nonzero);
  
  // Now iterate across both matrices.
  typename SpProxy<T1>::const_iterator_type x_it  = pa.begin();
  typename SpProxy<T1>::const_iterator_type x_end = pa.end();
  
  typename SpProxy<T2>::const_iterator_type y_it  = pb.begin();
  typename SpProxy<T2>::const_iterator_type y_end = pb.end();
  
  uword count = 0;
  
  while( (x_it != x_end) || (y_it != y_end) )
    {
    eT out_val;
    
    const uword x_it_col = x_it.col();
    const uword x_it_row = x_it.row();
    
    const uword y_it_col = y_it.col();
    const uword y_it_row = y_it.row();
    
    bool use_y_loc = false;
    
    if(x_it == y_it)
      {
      out_val = elem_max(eT(*x_it), eT(*y_it));
      
      ++x_it;
      ++y_it;
      }
    else
      {
      if((x_it_col < y_it_col) || ((x_it_col == y_it_col) && (x_it_row < y_it_row))) // if y is closer to the end
        {
        out_val = elem_max(eT(*x_it), eT(0));
        
        ++x_it;
        }
      else
        {
        out_val = elem_max(eT(*y_it), eT(0));
        
        ++y_it;
        
        use_y_loc = true;
        }
      }
    
    if(out_val != eT(0))
      {
      access::rw(out.values[count]) = out_val;
      
      const uword out_row = (use_y_loc == false) ? x_it_row : y_it_row;
      const uword out_col = (use_y_loc == false) ? x_it_col : y_it_col;
      
      access::rw(out.row_indices[count]) = out_row;
      access::rw(out.col_ptrs[out_col + 1])++;
      ++count;
      }
    
    arma_check( (count > max_n_nonzero), "internal error: spglue_max::apply_noalias(): count > max_n_nonzero" );
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



template<typename eT>
inline
void
spglue_max::apply_noalias(SpMat<eT>& out, const SpMat<eT>& A, const SpMat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  const SpProxy< SpMat<eT> > pa(A);
  const SpProxy< SpMat<eT> > pb(B);
  
  spglue_max::apply_noalias(out, pa, pb);
  }



template<typename eT, typename T1, typename T2>
inline
void
spglue_max::dense_sparse_max(Mat<eT>& out, const Base<eT,T1>& X, const SpBase<eT,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: this function assumes there is no aliasing between matrix 'out' and X
  
  const   Proxy<T1> pa(X.get_ref());
  const SpProxy<T2> pb(Y.get_ref());
  
  const uword n_rows = pa.get_n_rows();
  const uword n_cols = pa.get_n_cols();
  
  arma_debug_assert_same_size( n_rows, n_cols, pb.get_n_rows(), pb.get_n_cols(), "element-wise max()" );
  
  out.set_size(n_rows, n_cols);
  
  for(uword c=0; c < n_cols; ++c)
  for(uword r=0; r < n_rows; ++r)
    {
    out.at(r,c) = elem_max(pa.at(r,c), pb.at(r,c));
    }
  }



//! max of non-complex elements
template<typename eT>
inline
typename enable_if2<is_cx<eT>::no, eT>::result
spglue_max::elem_max(const eT& a, const eT& b)
  {
  return (std::max)(a, b);
  }



//! max of complex elements
template<typename eT>
inline
typename enable_if2<is_cx<eT>::yes, eT>::result
spglue_max::elem_max(const eT& a, const eT& b)
  {
  return (std::abs(a) > std::abs(b)) ? a : b;
  }



//! @}
