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


//! \addtogroup spglue_schur
//! @{



template<typename T1, typename T2>
arma_hot
inline
void
spglue_schur::apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_schur>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> pa(X.A);
  const SpProxy<T2> pb(X.B);
  
  const bool is_alias = pa.is_alias(out) || pb.is_alias(out);
  
  if(is_alias == false)
    {
    spglue_schur::apply_noalias(out, pa, pb);
    }
  else
    {
    SpMat<eT> tmp;
    
    spglue_schur::apply_noalias(tmp, pa, pb);
    
    out.steal_mem(tmp);
    }
  }



template<typename eT, typename T1, typename T2>
arma_hot
inline
void
spglue_schur::apply_noalias(SpMat<eT>& out, const SpProxy<T1>& pa, const SpProxy<T2>& pb)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(pa.get_n_rows(), pa.get_n_cols(), pb.get_n_rows(), pb.get_n_cols(), "element-wise multiplication");
  
  if( (pa.get_n_nonzero() == 0) || (pb.get_n_nonzero() == 0) )
    {
    out.zeros(pa.get_n_rows(), pa.get_n_cols());
    return;
    }
  
  const uword max_n_nonzero = (std::min)(pa.get_n_nonzero(), pb.get_n_nonzero());
  
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
    const uword x_it_row = x_it.row();
    const uword x_it_col = x_it.col();
    
    const uword y_it_row = y_it.row();
    const uword y_it_col = y_it.col();
    
    if(x_it == y_it)
      {
      const eT out_val = (*x_it) * (*y_it);
      
      if(out_val != eT(0))
        {
        access::rw(out.values[count]) = out_val;
        
        access::rw(out.row_indices[count]) = x_it_row;
        access::rw(out.col_ptrs[x_it_col + 1])++;
        ++count;
        }
      
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
    
    arma_check( (count > max_n_nonzero), "internal error: spglue_schur::apply_noalias(): count > max_n_nonzero" );
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
arma_hot
inline
void
spglue_schur::apply_noalias(SpMat<eT>& out, const SpMat<eT>& A, const SpMat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  const SpProxy< SpMat<eT> > pa(A);
  const SpProxy< SpMat<eT> > pb(B);
  
  spglue_schur::apply_noalias(out, pa, pb);
  }



//
//
//



template<typename T1, typename T2>
inline
void
spglue_schur_misc::dense_schur_sparse(SpMat<typename T1::elem_type>& out, const T1& x, const T2& y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const   Proxy<T1> pa(x);
  const SpProxy<T2> pb(y);
  
  arma_debug_assert_same_size(pa.get_n_rows(), pa.get_n_cols(), pb.get_n_rows(), pb.get_n_cols(), "element-wise multiplication");
  
  const uword max_n_nonzero = pb.get_n_nonzero();
  
  // Resize memory to upper bound.
  out.reserve(pa.get_n_rows(), pa.get_n_cols(), max_n_nonzero);
  
  uword count = 0;
  
  typename SpProxy<T2>::const_iterator_type it     = pb.begin();
  typename SpProxy<T2>::const_iterator_type it_end = pb.end();
  
  while(it != it_end)
    {
    const uword it_row = it.row();
    const uword it_col = it.col();
    
    const eT val = (*it) * pa.at(it_row, it_col);
    
    if(val != eT(0))
      {
      access::rw(        out.values[count]) = val;
      access::rw(   out.row_indices[count]) = it_row;
      access::rw(out.col_ptrs[it_col + 1])++;
      ++count;
      }
    
    ++it;
    
    arma_check( (count > max_n_nonzero), "internal error: spglue_schur_misc::dense_schur_sparse(): count > max_n_nonzero" );
    }
  
  // Fix column pointers.
  for(uword c = 1; c <= out.n_cols; ++c)
    {
    access::rw(out.col_ptrs[c]) += out.col_ptrs[c - 1];
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
spglue_schur_mixed::apply(SpMat<typename eT_promoter<T1,T2>::eT>& out, const mtSpGlue<typename eT_promoter<T1,T2>::eT, T1, T2, spglue_schur_mixed>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  if( (is_same_type<eT1,out_eT>::no) && (is_same_type<eT2,out_eT>::yes) )
    {
    // upgrade T1
    
    const unwrap_spmat<T1> UA(expr.A);
    const unwrap_spmat<T2> UB(expr.B);
    
    const SpMat<eT1>& A = UA.M;
    const SpMat<eT2>& B = UB.M;
    
    SpMat<out_eT> AA(arma_layout_indicator(), A);
    
    for(uword i=0; i < A.n_nonzero; ++i)  { access::rw(AA.values[i]) = out_eT(A.values[i]); }
    
    const SpMat<out_eT>& BB = reinterpret_cast< const SpMat<out_eT>& >(B);
    
    out = AA % BB;
    }
  else
  if( (is_same_type<eT1,out_eT>::yes) && (is_same_type<eT2,out_eT>::no) )
    {
    // upgrade T2 
    
    const unwrap_spmat<T1> UA(expr.A);
    const unwrap_spmat<T2> UB(expr.B);
    
    const SpMat<eT1>& A = UA.M;
    const SpMat<eT2>& B = UB.M;
    
    const SpMat<out_eT>& AA = reinterpret_cast< const SpMat<out_eT>& >(A);
    
    SpMat<out_eT> BB(arma_layout_indicator(), B);
    
    for(uword i=0; i < B.n_nonzero; ++i)  { access::rw(BB.values[i]) = out_eT(B.values[i]); }
    
    out = AA % BB;
    }
  else
    {
    // upgrade T1 and T2
    
    const unwrap_spmat<T1> UA(expr.A);
    const unwrap_spmat<T2> UB(expr.B);
    
    const SpMat<eT1>& A = UA.M;
    const SpMat<eT2>& B = UB.M;
    
    SpMat<out_eT> AA(arma_layout_indicator(), A);
    SpMat<out_eT> BB(arma_layout_indicator(), B);
    
    for(uword i=0; i < A.n_nonzero; ++i)  { access::rw(AA.values[i]) = out_eT(A.values[i]); }
    for(uword i=0; i < B.n_nonzero; ++i)  { access::rw(BB.values[i]) = out_eT(B.values[i]); }
    
    out = AA % BB;
    }
  }



template<typename T1, typename T2>
inline
void
spglue_schur_mixed::dense_schur_sparse(SpMat< typename promote_type<typename T1::elem_type, typename T2::elem_type >::result>& out, const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const   Proxy<T1> pa(X);
  const SpProxy<T2> pb(Y);
  
  arma_debug_assert_same_size(pa.get_n_rows(), pa.get_n_cols(), pb.get_n_rows(), pb.get_n_cols(), "element-wise multiplication");
  
  // count new size
  uword new_n_nonzero = 0;
  
  typename SpProxy<T2>::const_iterator_type it     = pb.begin();
  typename SpProxy<T2>::const_iterator_type it_end = pb.end();
  
  while(it != it_end)
    {
    if( (out_eT(*it) * out_eT(pa.at(it.row(), it.col()))) != out_eT(0) )  { ++new_n_nonzero; }
    
    ++it;
    }
  
  // Resize memory accordingly.
  out.reserve(pa.get_n_rows(), pa.get_n_cols(), new_n_nonzero);
  
  uword count = 0;
  
  typename SpProxy<T2>::const_iterator_type it2 = pb.begin();
  
  while(it2 != it_end)
    {
    const uword it2_row = it2.row();
    const uword it2_col = it2.col();
    
    const out_eT val = out_eT(*it2) * out_eT(pa.at(it2_row, it2_col));
    
    if(val != out_eT(0))
      {
      access::rw(        out.values[count]) = val;
      access::rw(   out.row_indices[count]) = it2_row;
      access::rw(out.col_ptrs[it2_col + 1])++;
      ++count;
      }
    
    ++it2;
    }
  
  // Fix column pointers.
  for(uword c = 1; c <= out.n_cols; ++c)
    {
    access::rw(out.col_ptrs[c]) += out.col_ptrs[c - 1];
    }
  }



//! @}
