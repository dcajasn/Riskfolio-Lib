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


//! \addtogroup glue_mixed
//! @{



//! matrix multiplication with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_times::apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type in_eT1;
  typedef typename T2::elem_type in_eT2;
  
  typedef typename eT_promoter<T1,T2>::eT out_eT;
  
  const partial_unwrap<T1> tmp1(X.A);
  const partial_unwrap<T2> tmp2(X.B);
  
  const typename partial_unwrap<T1>::stored_type& A = tmp1.M;
  const typename partial_unwrap<T2>::stored_type& B = tmp2.M;
  
  const bool   use_alpha = partial_unwrap<T1>::do_times || partial_unwrap<T2>::do_times;
  const out_eT     alpha = use_alpha ? (upgrade_val<in_eT1,in_eT2>::apply(tmp1.get_val()) * upgrade_val<in_eT1,in_eT2>::apply(tmp2.get_val())) : out_eT(0);
  
  const bool do_trans_A = partial_unwrap<T1>::do_trans;
  const bool do_trans_B = partial_unwrap<T2>::do_trans;
  
  arma_debug_assert_trans_mul_size<do_trans_A, do_trans_B>(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  const uword out_n_rows = (do_trans_A == false) ? A.n_rows : A.n_cols;
  const uword out_n_cols = (do_trans_B == false) ? B.n_cols : B.n_rows;
  
  const bool alias = tmp1.is_alias(out) || tmp2.is_alias(out);
  
  if(alias == false)
    {
    out.set_size(out_n_rows, out_n_cols);
    
    gemm_mixed<do_trans_A, do_trans_B, use_alpha, false>::apply(out, A, B, alpha);
    }
  else
    {
    Mat<out_eT> tmp(out_n_rows, out_n_cols, arma_nozeros_indicator());
    
    gemm_mixed<do_trans_A, do_trans_B, use_alpha, false>::apply(tmp, A, B, alpha);
    
    out.steal_mem(tmp);
    }
  }



//! matrix addition with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_plus::apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_plus>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "addition");
  
  const uword n_rows = A.get_n_rows();
  const uword n_cols = A.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
        out_eT* out_mem = out.memptr();
  const uword   n_elem  = out.n_elem;
    
  const bool use_at = (Proxy<T1>::use_at || Proxy<T2>::use_at);
  
  if(use_at == false)
    {
    typename Proxy<T1>::ea_type AA = A.get_ea();
    typename Proxy<T2>::ea_type BB = B.get_ea();
    
    if(memory::is_aligned(out_mem))
      {
      memory::mark_as_aligned(out_mem);
      
      for(uword i=0; i<n_elem; ++i)
        {
        out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) + upgrade_val<eT1,eT2>::apply(BB[i]);
        }
      }
    else
      {
      for(uword i=0; i<n_elem; ++i)
        {
        out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) + upgrade_val<eT1,eT2>::apply(BB[i]);
        }
      }
    }
  else
    {
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      (*out_mem) = upgrade_val<eT1,eT2>::apply(A.at(row,col)) + upgrade_val<eT1,eT2>::apply(B.at(row,col));
      out_mem++;
      }
    }
  }



//! matrix subtraction with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_minus::apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_minus>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "subtraction");
  
  const uword n_rows = A.get_n_rows();
  const uword n_cols = A.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
        out_eT* out_mem = out.memptr();
  const uword   n_elem  = out.n_elem;
    
  const bool use_at = (Proxy<T1>::use_at || Proxy<T2>::use_at);
  
  if(use_at == false)
    {
    typename Proxy<T1>::ea_type AA = A.get_ea();
    typename Proxy<T2>::ea_type BB = B.get_ea();
    
    if(memory::is_aligned(out_mem))
      {
      memory::mark_as_aligned(out_mem);
      
      for(uword i=0; i<n_elem; ++i)
        {
        out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) - upgrade_val<eT1,eT2>::apply(BB[i]);
        }
      }
    else
      {
      for(uword i=0; i<n_elem; ++i)
        {
        out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) - upgrade_val<eT1,eT2>::apply(BB[i]);
        }
      }
    }
  else
    {
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      (*out_mem) = upgrade_val<eT1,eT2>::apply(A.at(row,col)) - upgrade_val<eT1,eT2>::apply(B.at(row,col));
      out_mem++;
      }
    }
  }



//! element-wise matrix division with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_div::apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_div>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "element-wise division");
  
  const uword n_rows = A.get_n_rows();
  const uword n_cols = A.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
        out_eT* out_mem = out.memptr();
  const uword   n_elem  = out.n_elem;
    
  const bool use_at = (Proxy<T1>::use_at || Proxy<T2>::use_at);
  
  if(use_at == false)
    {
    typename Proxy<T1>::ea_type AA = A.get_ea();
    typename Proxy<T2>::ea_type BB = B.get_ea();
    
    if(memory::is_aligned(out_mem))
      {
      memory::mark_as_aligned(out_mem);
      
      for(uword i=0; i<n_elem; ++i)
        {
        out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) / upgrade_val<eT1,eT2>::apply(BB[i]);
        }
      }
    else
      {
      for(uword i=0; i<n_elem; ++i)
        {
        out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) / upgrade_val<eT1,eT2>::apply(BB[i]);
        }
      }
    }
  else
    {
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      (*out_mem) = upgrade_val<eT1,eT2>::apply(A.at(row,col)) / upgrade_val<eT1,eT2>::apply(B.at(row,col));
      out_mem++;
      }
    }
  }



//! element-wise matrix multiplication with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_schur::apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_schur>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "element-wise multiplication");
  
  const uword n_rows = A.get_n_rows();
  const uword n_cols = A.get_n_cols();
  
  out.set_size(n_rows, n_cols);
  
        out_eT* out_mem = out.memptr();
  const uword   n_elem  = out.n_elem;
    
  const bool use_at = (Proxy<T1>::use_at || Proxy<T2>::use_at);
  
  if(use_at == false)
    {
    typename Proxy<T1>::ea_type AA = A.get_ea();
    typename Proxy<T2>::ea_type BB = B.get_ea();
    
    if(memory::is_aligned(out_mem))
      {
      memory::mark_as_aligned(out_mem);
      
      for(uword i=0; i<n_elem; ++i)
        {
        out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) * upgrade_val<eT1,eT2>::apply(BB[i]);
        }
      }
    else
      {
      for(uword i=0; i<n_elem; ++i)
        {
        out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) * upgrade_val<eT1,eT2>::apply(BB[i]);
        }
      }
    }
  else
    {
    for(uword col=0; col < n_cols; ++col)
    for(uword row=0; row < n_rows; ++row)
      {
      (*out_mem) = upgrade_val<eT1,eT2>::apply(A.at(row,col)) * upgrade_val<eT1,eT2>::apply(B.at(row,col));
      out_mem++;
      }
    }
  }



//
//
//



//! cube addition with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_plus::apply(Cube<typename eT_promoter<T1,T2>::eT>& out, const mtGlueCube<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_plus>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const ProxyCube<T1> A(X.A);
  const ProxyCube<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "addition");
  
  const uword n_rows   = A.get_n_rows();
  const uword n_cols   = A.get_n_cols();
  const uword n_slices = A.get_n_slices();

  out.set_size(n_rows, n_cols, n_slices);
  
        out_eT* out_mem = out.memptr();
  const uword    n_elem = out.n_elem;
  
  const bool use_at = (ProxyCube<T1>::use_at || ProxyCube<T2>::use_at);
  
  if(use_at == false)
    {
    typename ProxyCube<T1>::ea_type AA = A.get_ea();
    typename ProxyCube<T2>::ea_type BB = B.get_ea();
    
    for(uword i=0; i<n_elem; ++i)
      {
      out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) + upgrade_val<eT1,eT2>::apply(BB[i]);
      }
    }
  else
    {
    for(uword slice = 0; slice < n_slices; ++slice)
    for(uword col   = 0; col   < n_cols;   ++col  )
    for(uword row   = 0; row   < n_rows;   ++row  )
      {
      (*out_mem) = upgrade_val<eT1,eT2>::apply(A.at(row,col,slice)) + upgrade_val<eT1,eT2>::apply(B.at(row,col,slice));
      out_mem++;
      }
    }
  }



//! cube subtraction with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_minus::apply(Cube<typename eT_promoter<T1,T2>::eT>& out, const mtGlueCube<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_minus>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const ProxyCube<T1> A(X.A);
  const ProxyCube<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "subtraction");
  
  const uword n_rows   = A.get_n_rows();
  const uword n_cols   = A.get_n_cols();
  const uword n_slices = A.get_n_slices();

  out.set_size(n_rows, n_cols, n_slices);
  
        out_eT* out_mem = out.memptr();
  const uword    n_elem = out.n_elem;
  
  const bool use_at = (ProxyCube<T1>::use_at || ProxyCube<T2>::use_at);
  
  if(use_at == false)
    {
    typename ProxyCube<T1>::ea_type AA = A.get_ea();
    typename ProxyCube<T2>::ea_type BB = B.get_ea();
    
    for(uword i=0; i<n_elem; ++i)
      {
      out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) - upgrade_val<eT1,eT2>::apply(BB[i]);
      }
    }
  else
    {
    for(uword slice = 0; slice < n_slices; ++slice)
    for(uword col   = 0; col   < n_cols;   ++col  )
    for(uword row   = 0; row   < n_rows;   ++row  )
      {
      (*out_mem) = upgrade_val<eT1,eT2>::apply(A.at(row,col,slice)) - upgrade_val<eT1,eT2>::apply(B.at(row,col,slice));
      out_mem++;
      }
    }
  }



//! element-wise cube division with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_div::apply(Cube<typename eT_promoter<T1,T2>::eT>& out, const mtGlueCube<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_div>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const ProxyCube<T1> A(X.A);
  const ProxyCube<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "element-wise division");
  
  const uword n_rows   = A.get_n_rows();
  const uword n_cols   = A.get_n_cols();
  const uword n_slices = A.get_n_slices();

  out.set_size(n_rows, n_cols, n_slices);
  
        out_eT* out_mem = out.memptr();
  const uword    n_elem = out.n_elem;
  
  const bool use_at = (ProxyCube<T1>::use_at || ProxyCube<T2>::use_at);
  
  if(use_at == false)
    {
    typename ProxyCube<T1>::ea_type AA = A.get_ea();
    typename ProxyCube<T2>::ea_type BB = B.get_ea();
    
    for(uword i=0; i<n_elem; ++i)
      {
      out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) / upgrade_val<eT1,eT2>::apply(BB[i]);
      }
    }
  else
    {
    for(uword slice = 0; slice < n_slices; ++slice)
    for(uword col   = 0; col   < n_cols;   ++col  )
    for(uword row   = 0; row   < n_rows;   ++row  )
      {
      (*out_mem) = upgrade_val<eT1,eT2>::apply(A.at(row,col,slice)) / upgrade_val<eT1,eT2>::apply(B.at(row,col,slice));
      out_mem++;
      }
    }
  }



//! element-wise cube multiplication with different element types
template<typename T1, typename T2>
inline
void
glue_mixed_schur::apply(Cube<typename eT_promoter<T1,T2>::eT>& out, const mtGlueCube<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_schur>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const ProxyCube<T1> A(X.A);
  const ProxyCube<T2> B(X.B);
  
  arma_debug_assert_same_size(A, B, "element-wise multiplication");
  
  const uword n_rows   = A.get_n_rows();
  const uword n_cols   = A.get_n_cols();
  const uword n_slices = A.get_n_slices();

  out.set_size(n_rows, n_cols, n_slices);
  
        out_eT* out_mem = out.memptr();
  const uword    n_elem = out.n_elem;
  
  const bool use_at = (ProxyCube<T1>::use_at || ProxyCube<T2>::use_at);
  
  if(use_at == false)
    {
    typename ProxyCube<T1>::ea_type AA = A.get_ea();
    typename ProxyCube<T2>::ea_type BB = B.get_ea();
    
    for(uword i=0; i<n_elem; ++i)
      {
      out_mem[i] = upgrade_val<eT1,eT2>::apply(AA[i]) * upgrade_val<eT1,eT2>::apply(BB[i]);
      }
    }
  else
    {
    for(uword slice = 0; slice < n_slices; ++slice)
    for(uword col   = 0; col   < n_cols;   ++col  )
    for(uword row   = 0; row   < n_rows;   ++row  )
      {
      (*out_mem) = upgrade_val<eT1,eT2>::apply(A.at(row,col,slice)) * upgrade_val<eT1,eT2>::apply(B.at(row,col,slice));
      out_mem++;
      }
    }
  }



//! @}
