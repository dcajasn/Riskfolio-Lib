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


//! \addtogroup spglue_times
//! @{



template<typename T1, typename T2>
inline
void
spglue_times::apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> UA(X.A);
  const unwrap_spmat<T2> UB(X.B);
  
  const bool is_alias = (UA.is_alias(out) || UB.is_alias(out));
  
  if(is_alias == false)
    {
    spglue_times::apply_noalias(out, UA.M, UB.M);
    }
  else
    {
    SpMat<eT> tmp;
    
    spglue_times::apply_noalias(tmp, UA.M, UB.M);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1, typename T2>
inline
void
spglue_times::apply(SpMat<typename T1::elem_type>& out, const SpGlue<SpOp<T1,spop_scalar_times>,T2,spglue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> UA(X.A.m);
  const unwrap_spmat<T2> UB(X.B);
  
  const bool is_alias = (UA.is_alias(out) || UB.is_alias(out));
  
  if(is_alias == false)
    {
    spglue_times::apply_noalias(out, UA.M, UB.M);
    }
  else
    {
    SpMat<eT> tmp;
    
    spglue_times::apply_noalias(tmp, UA.M, UB.M);
    
    out.steal_mem(tmp);
    }
  
  out *= X.A.aux;
  }



template<typename eT>
arma_hot
inline
void
spglue_times::apply_noalias(SpMat<eT>& c, const SpMat<eT>& x, const SpMat<eT>& y)
  {
  arma_extra_debug_sigprint();
  
  const uword x_n_rows = x.n_rows;
  const uword x_n_cols = x.n_cols;
  const uword y_n_rows = y.n_rows;
  const uword y_n_cols = y.n_cols;

  arma_debug_assert_mul_size(x_n_rows, x_n_cols, y_n_rows, y_n_cols, "matrix multiplication");

  // First we must determine the structure of the new matrix (column pointers).
  // This follows the algorithm described in 'Sparse Matrix Multiplication
  // Package (SMMP)' (R.E. Bank and C.C. Douglas, 2001).  Their description of
  // "SYMBMM" does not include anything about memory allocation.  In addition it
  // does not consider that there may be elements which space may be allocated
  // for but which evaluate to zero anyway.  So we have to modify the algorithm
  // to work that way.  For the "SYMBMM" implementation we will not determine
  // the row indices but instead just the column pointers.
  
  //SpMat<typename T1::elem_type> c(x_n_rows, y_n_cols); // Initializes col_ptrs to 0.
  c.zeros(x_n_rows, y_n_cols);
  
  //if( (x.n_elem == 0) || (y.n_elem == 0) )  { return; }
  if( (x.n_nonzero == 0) || (y.n_nonzero == 0) )  { return; }
  
  // Auxiliary storage which denotes when items have been found.
  podarray<uword> index(x_n_rows);
  index.fill(x_n_rows); // Fill with invalid links.
  
  typename SpMat<eT>::const_iterator y_it  = y.begin();
  typename SpMat<eT>::const_iterator y_end = y.end();

  // SYMBMM: calculate column pointers for resultant matrix to obtain a good
  // upper bound on the number of nonzero elements.
  uword cur_col_length = 0;
  uword last_ind = x_n_rows + 1;
  do
    {
    const uword y_it_row = y_it.row();
    
    // Look through the column that this point (*y_it) could affect.
    typename SpMat<eT>::const_iterator x_it = x.begin_col_no_sync(y_it_row);
    
    while(x_it.col() == y_it_row)
      {
      const uword x_it_row = x_it.row();
      
      // A point at x(i, j) and y(j, k) implies a point at c(i, k).
      if(index[x_it_row] == x_n_rows)
        {
        index[x_it_row] = last_ind;
        last_ind = x_it_row;
        ++cur_col_length;
        }

      ++x_it;
      }

    const uword old_col = y_it.col();
    ++y_it;

    // See if column incremented.
    if(old_col != y_it.col())
      {
      // Set column pointer (this is not a cumulative count; that is done later).
      access::rw(c.col_ptrs[old_col + 1]) = cur_col_length;
      cur_col_length = 0;

      // Return index markers to zero.  Use last_ind for traversal.
      while(last_ind != x_n_rows + 1)
        {
        const uword tmp = index[last_ind];
        index[last_ind] = x_n_rows;
        last_ind = tmp;
        }
      }
    }
  while(y_it != y_end);

  // Accumulate column pointers.
  for(uword i = 0; i < c.n_cols; ++i)
    {
    access::rw(c.col_ptrs[i + 1]) += c.col_ptrs[i];
    }

  // Now that we know a decent bound on the number of nonzero elements, allocate
  // the memory and fill it.
  c.mem_resize(c.col_ptrs[c.n_cols]);

  // Now the implementation of the NUMBMM algorithm.
  uword cur_pos = 0; // Current position in c matrix.
  podarray<eT> sums(x_n_rows); // Partial sums.
  sums.zeros();
  
  podarray<uword> sorted_indices(x_n_rows);  // upper bound
  
  // last_ind is already set to x_n_rows, and cur_col_length is already set to 0.
  // We will loop through all columns as necessary.
  uword cur_col = 0;
  while(cur_col < c.n_cols)
    {
    // Skip to next column with elements in it.
    while((cur_col < c.n_cols) && (c.col_ptrs[cur_col] == c.col_ptrs[cur_col + 1]))
      {
      // Update current column pointer to actual number of nonzero elements up
      // to this point.
      access::rw(c.col_ptrs[cur_col]) = cur_pos;
      ++cur_col;
      }

    if(cur_col == c.n_cols)
      {
      break;
      }

    // Update current column pointer.
    access::rw(c.col_ptrs[cur_col]) = cur_pos;

    // Check all elements in this column.
    typename SpMat<eT>::const_iterator y_col_it = y.begin_col_no_sync(cur_col);
    
    while(y_col_it.col() == cur_col)
      {
      const uword y_col_it_row = y_col_it.row();
      
      // Check all elements in the column of the other matrix corresponding to
      // the row of this column.
      typename SpMat<eT>::const_iterator x_col_it = x.begin_col_no_sync(y_col_it_row);

      const eT y_value = (*y_col_it);

      while(x_col_it.col() == y_col_it_row)
        {
        const uword x_col_it_row = x_col_it.row();
        
        // A point at x(i, j) and y(j, k) implies a point at c(i, k).
        // Add to partial sum.
        const eT x_value = (*x_col_it);
        sums[x_col_it_row] += (x_value * y_value);

        // Add point if it hasn't already been marked.
        if(index[x_col_it_row] == x_n_rows)
          {
          index[x_col_it_row] = last_ind;
          last_ind = x_col_it_row;
          }

        ++x_col_it;
        }

      ++y_col_it;
      }

    // Now sort the indices that were used in this column.
    uword cur_index = 0;
    while(last_ind != x_n_rows + 1)
      {
      const uword tmp = last_ind;

      // Check that it wasn't a "fake" nonzero element.
      if(sums[tmp] != eT(0))
        {
        // Assign to next open position.
        sorted_indices[cur_index] = tmp;
        ++cur_index;
        }

      last_ind = index[tmp];
      index[tmp] = x_n_rows;
      }

    // Now sort the indices.
    if(cur_index != 0)
      {
      op_sort::direct_sort_ascending(sorted_indices.memptr(), cur_index);

      for(uword k = 0; k < cur_index; ++k)
        {
        const uword row = sorted_indices[k];
        access::rw(c.row_indices[cur_pos]) = row;
        access::rw(c.values[cur_pos]) = sums[row];
        sums[row] = eT(0);
        ++cur_pos;
        }
      }

    // Move to next column.
    ++cur_col;
    }

  // Update last column pointer and resize to actual memory size.
  access::rw(c.col_ptrs[c.n_cols]) = cur_pos;
  c.mem_resize(cur_pos);
  }



//
//
//



template<typename T1, typename T2>
inline
void
spglue_times_misc::sparse_times_dense(Mat<typename T1::elem_type>& out, const T1& x, const T2& y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_op_diagmat<T2>::value)
    {
    const SpMat<eT> tmp(y);
    
    out = x * tmp;
    }
  else
    {
    const unwrap_spmat<T1> UA(x);
    const quasi_unwrap<T2> UB(y);
    
    const SpMat<eT>& A = UA.M;
    const   Mat<eT>& B = UB.M;
    
    if( (resolves_to_vector<T2>::no) && (B.is_vec() == false) && B.is_diagmat() )
      {
      const SpMat<eT> tmp(diagmat(B));
      
      out = A * tmp;
      
      return;
      }
    
    const uword A_n_rows = A.n_rows;
    const uword A_n_cols = A.n_cols;
    
    const uword B_n_rows = B.n_rows;
    const uword B_n_cols = B.n_cols;
    
    arma_debug_assert_mul_size(A_n_rows, A_n_cols, B_n_rows, B_n_cols, "matrix multiplication");
    
    if(B_n_cols >= (B_n_rows / uword(100)))
      {
      arma_extra_debug_print("using transpose-based multiplication");
      
      const SpMat<eT> At = A.st();
      const   Mat<eT> Bt = B.st();
      
      if(A_n_rows == B_n_cols)
        {
        spglue_times_misc::dense_times_sparse(out, Bt, At);
        
        op_strans::apply_mat(out, out);  // since 'out' is square-sized, this will do an inplace transpose
        }
      else
        {
        Mat<eT> tmp;
        
        spglue_times_misc::dense_times_sparse(tmp, Bt, At);
        
        op_strans::apply_mat(out, tmp);
        }
      }
    else
      {
      arma_extra_debug_print("using standard multiplication");
      
      out.zeros(A_n_rows, B_n_cols);
      
      typename SpMat<eT>::const_iterator A_it     = A.begin();
      typename SpMat<eT>::const_iterator A_it_end = A.end();
      
      while(A_it != A_it_end)
        {
        const eT    A_it_val = (*A_it);
        const uword A_it_row = A_it.row();
        const uword A_it_col = A_it.col();
        
        for(uword col = 0; col < B_n_cols; ++col)
          {
          out.at(A_it_row, col) += A_it_val * B.at(A_it_col, col);
          }
        
        ++A_it;
        }
      }
    }
  }



template<typename T1, typename T2>
inline
void
spglue_times_misc::dense_times_sparse(Mat<typename T1::elem_type>& out, const T1& x, const T2& y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_op_diagmat<T1>::value)
    {
    const SpMat<eT> tmp(x);
    
    out = tmp * y;
    }
  else
    {
    const quasi_unwrap<T1> UA(x);
    const unwrap_spmat<T2> UB(y);
    
    const   Mat<eT>& A = UA.M;
    const SpMat<eT>& B = UB.M;
    
    if( (resolves_to_vector<T1>::no) && (A.is_vec() == false) && A.is_diagmat() )
      {
      const SpMat<eT> tmp(diagmat(A));
      
      out = tmp * B;
      
      return;
      }
    
    arma_debug_assert_mul_size(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
    
    out.zeros(A.n_rows, B.n_cols);
    
    if( (A.n_elem > 0) && (B.n_nonzero > 0) )
      {
      if( (arma_config::openmp) && (mp_thread_limit::in_parallel() == false) && (A.n_rows <= (A.n_cols / uword(100))) )
        {
        #if defined(ARMA_USE_OPENMP)
          {
          arma_extra_debug_print("using parallelised multiplication");
          
          const uword B_n_cols  = B.n_cols;
          const int   n_threads = mp_thread_limit::get();
          
          #pragma omp parallel for schedule(static) num_threads(n_threads)
          for(uword i=0; i < B_n_cols; ++i)
            {
            const uword col_offset_1 = B.col_ptrs[i  ];
            const uword col_offset_2 = B.col_ptrs[i+1];
            
            const uword col_offset_delta = col_offset_2 - col_offset_1;
            
            const uvec    indices(const_cast<uword*>(&(B.row_indices[col_offset_1])), col_offset_delta, false, false);
            const Col<eT>   B_col(const_cast<   eT*>(&(     B.values[col_offset_1])), col_offset_delta, false, false);
            
            out.col(i) = A.cols(indices) * B_col;
            }
          }
        #endif
        }
      else
        {
        arma_extra_debug_print("using standard multiplication");
        
        typename SpMat<eT>::const_iterator B_it     = B.begin();
        typename SpMat<eT>::const_iterator B_it_end = B.end();
        
        const uword out_n_rows = out.n_rows;
        
        while(B_it != B_it_end)
          {
          const eT    B_it_val = (*B_it);
          const uword B_it_col = B_it.col();
          const uword B_it_row = B_it.row();
          
          eT* out_col = out.colptr(B_it_col);
          
          for(uword row = 0; row < out_n_rows; ++row)
            {
            out_col[row] += A.at(row, B_it_row) * B_it_val;
            }
          
          ++B_it;
          }
        }
      }
    }
  }



//



template<typename T1, typename T2>
inline
void
spglue_times_mixed::apply(SpMat<typename eT_promoter<T1,T2>::eT>& out, const mtSpGlue<typename eT_promoter<T1,T2>::eT, T1, T2, spglue_times_mixed>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename eT_promoter<T1,T2>::eT out_eT;
  
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
    
    out = AA * BB;
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
    
    out = AA * BB;
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
    
    out = AA * BB;
    }
  }



template<typename T1, typename T2>
inline
void
spglue_times_mixed::sparse_times_dense(Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >& out, const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  if( (is_same_type<eT1,out_eT>::no) && (is_same_type<eT2,out_eT>::yes) )
    {
    // upgrade T1
    
    const unwrap_spmat<T1> UA(X);
    const quasi_unwrap<T2> UB(Y);
    
    const SpMat<eT1>& A = UA.M;
    const   Mat<eT2>& B = UB.M;
    
    SpMat<out_eT> AA(arma_layout_indicator(), A);
    
    for(uword i=0; i < A.n_nonzero; ++i)  { access::rw(AA.values[i]) = out_eT(A.values[i]); }
    
    const Mat<out_eT>& BB = reinterpret_cast< const Mat<out_eT>& >(B);
    
    spglue_times_misc::sparse_times_dense(out, AA, BB);
    }
  else
  if( (is_same_type<eT1,out_eT>::yes) && (is_same_type<eT2,out_eT>::no) )
    {
    // upgrade T2 
    
    const unwrap_spmat<T1> UA(X);
    const quasi_unwrap<T2> UB(Y);
    
    const SpMat<eT1>& A = UA.M;
    const   Mat<eT2>& B = UB.M;
    
    const SpMat<out_eT>& AA = reinterpret_cast< const SpMat<out_eT>& >(A);
    
    const Mat<out_eT> BB = conv_to< Mat<out_eT> >::from(B);
    
    spglue_times_misc::sparse_times_dense(out, AA, BB);
    }
  else
    {
    // upgrade T1 and T2
    
    const unwrap_spmat<T1> UA(X);
    const quasi_unwrap<T2> UB(Y);
    
    const SpMat<eT1>& A = UA.M;
    const   Mat<eT2>& B = UB.M;
    
    SpMat<out_eT> AA(arma_layout_indicator(), A);
    
    for(uword i=0; i < A.n_nonzero; ++i)  { access::rw(AA.values[i]) = out_eT(A.values[i]); }
    
    const Mat<out_eT> BB = conv_to< Mat<out_eT> >::from(B);
    
    spglue_times_misc::sparse_times_dense(out, AA, BB);
    }
  }



template<typename T1, typename T2>
inline
void
spglue_times_mixed::dense_times_sparse(Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >& out, const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  if( (is_same_type<eT1,out_eT>::no) && (is_same_type<eT2,out_eT>::yes) )
    {
    // upgrade T1
    
    const quasi_unwrap<T1> UA(X);
    const unwrap_spmat<T2> UB(Y);
    
    const   Mat<eT1>& A = UA.M;
    const SpMat<eT2>& B = UB.M;
    
    const Mat<out_eT> AA = conv_to< Mat<out_eT> >::from(A);
    
    const SpMat<out_eT>& BB = reinterpret_cast< const SpMat<out_eT>& >(B);
    
    spglue_times_misc::dense_times_sparse(out, AA, BB);
    }
  else
  if( (is_same_type<eT1,out_eT>::yes) && (is_same_type<eT2,out_eT>::no) )
    {
    // upgrade T2 
    
    const quasi_unwrap<T1> UA(X);
    const unwrap_spmat<T2> UB(Y);
    
    const   Mat<eT1>& A = UA.M;
    const SpMat<eT2>& B = UB.M;
    
    const Mat<out_eT>& AA = reinterpret_cast< const Mat<out_eT>& >(A);
    
    SpMat<out_eT> BB(arma_layout_indicator(), B);
    
    for(uword i=0; i < B.n_nonzero; ++i)  { access::rw(BB.values[i]) = out_eT(B.values[i]); }
    
    spglue_times_misc::dense_times_sparse(out, AA, BB);
    }
  else
    {
    // upgrade T1 and T2
    
    const quasi_unwrap<T1> UA(X);
    const unwrap_spmat<T2> UB(Y);
    
    const   Mat<eT1>& A = UA.M;
    const SpMat<eT2>& B = UB.M;
    
    const Mat<out_eT> AA = conv_to< Mat<out_eT> >::from(A);
    
    SpMat<out_eT> BB(arma_layout_indicator(), B);
    
    for(uword i=0; i < B.n_nonzero; ++i)  { access::rw(BB.values[i]) = out_eT(B.values[i]); }
    
    spglue_times_misc::dense_times_sparse(out, AA, BB);
    }
  }



//! @}
