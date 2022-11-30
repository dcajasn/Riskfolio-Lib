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


//! \addtogroup SpMat
//! @{


/**
 * Initialize a sparse matrix with size 0x0 (empty).
 */
template<typename eT>
inline
SpMat<eT>::SpMat()
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold(0,0);
  }



/**
 * Clean up the memory of a sparse matrix and destruct it.
 */
template<typename eT>
inline
SpMat<eT>::~SpMat()
  {
  arma_extra_debug_sigprint_this(this);
  
  if(values     )  { memory::release(access::rw(values));      }
  if(row_indices)  { memory::release(access::rw(row_indices)); }
  if(col_ptrs   )  { memory::release(access::rw(col_ptrs));    }
  }



/**
 * Constructor with size given.
 */
template<typename eT>
inline
SpMat<eT>::SpMat(const uword in_rows, const uword in_cols)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold(in_rows, in_cols);
  }



template<typename eT>
inline
SpMat<eT>::SpMat(const SizeMat& s)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold(s.n_rows, s.n_cols);
  }



template<typename eT>
inline
SpMat<eT>::SpMat(const arma_reserve_indicator&, const uword in_rows, const uword in_cols, const uword new_n_nonzero)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold(in_rows, in_cols, new_n_nonzero);
  }



template<typename eT>
template<typename eT2>
inline
SpMat<eT>::SpMat(const arma_layout_indicator&, const SpMat<eT2>& x)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold(x.n_rows, x.n_cols, x.n_nonzero);
  
  if(x.n_nonzero == 0)  { return; }
  
  if(x.row_indices)  { arrayops::copy(access::rwp(row_indices), x.row_indices, x.n_nonzero + 1); }
  if(x.col_ptrs   )  { arrayops::copy(access::rwp(col_ptrs),    x.col_ptrs,    x.n_cols    + 1); }
  
  // NOTE: 'values' array is not initialised
  }



/**
 * Assemble from text.
 */
template<typename eT>
inline
SpMat<eT>::SpMat(const char* text)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(std::string(text));
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();
  
  init(std::string(text));
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>::SpMat(const std::string& text)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint();
  
  init(text);
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator=(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  init(text);
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>::SpMat(const SpMat<eT>& x)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(x);
  }



template<typename eT>
inline
SpMat<eT>::SpMat(SpMat<eT>&& in_mat)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  arma_extra_debug_sigprint(arma_str::format("this = %x   in_mat = %x") % this % &in_mat);
  
  (*this).steal_mem(in_mat);
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator=(SpMat<eT>&& in_mat)
  {
  arma_extra_debug_sigprint(arma_str::format("this = %x   in_mat = %x") % this % &in_mat);
  
  (*this).steal_mem(in_mat);
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>::SpMat(const MapMat<eT>& x)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(x);
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator=(const MapMat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  init(x);
  
  return *this;
  }



//! Insert a large number of values at once.
//! locations.row[0] should be row indices, locations.row[1] should be column indices,
//! and values should be the corresponding values.
//! If sort_locations is false, then it is assumed that the locations and values
//! are already sorted in column-major ordering.
template<typename eT>
template<typename T1, typename T2>
inline
SpMat<eT>::SpMat(const Base<uword,T1>& locations_expr, const Base<eT,T2>& vals_expr, const bool sort_locations)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  const unwrap<T1> locs_tmp( locations_expr.get_ref() );
  const unwrap<T2> vals_tmp(      vals_expr.get_ref() );
  
  const Mat<uword>& locs = locs_tmp.M;
  const Mat<eT>&    vals = vals_tmp.M;
  
  arma_debug_check( (vals.is_vec() == false),     "SpMat::SpMat(): given 'values' object must be a vector"                 );
  arma_debug_check( (locs.n_rows != 2),           "SpMat::SpMat(): locations matrix must have two rows"                    );
  arma_debug_check( (locs.n_cols != vals.n_elem), "SpMat::SpMat(): number of locations is different than number of values" );

  // If there are no elements in the list, max() will fail.
  if(locs.n_cols == 0)  { init_cold(0, 0); return; }
  
  // Automatically determine size before pruning zeros.
  uvec bounds = arma::max(locs, 1);
  init_cold(bounds[0] + 1, bounds[1] + 1);
  
  // Ensure that there are no zeros
  const uword N_old = vals.n_elem;
        uword N_new = 0;
  
  for(uword i=0; i < N_old; ++i)  { N_new += (vals[i] != eT(0)) ? uword(1) : uword(0); }
  
  if(N_new != N_old)
    {
    Col<eT>    filtered_vals(   N_new, arma_nozeros_indicator());
    Mat<uword> filtered_locs(2, N_new, arma_nozeros_indicator());
    
    uword index = 0;
    for(uword i = 0; i < N_old; ++i)
      {
      if(vals[i] != eT(0))
        {
        filtered_vals[index] = vals[i];
        
        filtered_locs.at(0, index) = locs.at(0, i);
        filtered_locs.at(1, index) = locs.at(1, i);
        
        ++index;
        }
      }
    
    init_batch_std(filtered_locs, filtered_vals, sort_locations);
    }
  else
    {
    init_batch_std(locs, vals, sort_locations);
    }
  }



//! Insert a large number of values at once.
//! locations.row[0] should be row indices, locations.row[1] should be column indices,
//! and values should be the corresponding values.
//! If sort_locations is false, then it is assumed that the locations and values
//! are already sorted in column-major ordering.
//! In this constructor the size is explicitly given.
template<typename eT>
template<typename T1, typename T2>
inline
SpMat<eT>::SpMat(const Base<uword,T1>& locations_expr, const Base<eT,T2>& vals_expr, const uword in_n_rows, const uword in_n_cols, const bool sort_locations, const bool check_for_zeros)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  const unwrap<T1> locs_tmp( locations_expr.get_ref() );
  const unwrap<T2> vals_tmp(      vals_expr.get_ref() );
  
  const Mat<uword>& locs = locs_tmp.M;
  const Mat<eT>&    vals = vals_tmp.M;
  
  arma_debug_check( (vals.is_vec() == false),     "SpMat::SpMat(): given 'values' object must be a vector"                 );
  arma_debug_check( (locs.n_rows != 2),           "SpMat::SpMat(): locations matrix must have two rows"                    );
  arma_debug_check( (locs.n_cols != vals.n_elem), "SpMat::SpMat(): number of locations is different than number of values" );
  
  init_cold(in_n_rows, in_n_cols);
  
  // Ensure that there are no zeros, unless the user asked not to.
  if(check_for_zeros)
    {
    const uword N_old = vals.n_elem;
          uword N_new = 0;
    
    for(uword i=0; i < N_old; ++i)  { N_new += (vals[i] != eT(0)) ? uword(1) : uword(0); }
    
    if(N_new != N_old)
      {
      Col<eT>    filtered_vals(   N_new, arma_nozeros_indicator());
      Mat<uword> filtered_locs(2, N_new, arma_nozeros_indicator());
      
      uword index = 0;
      for(uword i = 0; i < N_old; ++i)
        {
        if(vals[i] != eT(0))
          {
          filtered_vals[index] = vals[i];
          
          filtered_locs.at(0, index) = locs.at(0, i);
          filtered_locs.at(1, index) = locs.at(1, i);
          
          ++index;
          }
        }
      
      init_batch_std(filtered_locs, filtered_vals, sort_locations);
      }
    else
      {
      init_batch_std(locs, vals, sort_locations);
      }
    }
  else
    {
    init_batch_std(locs, vals, sort_locations);
    }
  }



template<typename eT>
template<typename T1, typename T2>
inline
SpMat<eT>::SpMat(const bool add_values, const Base<uword,T1>& locations_expr, const Base<eT,T2>& vals_expr, const uword in_n_rows, const uword in_n_cols, const bool sort_locations, const bool check_for_zeros)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  const unwrap<T1> locs_tmp( locations_expr.get_ref() );
  const unwrap<T2> vals_tmp(      vals_expr.get_ref() );
  
  const Mat<uword>& locs = locs_tmp.M;
  const Mat<eT>&    vals = vals_tmp.M;
  
  arma_debug_check( (vals.is_vec() == false),     "SpMat::SpMat(): given 'values' object must be a vector"                 );
  arma_debug_check( (locs.n_rows != 2),           "SpMat::SpMat(): locations matrix must have two rows"                    );
  arma_debug_check( (locs.n_cols != vals.n_elem), "SpMat::SpMat(): number of locations is different than number of values" );
  
  init_cold(in_n_rows, in_n_cols);
  
  // Ensure that there are no zeros, unless the user asked not to.
  if(check_for_zeros)
    {
    const uword N_old = vals.n_elem;
          uword N_new = 0;
    
    for(uword i=0; i < N_old; ++i)  { N_new += (vals[i] != eT(0)) ? uword(1) : uword(0); }
    
    if(N_new != N_old)
      {
      Col<eT>    filtered_vals(   N_new, arma_nozeros_indicator());
      Mat<uword> filtered_locs(2, N_new, arma_nozeros_indicator());
      
      uword index = 0;
      for(uword i = 0; i < N_old; ++i)
        {
        if(vals[i] != eT(0))
          {
          filtered_vals[index] = vals[i];
          
          filtered_locs.at(0, index) = locs.at(0, i);
          filtered_locs.at(1, index) = locs.at(1, i);
          
          ++index;
          }
        }
      
      add_values ? init_batch_add(filtered_locs, filtered_vals, sort_locations) : init_batch_std(filtered_locs, filtered_vals, sort_locations);
      }
    else
      {
      add_values ? init_batch_add(locs, vals, sort_locations) : init_batch_std(locs, vals, sort_locations);
      }
    }
  else
    {
    add_values ? init_batch_add(locs, vals, sort_locations) : init_batch_std(locs, vals, sort_locations);
    }
  }



//! Insert a large number of values at once.
//! Per CSC format, rowind_expr should be row indices, 
//! colptr_expr should column ptr indices locations,
//! and values should be the corresponding values.
//! In this constructor the size is explicitly given.
//! Values are assumed to be sorted, and the size 
//! information is trusted
template<typename eT>
template<typename T1, typename T2, typename T3>
inline
SpMat<eT>::SpMat
  (
  const Base<uword,T1>& rowind_expr, 
  const Base<uword,T2>& colptr_expr, 
  const Base<eT,   T3>& values_expr, 
  const uword           in_n_rows, 
  const uword           in_n_cols
  )
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  const unwrap<T1> rowind_tmp( rowind_expr.get_ref() );
  const unwrap<T2> colptr_tmp( colptr_expr.get_ref() );
  const unwrap<T3>   vals_tmp( values_expr.get_ref() );
  
  const Mat<uword>& rowind = rowind_tmp.M;
  const Mat<uword>& colptr = colptr_tmp.M;
  const Mat<eT>&      vals = vals_tmp.M;
  
  arma_debug_check( (rowind.is_vec() == false), "SpMat::SpMat(): given 'rowind' object must be a vector" );
  arma_debug_check( (colptr.is_vec() == false), "SpMat::SpMat(): given 'colptr' object must be a vector" );
  arma_debug_check( (vals.is_vec()   == false), "SpMat::SpMat(): given 'values' object must be a vector" );
  
  // Resize to correct number of elements (this also sets n_nonzero)
  init_cold(in_n_rows, in_n_cols, vals.n_elem);
  
  arma_debug_check( (rowind.n_elem != vals.n_elem), "SpMat::SpMat(): number of row indices is not equal to number of values" );
  arma_debug_check( (colptr.n_elem != (n_cols+1) ), "SpMat::SpMat(): number of column pointers is not equal to n_cols+1" );
  
  // copy supplied values into sparse matrix -- not checked for consistency
  arrayops::copy(access::rwp(row_indices), rowind.memptr(), rowind.n_elem );
  arrayops::copy(access::rwp(col_ptrs),    colptr.memptr(), colptr.n_elem );
  arrayops::copy(access::rwp(values),      vals.memptr(),   vals.n_elem   );
  
  // important: set the sentinel as well
  access::rw(col_ptrs[n_cols + 1]) = std::numeric_limits<uword>::max();
  
  // make sure no zeros are stored
  remove_zeros();
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(val != eT(0))
    {
    // Resize to 1x1 then set that to the right value.
    init(1, 1, 1); // Sets col_ptrs to 0.
    
    // Manually set element.
    access::rw(values[0])      = val;
    access::rw(row_indices[0]) = 0;
    access::rw(col_ptrs[1])    = 1;
    }
  else
    {
    init(0, 0);
    }
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator*=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  if(val != eT(0))
    {
    sync_csc();
    invalidate_cache();
    
    const uword n_nz = n_nonzero;
    
    eT* vals = access::rwp(values);
    
    bool has_zero = false;
    
    for(uword i=0; i<n_nz; ++i)
      {
      eT& vals_i = vals[i];
      
      vals_i *= val;
      
      if(vals_i == eT(0))  { has_zero = true; }
      }
    
    if(has_zero)  { remove_zeros(); }
    }
  else
    {
    (*this).zeros();
    }
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator/=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (val == eT(0)), "element-wise division: division by zero" );
  
  sync_csc();
  invalidate_cache();
  
  const uword n_nz = n_nonzero;
  
  eT* vals = access::rwp(values);
  
  bool has_zero = false;
  
  for(uword i=0; i<n_nz; ++i)
    {
    eT& vals_i = vals[i];
    
    vals_i /= val;
    
    if(vals_i == eT(0))  { has_zero = true; }
    }
  
  if(has_zero)  { remove_zeros(); }
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator=(const SpMat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  init(x);
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator+=(const SpMat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  SpMat<eT> out = (*this) + x;
  
  steal_mem(out);
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator-=(const SpMat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  SpMat<eT> out = (*this) - x;
  
  steal_mem(out);
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator*=(const SpMat<eT>& y)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  SpMat<eT> z = (*this) * y;
  
  steal_mem(z);
  
  return *this;
  }



// This is in-place element-wise matrix multiplication.
template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator%=(const SpMat<eT>& y)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  SpMat<eT> z = (*this) % y;
  
  steal_mem(z);
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator/=(const SpMat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: use of this function is not advised; it is implemented only for completeness
  
  arma_debug_assert_same_size(n_rows, n_cols, x.n_rows, x.n_cols, "element-wise division");
  
  for(uword c = 0; c < n_cols; ++c)
  for(uword r = 0; r < n_rows; ++r)
    {
    at(r, c) /= x.at(r, c);
    }
  
  return *this;
  }



template<typename eT>
template<typename T1, typename op_type>
inline
SpMat<eT>::SpMat(const SpToDOp<T1, op_type>& expr)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  typedef typename T1::elem_type T;
  
  // Make sure the type is compatible.
  arma_type_check(( is_same_type< eT, T >::no ));
  
  op_type::apply(*this, expr);
  }



// Construct a complex matrix out of two non-complex matrices
template<typename eT>
template<typename T1, typename T2>
inline
SpMat<eT>::SpMat
  (
  const SpBase<typename SpMat<eT>::pod_type, T1>& A,
  const SpBase<typename SpMat<eT>::pod_type, T2>& B
  )
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type T;
  
  // Make sure eT is complex and T is not (compile-time check).
  arma_type_check(( is_cx<eT>::no  ));
  arma_type_check(( is_cx< T>::yes ));
  
  // Compile-time abort if types are not compatible.
  arma_type_check(( is_same_type< std::complex<T>, eT >::no ));
  
  const unwrap_spmat<T1> tmp1(A.get_ref());
  const unwrap_spmat<T2> tmp2(B.get_ref());
  
  const SpMat<T>& X = tmp1.M;
  const SpMat<T>& Y = tmp2.M;
  
  arma_debug_assert_same_size(X.n_rows, X.n_cols, Y.n_rows, Y.n_cols, "SpMat()");
  
  const uword l_n_rows = X.n_rows;
  const uword l_n_cols = X.n_cols;
  
  // Set size of matrix correctly.
  init_cold(l_n_rows, l_n_cols, n_unique(X, Y, op_n_unique_count()));
  
  // Now on a second iteration, fill it.
  typename SpMat<T>::const_iterator x_it  = X.begin();
  typename SpMat<T>::const_iterator x_end = X.end();
  
  typename SpMat<T>::const_iterator y_it  = Y.begin();
  typename SpMat<T>::const_iterator y_end = Y.end();
  
  uword cur_pos = 0;
  
  while((x_it != x_end) || (y_it != y_end))
    {
    if(x_it == y_it) // if we are at the same place
      {
      access::rw(values[cur_pos]) = std::complex<T>((T) *x_it, (T) *y_it);
      access::rw(row_indices[cur_pos]) = x_it.row();
      ++access::rw(col_ptrs[x_it.col() + 1]);
      
      ++x_it;
      ++y_it;
      }
    else
      {
      if((x_it.col() < y_it.col()) || ((x_it.col() == y_it.col()) && (x_it.row() < y_it.row()))) // if y is closer to the end
        {
        access::rw(values[cur_pos]) = std::complex<T>((T) *x_it, T(0));
        access::rw(row_indices[cur_pos]) = x_it.row();
        ++access::rw(col_ptrs[x_it.col() + 1]);
        
        ++x_it;
        }
      else // x is closer to the end
        {
        access::rw(values[cur_pos]) = std::complex<T>(T(0), (T) *y_it);
        access::rw(row_indices[cur_pos]) = y_it.row();
        ++access::rw(col_ptrs[y_it.col() + 1]);
        
        ++y_it;
        }
      }
    
    ++cur_pos;
    }
  
  // Now fix the column pointers; they are supposed to be a sum.
  for(uword c = 1; c <= n_cols; ++c)
    {
    access::rw(col_ptrs[c]) += col_ptrs[c - 1];
    }
  
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>::SpMat(const Base<eT, T1>& x)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  (*this).operator=(x);
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator=(const Base<eT, T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  if(is_same_type< T1, Gen<Mat<eT>, gen_zeros> >::yes)
    {
    const Proxy<T1> P(expr.get_ref());
    
    (*this).zeros( P.get_n_rows(), P.get_n_cols() );
    
    return *this;
    }
  
  if(is_same_type< T1, Gen<Mat<eT>, gen_eye> >::yes)
    {
    const Proxy<T1> P(expr.get_ref());
    
    (*this).eye( P.get_n_rows(), P.get_n_cols() );
    
    return *this;
    }
  
  const quasi_unwrap<T1> tmp(expr.get_ref());
  const Mat<eT>& x     = tmp.M;
  
  const uword x_n_rows = x.n_rows;
  const uword x_n_cols = x.n_cols;
  const uword x_n_elem = x.n_elem;
  
  // Count number of nonzero elements in base object.
  uword n = 0;
  
  const eT* x_mem = x.memptr();
  
  for(uword i=0; i < x_n_elem; ++i)  { n += (x_mem[i] != eT(0)) ? uword(1) : uword(0); }
  
  init(x_n_rows, x_n_cols, n);
  
  if(n == 0)  { return *this; }
  
  // Now the memory is resized correctly; set nonzero elements.
  n = 0;
  for(uword j = 0; j < x_n_cols; ++j)
  for(uword i = 0; i < x_n_rows; ++i)
    {
    const eT val = (*x_mem);  x_mem++;
    
    if(val != eT(0))
      {
      access::rw(values[n])      = val;
      access::rw(row_indices[n]) = i;
      access::rw(col_ptrs[j + 1])++;
      ++n;
      }
    }
  
  // Sum column counts to be column pointers.
  for(uword c = 1; c <= n_cols; ++c)
    {
    access::rw(col_ptrs[c]) += col_ptrs[c - 1];
    }
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator+=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  return (*this).operator=( (*this) + x.get_ref() );
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator-=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  return (*this).operator=( (*this) - x.get_ref() );
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator*=(const Base<eT, T1>& y)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  const Proxy<T1> p(y.get_ref());
  
  arma_debug_assert_mul_size(n_rows, n_cols, p.get_n_rows(), p.get_n_cols(), "matrix multiplication");
  
  // We assume the matrix structure is such that we will end up with a sparse
  // matrix.  Assuming that every entry in the dense matrix is nonzero (which is
  // a fairly valid assumption), each row with any nonzero elements in it (in this
  // matrix) implies an entire nonzero column.  Therefore, we iterate over all
  // the row_indices and count the number of rows with any elements in them
  // (using the quasi-linked-list idea from SYMBMM -- see spglue_times_meat.hpp).
  podarray<uword> index(n_rows);
  index.fill(n_rows); // Fill with invalid links.
  
  uword last_index = n_rows + 1;
  for(uword i = 0; i < n_nonzero; ++i)
    {
    if(index[row_indices[i]] == n_rows)
      {
      index[row_indices[i]] = last_index;
      last_index = row_indices[i];
      }
    }
  
  // Now count the number of rows which have nonzero elements.
  uword nonzero_rows = 0;
  while(last_index != n_rows + 1)
    {
    ++nonzero_rows;
    last_index = index[last_index];
    }
  
  SpMat<eT> z(arma_reserve_indicator(), n_rows, p.get_n_cols(), (nonzero_rows * p.get_n_cols())); // upper bound on size
  
  // Now we have to fill all the elements using a modification of the NUMBMM algorithm.
  uword cur_pos = 0;
  
  podarray<eT> partial_sums(n_rows);
  partial_sums.zeros();
  
  for(uword lcol = 0; lcol < n_cols; ++lcol)
    {
    const_iterator it     = begin();
    const_iterator it_end = end();
    
    while(it != it_end)
      {
      const eT value = (*it);
      
      partial_sums[it.row()] += (value * p.at(it.col(), lcol));
      
      ++it;
      }
    
    // Now add all partial sums to the matrix.
    for(uword i = 0; i < n_rows; ++i)
      {
      if(partial_sums[i] != eT(0))
        {
        access::rw(z.values[cur_pos]) = partial_sums[i];
        access::rw(z.row_indices[cur_pos]) = i;
        ++access::rw(z.col_ptrs[lcol + 1]);
        //printf("colptr %d now %d\n", lcol + 1, z.col_ptrs[lcol + 1]);
        ++cur_pos;
        partial_sums[i] = 0; // Would it be faster to do this in batch later?
        }
      }
    }
  
  // Now fix the column pointers.
  for(uword c = 1; c <= z.n_cols; ++c)
    {
    access::rw(z.col_ptrs[c]) += z.col_ptrs[c - 1];
    }
  
  // Resize to final correct size.
  z.mem_resize(z.col_ptrs[z.n_cols]);
  
  // Now take the memory of the temporary matrix.
  steal_mem(z);
  
  return *this;
  }



// NOTE: use of this function is not advised; it is implemented only for completeness
template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator/=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  SpMat<eT> tmp = (*this) / x.get_ref();
  
  steal_mem(tmp);
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator%=(const Base<eT, T1>& x)
  {
  arma_extra_debug_sigprint();
  
  SpMat<eT> tmp;
  
  // Just call the other order (these operations are commutative)
  // TODO: if there is a matrix size mismatch, the debug assert will print the matrix sizes in wrong order
  spglue_schur_misc::dense_schur_sparse(tmp, x.get_ref(), (*this));
  
  steal_mem(tmp);
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>::SpMat(const Op<T1, op_diagmat>& expr)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);

  (*this).operator=(expr);
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator=(const Op<T1, op_diagmat>& expr)
  {
  arma_extra_debug_sigprint();
  
  const diagmat_proxy<T1> P(expr.m);
  
  const uword max_n_nonzero = (std::min)(P.n_rows, P.n_cols);
  
  // resize memory to upper bound
  init(P.n_rows, P.n_cols, max_n_nonzero);
  
  uword count = 0;
  
  for(uword i=0; i < max_n_nonzero; ++i)
    {
    const eT val = P[i];
    
    if(val != eT(0))
      {
      access::rw(values[count])      = val;
      access::rw(row_indices[count]) = i;
      access::rw(col_ptrs[i + 1])++;
      ++count;
      }
    }
  
  // fix column pointers to be cumulative
  for(uword i = 1; i < n_cols + 1; ++i)
    {
    access::rw(col_ptrs[i]) += col_ptrs[i - 1];
    }
  
  // quick resize without reallocating memory and copying data
  access::rw(         n_nonzero) = count;
  access::rw(     values[count]) = eT(0);
  access::rw(row_indices[count]) = uword(0);
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator+=(const Op<T1, op_diagmat>& expr)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(expr);
  
  return (*this).operator+=(tmp);
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator-=(const Op<T1, op_diagmat>& expr)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(expr);
  
  return (*this).operator-=(tmp);
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator*=(const Op<T1, op_diagmat>& expr)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(expr);
  
  return (*this).operator*=(tmp);
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator/=(const Op<T1, op_diagmat>& expr)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(expr);
  
  return (*this).operator/=(tmp);
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator%=(const Op<T1, op_diagmat>& expr)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(expr);
  
  return (*this).operator%=(tmp);
  }



/**
 * Functions on subviews.
 */
template<typename eT>
inline
SpMat<eT>::SpMat(const SpSubview<eT>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  (*this).operator=(X);
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator=(const SpSubview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  if(X.n_nonzero == 0)  { zeros(X.n_rows, X.n_cols); return *this; }
  
  X.m.sync_csc();
  
  const bool alias = (this == &(X.m));
  
  if(alias)
    {
    SpMat<eT> tmp(X);
    
    steal_mem(tmp);
    }
  else
    {
    init(X.n_rows, X.n_cols, X.n_nonzero);
    
    if(X.n_rows == X.m.n_rows)
      {
      const uword sv_col_start = X.aux_col1;
      const uword sv_col_end   = X.aux_col1 + X.n_cols - 1;
      
      typename SpMat<eT>::const_col_iterator m_it     = X.m.begin_col_no_sync(sv_col_start);
      typename SpMat<eT>::const_col_iterator m_it_end = X.m.end_col_no_sync(sv_col_end);
      
      uword count = 0;
      
      while(m_it != m_it_end)
        {
        const uword m_it_col_adjusted = m_it.col() - sv_col_start;
        
        access::rw(row_indices[count]) = m_it.row();
        access::rw(values[count]) = (*m_it);
        ++access::rw(col_ptrs[m_it_col_adjusted + 1]);
        
        count++;
        
        ++m_it;
        }
      }
    else
      {
      typename SpSubview<eT>::const_iterator it     = X.begin();
      typename SpSubview<eT>::const_iterator it_end = X.end();
      
      while(it != it_end)
        {
        const uword it_pos = it.pos();
        
        access::rw(row_indices[it_pos]) = it.row();
        access::rw(values[it_pos]) = (*it);
        ++access::rw(col_ptrs[it.col() + 1]);
        ++it;
        }
      }
    
    // Now sum column pointers.
    for(uword c = 1; c <= n_cols; ++c)
      {
      access::rw(col_ptrs[c]) += col_ptrs[c - 1];
      }
    }
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator+=(const SpSubview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  SpMat<eT> tmp = (*this) + X;
  
  steal_mem(tmp);
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator-=(const SpSubview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  SpMat<eT> tmp = (*this) - X;
  
  steal_mem(tmp);
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator*=(const SpSubview<eT>& y)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  SpMat<eT> z = (*this) * y;
  
  steal_mem(z);
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator%=(const SpSubview<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  SpMat<eT> tmp = (*this) % x;
  
  steal_mem(tmp);
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator/=(const SpSubview<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(n_rows, n_cols, x.n_rows, x.n_cols, "element-wise division");
  
  // There is no pretty way to do this.
  for(uword elem = 0; elem < n_elem; elem++)
    {
    at(elem) /= x(elem);
    }
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>::SpMat(const SpSubview_col_list<eT,T1>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  SpSubview_col_list<eT,T1>::extract(*this, X);
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator=(const SpSubview_col_list<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool alias = (this == &(X.m));
  
  if(alias == false)
    {
    SpSubview_col_list<eT,T1>::extract(*this, X);
    }
  else
    {
    SpMat<eT> tmp(X);
    
    steal_mem(tmp);
    }
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator+=(const SpSubview_col_list<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  SpSubview_col_list<eT,T1>::plus_inplace(*this, X);
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator-=(const SpSubview_col_list<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  SpSubview_col_list<eT,T1>::minus_inplace(*this, X);
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator*=(const SpSubview_col_list<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  SpMat<eT> z = (*this) * X;
  
  steal_mem(z);
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator%=(const SpSubview_col_list<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  SpSubview_col_list<eT,T1>::schur_inplace(*this, X);
  
  return *this;
  }



template<typename eT>
template<typename T1>
inline
SpMat<eT>&
SpMat<eT>::operator/=(const SpSubview_col_list<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  SpSubview_col_list<eT,T1>::div_inplace(*this, X);
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>::SpMat(const spdiagview<eT>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  spdiagview<eT>::extract(*this, X);
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator=(const spdiagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  spdiagview<eT>::extract(*this, X);
  
  return *this;
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator+=(const spdiagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(X);
  
  return (*this).operator+=(tmp);
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator-=(const spdiagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(X);
  
  return (*this).operator-=(tmp);
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator*=(const spdiagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(X);
  
  return (*this).operator*=(tmp);
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator%=(const spdiagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(X);
  
  return (*this).operator%=(tmp);
  }



template<typename eT>
inline
SpMat<eT>&
SpMat<eT>::operator/=(const spdiagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT> tmp(X);
  
  return (*this).operator/=(tmp);
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>::SpMat(const SpOp<T1, spop_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr) // set in application of sparse operation
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  spop_type::apply(*this, X);
  
  sync_csc();          // in case apply() used element accessors
  invalidate_cache();  // in case apply() modified the CSC representation
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>&
SpMat<eT>::operator=(const SpOp<T1, spop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  spop_type::apply(*this, X);
  
  sync_csc();          // in case apply() used element accessors
  invalidate_cache();  // in case apply() modified the CSC representation
  
  return *this;
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>&
SpMat<eT>::operator+=(const SpOp<T1, spop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator+=(m);
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>&
SpMat<eT>::operator-=(const SpOp<T1, spop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator-=(m);
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>&
SpMat<eT>::operator*=(const SpOp<T1, spop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator*=(m);
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>&
SpMat<eT>::operator%=(const SpOp<T1, spop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator%=(m);
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>&
SpMat<eT>::operator/=(const SpOp<T1, spop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator/=(m);
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>::SpMat(const SpGlue<T1, T2, spglue_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  spglue_type::apply(*this, X);
  
  sync_csc();          // in case apply() used element accessors
  invalidate_cache();  // in case apply() modified the CSC representation
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>&
SpMat<eT>::operator=(const SpGlue<T1, T2, spglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  spglue_type::apply(*this, X);
  
  sync_csc();          // in case apply() used element accessors
  invalidate_cache();  // in case apply() modified the CSC representation
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>&
SpMat<eT>::operator+=(const SpGlue<T1, T2, spglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator+=(m);
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>&
SpMat<eT>::operator-=(const SpGlue<T1, T2, spglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator-=(m);
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>&
SpMat<eT>::operator*=(const SpGlue<T1, T2, spglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator*=(m);
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>&
SpMat<eT>::operator%=(const SpGlue<T1, T2, spglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator%=(m);
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>&
SpMat<eT>::operator/=(const SpGlue<T1, T2, spglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_same_type< eT, typename T1::elem_type >::no ));
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator/=(m);
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>::SpMat(const mtSpOp<eT, T1, spop_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  spop_type::apply(*this, X);
  
  sync_csc();          // in case apply() used element accessors
  invalidate_cache();  // in case apply() modified the CSC representation
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>&
SpMat<eT>::operator=(const mtSpOp<eT, T1, spop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  spop_type::apply(*this, X);
  
  sync_csc();          // in case apply() used element accessors
  invalidate_cache();  // in case apply() modified the CSC representation
  
  return *this;
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>&
SpMat<eT>::operator+=(const mtSpOp<eT, T1, spop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator+=(m);
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>&
SpMat<eT>::operator-=(const mtSpOp<eT, T1, spop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator-=(m);
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>&
SpMat<eT>::operator*=(const mtSpOp<eT, T1, spop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator*=(m);
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>&
SpMat<eT>::operator%=(const mtSpOp<eT, T1, spop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator%=(m);
  }



template<typename eT>
template<typename T1, typename spop_type>
inline
SpMat<eT>&
SpMat<eT>::operator/=(const mtSpOp<eT, T1, spop_type>& X)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator/=(m);
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>::SpMat(const mtSpGlue<eT, T1, T2, spglue_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(0)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  spglue_type::apply(*this, X);
  
  sync_csc();          // in case apply() used element accessors
  invalidate_cache();  // in case apply() modified the CSC representation
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>&
SpMat<eT>::operator=(const mtSpGlue<eT, T1, T2, spglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  spglue_type::apply(*this, X);
  
  sync_csc();          // in case apply() used element accessors
  invalidate_cache();  // in case apply() modified the CSC representation
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>&
SpMat<eT>::operator+=(const mtSpGlue<eT, T1, T2, spglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator+=(m);
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>&
SpMat<eT>::operator-=(const mtSpGlue<eT, T1, T2, spglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator-=(m);
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>&
SpMat<eT>::operator*=(const mtSpGlue<eT, T1, T2, spglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator*=(m);
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>&
SpMat<eT>::operator%=(const mtSpGlue<eT, T1, T2, spglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator%=(m);
  }



template<typename eT>
template<typename T1, typename T2, typename spglue_type>
inline
SpMat<eT>&
SpMat<eT>::operator/=(const mtSpGlue<eT, T1, T2, spglue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  const SpMat<eT> m(X);
  
  return (*this).operator/=(m);
  }



template<typename eT>
arma_inline
SpSubview_row<eT>
SpMat<eT>::row(const uword row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds(row_num >= n_rows, "SpMat::row(): out of bounds");
  
  return SpSubview_row<eT>(*this, row_num);
  }



template<typename eT>
arma_inline
const SpSubview_row<eT>
SpMat<eT>::row(const uword row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds(row_num >= n_rows, "SpMat::row(): out of bounds");
  
  return SpSubview_row<eT>(*this, row_num);
  }



template<typename eT>
inline
SpSubview_row<eT>
SpMat<eT>::operator()(const uword row_num, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  const bool col_all = col_span.whole;
  
  const uword local_n_cols = n_cols;
  
  const uword in_col1       = col_all ? 0            : col_span.a;
  const uword in_col2       =                          col_span.b;
  const uword submat_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;
  
  arma_debug_check_bounds
    (
    (row_num >= n_rows)
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) )
    ,
    "SpMat::operator(): indices out of bounds or incorrectly used"
    );
  
  return SpSubview_row<eT>(*this, row_num, in_col1, submat_n_cols);
  }



template<typename eT>
inline
const SpSubview_row<eT>
SpMat<eT>::operator()(const uword row_num, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  const bool col_all = col_span.whole;
  
  const uword local_n_cols = n_cols;
  
  const uword in_col1       = col_all ? 0            : col_span.a;
  const uword in_col2       =                          col_span.b;
  const uword submat_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1;
  
  arma_debug_check_bounds
    (
    (row_num >= n_rows)
    ||
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) )
    ,
    "SpMat::operator(): indices out of bounds or incorrectly used"
    );
  
  return SpSubview_row<eT>(*this, row_num, in_col1, submat_n_cols);
  }



template<typename eT>
arma_inline
SpSubview_col<eT>
SpMat<eT>::col(const uword col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds(col_num >= n_cols, "SpMat::col(): out of bounds");
  
  return SpSubview_col<eT>(*this, col_num);
  }



template<typename eT>
arma_inline
const SpSubview_col<eT>
SpMat<eT>::col(const uword col_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds(col_num >= n_cols, "SpMat::col(): out of bounds");
  
  return SpSubview_col<eT>(*this, col_num);
  }



template<typename eT>
inline
SpSubview_col<eT>
SpMat<eT>::operator()(const span& row_span, const uword col_num)
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;
  
  const uword local_n_rows = n_rows;
  
  const uword in_row1       = row_all ? 0            : row_span.a;
  const uword in_row2       =                          row_span.b;
  const uword submat_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;
  
  arma_debug_check_bounds
    (
    (col_num >= n_cols)
    ||
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) )
    ,
    "SpMat::operator(): indices out of bounds or incorrectly used"
    );
  
  return SpSubview_col<eT>(*this, col_num, in_row1, submat_n_rows);
  }



template<typename eT>
inline
const SpSubview_col<eT>
SpMat<eT>::operator()(const span& row_span, const uword col_num) const
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;
  
  const uword local_n_rows = n_rows;
  
  const uword in_row1       = row_all ? 0            : row_span.a;
  const uword in_row2       =                          row_span.b;
  const uword submat_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1;
  
  arma_debug_check_bounds
    (
    (col_num >= n_cols)
    ||
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) )
    ,
    "SpMat::operator(): indices out of bounds or incorrectly used"
    );
  
  return SpSubview_col<eT>(*this, col_num, in_row1, submat_n_rows);
  }



template<typename eT>
arma_inline
SpSubview<eT>
SpMat<eT>::rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "SpMat::rows(): indices out of bounds or incorrectly used"
    );
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  return SpSubview<eT>(*this, in_row1, 0, subview_n_rows, n_cols);
  }



template<typename eT>
arma_inline
const SpSubview<eT>
SpMat<eT>::rows(const uword in_row1, const uword in_row2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "SpMat::rows(): indices out of bounds or incorrectly used"
    );
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  
  return SpSubview<eT>(*this, in_row1, 0, subview_n_rows, n_cols);
  }



template<typename eT>
arma_inline
SpSubview<eT>
SpMat<eT>::cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "SpMat::cols(): indices out of bounds or incorrectly used"
    );
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  return SpSubview<eT>(*this, 0, in_col1, n_rows, subview_n_cols);
  }



template<typename eT>
arma_inline
const SpSubview<eT>
SpMat<eT>::cols(const uword in_col1, const uword in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "SpMat::cols(): indices out of bounds or incorrectly used"
    );
  
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  return SpSubview<eT>(*this, 0, in_col1, n_rows, subview_n_cols);
  }



template<typename eT>
arma_inline
SpSubview<eT>
SpMat<eT>::submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_row1 > in_row2) || (in_col1 >  in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "SpMat::submat(): indices out of bounds or incorrectly used"
    );
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  return SpSubview<eT>(*this, in_row1, in_col1, subview_n_rows, subview_n_cols);
  }



template<typename eT>
arma_inline
const SpSubview<eT>
SpMat<eT>::submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_row1 > in_row2) || (in_col1 >  in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "SpMat::submat(): indices out of bounds or incorrectly used"
    );
  
  const uword subview_n_rows = in_row2 - in_row1 + 1;
  const uword subview_n_cols = in_col2 - in_col1 + 1;
  
  return SpSubview<eT>(*this, in_row1, in_col1, subview_n_rows, subview_n_cols);
  }



template<typename eT>
arma_inline
SpSubview<eT>
SpMat<eT>::submat(const uword in_row1, const uword in_col1, const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  const uword l_n_rows = n_rows;
  const uword l_n_cols = n_cols;
  
  const uword s_n_rows = s.n_rows;
  const uword s_n_cols = s.n_cols;
  
  arma_debug_check_bounds
    (
    ((in_row1 >= l_n_rows) || (in_col1 >= l_n_cols) || ((in_row1 + s_n_rows) > l_n_rows) || ((in_col1 + s_n_cols) > l_n_cols)),
    "SpMat::submat(): indices or size out of bounds"
    );
  
  return SpSubview<eT>(*this, in_row1, in_col1, s_n_rows, s_n_cols);
  }



template<typename eT>
arma_inline
const SpSubview<eT>
SpMat<eT>::submat(const uword in_row1, const uword in_col1, const SizeMat& s) const
  {
  arma_extra_debug_sigprint();
  
  const uword l_n_rows = n_rows;
  const uword l_n_cols = n_cols;
  
  const uword s_n_rows = s.n_rows;
  const uword s_n_cols = s.n_cols;
  
  arma_debug_check_bounds
    (
    ((in_row1 >= l_n_rows) || (in_col1 >= l_n_cols) || ((in_row1 + s_n_rows) > l_n_rows) || ((in_col1 + s_n_cols) > l_n_cols)),
    "SpMat::submat(): indices or size out of bounds"
    );
  
  return SpSubview<eT>(*this, in_row1, in_col1, s_n_rows, s_n_cols);
  }



template<typename eT>
inline
SpSubview<eT>
SpMat<eT>::submat(const span& row_span, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;
  const bool col_all = col_span.whole;
  
  const uword local_n_rows = n_rows;
  const uword local_n_cols = n_cols;
  
  const uword in_row1       = row_all ? 0            : row_span.a;
  const uword in_row2       =                          row_span.b;
  const uword submat_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1; 
  
  const uword in_col1       = col_all ? 0            : col_span.a;
  const uword in_col2       =                          col_span.b;
  const uword submat_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1; 
  
  arma_debug_check_bounds
    (    
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) )
    ||   
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) )
    ,    
    "SpMat::submat(): indices out of bounds or incorrectly used"
    );   
  
  return SpSubview<eT>(*this, in_row1, in_col1, submat_n_rows, submat_n_cols);
  }



template<typename eT>
inline
const SpSubview<eT>
SpMat<eT>::submat(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  const bool row_all = row_span.whole;
  const bool col_all = col_span.whole;
  
  const uword local_n_rows = n_rows;
  const uword local_n_cols = n_cols;
  
  const uword in_row1       = row_all ? 0            : row_span.a;
  const uword in_row2       =                          row_span.b;
  const uword submat_n_rows = row_all ? local_n_rows : in_row2 - in_row1 + 1; 
  
  const uword in_col1       = col_all ? 0            : col_span.a;
  const uword in_col2       =                          col_span.b;
  const uword submat_n_cols = col_all ? local_n_cols : in_col2 - in_col1 + 1; 
  
  arma_debug_check_bounds
    (    
    ( row_all ? false : ((in_row1 > in_row2) || (in_row2 >= local_n_rows)) )
    ||   
    ( col_all ? false : ((in_col1 > in_col2) || (in_col2 >= local_n_cols)) )
    ,    
    "SpMat::submat(): indices out of bounds or incorrectly used"
    );   
  
  return SpSubview<eT>(*this, in_row1, in_col1, submat_n_rows, submat_n_cols);
  }



template<typename eT>
inline
SpSubview<eT>
SpMat<eT>::operator()(const span& row_span, const span& col_span)
  {
  arma_extra_debug_sigprint();
  
  return submat(row_span, col_span);
  }



template<typename eT>
inline
const SpSubview<eT>
SpMat<eT>::operator()(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  return submat(row_span, col_span);
  }



template<typename eT>
arma_inline
SpSubview<eT>
SpMat<eT>::operator()(const uword in_row1, const uword in_col1, const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  return (*this).submat(in_row1, in_col1, s);
  }



template<typename eT>
arma_inline
const SpSubview<eT>
SpMat<eT>::operator()(const uword in_row1, const uword in_col1, const SizeMat& s) const
  {
  arma_extra_debug_sigprint();
  
  return (*this).submat(in_row1, in_col1, s);
  }



template<typename eT>
inline
SpSubview<eT>
SpMat<eT>::head_rows(const uword N)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds( (N > n_rows), "SpMat::head_rows(): size out of bounds" );
  
  return SpSubview<eT>(*this, 0, 0, N, n_cols);
  }



template<typename eT>
inline
const SpSubview<eT>
SpMat<eT>::head_rows(const uword N) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds( (N > n_rows), "SpMat::head_rows(): size out of bounds" );
  
  return SpSubview<eT>(*this, 0, 0, N, n_cols);
  }



template<typename eT>
inline
SpSubview<eT>
SpMat<eT>::tail_rows(const uword N)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds( (N > n_rows), "SpMat::tail_rows(): size out of bounds" );
  
  const uword start_row = n_rows - N;
  
  return SpSubview<eT>(*this, start_row, 0, N, n_cols);
  }



template<typename eT>
inline
const SpSubview<eT>
SpMat<eT>::tail_rows(const uword N) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds( (N > n_rows), "SpMat::tail_rows(): size out of bounds" );
  
  const uword start_row = n_rows - N;
  
  return SpSubview<eT>(*this, start_row, 0, N, n_cols);
  }



template<typename eT>
inline
SpSubview<eT>
SpMat<eT>::head_cols(const uword N)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds( (N > n_cols), "SpMat::head_cols(): size out of bounds" );
  
  return SpSubview<eT>(*this, 0, 0, n_rows, N);
  }



template<typename eT>
inline
const SpSubview<eT>
SpMat<eT>::head_cols(const uword N) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds( (N > n_cols), "SpMat::head_cols(): size out of bounds" );
  
  return SpSubview<eT>(*this, 0, 0, n_rows, N);
  }



template<typename eT>
inline
SpSubview<eT>
SpMat<eT>::tail_cols(const uword N)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds( (N > n_cols), "SpMat::tail_cols(): size out of bounds" );
  
  const uword start_col = n_cols - N;
  
  return SpSubview<eT>(*this, 0, start_col, n_rows, N);
  }



template<typename eT>
inline
const SpSubview<eT>
SpMat<eT>::tail_cols(const uword N) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds( (N > n_cols), "SpMat::tail_cols(): size out of bounds" );
  
  const uword start_col = n_cols - N;
  
  return SpSubview<eT>(*this, 0, start_col, n_rows, N);
  }



template<typename eT>
template<typename T1>
arma_inline
SpSubview_col_list<eT, T1>
SpMat<eT>::cols(const Base<uword, T1>& indices)
  {
  arma_extra_debug_sigprint();
  
  return SpSubview_col_list<eT, T1>(*this, indices);
  }



template<typename eT>
template<typename T1>
arma_inline
const SpSubview_col_list<eT, T1>
SpMat<eT>::cols(const Base<uword, T1>& indices) const
  {
  arma_extra_debug_sigprint();
  
  return SpSubview_col_list<eT, T1>(*this, indices);
  }



//! creation of spdiagview (diagonal)
template<typename eT>
inline
spdiagview<eT>
SpMat<eT>::diag(const sword in_id)
  {
  arma_extra_debug_sigprint();
  
  const uword row_offset = (in_id < 0) ? uword(-in_id) : 0;
  const uword col_offset = (in_id > 0) ? uword( in_id) : 0;
  
  arma_debug_check_bounds
    (
    ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)),
    "SpMat::diag(): requested diagonal out of bounds"
    );
  
  const uword len = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  return spdiagview<eT>(*this, row_offset, col_offset, len);
  }



//! creation of spdiagview (diagonal)
template<typename eT>
inline
const spdiagview<eT>
SpMat<eT>::diag(const sword in_id) const
  {
  arma_extra_debug_sigprint();
  
  const uword row_offset = uword( (in_id < 0) ? -in_id : 0 );
  const uword col_offset = uword( (in_id > 0) ?  in_id : 0 );
  
  arma_debug_check_bounds
    (
    ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)),
    "SpMat::diag(): requested diagonal out of bounds"
    );
  
  const uword len = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  return spdiagview<eT>(*this, row_offset, col_offset, len);
  }



template<typename eT>
inline
void
SpMat<eT>::swap_rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds( ((in_row1 >= n_rows) || (in_row2 >= n_rows)), "SpMat::swap_rows(): out of bounds" );
  
  if(in_row1 == in_row2)  { return; }
  
  sync_csc();
  invalidate_cache();
  
  // The easier way to do this, instead of collecting all the elements in one row and then swapping with the other, will be
  // to iterate over each column of the matrix (since we store in column-major format) and then swap the two elements in the two rows at that time.
  // We will try to avoid using the at() call since it is expensive, instead preferring to use an iterator to track our position.
  uword col1 = (in_row1 < in_row2) ? in_row1 : in_row2;
  uword col2 = (in_row1 < in_row2) ? in_row2 : in_row1;
  
  for(uword lcol = 0; lcol < n_cols; lcol++)
    {
    // If there is nothing in this column we can ignore it.
    if(col_ptrs[lcol] == col_ptrs[lcol + 1])
      {
      continue;
      }
    
    // These will represent the positions of the items themselves.
    uword loc1 = n_nonzero + 1;
    uword loc2 = n_nonzero + 1;
    
    for(uword search_pos = col_ptrs[lcol]; search_pos < col_ptrs[lcol + 1]; search_pos++)
      {
      if(row_indices[search_pos] == col1)
        {
        loc1 = search_pos;
        }
      
      if(row_indices[search_pos] == col2)
        {
        loc2 = search_pos;
        break; // No need to look any further.
        }
      }
    
    // There are four cases: we found both elements; we found one element (loc1); we found one element (loc2); we found zero elements.
    // If we found zero elements no work needs to be done and we can continue to the next column.
    if((loc1 != (n_nonzero + 1)) && (loc2 != (n_nonzero + 1)))
      {
      // This is an easy case: just swap the values.  No index modifying necessary.
      eT tmp = values[loc1];
      access::rw(values[loc1]) = values[loc2];
      access::rw(values[loc2]) = tmp;
      }
    else if(loc1 != (n_nonzero + 1)) // We only found loc1 and not loc2.
      {
      // We need to find the correct place to move our value to.  It will be forward (not backwards) because in_row2 > in_row1.
      // Each iteration of the loop swaps the current value (loc1) with (loc1 + 1); in this manner we move our value down to where it should be.
      while(((loc1 + 1) < col_ptrs[lcol + 1]) && (row_indices[loc1 + 1] < in_row2))
        {
        // Swap both the values and the indices.  The column should not change.
        eT tmp = values[loc1];
        access::rw(values[loc1]) = values[loc1 + 1];
        access::rw(values[loc1 + 1]) = tmp;
        
        uword tmp_index = row_indices[loc1];
        access::rw(row_indices[loc1]) = row_indices[loc1 + 1];
        access::rw(row_indices[loc1 + 1]) = tmp_index;
        
        loc1++; // And increment the counter.
        }
      
      // Now set the row index correctly.
      access::rw(row_indices[loc1]) = in_row2;
      
      }
    else if(loc2 != (n_nonzero + 1))
      {
      // We need to find the correct place to move our value to.  It will be backwards (not forwards) because in_row1 < in_row2.
      // Each iteration of the loop swaps the current value (loc2) with (loc2 - 1); in this manner we move our value up to where it should be.
      while(((loc2 - 1) >= col_ptrs[lcol]) && (row_indices[loc2 - 1] > in_row1))
        {
        // Swap both the values and the indices.  The column should not change.
        eT tmp = values[loc2];
        access::rw(values[loc2]) = values[loc2 - 1];
        access::rw(values[loc2 - 1]) = tmp;
        
        uword tmp_index = row_indices[loc2];
        access::rw(row_indices[loc2]) = row_indices[loc2 - 1];
        access::rw(row_indices[loc2 - 1]) = tmp_index;
        
        loc2--; // And decrement the counter.
        }
      
      // Now set the row index correctly.
      access::rw(row_indices[loc2]) = in_row1;
      
      }
    /* else: no need to swap anything; both values are zero */
    }
  }



template<typename eT>
inline
void
SpMat<eT>::swap_cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds( ((in_col1 >= n_cols) || (in_col2 >= n_cols)), "SpMat::swap_cols(): out of bounds" );
  
  if(in_col1 == in_col2)  { return; }
  
  // TODO: this is a rudimentary implementation
  
  SpMat<eT> tmp = (*this);
  
  tmp.col(in_col1) = (*this).col(in_col2);
  tmp.col(in_col2) = (*this).col(in_col1);
  
  steal_mem(tmp);
  
  // for(uword lrow = 0; lrow < n_rows; ++lrow)
  //   {
  //   const eT tmp = at(lrow, in_col1);
  //   at(lrow, in_col1) = eT( at(lrow, in_col2) );
  //   at(lrow, in_col2) = tmp;
  //   }
  }



template<typename eT>
inline
void
SpMat<eT>::shed_row(const uword row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds(row_num >= n_rows, "SpMat::shed_row(): out of bounds");
  
  shed_rows (row_num, row_num);
  }



template<typename eT>
inline
void
SpMat<eT>::shed_col(const uword col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds(col_num >= n_cols, "SpMat::shed_col(): out of bounds");
  
  shed_cols(col_num, col_num);
  }



template<typename eT>
inline
void
SpMat<eT>::shed_rows(const uword in_row1, const uword in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "SpMat::shed_rows(): indices out of bounds or incorectly used"
    );
  
  sync_csc();
  
  SpMat<eT> newmat(n_rows - (in_row2 - in_row1 + 1), n_cols);
  
  // First, count the number of elements we will be removing.
  uword removing = 0;
  for(uword i = 0; i < n_nonzero; ++i)
    {
    const uword lrow = row_indices[i];
    if(lrow >= in_row1 && lrow <= in_row2)
      {
      ++removing;
      }
    }
  
  // Obtain counts of the number of points in each column and store them as the
  // (invalid) column pointers of the new matrix.
  for(uword i = 1; i < n_cols + 1; ++i)
    {
    access::rw(newmat.col_ptrs[i]) = col_ptrs[i] - col_ptrs[i - 1];
    }
  
  // Now initialize memory for the new matrix.
  newmat.mem_resize(n_nonzero - removing);
  
  // Now, copy over the elements.
  // i is the index in the old matrix; j is the index in the new matrix.
  const_iterator it     = begin();
  const_iterator it_end = end();
  
  uword j = 0; // The index in the new matrix.
  while(it != it_end)
    {
    const uword lrow = it.row();
    const uword lcol = it.col();
    
    if(lrow >= in_row1 && lrow <= in_row2)
      {
      // This element is being removed.  Subtract it from the column counts.
      --access::rw(newmat.col_ptrs[lcol + 1]);
      }
    else
      {
      // This element is being kept.  We may need to map the row index,
      // if it is past the section of rows we are removing.
      if(lrow > in_row2)
        {
        access::rw(newmat.row_indices[j]) = lrow - (in_row2 - in_row1 + 1);
        }
      else
        {
        access::rw(newmat.row_indices[j]) = lrow;
        }

      access::rw(newmat.values[j]) = (*it);
      ++j; // Increment index in new matrix.
      }
    
    ++it;
    }
  
  // Finally, sum the column counts so they are correct column pointers.
  for(uword i = 1; i < n_cols + 1; ++i)
    {
    access::rw(newmat.col_ptrs[i]) += newmat.col_ptrs[i - 1];
    }
  
  // Now steal the memory of the new matrix.
  steal_mem(newmat);
  }



template<typename eT>
inline
void
SpMat<eT>::shed_cols(const uword in_col1, const uword in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "SpMat::shed_cols(): indices out of bounds or incorrectly used"
    );
  
  sync_csc();
  invalidate_cache();
  
  // First we find the locations in values and row_indices for the column entries.
  uword col_beg = col_ptrs[in_col1];
  uword col_end = col_ptrs[in_col2 + 1];
  
  // Then we find the number of entries in the column.
  uword diff = col_end - col_beg;
  
  if(diff > 0)
    {
    eT*    new_values      = memory::acquire<eT>   (n_nonzero - diff);
    uword* new_row_indices = memory::acquire<uword>(n_nonzero - diff);
    
    // Copy first part.
    if(col_beg != 0)
      {
      arrayops::copy(new_values, values, col_beg);
      arrayops::copy(new_row_indices, row_indices, col_beg);
      }
    
    // Copy second part.
    if(col_end != n_nonzero)
      {
      arrayops::copy(new_values + col_beg, values + col_end, n_nonzero - col_end);
      arrayops::copy(new_row_indices + col_beg, row_indices + col_end, n_nonzero - col_end);
      }
    
    if(values)       { memory::release(access::rw(values));      }
    if(row_indices)  { memory::release(access::rw(row_indices)); }
    
    access::rw(values)      = new_values;
    access::rw(row_indices) = new_row_indices;
    
    // Update counts and such.
    access::rw(n_nonzero) -= diff;
    }
  
  // Update column pointers.
  const uword new_n_cols = n_cols - ((in_col2 - in_col1) + 1);
  
  uword* new_col_ptrs = memory::acquire<uword>(new_n_cols + 2);
  new_col_ptrs[new_n_cols + 1] = std::numeric_limits<uword>::max();
  
  // Copy first set of columns (no manipulation required).
  if(in_col1 != 0)
    {
    arrayops::copy(new_col_ptrs, col_ptrs, in_col1);
    }
  
  // Copy second set of columns (manipulation required).
  uword cur_col = in_col1;
  for(uword i = in_col2 + 1; i <= n_cols; ++i, ++cur_col)
    {
    new_col_ptrs[cur_col] = col_ptrs[i] - diff;
    }
  
  if(col_ptrs)  { memory::release(access::rw(col_ptrs)); }
  access::rw(col_ptrs) = new_col_ptrs;
  
  // We update the element and column counts, and we're done.
  access::rw(n_cols) = new_n_cols;
  access::rw(n_elem) = n_cols * n_rows;
  }



/**
 * Element access; acces the i'th element (works identically to the Mat accessors).
 * If there is nothing at element i, 0 is returned.
 */

template<typename eT>
arma_inline
arma_warn_unused
SpMat_MapMat_val<eT>
SpMat<eT>::operator[](const uword i)
  {
  const uword in_col = i / n_rows;
  const uword in_row = i % n_rows;
  
  return SpMat_MapMat_val<eT>((*this), cache, in_row, in_col);
  }



template<typename eT>
arma_inline
arma_warn_unused
eT
SpMat<eT>::operator[](const uword i) const
  {
  return get_value(i);
  }



template<typename eT>
arma_inline
arma_warn_unused
SpMat_MapMat_val<eT>
SpMat<eT>::at(const uword i)
  {
  const uword in_col = i / n_rows;
  const uword in_row = i % n_rows;
  
  return SpMat_MapMat_val<eT>((*this), cache, in_row, in_col);
  }



template<typename eT>
arma_inline
arma_warn_unused
eT
SpMat<eT>::at(const uword i) const
  {
  return get_value(i);
  }



template<typename eT>
arma_inline
arma_warn_unused
SpMat_MapMat_val<eT>
SpMat<eT>::operator()(const uword i)
  {
  arma_debug_check_bounds( (i >= n_elem), "SpMat::operator(): out of bounds" );
  
  const uword in_col = i / n_rows;
  const uword in_row = i % n_rows;
  
  return SpMat_MapMat_val<eT>((*this), cache, in_row, in_col);
  }



template<typename eT>
arma_inline
arma_warn_unused
eT
SpMat<eT>::operator()(const uword i) const
  {
  arma_debug_check_bounds( (i >= n_elem), "SpMat::operator(): out of bounds" );
  
  return get_value(i);
  }



/**
 * Element access; access the element at row in_rows and column in_col.
 * If there is nothing at that position, 0 is returned.
 */

#if defined(__cpp_multidimensional_subscript)
  
  template<typename eT>
  arma_inline
  arma_warn_unused
  SpMat_MapMat_val<eT>
  SpMat<eT>::operator[] (const uword in_row, const uword in_col)
    {
    return SpMat_MapMat_val<eT>((*this), cache, in_row, in_col);
    }
  
  
  
  template<typename eT>
  arma_inline
  arma_warn_unused
  eT
  SpMat<eT>::operator[] (const uword in_row, const uword in_col) const
    {
    return get_value(in_row, in_col);
    }
  
#endif



template<typename eT>
arma_inline
arma_warn_unused
SpMat_MapMat_val<eT>
SpMat<eT>::at(const uword in_row, const uword in_col)
  {
  return SpMat_MapMat_val<eT>((*this), cache, in_row, in_col);
  }



template<typename eT>
arma_inline
arma_warn_unused
eT
SpMat<eT>::at(const uword in_row, const uword in_col) const
  {
  return get_value(in_row, in_col);
  }



template<typename eT>
arma_inline
arma_warn_unused
SpMat_MapMat_val<eT>
SpMat<eT>::operator()(const uword in_row, const uword in_col)
  {
  arma_debug_check_bounds( ((in_row >= n_rows) || (in_col >= n_cols)), "SpMat::operator(): out of bounds" );
  
  return SpMat_MapMat_val<eT>((*this), cache, in_row, in_col);
  }



template<typename eT>
arma_inline
arma_warn_unused
eT
SpMat<eT>::operator()(const uword in_row, const uword in_col) const
  {
  arma_debug_check_bounds( ((in_row >= n_rows) || (in_col >= n_cols)), "SpMat::operator(): out of bounds" );
  
  return get_value(in_row, in_col);
  }



/**
 * Check if matrix is empty (no size, no values).
 */
template<typename eT>
arma_inline
arma_warn_unused
bool
SpMat<eT>::is_empty() const
  {
  return (n_elem == 0);
  }



//! returns true if the object can be interpreted as a column or row vector
template<typename eT>
arma_inline
arma_warn_unused
bool
SpMat<eT>::is_vec() const
  {
  return ( (n_rows == 1) || (n_cols == 1) );
  }



//! returns true if the object can be interpreted as a row vector
template<typename eT>
arma_inline
arma_warn_unused
bool
SpMat<eT>::is_rowvec() const
  {
  return (n_rows == 1);
  }



//! returns true if the object can be interpreted as a column vector
template<typename eT>
arma_inline
arma_warn_unused
bool
SpMat<eT>::is_colvec() const
  {
  return (n_cols == 1);
  }



//! returns true if the object has the same number of non-zero rows and columnns
template<typename eT>
arma_inline
arma_warn_unused
bool
SpMat<eT>::is_square() const
  {
  return (n_rows == n_cols);
  }



//! returns true if all of the elements are finite
template<typename eT>
inline
arma_warn_unused
bool
SpMat<eT>::is_finite() const
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  return arrayops::is_finite(values, n_nonzero);
  }



template<typename eT>
inline
arma_warn_unused
bool
SpMat<eT>::is_symmetric() const
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT>& A = (*this);
  
  if(A.n_rows != A.n_cols)  { return false; }
  
  const SpMat<eT> tmp = A - A.st();
  
  return (tmp.n_nonzero == uword(0));
  }



template<typename eT>
inline
arma_warn_unused
bool
SpMat<eT>::is_symmetric(const typename get_pod_type<elem_type>::result tol) const
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  if(tol == T(0))  { return (*this).is_symmetric(); }
  
  arma_debug_check( (tol < T(0)), "is_symmetric(): parameter 'tol' must be >= 0" );
  
  const SpMat<eT>& A = (*this);
  
  if(A.n_rows != A.n_cols)  { return false; }
  
  const T norm_A = as_scalar( arma::max(sum(abs(A), 1), 0) );
  
  if(norm_A == T(0))  { return true; }
  
  const T norm_A_Ast = as_scalar( arma::max(sum(abs(A - A.st()), 1), 0) );
  
  return ( (norm_A_Ast / norm_A) <= tol );
  }



template<typename eT>
inline
arma_warn_unused
bool
SpMat<eT>::is_hermitian() const
  {
  arma_extra_debug_sigprint();
  
  const SpMat<eT>& A = (*this);
  
  if(A.n_rows != A.n_cols)  { return false; }
  
  const SpMat<eT> tmp = A - A.t();
  
  return (tmp.n_nonzero == uword(0));
  }



template<typename eT>
inline
arma_warn_unused
bool
SpMat<eT>::is_hermitian(const typename get_pod_type<elem_type>::result tol) const
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  if(tol == T(0))  { return (*this).is_hermitian(); }
  
  arma_debug_check( (tol < T(0)), "is_hermitian(): parameter 'tol' must be >= 0" );
  
  const SpMat<eT>& A = (*this);
  
  if(A.n_rows != A.n_cols)  { return false; }
  
  const T norm_A = as_scalar( arma::max(sum(abs(A), 1), 0) );
  
  if(norm_A == T(0))  { return true; }
  
  const T norm_A_At = as_scalar( arma::max(sum(abs(A - A.t()), 1), 0) );
  
  return ( (norm_A_At / norm_A) <= tol );
  }



template<typename eT>
inline
arma_warn_unused
bool
SpMat<eT>::has_inf() const
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  return arrayops::has_inf(values, n_nonzero);
  }



template<typename eT>
inline
arma_warn_unused
bool
SpMat<eT>::has_nan() const
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  return arrayops::has_nan(values, n_nonzero);
  }



template<typename eT>
inline
arma_warn_unused
bool
SpMat<eT>::has_nonfinite() const
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  return (arrayops::is_finite(values, n_nonzero) == false);
  }



//! returns true if the given index is currently in range
template<typename eT>
arma_inline
arma_warn_unused
bool
SpMat<eT>::in_range(const uword i) const
  {
  return (i < n_elem);
  }


//! returns true if the given start and end indices are currently in range
template<typename eT>
arma_inline
arma_warn_unused
bool
SpMat<eT>::in_range(const span& x) const
  {
  arma_extra_debug_sigprint();
  
  if(x.whole)
    {
    return true;
    }
  else
    {
    const uword a = x.a;
    const uword b = x.b;
    
    return ( (a <= b) && (b < n_elem) );
    }
  }



//! returns true if the given location is currently in range
template<typename eT>
arma_inline
arma_warn_unused
bool
SpMat<eT>::in_range(const uword in_row, const uword in_col) const
  {
  return ( (in_row < n_rows) && (in_col < n_cols) );
  }



template<typename eT>
arma_inline
arma_warn_unused
bool
SpMat<eT>::in_range(const span& row_span, const uword in_col) const
  {
  arma_extra_debug_sigprint();
  
  if(row_span.whole)
    {
    return (in_col < n_cols);
    }
  else
    {
    const uword in_row1 = row_span.a;
    const uword in_row2 = row_span.b;
    
    return ( (in_row1 <= in_row2) && (in_row2 < n_rows) && (in_col < n_cols) );
    }
  }



template<typename eT>
arma_inline
arma_warn_unused
bool
SpMat<eT>::in_range(const uword in_row, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  if(col_span.whole)
    {
    return (in_row < n_rows);
    }
  else
    {
    const uword in_col1 = col_span.a;
    const uword in_col2 = col_span.b;
    
    return ( (in_row < n_rows) && (in_col1 <= in_col2) && (in_col2 < n_cols) );
    }
  }



template<typename eT>
arma_inline
arma_warn_unused
bool
SpMat<eT>::in_range(const span& row_span, const span& col_span) const
  {
  arma_extra_debug_sigprint();
  
  const uword in_row1 = row_span.a;
  const uword in_row2 = row_span.b;
  
  const uword in_col1 = col_span.a;
  const uword in_col2 = col_span.b;
  
  const bool rows_ok = row_span.whole ? true : ( (in_row1 <= in_row2) && (in_row2 < n_rows) );
  const bool cols_ok = col_span.whole ? true : ( (in_col1 <= in_col2) && (in_col2 < n_cols) );
  
  return ( rows_ok && cols_ok );
  }



template<typename eT>
arma_inline
arma_warn_unused
bool
SpMat<eT>::in_range(const uword in_row, const uword in_col, const SizeMat& s) const
  {
  const uword l_n_rows = n_rows;
  const uword l_n_cols = n_cols;
  
  if( (in_row >= l_n_rows) || (in_col >= l_n_cols) || ((in_row + s.n_rows) > l_n_rows) || ((in_col + s.n_cols) > l_n_cols) )
    {
    return false;
    }
  else
    {
    return true;
    }
  }



//! Set the size to the size of another matrix.
template<typename eT>
template<typename eT2>
inline
void
SpMat<eT>::copy_size(const SpMat<eT2>& m)
  {
  arma_extra_debug_sigprint();
  
  set_size(m.n_rows, m.n_cols);
  }



template<typename eT>
template<typename eT2>
inline
void
SpMat<eT>::copy_size(const Mat<eT2>& m)
  {
  arma_extra_debug_sigprint();
  
  set_size(m.n_rows, m.n_cols);
  }



template<typename eT>
inline
void
SpMat<eT>::set_size(const uword in_elem)
  {
  arma_extra_debug_sigprint();
  
  // If this is a row vector, we resize to a row vector.
  if(vec_state == 2)
    {
    set_size(1, in_elem);
    }
  else
    {
    set_size(in_elem, 1);
    }
  }



template<typename eT>
inline
void
SpMat<eT>::set_size(const uword in_rows, const uword in_cols)
  {
  arma_extra_debug_sigprint();
  
  invalidate_cache(); // placed here, as set_size() is used during matrix modification
  
  if( (n_rows == in_rows) && (n_cols == in_cols) )
    {
    return;
    }
  else
    {
    init(in_rows, in_cols);
    }
  }



template<typename eT>
inline
void
SpMat<eT>::set_size(const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  (*this).set_size(s.n_rows, s.n_cols);
  }



template<typename eT>
inline
void
SpMat<eT>::resize(const uword in_rows, const uword in_cols)
  {
  arma_extra_debug_sigprint();
  
  if( (n_rows == in_rows) && (n_cols == in_cols) )  { return; }
  
  if( (n_elem == 0) || (n_nonzero == 0) )
    {
    set_size(in_rows, in_cols);
    return;
    }
  
  SpMat<eT> tmp(in_rows, in_cols);
  
  if(tmp.n_elem > 0)
    {
    sync_csc();
    
    const uword last_row = (std::min)(in_rows, n_rows) - 1;
    const uword last_col = (std::min)(in_cols, n_cols) - 1;
    
    tmp.submat(0, 0, last_row, last_col) = (*this).submat(0, 0, last_row, last_col);
    }
  
  steal_mem(tmp);
  }



template<typename eT>
inline
void
SpMat<eT>::resize(const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  (*this).resize(s.n_rows, s.n_cols);
  }



template<typename eT>
inline
void
SpMat<eT>::reshape(const uword in_rows, const uword in_cols)
  {
  arma_extra_debug_sigprint();
  
  arma_check( ((in_rows*in_cols) != n_elem), "SpMat::reshape(): changing the number of elements in a sparse matrix is currently not supported" );
  
  if( (n_rows == in_rows) && (n_cols == in_cols) )  { return; }
  
  if(vec_state == 1)  { arma_debug_check( (in_cols != 1), "SpMat::reshape(): object is a column vector; requested size is not compatible" ); }
  if(vec_state == 2)  { arma_debug_check( (in_rows != 1), "SpMat::reshape(): object is a row vector; requested size is not compatible"    ); }
  
  if(n_nonzero == 0)
    {
    (*this).zeros(in_rows, in_cols);
    return;
    }
  
  if(in_cols == 1)
    {
    (*this).reshape_helper_intovec();
    }
  else
    {
    (*this).reshape_helper_generic(in_rows, in_cols);
    }
  }



template<typename eT>
inline
void
SpMat<eT>::reshape(const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  (*this).reshape(s.n_rows, s.n_cols);
  }



template<typename eT>
inline
void
SpMat<eT>::reshape_helper_generic(const uword in_rows, const uword in_cols)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  invalidate_cache();
  
  // We have to modify all of the relevant row indices and the relevant column pointers.
  // Iterate over all the points to do this.  We won't be deleting any points, but we will be modifying
  // columns and rows. We'll have to store a new set of column vectors.
  uword* new_col_ptrs    = memory::acquire<uword>(in_cols + 2);
  new_col_ptrs[in_cols + 1] = std::numeric_limits<uword>::max();
  
  uword* new_row_indices = memory::acquire<uword>(n_nonzero + 1);
  access::rw(new_row_indices[n_nonzero]) = 0;
  
  arrayops::fill_zeros(new_col_ptrs, in_cols + 1);
  
  const_iterator it     = begin();
  const_iterator it_end = end();
  
  for(; it != it_end; ++it)
    {
    uword vector_position = (it.col() * n_rows) + it.row();
    new_row_indices[it.pos()] = vector_position % in_rows;
    ++new_col_ptrs[vector_position / in_rows + 1];
    }
  
  // Now sum the column counts to get the new column pointers.
  for(uword i = 1; i <= in_cols; i++)
    {
    access::rw(new_col_ptrs[i]) += new_col_ptrs[i - 1];
    }
  
  // Copy the new row indices.
  if(row_indices)  { memory::release(access::rw(row_indices)); }
  if(col_ptrs)     { memory::release(access::rw(col_ptrs));    } 
  
  access::rw(row_indices) = new_row_indices;
  access::rw(col_ptrs)    = new_col_ptrs;
  
  // Now set the size.
  access::rw(n_rows) = in_rows;
  access::rw(n_cols) = in_cols;
  }



template<typename eT>
inline
void
SpMat<eT>::reshape_helper_intovec()
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  invalidate_cache();
  
  const_iterator it = begin();
  
  const uword t_n_rows    = n_rows;
  const uword t_n_nonzero = n_nonzero;
  
  for(uword i=0; i < t_n_nonzero; ++i)
    {
    const uword t_index = (it.col() * t_n_rows) + it.row();
    
    // ensure the iterator is pointing to the next element
    // before we overwrite the row index of the current element
    ++it;
    
    access::rw(row_indices[i]) = t_index;
    }
  
  access::rw(row_indices[n_nonzero]) = 0;
  
  access::rw(col_ptrs[0]) = 0;
  access::rw(col_ptrs[1]) = n_nonzero;
  access::rw(col_ptrs[2]) = std::numeric_limits<uword>::max();
  
  access::rw(n_rows) = (n_rows * n_cols);
  access::rw(n_cols) = 1;
  }



//! apply a functor to each non-zero element
template<typename eT>
template<typename functor>
inline
const SpMat<eT>&
SpMat<eT>::for_each(functor F)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  const uword N = (*this).n_nonzero;
  
  eT* rw_values = access::rwp(values);
  
  bool modified = false;
  bool has_zero = false;
  
  for(uword i=0; i < N; ++i)
    {
          eT& new_value = rw_values[i];
    const eT  old_value = new_value;
    
    F(new_value);
    
    if(new_value != old_value)  { modified = true; }
    if(new_value == eT(0)    )  { has_zero = true; }
    }
  
  if(modified)  { invalidate_cache(); }
  if(has_zero)  { remove_zeros();     }
  
  return *this;
  }



template<typename eT>
template<typename functor>
inline
const SpMat<eT>&
SpMat<eT>::for_each(functor F) const
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  const uword N = (*this).n_nonzero;
  
  for(uword i=0; i < N; ++i)
    {
    F(values[i]);
    }
  
  return *this;
  }



//! transform each non-zero element using a functor
template<typename eT>
template<typename functor>
inline
const SpMat<eT>&
SpMat<eT>::transform(functor F)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  invalidate_cache();
  
  const uword N = (*this).n_nonzero;
  
  eT* rw_values = access::rwp(values);
  
  bool has_zero = false;
  
  for(uword i=0; i < N; ++i)
    {
    eT& rw_values_i = rw_values[i];
    
    rw_values_i = eT( F(rw_values_i) );
    
    if(rw_values_i == eT(0))  { has_zero = true; }
    }
  
  if(has_zero)  { remove_zeros(); }
  
  return *this;
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::replace(const eT old_val, const eT new_val)
  {
  arma_extra_debug_sigprint();
  
  if(old_val == eT(0))
    {
    arma_debug_warn_level(1, "SpMat::replace(): replacement not done, as old_val = 0");
    }
  else
    {
    sync_csc();
    invalidate_cache();
    
    arrayops::replace(access::rwp(values), n_nonzero, old_val, new_val);
    
    if(new_val == eT(0))  { remove_zeros(); }
    }
  
  return *this;
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::clean(const typename get_pod_type<eT>::result threshold)
  {
  arma_extra_debug_sigprint();
  
  if(n_nonzero == 0)  { return *this; }
  
  sync_csc();
  invalidate_cache();
  
  arrayops::clean(access::rwp(values), n_nonzero, threshold);
  
  remove_zeros();
  
  return *this;
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::clamp(const eT min_val, const eT max_val)
  {
  arma_extra_debug_sigprint();
  
  if(is_cx<eT>::no)
    {
    arma_debug_check( (access::tmp_real(min_val) > access::tmp_real(max_val)), "SpMat::clamp(): min_val must be less than max_val" );
    }
  else
    {
    arma_debug_check( (access::tmp_real(min_val) > access::tmp_real(max_val)), "SpMat::clamp(): real(min_val) must be less than real(max_val)" );
    arma_debug_check( (access::tmp_imag(min_val) > access::tmp_imag(max_val)), "SpMat::clamp(): imag(min_val) must be less than imag(max_val)" );
    }
  
  if(n_nonzero == 0)  { return *this; }
  
  sync_csc();
  invalidate_cache();
  
  arrayops::clamp(access::rwp(values), n_nonzero, min_val, max_val);
  
  if( (min_val == eT(0)) || (max_val == eT(0)) )  { remove_zeros(); }
  
  return *this;
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  const bool already_done = ( (sync_state != 1) && (n_nonzero == 0) );
  
  if(already_done == false)
    {
    init(n_rows, n_cols);
    }
  
  return *this;
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::zeros(const uword in_elem)
  {
  arma_extra_debug_sigprint();
  
  if(vec_state == 2)
    {
    zeros(1, in_elem); // Row vector
    }
  else
    {
    zeros(in_elem, 1);
    }
  
  return *this;
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::zeros(const uword in_rows, const uword in_cols)
  {
  arma_extra_debug_sigprint();
  
  const bool already_done = ( (sync_state != 1) && (n_nonzero == 0) && (n_rows == in_rows) && (n_cols == in_cols) );
  
  if(already_done == false)
    {
    init(in_rows, in_cols);
    }
  
  return *this;
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::zeros(const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  return (*this).zeros(s.n_rows, s.n_cols);
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::eye()
  {
  arma_extra_debug_sigprint();
  
  return (*this).eye(n_rows, n_cols);
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::eye(const uword in_rows, const uword in_cols)
  {
  arma_extra_debug_sigprint();
  
  const uword N = (std::min)(in_rows, in_cols);
  
  init(in_rows, in_cols, N);
  
  arrayops::inplace_set(access::rwp(values), eT(1), N);
  
  for(uword i = 0; i <  N; ++i) { access::rw(row_indices[i]) = i; }
  
  for(uword i = 0; i <= N; ++i) { access::rw(col_ptrs[i])    = i; }
  
  // take into account non-square matrices
  for(uword i = (N+1); i <= in_cols; ++i)  { access::rw(col_ptrs[i]) = N; }
  
  access::rw(n_nonzero) = N;
  
  return *this;
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::eye(const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  return (*this).eye(s.n_rows, s.n_cols);
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::speye()
  {
  arma_extra_debug_sigprint();
  
  return (*this).eye(n_rows, n_cols);
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::speye(const uword in_n_rows, const uword in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  return (*this).eye(in_n_rows, in_n_cols);
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::speye(const SizeMat& s)
  {
  arma_extra_debug_sigprint();
  
  return (*this).eye(s.n_rows, s.n_cols);
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::sprandu(const uword in_rows, const uword in_cols, const double density)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (density < double(0)) || (density > double(1)) ), "sprandu(): density must be in the [0,1] interval" );
  
  const uword new_n_nonzero = uword(density * double(in_rows) * double(in_cols) + 0.5);
  
  init(in_rows, in_cols, new_n_nonzero);
  
  if(new_n_nonzero == 0)  { return *this; }
  
  arma_rng::randu<eT>::fill( access::rwp(values), new_n_nonzero );
  
  uvec indices = linspace<uvec>( 0u, in_rows*in_cols-1, new_n_nonzero );
  
  // perturb the indices
  for(uword i=1; i < new_n_nonzero-1; ++i)
    {
    const uword index_left  = indices[i-1];
    const uword index_right = indices[i+1];
    
    const uword center = (index_left + index_right) / 2;
    
    const uword delta1 = center      - index_left - 1;
    const uword delta2 = index_right - center     - 1;
    
    const uword min_delta = (std::min)(delta1, delta2);
    
    uword index_new = uword( double(center) + double(min_delta) * (2.0*randu()-1.0) );
    
    // paranoia, but better be safe than sorry
    if( (index_left < index_new) && (index_new < index_right) )
      {
      indices[i] = index_new;
      }
    }
  
  uword cur_index = 0;
  uword count     = 0;  
  
  for(uword lcol = 0; lcol < in_cols; ++lcol)
  for(uword lrow = 0; lrow < in_rows; ++lrow)
    {
    if(count == indices[cur_index])
      {
      access::rw(row_indices[cur_index]) = lrow;
      access::rw(col_ptrs[lcol + 1])++;
      ++cur_index;
      }
    
    ++count;
    }
  
  if(cur_index != new_n_nonzero)
    {
    // Fix size to correct size.
    mem_resize(cur_index);
    }
  
  // Sum column pointers.
  for(uword lcol = 1; lcol <= in_cols; ++lcol)
    {
    access::rw(col_ptrs[lcol]) += col_ptrs[lcol - 1];
    }
  
  return *this;
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::sprandu(const SizeMat& s, const double density)
  {
  arma_extra_debug_sigprint();
  
  return (*this).sprandu(s.n_rows, s.n_cols, density);
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::sprandn(const uword in_rows, const uword in_cols, const double density)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (density < double(0)) || (density > double(1)) ), "sprandn(): density must be in the [0,1] interval" );
  
  const uword new_n_nonzero = uword(density * double(in_rows) * double(in_cols) + 0.5);
  
  init(in_rows, in_cols, new_n_nonzero);
  
  if(new_n_nonzero == 0)  { return *this; }
  
  arma_rng::randn<eT>::fill( access::rwp(values), new_n_nonzero );
  
  uvec indices = linspace<uvec>( 0u, in_rows*in_cols-1, new_n_nonzero );
  
  // perturb the indices
  for(uword i=1; i < new_n_nonzero-1; ++i)
    {
    const uword index_left  = indices[i-1];
    const uword index_right = indices[i+1];
    
    const uword center = (index_left + index_right) / 2;
    
    const uword delta1 = center      - index_left - 1;
    const uword delta2 = index_right - center     - 1;
    
    const uword min_delta = (std::min)(delta1, delta2);
    
    uword index_new = uword( double(center) + double(min_delta) * (2.0*randu()-1.0) );
    
    // paranoia, but better be safe than sorry
    if( (index_left < index_new) && (index_new < index_right) )
      {
      indices[i] = index_new;
      }
    }
  
  uword cur_index = 0;
  uword count     = 0;  
  
  for(uword lcol = 0; lcol < in_cols; ++lcol)
  for(uword lrow = 0; lrow < in_rows; ++lrow)
    {
    if(count == indices[cur_index])
      {
      access::rw(row_indices[cur_index]) = lrow;
      access::rw(col_ptrs[lcol + 1])++;
      ++cur_index;
      }
    
    ++count;
    }
  
  if(cur_index != new_n_nonzero)
    {
    // Fix size to correct size.
    mem_resize(cur_index);
    }
  
  // Sum column pointers.
  for(uword lcol = 1; lcol <= in_cols; ++lcol)
    {
    access::rw(col_ptrs[lcol]) += col_ptrs[lcol - 1];
    }
  
  return *this;
  }



template<typename eT>
inline
const SpMat<eT>&
SpMat<eT>::sprandn(const SizeMat& s, const double density)
  {
  arma_extra_debug_sigprint();
  
  return (*this).sprandn(s.n_rows, s.n_cols, density);
  }



template<typename eT>
inline
void
SpMat<eT>::reset()
  {
  arma_extra_debug_sigprint();

  switch(vec_state)
    {
    default:
      init(0, 0);
      break;
      
    case 1:
      init(0, 1);
      break;
    
    case 2:
      init(1, 0);
      break;
    }
  }



template<typename eT>
inline
void
SpMat<eT>::reset_cache()
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  #if defined(ARMA_USE_OPENMP)
    {
    #pragma omp critical (arma_SpMat_cache)
      {
      cache.reset();
      
      sync_state = 0;
      }
    }
  #elif (!defined(ARMA_DONT_USE_STD_MUTEX))
    {
    cache_mutex.lock();
    
    cache.reset();
    
    sync_state = 0;
    
    cache_mutex.unlock();
    }
  #else
    {
    cache.reset();
    
    sync_state = 0;
    }
  #endif
  }



template<typename eT>
inline
void
SpMat<eT>::reserve(const uword in_rows, const uword in_cols, const uword new_n_nonzero)
  {
  arma_extra_debug_sigprint();
  
  init(in_rows, in_cols, new_n_nonzero);
  }



template<typename eT>
template<typename T1>
inline
void
SpMat<eT>::set_real(const SpBase<typename SpMat<eT>::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  SpMat_aux::set_real(*this, X);
  }



template<typename eT>
template<typename T1>
inline
void
SpMat<eT>::set_imag(const SpBase<typename SpMat<eT>::pod_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  SpMat_aux::set_imag(*this, X);
  }



//! save the matrix to a file
template<typename eT>
inline
arma_cold
bool
SpMat<eT>::save(const std::string name, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  bool save_okay;
  
  switch(type)
    {
    case csv_ascii:
      return (*this).save(csv_name(name), type);
      break;
    
    case ssv_ascii:
      return (*this).save(csv_name(name), type);
      break;
    
    case arma_binary:
      save_okay = diskio::save_arma_binary(*this, name);
      break;
    
    case coord_ascii:
      save_okay = diskio::save_coord_ascii(*this, name);
      break;
    
    default:
      arma_debug_warn_level(1, "SpMat::save(): unsupported file type");
      save_okay = false;
    }
  
  if(save_okay == false)  { arma_debug_warn_level(3, "SpMat::save(): couldn't write; file: ", name); }
  
  return save_okay;
  }



template<typename eT>
inline
arma_cold
bool
SpMat<eT>::save(const csv_name& spec, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  if( (type != csv_ascii) && (type != ssv_ascii) ) 
    {
    arma_stop_runtime_error("SpMat::save(): unsupported file type for csv_name()");
    return false;
    }
  
  const bool   do_trans     = bool(spec.opts.flags & csv_opts::flag_trans      );
  const bool   no_header    = bool(spec.opts.flags & csv_opts::flag_no_header  );
        bool with_header    = bool(spec.opts.flags & csv_opts::flag_with_header);
  const bool  use_semicolon = bool(spec.opts.flags & csv_opts::flag_semicolon  ) || (type == ssv_ascii);
  
  arma_extra_debug_print("SpMat::save(csv_name): enabled flags:");
  
  if(do_trans     )  { arma_extra_debug_print("trans");       }
  if(no_header    )  { arma_extra_debug_print("no_header");   }
  if(with_header  )  { arma_extra_debug_print("with_header"); }
  if(use_semicolon)  { arma_extra_debug_print("semicolon");   }
  
  const char separator = (use_semicolon) ? char(';') : char(',');
  
  if(no_header)  { with_header = false; }
  
  if(with_header)
    {
    if( (spec.header_ro.n_cols != 1) && (spec.header_ro.n_rows != 1) )
      {
      arma_debug_warn_level(1, "SpMat::save(): given header must have a vector layout");
      return false;
      }
    
    for(uword i=0; i < spec.header_ro.n_elem; ++i)
      {
      const std::string& token = spec.header_ro.at(i);
      
      if(token.find(separator) != std::string::npos)
        {
        arma_debug_warn_level(1, "SpMat::save(): token within the header contains the separator character: '", token, "'");
        return false;
        }
      }
    
    const uword save_n_cols = (do_trans) ? (*this).n_rows : (*this).n_cols;
    
    if(spec.header_ro.n_elem != save_n_cols)
      {
      arma_debug_warn_level(1, "SpMat::save(): size mistmach between header and matrix");
      return false;
      }
    }
  
  bool save_okay = false;
  
  if(do_trans)
    {
    const SpMat<eT> tmp = (*this).st();
    
    save_okay = diskio::save_csv_ascii(tmp, spec.filename, spec.header_ro, with_header, separator);
    }
  else
    {
    save_okay = diskio::save_csv_ascii(*this, spec.filename, spec.header_ro, with_header, separator);
    }
  
  if(save_okay == false)  { arma_debug_warn_level(3, "SpMat::save(): couldn't write; file: ", spec.filename); }
  
  return save_okay;
  }



//! save the matrix to a stream
template<typename eT>
inline
arma_cold
bool
SpMat<eT>::save(std::ostream& os, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  bool save_okay;
  
  switch(type)
    {
    case csv_ascii:
      save_okay = diskio::save_csv_ascii(*this, os, char(','));
      break;
    
    case ssv_ascii:
      save_okay = diskio::save_csv_ascii(*this, os, char(';'));
      break;
    
    case arma_binary:
      save_okay = diskio::save_arma_binary(*this, os);
      break;
    
    case coord_ascii:
      save_okay = diskio::save_coord_ascii(*this, os);
      break;
    
    default:
      arma_debug_warn_level(1, "SpMat::save(): unsupported file type");
      save_okay = false;
    }
  
  if(save_okay == false)  { arma_debug_warn_level(3, "SpMat::save(): couldn't write to stream"); }
  
  return save_okay;
  }



//! load a matrix from a file
template<typename eT>
inline
arma_cold
bool
SpMat<eT>::load(const std::string name, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  invalidate_cache();
  
  bool load_okay;
  std::string err_msg;
  
  switch(type)
    {
    // case auto_detect:
    //   load_okay = diskio::load_auto_detect(*this, name, err_msg);
    //   break;
    
    case csv_ascii:
      return (*this).load(csv_name(name), type);
      break;
    
    case ssv_ascii:
      return (*this).load(csv_name(name), type);
      break;
    
    case arma_binary:
      load_okay = diskio::load_arma_binary(*this, name, err_msg);
      break;
    
    case coord_ascii:
      load_okay = diskio::load_coord_ascii(*this, name, err_msg);
      break;
    
    default:
      arma_debug_warn_level(1, "SpMat::load(): unsupported file type");
      load_okay = false;
    }
  
  if(load_okay == false)
    {
    if(err_msg.length() > 0)
      {
      arma_debug_warn_level(3, "SpMat::load(): ", err_msg, "; file: ", name);
      }
    else
      {
      arma_debug_warn_level(3, "SpMat::load(): couldn't read; file: ", name);
      }
    }
  
  if(load_okay == false)  { (*this).reset(); }
  
  return load_okay;
  }



template<typename eT>
inline
arma_cold
bool
SpMat<eT>::load(const csv_name& spec, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  if( (type != csv_ascii) && (type != ssv_ascii) ) 
    {
    arma_stop_runtime_error("SpMat::load(): unsupported file type for csv_name()");
    return false;
    }
  
  const bool   do_trans     = bool(spec.opts.flags & csv_opts::flag_trans      );
  const bool   no_header    = bool(spec.opts.flags & csv_opts::flag_no_header  );
        bool with_header    = bool(spec.opts.flags & csv_opts::flag_with_header);
  const bool  use_semicolon = bool(spec.opts.flags & csv_opts::flag_semicolon  ) || (type == ssv_ascii);
  
  arma_extra_debug_print("SpMat::load(csv_name): enabled flags:");
  
  if(do_trans     )  { arma_extra_debug_print("trans");       }
  if(no_header    )  { arma_extra_debug_print("no_header");   }
  if(with_header  )  { arma_extra_debug_print("with_header"); }
  if(use_semicolon)  { arma_extra_debug_print("semicolon");   }
  
  const char separator = (use_semicolon) ? char(';') : char(',');
  
  if(no_header)  { with_header = false; }
  
  bool load_okay = false;
  std::string err_msg;
  
  if(do_trans)
    {
    SpMat<eT> tmp_mat;
    
    load_okay = diskio::load_csv_ascii(tmp_mat, spec.filename, err_msg, spec.header_rw, with_header, separator);
    
    if(load_okay)
      {
      (*this) = tmp_mat.st();
      
      if(with_header)
        {
        // field::set_size() preserves data if the number of elements hasn't changed
        spec.header_rw.set_size(spec.header_rw.n_elem, 1);
        }
      }
    }
  else
    {
    load_okay = diskio::load_csv_ascii(*this, spec.filename, err_msg, spec.header_rw, with_header, separator);
    }
  
  if(load_okay == false)
    {
    if(err_msg.length() > 0)
      {
      arma_debug_warn_level(3, "SpMat::load(): ", err_msg, "; file: ", spec.filename);
      }
    else
      {
      arma_debug_warn_level(3, "SpMat::load(): couldn't read; file: ", spec.filename);
      }
    }
  else
    {
    const uword load_n_cols = (do_trans) ? (*this).n_rows : (*this).n_cols;
    
    if(with_header && (spec.header_rw.n_elem != load_n_cols))
      {
      arma_debug_warn_level(3, "SpMat::load(): size mistmach between header and matrix");
      }
    }
  
  if(load_okay == false)
    {
    (*this).reset();
    
    if(with_header)  { spec.header_rw.reset(); }
    }
  
  return load_okay;
  }



//! load a matrix from a stream
template<typename eT>
inline
arma_cold
bool
SpMat<eT>::load(std::istream& is, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  invalidate_cache();
  
  bool load_okay;
  std::string err_msg;
  
  switch(type)
    {
    // case auto_detect:
    //   load_okay = diskio::load_auto_detect(*this, is, err_msg);
    //   break;
    
    case csv_ascii:
      load_okay = diskio::load_csv_ascii(*this, is, err_msg, char(','));
      break;
    
    case ssv_ascii:
      load_okay = diskio::load_csv_ascii(*this, is, err_msg, char(';'));
      break;
    
    case arma_binary:
      load_okay = diskio::load_arma_binary(*this, is, err_msg);
      break;
    
    case coord_ascii:
      load_okay = diskio::load_coord_ascii(*this, is, err_msg);
      break;
    
    default:
      arma_debug_warn_level(1, "SpMat::load(): unsupported file type");
      load_okay = false;
    }
  
  if(load_okay == false)
    {
    if(err_msg.length() > 0)
      {
      arma_debug_warn_level(3, "SpMat::load(): ", err_msg);
      }
    else
      {
      arma_debug_warn_level(3, "SpMat::load(): couldn't load from stream");
      }
    }
  
  if(load_okay == false)  { (*this).reset(); }
  
  return load_okay;
  }



//! save the matrix to a file, without printing any error messages
template<typename eT>
inline
arma_cold
bool
SpMat<eT>::quiet_save(const std::string name, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  return (*this).save(name, type);
  }



//! save the matrix to a stream, without printing any error messages
template<typename eT>
inline
arma_cold
bool
SpMat<eT>::quiet_save(std::ostream& os, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  return (*this).save(os, type);
  }



//! load a matrix from a file, without printing any error messages
template<typename eT>
inline
arma_cold
bool
SpMat<eT>::quiet_load(const std::string name, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  return (*this).load(name, type);
  }



//! load a matrix from a stream, without printing any error messages
template<typename eT>
inline
arma_cold
bool
SpMat<eT>::quiet_load(std::istream& is, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  return (*this).load(is, type);
  }



/**
 * Initialize the matrix to the specified size.  Data is not preserved, so the matrix is assumed to be entirely sparse (empty).
 */
template<typename eT>
inline
void
SpMat<eT>::init(uword in_rows, uword in_cols, const uword new_n_nonzero)
  {
  arma_extra_debug_sigprint();
  
  invalidate_cache(); // placed here, as init() is used during matrix modification
  
  // Clean out the existing memory.
  if(values     )  { memory::release(access::rw(values));      }
  if(row_indices)  { memory::release(access::rw(row_indices)); }
  if(col_ptrs   )  { memory::release(access::rw(col_ptrs));    }
  
  // in case init_cold() throws an exception
  access::rw(n_rows)      = 0;
  access::rw(n_cols)      = 0;
  access::rw(n_elem)      = 0;
  access::rw(n_nonzero)   = 0;
  access::rw(values)      = nullptr;
  access::rw(row_indices) = nullptr;
  access::rw(col_ptrs)    = nullptr;
  
  init_cold(in_rows, in_cols, new_n_nonzero);
  }



template<typename eT>
inline
void
arma_cold
SpMat<eT>::init_cold(uword in_rows, uword in_cols, const uword new_n_nonzero)
  {
  arma_extra_debug_sigprint();
  
  // Verify that we are allowed to do this.
  if(vec_state > 0)
    {
    if((in_rows == 0) && (in_cols == 0))
      {
      if(vec_state == 1)  { in_cols = 1; }
      if(vec_state == 2)  { in_rows = 1; }
      }
    else
      {
      if(vec_state == 1)  { arma_debug_check( (in_cols != 1), "SpMat::init(): object is a column vector; requested size is not compatible" ); }
      if(vec_state == 2)  { arma_debug_check( (in_rows != 1), "SpMat::init(): object is a row vector; requested size is not compatible"    ); }
      }
    }
  
  #if defined(ARMA_64BIT_WORD)
    const char* error_message = "SpMat::init(): requested size is too large";
  #else
    const char* error_message = "SpMat::init(): requested size is too large; suggest to enable ARMA_64BIT_WORD";
  #endif
  
  // Ensure that n_elem can hold the result of (n_rows * n_cols)
  arma_debug_check
    (
      (
      ( (in_rows > ARMA_MAX_UHWORD) || (in_cols > ARMA_MAX_UHWORD) )
        ? ( (double(in_rows) * double(in_cols)) > double(ARMA_MAX_UWORD) )
        : false
      ),
      error_message
    );
  
  access::rw(col_ptrs)    = memory::acquire<uword>(in_cols + 2);
  access::rw(values)      = memory::acquire<eT>   (new_n_nonzero + 1);
  access::rw(row_indices) = memory::acquire<uword>(new_n_nonzero + 1);
  
  // fill column pointers with 0,
  // except for the last element which contains the maximum possible element
  // (so iterators terminate correctly).
  arrayops::fill_zeros(access::rwp(col_ptrs), in_cols + 1);
  
  access::rw(col_ptrs[in_cols + 1]) = std::numeric_limits<uword>::max();
  
  access::rw(     values[new_n_nonzero]) = 0;
  access::rw(row_indices[new_n_nonzero]) = 0;
  
  // Set the new size accordingly.
  access::rw(n_rows)    = in_rows;
  access::rw(n_cols)    = in_cols;
  access::rw(n_elem)    = (in_rows * in_cols);
  access::rw(n_nonzero) = new_n_nonzero;
  }



template<typename eT>
inline
void
SpMat<eT>::init(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp(text);
  
  if(vec_state == 1)
    {
    if((tmp.n_elem > 0) && tmp.is_vec())
      {
      access::rw(tmp.n_rows) = tmp.n_elem;
      access::rw(tmp.n_cols) = 1;
      }
    }
  
  if(vec_state == 2)
    {
    if((tmp.n_elem > 0) && tmp.is_vec())
      {
      access::rw(tmp.n_rows) = 1;
      access::rw(tmp.n_cols) = tmp.n_elem;
      }
    }
  
  (*this).operator=(tmp);
  }



template<typename eT>
inline
void
SpMat<eT>::init(const SpMat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  if(this == &x)  { return; }
  
  bool init_done = false;
  
  #if defined(ARMA_USE_OPENMP)
    if(x.sync_state == 1)
      {
      #pragma omp critical (arma_SpMat_init)
      if(x.sync_state == 1)
        {
        (*this).init(x.cache);
        init_done = true;
        }
      }
  #elif (!defined(ARMA_DONT_USE_STD_MUTEX))
    if(x.sync_state == 1)
      {
      x.cache_mutex.lock();
      if(x.sync_state == 1)
        {
        (*this).init(x.cache);
        init_done = true;
        }
      x.cache_mutex.unlock();
      }
  #else
    if(x.sync_state == 1)
      {
      (*this).init(x.cache);
      init_done = true;
      }
  #endif
  
  if(init_done == false)
    {
    (*this).init_simple(x);
    }
  }



template<typename eT>
inline
void
SpMat<eT>::init(const MapMat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  const uword x_n_rows = x.n_rows;
  const uword x_n_cols = x.n_cols;
  const uword x_n_nz   = x.get_n_nonzero();
  
  init(x_n_rows, x_n_cols, x_n_nz);
  
  if(x_n_nz == 0)  { return; }
  
  typename MapMat<eT>::map_type& x_map_ref = *(x.map_ptr);
  
  typename MapMat<eT>::map_type::const_iterator x_it = x_map_ref.begin();
  
  uword x_col             = 0;
  uword x_col_index_start = 0;
  uword x_col_index_endp1 = x_n_rows;
  
  for(uword i=0; i < x_n_nz; ++i)
    {
    const std::pair<uword, eT>& x_entry = (*x_it);
    
    const uword x_index = x_entry.first;
    const eT    x_val   = x_entry.second;
    
    // have we gone past the curent column?
    if(x_index >= x_col_index_endp1)
      {
      x_col = x_index / x_n_rows;
      
      x_col_index_start = x_col * x_n_rows;
      x_col_index_endp1 = x_col_index_start + x_n_rows;
      }
    
    const uword x_row = x_index - x_col_index_start;
    
    // // sanity check
    // 
    // const uword tmp_x_row = x_index % x_n_rows;
    // const uword tmp_x_col = x_index / x_n_rows;
    // 
    // if(x_row != tmp_x_row)  { cout << "x_row != tmp_x_row" << endl; exit(-1); }
    // if(x_col != tmp_x_col)  { cout << "x_col != tmp_x_col" << endl; exit(-1); }
    
    access::rw(values[i])      = x_val;
    access::rw(row_indices[i]) = x_row;
    
    access::rw(col_ptrs[ x_col + 1 ])++;
    
    ++x_it;
    }
  
  
  for(uword i = 0; i < x_n_cols; ++i)
    {
    access::rw(col_ptrs[i + 1]) += col_ptrs[i];
    }
  
  
  // // OLD METHOD
  // 
  // for(uword i=0; i < x_n_nz; ++i)
  //   {
  //   const std::pair<uword, eT>& x_entry = (*x_it);
  //   
  //   const uword x_index = x_entry.first;
  //   const eT    x_val   = x_entry.second;
  // 
  // const uword x_row = x_index % x_n_rows;
  // const uword x_col = x_index / x_n_rows;
  // 
  // access::rw(values[i])      = x_val;
  // access::rw(row_indices[i]) = x_row;
  // 
  // access::rw(col_ptrs[ x_col + 1 ])++;
  // 
  // ++x_it;
  // }
  // 
  // 
  // for(uword i = 0; i < x_n_cols; ++i)
  // {
  // access::rw(col_ptrs[i + 1]) += col_ptrs[i];
  // }
  }



template<typename eT>
inline
void
SpMat<eT>::init_simple(const SpMat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  if(this == &x)  { return; }
  
  init(x.n_rows, x.n_cols, x.n_nonzero);
  
  if(x.values     )  { arrayops::copy(access::rwp(values),      x.values,      x.n_nonzero + 1); }
  if(x.row_indices)  { arrayops::copy(access::rwp(row_indices), x.row_indices, x.n_nonzero + 1); }
  if(x.col_ptrs   )  { arrayops::copy(access::rwp(col_ptrs),    x.col_ptrs,    x.n_cols    + 1); }
  }



template<typename eT>
inline
void
SpMat<eT>::init_batch_std(const Mat<uword>& locs, const Mat<eT>& vals, const bool sort_locations)
  {
  arma_extra_debug_sigprint();
  
  // Resize to correct number of elements.
  mem_resize(vals.n_elem);
  
  // Reset column pointers to zero.
  arrayops::fill_zeros(access::rwp(col_ptrs), n_cols + 1);
  
  bool actually_sorted = true;
  
  if(sort_locations)
    {
    // check if we really need a time consuming sort
    
    const uword locs_n_cols = locs.n_cols;
    
    for(uword i = 1; i < locs_n_cols; ++i)
      {
      const uword* locs_i   = locs.colptr(i  );
      const uword* locs_im1 = locs.colptr(i-1);
      
      const uword row_i = locs_i[0];
      const uword col_i = locs_i[1];
      
      const uword row_im1 = locs_im1[0];
      const uword col_im1 = locs_im1[1];
      
      if( (col_i < col_im1) || ((col_i == col_im1) && (row_i <= row_im1)) )
        {
        actually_sorted = false;
        break;
        }
      }
    
    if(actually_sorted == false)
      {
      // see op_sort_index_bones.hpp for the definition of arma_sort_index_packet and arma_sort_index_helper_ascend
      
      std::vector< arma_sort_index_packet<uword> > packet_vec(locs_n_cols);
      
      const uword* locs_mem = locs.memptr();
      
      for(uword i = 0; i < locs_n_cols; ++i)
        {
        const uword row = (*locs_mem);  locs_mem++;
        const uword col = (*locs_mem);  locs_mem++;
        
        packet_vec[i].val   = (col * n_rows) + row;
        packet_vec[i].index = i;
        }
      
      arma_sort_index_helper_ascend<uword> comparator;
      
      std::sort( packet_vec.begin(), packet_vec.end(), comparator );
      
      // insert the elements in the sorted order
      for(uword i = 0; i < locs_n_cols; ++i)
        {
        const uword index = packet_vec[i].index;
        
        const uword* locs_i = locs.colptr(index);
        
        const uword row_i = locs_i[0];
        const uword col_i = locs_i[1];
        
        arma_debug_check( ( (row_i >= n_rows) || (col_i >= n_cols) ), "SpMat::SpMat(): invalid row or column index" );
        
        if(i > 0)
          {
          const uword prev_index = packet_vec[i-1].index;
          
          const uword* locs_im1 = locs.colptr(prev_index);
          
          const uword row_im1 = locs_im1[0];
          const uword col_im1 = locs_im1[1];
          
          arma_debug_check( ( (row_i == row_im1) && (col_i == col_im1) ), "SpMat::SpMat(): detected identical locations" );
          }
        
        access::rw(values[i])      = vals[index];
        access::rw(row_indices[i]) = row_i;
        
        access::rw(col_ptrs[ col_i + 1 ])++;
        }
      }
    }
  
  if( (sort_locations == false) || (actually_sorted == true) )
    {
    // Now set the values and row indices correctly.
    // Increment the column pointers in each column (so they are column "counts").
    
    const uword locs_n_cols = locs.n_cols;
    
    for(uword i=0; i < locs_n_cols; ++i)
      {
      const uword* locs_i = locs.colptr(i);
      
      const uword row_i = locs_i[0];
      const uword col_i = locs_i[1];
      
      arma_debug_check( ( (row_i >= n_rows) || (col_i >= n_cols) ), "SpMat::SpMat(): invalid row or column index" );
      
      if(i > 0)
        {
        const uword* locs_im1 = locs.colptr(i-1);
        
        const uword row_im1 = locs_im1[0];
        const uword col_im1 = locs_im1[1];
        
        arma_debug_check
          (
          ( (col_i < col_im1) || ((col_i == col_im1) && (row_i < row_im1)) ),
          "SpMat::SpMat(): out of order points; either pass sort_locations = true, or sort points in column-major ordering"
          );
        
        arma_debug_check( ( (col_i == col_im1) && (row_i == row_im1) ), "SpMat::SpMat(): detected identical locations" );
        }
      
      access::rw(values[i])      = vals[i];
      access::rw(row_indices[i]) = row_i;
      
      access::rw(col_ptrs[ col_i + 1 ])++;
      }
    }
  
  // Now fix the column pointers.
  for(uword i = 0; i < n_cols; ++i)
    {
    access::rw(col_ptrs[i + 1]) += col_ptrs[i];
    }
  }



template<typename eT>
inline
void
SpMat<eT>::init_batch_add(const Mat<uword>& locs, const Mat<eT>& vals, const bool sort_locations)
  {
  arma_extra_debug_sigprint();
  
  if(locs.n_cols < 2)
    {
    init_batch_std(locs, vals, false);
    return;
    }
  
  // Reset column pointers to zero.
  arrayops::fill_zeros(access::rwp(col_ptrs), n_cols + 1);
  
  bool actually_sorted = true;
  
  if(sort_locations)
    {
    // sort_index() uses std::sort() which may use quicksort... so we better
    // make sure it's not already sorted before taking an O(N^2) sort penalty.
    for(uword i = 1; i < locs.n_cols; ++i)
      {
      const uword* locs_i   = locs.colptr(i  );
      const uword* locs_im1 = locs.colptr(i-1);
      
      if( (locs_i[1] < locs_im1[1]) || (locs_i[1] == locs_im1[1]  &&  locs_i[0] <= locs_im1[0]) )
        {
        actually_sorted = false;
        break;
        }
      }
    
    if(actually_sorted == false)
      {
      // This may not be the fastest possible implementation but it maximizes code reuse.
      Col<uword> abslocs(locs.n_cols, arma_nozeros_indicator());
      
      for(uword i = 0; i < locs.n_cols; ++i)
        {
        const uword* locs_i = locs.colptr(i);
        
        abslocs[i] = locs_i[1] * n_rows + locs_i[0];
        }
      
      uvec sorted_indices = sort_index(abslocs); // Ascending sort.
      
      // work out the number of unique elments 
      uword n_unique = 1;  // first element is unique
      
      for(uword i=1; i < sorted_indices.n_elem; ++i)
        {
        const uword* locs_i   = locs.colptr( sorted_indices[i  ] );
        const uword* locs_im1 = locs.colptr( sorted_indices[i-1] );
        
        if( (locs_i[1] != locs_im1[1]) || (locs_i[0] != locs_im1[0]) )  { ++n_unique; }
        }
      
      // resize to correct number of elements
      mem_resize(n_unique);
      
      // Now we add the elements in this sorted order.
      uword count = 0;
      
      // first element
        {
        const uword  i      = 0;
        const uword* locs_i = locs.colptr( sorted_indices[i] );
        
        arma_debug_check( ( (locs_i[0] >= n_rows) || (locs_i[1] >= n_cols) ), "SpMat::SpMat(): invalid row or column index" );
        
        access::rw(values[count])      = vals[ sorted_indices[i] ];
        access::rw(row_indices[count]) = locs_i[0];
        
        access::rw(col_ptrs[ locs_i[1] + 1 ])++;
        }
      
      for(uword i=1; i < sorted_indices.n_elem; ++i)
        {
        const uword* locs_i   = locs.colptr( sorted_indices[i  ] );
        const uword* locs_im1 = locs.colptr( sorted_indices[i-1] );
        
        arma_debug_check( ( (locs_i[0] >= n_rows) || (locs_i[1] >= n_cols) ), "SpMat::SpMat(): invalid row or column index" );
        
        if( (locs_i[1] == locs_im1[1]) && (locs_i[0] == locs_im1[0]) )
          {
          access::rw(values[count]) += vals[ sorted_indices[i] ];
          }
        else
          {
          count++;
          access::rw(values[count])      = vals[ sorted_indices[i] ];
          access::rw(row_indices[count]) = locs_i[0];
          
          access::rw(col_ptrs[ locs_i[1] + 1 ])++;
          }
        }
      }
    }
  
  if( (sort_locations == false) || (actually_sorted == true) )
    {
    // work out the number of unique elments 
    uword n_unique = 1;  // first element is unique
    
    for(uword i=1; i < locs.n_cols; ++i)
      {
      const uword* locs_i   = locs.colptr(i  );
      const uword* locs_im1 = locs.colptr(i-1);
      
      if( (locs_i[1] != locs_im1[1]) || (locs_i[0] != locs_im1[0]) )  { ++n_unique; }
      }
    
    // resize to correct number of elements
    mem_resize(n_unique);
    
    // Now set the values and row indices correctly.
    // Increment the column pointers in each column (so they are column "counts").
    
    uword count = 0;
    
    // first element
      {
      const uword  i      = 0;
      const uword* locs_i = locs.colptr(i);
      
      arma_debug_check( ( (locs_i[0] >= n_rows) || (locs_i[1] >= n_cols) ), "SpMat::SpMat(): invalid row or column index" );
      
      access::rw(values[count])      = vals[i];
      access::rw(row_indices[count]) = locs_i[0];
      
      access::rw(col_ptrs[ locs_i[1] + 1 ])++;
      }
    
    for(uword i=1; i < locs.n_cols; ++i)
      {
      const uword* locs_i   = locs.colptr(i  );
      const uword* locs_im1 = locs.colptr(i-1);
      
      arma_debug_check( ( (locs_i[0] >= n_rows) || (locs_i[1] >= n_cols) ), "SpMat::SpMat(): invalid row or column index" );
      
      arma_debug_check
        (
        ( (locs_i[1] < locs_im1[1]) || (locs_i[1] == locs_im1[1]  &&  locs_i[0] < locs_im1[0]) ),
        "SpMat::SpMat(): out of order points; either pass sort_locations = true, or sort points in column-major ordering"
        );
      
      if( (locs_i[1] == locs_im1[1]) && (locs_i[0] == locs_im1[0]) )
        {
        access::rw(values[count]) += vals[i];
        }
      else
        {
        count++;
        
        access::rw(values[count])      = vals[i];
        access::rw(row_indices[count]) = locs_i[0];
        
        access::rw(col_ptrs[ locs_i[1] + 1 ])++;
        }
      }
    }
  
  // Now fix the column pointers.
  for(uword i = 0; i < n_cols; ++i)
    {
    access::rw(col_ptrs[i + 1]) += col_ptrs[i];
    }
  }



//! constructor used by SpRow and SpCol classes
template<typename eT>
inline
SpMat<eT>::SpMat(const arma_vec_indicator&, const uword in_vec_state)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(in_vec_state)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  const uword in_n_rows = (in_vec_state == 2) ? 1 : 0;
  const uword in_n_cols = (in_vec_state == 1) ? 1 : 0;
  
  init_cold(in_n_rows, in_n_cols);
  }



//! constructor used by SpRow and SpCol classes
template<typename eT>
inline
SpMat<eT>::SpMat(const arma_vec_indicator&, const uword in_n_rows, const uword in_n_cols, const uword in_vec_state)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  , n_nonzero(0)
  , vec_state(in_vec_state)
  , values(nullptr)
  , row_indices(nullptr)
  , col_ptrs(nullptr)
  {
  arma_extra_debug_sigprint_this(this);
  
  init_cold(in_n_rows, in_n_cols);
  }



template<typename eT>
inline
void
SpMat<eT>::mem_resize(const uword new_n_nonzero)
  {
  arma_extra_debug_sigprint();
  
  invalidate_cache();  // placed here, as mem_resize() is used during matrix modification
  
  if(n_nonzero == new_n_nonzero)  { return; }
  
  eT*    new_values      = memory::acquire<eT>   (new_n_nonzero + 1);
  uword* new_row_indices = memory::acquire<uword>(new_n_nonzero + 1);
  
  if( (n_nonzero > 0 ) && (new_n_nonzero > 0) )
    {
    // Copy old elements.
    uword copy_len = (std::min)(n_nonzero, new_n_nonzero);
    
    arrayops::copy(new_values,      values,      copy_len);
    arrayops::copy(new_row_indices, row_indices, copy_len);
    }
  
  if(values)       { memory::release(access::rw(values));      }
  if(row_indices)  { memory::release(access::rw(row_indices)); }
  
  access::rw(values)      = new_values;
  access::rw(row_indices) = new_row_indices;
  
  // Set the "fake end" of the matrix by setting the last value and row index to 0.
  // This helps the iterators work correctly.
  access::rw(     values[new_n_nonzero]) = 0;
  access::rw(row_indices[new_n_nonzero]) = 0;
  
  access::rw(n_nonzero) = new_n_nonzero;
  }



template<typename eT>
inline
void
SpMat<eT>::sync() const
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  }



template<typename eT>
inline
void
SpMat<eT>::remove_zeros()
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  invalidate_cache();  // placed here, as remove_zeros() is used during matrix modification
  
  const uword old_n_nonzero = n_nonzero;
        uword new_n_nonzero = 0;
  
  const eT* old_values = values;
  
  for(uword i=0; i < old_n_nonzero; ++i)
    {
    new_n_nonzero += (old_values[i] != eT(0)) ? uword(1) : uword(0);
    }
  
  if(new_n_nonzero != old_n_nonzero)
    {
    if(new_n_nonzero == 0)  { init(n_rows, n_cols); return; }
    
    SpMat<eT> tmp(arma_reserve_indicator(), n_rows, n_cols, new_n_nonzero);
    
    uword new_index = 0;
    
    const_iterator it     = begin();
    const_iterator it_end = end();
    
    for(; it != it_end; ++it)
      {
      const eT val = eT(*it);
      
      if(val != eT(0))
        {
        access::rw(tmp.values[new_index])      = val;
        access::rw(tmp.row_indices[new_index]) = it.row();
        access::rw(tmp.col_ptrs[it.col() + 1])++;
        ++new_index;
        }
      }
    
    for(uword i=0; i < n_cols; ++i)
      {
      access::rw(tmp.col_ptrs[i + 1]) += tmp.col_ptrs[i];
      }
    
    steal_mem(tmp);
    }
  }



// Steal memory from another matrix.
template<typename eT>
inline
void
SpMat<eT>::steal_mem(SpMat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  if(this == &x)  { return; }
  
  bool layout_ok = false;
  
  if((*this).vec_state == x.vec_state)
    {
    layout_ok = true;
    }
  else
    {
    if( ((*this).vec_state == 1) && (x.n_cols == 1) )  { layout_ok = true; }
    if( ((*this).vec_state == 2) && (x.n_rows == 1) )  { layout_ok = true; }
    }
  
  if(layout_ok)
    {
    arma_extra_debug_print("SpMat::steal_mem(): stealing memory");
    
    x.sync_csc();
    
    steal_mem_simple(x);
    
    x.invalidate_cache();
    
    invalidate_cache();
    }
  else
    {
    arma_extra_debug_print("SpMat::steal_mem(): copying memory");
    
    (*this).operator=(x);
    }
  }



template<typename eT>
inline
void
SpMat<eT>::steal_mem_simple(SpMat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  if(this == &x)  { return; }
  
  if(values     )  { memory::release(access::rw(values));      }
  if(row_indices)  { memory::release(access::rw(row_indices)); }
  if(col_ptrs   )  { memory::release(access::rw(col_ptrs));    }
  
  access::rw(n_rows)    = x.n_rows;
  access::rw(n_cols)    = x.n_cols;
  access::rw(n_elem)    = x.n_elem;
  access::rw(n_nonzero) = x.n_nonzero;
  
  access::rw(values)      = x.values;
  access::rw(row_indices) = x.row_indices;
  access::rw(col_ptrs)    = x.col_ptrs;
  
  // Set other matrix to empty.
  access::rw(x.n_rows)    = 0;
  access::rw(x.n_cols)    = 0;
  access::rw(x.n_elem)    = 0;
  access::rw(x.n_nonzero) = 0;
  
  access::rw(x.values)      = nullptr;
  access::rw(x.row_indices) = nullptr;
  access::rw(x.col_ptrs)    = nullptr;
  }



template<typename eT>
template<typename T1, typename Functor>
inline
void
SpMat<eT>::init_xform(const SpBase<eT,T1>& A, const Functor& func)
  {
  arma_extra_debug_sigprint();
  
  // if possible, avoid doing a copy and instead apply func to the generated elements
  if(SpProxy<T1>::Q_is_generated)
    {
    (*this) = A.get_ref();
    
    const uword nnz = n_nonzero;
    
    eT* t_values = access::rwp(values);
    
    bool has_zero = false;
    
    for(uword i=0; i < nnz; ++i)
      {
      eT& t_values_i = t_values[i];
      
      t_values_i = func(t_values_i);
      
      if(t_values_i == eT(0))  { has_zero = true; }
      }
    
    if(has_zero)  { remove_zeros(); }
    }
  else
    {
    init_xform_mt(A.get_ref(), func);
    }
  }



template<typename eT>
template<typename eT2, typename T1, typename Functor>
inline
void
SpMat<eT>::init_xform_mt(const SpBase<eT2,T1>& A, const Functor& func)
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<T1> P(A.get_ref());
  
  if( P.is_alias(*this) || (is_SpMat<typename SpProxy<T1>::stored_type>::value) )
    {
    // NOTE: unwrap_spmat will convert a submatrix to a matrix, which in effect takes care of aliasing with submatrices;
    // NOTE: however, when more delayed ops are implemented, more elaborate handling of aliasing will be necessary
    const unwrap_spmat<typename SpProxy<T1>::stored_type> tmp(P.Q);
    
    const SpMat<eT2>& x = tmp.M;
    
    if(void_ptr(this) != void_ptr(&x))
      {
      init(x.n_rows, x.n_cols, x.n_nonzero);
      
      arrayops::copy(access::rwp(row_indices), x.row_indices, x.n_nonzero + 1);
      arrayops::copy(access::rwp(col_ptrs),    x.col_ptrs,    x.n_cols    + 1);
      }
    
    
    // initialise the elements array with a transformed version of the elements from x
    
    const uword nnz = n_nonzero;
    
    const eT2* x_values = x.values;
          eT*  t_values = access::rwp(values);
    
    bool has_zero = false;
    
    for(uword i=0; i < nnz; ++i)
      {
      eT& t_values_i = t_values[i];
      
      t_values_i = func(x_values[i]);   // NOTE: func() must produce a value of type eT (ie. act as a convertor between eT2 and eT)
      
      if(t_values_i == eT(0))  { has_zero = true; } 
      }
    
    if(has_zero)  { remove_zeros(); }
    }
  else
    {
    init(P.get_n_rows(), P.get_n_cols(), P.get_n_nonzero());
    
    typename SpProxy<T1>::const_iterator_type it     = P.begin();
    typename SpProxy<T1>::const_iterator_type it_end = P.end();
    
    bool has_zero = false;
    
    while(it != it_end)
      {
      const eT val = func(*it);   // NOTE: func() must produce a value of type eT (ie. act as a convertor between eT2 and eT)
      
      if(val == eT(0))  { has_zero = true; }
      
      const uword it_pos = it.pos();
      
      access::rw(row_indices[it_pos]) = it.row();
      access::rw(values[it_pos]) = val;
      ++access::rw(col_ptrs[it.col() + 1]);
      ++it;
      }
    
    // Now sum column pointers.
    for(uword c = 1; c <= n_cols; ++c)
      {
      access::rw(col_ptrs[c]) += col_ptrs[c - 1];
      }
    
    if(has_zero)  { remove_zeros(); }
    }
  }



template<typename eT>
arma_inline
bool
SpMat<eT>::is_alias(const SpMat<eT>& X) const
  {
  return (&X == this);
  }



template<typename eT>
inline
typename SpMat<eT>::iterator
SpMat<eT>::begin()
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  return iterator(*this);
  }



template<typename eT>
inline
typename SpMat<eT>::const_iterator
SpMat<eT>::begin() const
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  return const_iterator(*this);
  }



template<typename eT>
inline
typename SpMat<eT>::const_iterator
SpMat<eT>::cbegin() const
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  
  return const_iterator(*this);
  }



template<typename eT>
inline
typename SpMat<eT>::iterator
SpMat<eT>::end()
  {
  sync_csc();
  
  return iterator(*this, 0, n_cols, n_nonzero);
  }



template<typename eT>
inline
typename SpMat<eT>::const_iterator
SpMat<eT>::end() const
  {
  sync_csc();
  
  return const_iterator(*this, 0, n_cols, n_nonzero);
  }



template<typename eT>
inline
typename SpMat<eT>::const_iterator
SpMat<eT>::cend() const
  {
  sync_csc();
  
  return const_iterator(*this, 0, n_cols, n_nonzero);
  }



template<typename eT>
inline
typename SpMat<eT>::col_iterator
SpMat<eT>::begin_col(const uword col_num)
  {
  sync_csc();
  
  return col_iterator(*this, 0, col_num);
  }



template<typename eT>
inline
typename SpMat<eT>::const_col_iterator
SpMat<eT>::begin_col(const uword col_num) const
  {
  sync_csc();
  
  return const_col_iterator(*this, 0, col_num);
  }



template<typename eT>
inline
typename SpMat<eT>::col_iterator
SpMat<eT>::begin_col_no_sync(const uword col_num)
  {
  return col_iterator(*this, 0, col_num);
  }



template<typename eT>
inline
typename SpMat<eT>::const_col_iterator
SpMat<eT>::begin_col_no_sync(const uword col_num) const
  {
  return const_col_iterator(*this, 0, col_num);
  }



template<typename eT>
inline
typename SpMat<eT>::col_iterator
SpMat<eT>::end_col(const uword col_num)
  {
  sync_csc();
  
  return col_iterator(*this, 0, col_num + 1);
  }



template<typename eT>
inline
typename SpMat<eT>::const_col_iterator
SpMat<eT>::end_col(const uword col_num) const
  {
  sync_csc();
  
  return const_col_iterator(*this, 0, col_num + 1);
  }



template<typename eT>
inline
typename SpMat<eT>::col_iterator
SpMat<eT>::end_col_no_sync(const uword col_num)
  {
  return col_iterator(*this, 0, col_num + 1);
  }



template<typename eT>
inline
typename SpMat<eT>::const_col_iterator
SpMat<eT>::end_col_no_sync(const uword col_num) const
  {
  return const_col_iterator(*this, 0, col_num + 1);
  }



template<typename eT>
inline
typename SpMat<eT>::row_iterator
SpMat<eT>::begin_row(const uword row_num)
  {
  sync_csc();
  
  return row_iterator(*this, row_num, 0);
  }



template<typename eT>
inline
typename SpMat<eT>::const_row_iterator
SpMat<eT>::begin_row(const uword row_num) const
  {
  sync_csc();
  
  return const_row_iterator(*this, row_num, 0);
  }



template<typename eT>
inline
typename SpMat<eT>::row_iterator
SpMat<eT>::end_row()
  {
  sync_csc();
  
  return row_iterator(*this, n_nonzero);
  }



template<typename eT>
inline
typename SpMat<eT>::const_row_iterator
SpMat<eT>::end_row() const
  {
  sync_csc();
  
  return const_row_iterator(*this, n_nonzero);
  }



template<typename eT>
inline
typename SpMat<eT>::row_iterator
SpMat<eT>::end_row(const uword row_num)
  {
  sync_csc();
  
  return row_iterator(*this, row_num + 1, 0);
  }



template<typename eT>
inline
typename SpMat<eT>::const_row_iterator
SpMat<eT>::end_row(const uword row_num) const
  {
  sync_csc();
  
  return const_row_iterator(*this, row_num + 1, 0);
  }



template<typename eT>
inline
typename SpMat<eT>::row_col_iterator
SpMat<eT>::begin_row_col()
  {
  sync_csc();
  
  return begin();
  }



template<typename eT>
inline
typename SpMat<eT>::const_row_col_iterator
SpMat<eT>::begin_row_col() const
  {
  sync_csc();
  
  return begin();
  }



template<typename eT>
inline typename SpMat<eT>::row_col_iterator
SpMat<eT>::end_row_col()
  {
  sync_csc();
  
  return end();
  }



template<typename eT>
inline
typename SpMat<eT>::const_row_col_iterator
SpMat<eT>::end_row_col() const
  {
  sync_csc();
  
  return end();
  }



template<typename eT>
inline
void
SpMat<eT>::clear()
  {
  (*this).reset();
  }



template<typename eT>
inline
bool
SpMat<eT>::empty() const
  {
  return (n_elem == 0);
  }



template<typename eT>
inline
uword
SpMat<eT>::size() const
  {
  return n_elem;
  }



template<typename eT>
arma_inline
arma_warn_unused
SpMat_MapMat_val<eT>
SpMat<eT>::front()
  {
  arma_debug_check( (n_elem == 0), "SpMat::front(): matrix is empty" );
  
  return SpMat_MapMat_val<eT>((*this), cache, 0, 0);
  }



template<typename eT>
arma_inline
arma_warn_unused
eT
SpMat<eT>::front() const
  {
  arma_debug_check( (n_elem == 0), "SpMat::front(): matrix is empty" );
  
  return get_value(0,0);
  }



template<typename eT>
arma_inline
arma_warn_unused
SpMat_MapMat_val<eT>
SpMat<eT>::back()
  {
  arma_debug_check( (n_elem == 0), "SpMat::back(): matrix is empty" );
  
  return SpMat_MapMat_val<eT>((*this), cache, n_rows-1, n_cols-1);
  }



template<typename eT>
arma_inline
arma_warn_unused
eT
SpMat<eT>::back() const
  {
  arma_debug_check( (n_elem == 0), "SpMat::back(): matrix is empty" );
  
  return get_value(n_rows-1, n_cols-1);
  }



template<typename eT>
inline
arma_hot
arma_warn_unused
eT
SpMat<eT>::get_value(const uword i) const
  {
  const MapMat<eT>& const_cache = cache;  // declare as const for clarity of intent
  
  // get the element from the cache if it has more recent data than CSC
  
  return (sync_state == 1) ? const_cache.operator[](i) : get_value_csc(i);
  }



template<typename eT>
inline
arma_hot
arma_warn_unused
eT
SpMat<eT>::get_value(const uword in_row, const uword in_col) const
  {
  const MapMat<eT>& const_cache = cache;  // declare as const for clarity of intent
  
  // get the element from the cache if it has more recent data than CSC
  
  return (sync_state == 1) ? const_cache.at(in_row, in_col) : get_value_csc(in_row, in_col);
  }



template<typename eT>
inline
arma_hot
arma_warn_unused
eT
SpMat<eT>::get_value_csc(const uword i) const
  {
  // First convert to the actual location.
  uword lcol = i / n_rows; // Integer division.
  uword lrow = i % n_rows;
  
  return get_value_csc(lrow, lcol);
  }



template<typename eT>
inline
arma_hot
arma_warn_unused
const eT*
SpMat<eT>::find_value_csc(const uword in_row, const uword in_col) const
  {
  const uword      col_offset = col_ptrs[in_col    ];
  const uword next_col_offset = col_ptrs[in_col + 1];
  
  const uword* start_ptr = &row_indices[     col_offset];
  const uword*   end_ptr = &row_indices[next_col_offset];
  
  const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, in_row);  // binary search
  
  if( (pos_ptr != end_ptr) && ((*pos_ptr) == in_row) )
    {
    const uword offset = uword(pos_ptr - start_ptr);
    const uword index  = offset + col_offset;
    
    return &(values[index]);
    }
  
  return nullptr;
  }



template<typename eT>
inline
arma_hot
arma_warn_unused
eT
SpMat<eT>::get_value_csc(const uword in_row, const uword in_col) const
  {
  const eT* val_ptr = find_value_csc(in_row, in_col);
  
  return (val_ptr != nullptr) ? eT(*val_ptr) : eT(0);
  }



template<typename eT>
inline
arma_hot
arma_warn_unused
bool
SpMat<eT>::try_set_value_csc(const uword in_row, const uword in_col, const eT in_val)
  {
  const eT* val_ptr = find_value_csc(in_row, in_col);
  
  // element not found, ie. it's zero; fail if trying to set it to non-zero value
  if(val_ptr == nullptr)  { return (in_val == eT(0)); }
  
  // fail if trying to erase an existing element
  if(in_val == eT(0))  { return false; }
  
  access::rw(*val_ptr) = in_val;
  
  invalidate_cache();
  
  return true;
  }



template<typename eT>
inline
arma_hot
arma_warn_unused
bool
SpMat<eT>::try_add_value_csc(const uword in_row, const uword in_col, const eT in_val)
  {
  const eT* val_ptr = find_value_csc(in_row, in_col);
  
  // element not found, ie. it's zero; fail if trying to add a non-zero value
  if(val_ptr == nullptr)  { return (in_val == eT(0)); }
  
  const eT new_val = eT(*val_ptr) + in_val;
  
  // fail if trying to erase an existing element
  if(new_val == eT(0))  { return false; }
  
  access::rw(*val_ptr) = new_val;
  
  invalidate_cache();
  
  return true;
  }



template<typename eT>
inline
arma_hot
arma_warn_unused
bool
SpMat<eT>::try_sub_value_csc(const uword in_row, const uword in_col, const eT in_val)
  {
  const eT* val_ptr = find_value_csc(in_row, in_col);
  
  // element not found, ie. it's zero; fail if trying to subtract a non-zero value
  if(val_ptr == nullptr)  { return (in_val == eT(0)); }
  
  const eT new_val = eT(*val_ptr) - in_val;
  
  // fail if trying to erase an existing element
  if(new_val == eT(0))  { return false; }
  
  access::rw(*val_ptr) = new_val;
  
  invalidate_cache();
  
  return true;
  }



template<typename eT>
inline
arma_hot
arma_warn_unused
bool
SpMat<eT>::try_mul_value_csc(const uword in_row, const uword in_col, const eT in_val)
  {
  const eT* val_ptr = find_value_csc(in_row, in_col);
  
  // element not found, ie. it's zero; succeed if given value is finite; zero multiplied by anything is zero, except for nan and inf
  if(val_ptr == nullptr)  { return arma_isfinite(in_val); }
  
  const eT new_val = eT(*val_ptr) * in_val;
  
  // fail if trying to erase an existing element
  if(new_val == eT(0))  { return false; }
  
  access::rw(*val_ptr) = new_val;
  
  invalidate_cache();
  
  return true;
  }



template<typename eT>
inline
arma_hot
arma_warn_unused
bool
SpMat<eT>::try_div_value_csc(const uword in_row, const uword in_col, const eT in_val)
  {
  const eT* val_ptr = find_value_csc(in_row, in_col);
  
  // element not found, ie. it's zero; succeed if given value is not zero and not nan; zero divided by anything is zero, except for zero and nan
  if(val_ptr == nullptr)  { return ((in_val != eT(0)) && (arma_isnan(in_val) == false)); }
  
  const eT new_val = eT(*val_ptr) / in_val;
  
  // fail if trying to erase an existing element
  if(new_val == eT(0))  { return false; }
  
  access::rw(*val_ptr) = new_val;
  
  invalidate_cache();
  
  return true;
  }



/**
 * Insert an element at the given position, and return a reference to it.  
 * The element will be set to 0, unless otherwise specified.
 * If the element already exists, its value will be overwritten.
 */
template<typename eT>
inline
arma_warn_unused
eT&
SpMat<eT>::insert_element(const uword in_row, const uword in_col, const eT val)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  invalidate_cache();
  
  // We will assume the new element does not exist and begin the search for
  // where to insert it.  If we find that it already exists, we will then
  // overwrite it.
  uword colptr      = col_ptrs[in_col    ];
  uword next_colptr = col_ptrs[in_col + 1];
  
  uword pos = colptr; // The position in the matrix of this value.
  
  if(colptr != next_colptr)
    {
    // There are other elements in this column, so we must find where this
    // element will fit as compared to those.
    while(pos < next_colptr && in_row > row_indices[pos])
      {
      pos++;
      }
    
    // We aren't inserting into the last position, so it is still possible
    // that the element may exist.
    if(pos != next_colptr && row_indices[pos] == in_row)
      {
      // It already exists.  Then, just overwrite it.
      access::rw(values[pos]) = val;
      
      return access::rw(values[pos]);
      }
    }
  
  
  // 
  // Element doesn't exist, so we have to insert it
  // 
  
  // We have to update the rest of the column pointers.
  for(uword i = in_col + 1; i < n_cols + 1; i++)
    {
    access::rw(col_ptrs[i])++; // We are only inserting one new element.
    }
  
  const uword old_n_nonzero = n_nonzero;
  
  access::rw(n_nonzero)++; // Add to count of nonzero elements.
  
  // Allocate larger memory.
  eT*    new_values      = memory::acquire<eT>   (n_nonzero + 1);
  uword* new_row_indices = memory::acquire<uword>(n_nonzero + 1);
  
  // Copy things over, before the new element.
  if(pos > 0)
    {
    arrayops::copy(new_values,      values,      pos);
    arrayops::copy(new_row_indices, row_indices, pos);
    }
  
  // Insert the new element.
  new_values[pos]      = val;
  new_row_indices[pos] = in_row;
  
  // Copy the rest of things over (including the extra element at the end).
  arrayops::copy(new_values      + pos + 1, values      + pos, (old_n_nonzero - pos) + 1);
  arrayops::copy(new_row_indices + pos + 1, row_indices + pos, (old_n_nonzero - pos) + 1);
  
  // Assign new pointers.
  if(values)       { memory::release(access::rw(values));      }
  if(row_indices)  { memory::release(access::rw(row_indices)); }
  
  access::rw(values)      = new_values;
  access::rw(row_indices) = new_row_indices;
  
  return access::rw(values[pos]);
  }



/**
 * Delete an element at the given position.
 */
template<typename eT>
inline
void
SpMat<eT>::delete_element(const uword in_row, const uword in_col)
  {
  arma_extra_debug_sigprint();
  
  sync_csc();
  invalidate_cache();
  
  // We assume the element exists (although... it may not) and look for its
  // exact position.  If it doesn't exist... well, we don't need to do anything.
  uword colptr      = col_ptrs[in_col];
  uword next_colptr = col_ptrs[in_col + 1];
  
  if(colptr != next_colptr)
    {
    // There's at least one element in this column.
    // Let's see if we are one of them.
    for(uword pos = colptr; pos < next_colptr; pos++)
      {
      if(in_row == row_indices[pos])
        {
        --access::rw(n_nonzero); // Remove one from the count of nonzero elements.
        
        // Found it.  Now remove it.
        
        // Make new arrays.
        eT*    new_values      = memory::acquire<eT>   (n_nonzero + 1);
        uword* new_row_indices = memory::acquire<uword>(n_nonzero + 1);
        
        if(pos > 0)
          {
          arrayops::copy(new_values,      values,      pos);
          arrayops::copy(new_row_indices, row_indices, pos);
          }
        
        arrayops::copy(new_values      + pos, values      + pos + 1, (n_nonzero - pos) + 1);
        arrayops::copy(new_row_indices + pos, row_indices + pos + 1, (n_nonzero - pos) + 1);
        
        if(values)       { memory::release(access::rw(values));      }
        if(row_indices)  { memory::release(access::rw(row_indices)); }
        
        access::rw(values)      = new_values;
        access::rw(row_indices) = new_row_indices;
        
        // And lastly, update all the column pointers (decrement by one).
        for(uword i = in_col + 1; i < n_cols + 1; i++)
          {
          --access::rw(col_ptrs[i]); // We only removed one element.
          }
        
        return; // There is nothing left to do.
        }
      }
    }
  
  return; // The element does not exist, so there's nothing for us to do.
  }



template<typename eT>
arma_inline
void
SpMat<eT>::invalidate_cache() const
  {
  arma_extra_debug_sigprint();
  
  if(sync_state == 0)  { return; }
  
  cache.reset();
  
  sync_state = 0;
  }



template<typename eT>
arma_inline
void
SpMat<eT>::invalidate_csc() const
  {
  arma_extra_debug_sigprint();
  
  sync_state = 1;
  }



template<typename eT>
inline
void
SpMat<eT>::sync_cache() const
  {
  arma_extra_debug_sigprint();
  
  // using approach adapted from http://preshing.com/20130930/double-checked-locking-is-fixed-in-cpp11/
  // 
  // OpenMP mode:
  // sync_state uses atomic read/write, which has an implied flush;
  // flush is also implicitly executed at the entrance and the exit of critical section;
  // data races are prevented by the 'critical' directive
  // 
  // C++11  mode:
  // underlying type for sync_state is std::atomic<int>;
  // reading and writing to sync_state uses std::memory_order_seq_cst which has an implied fence;
  // data races are prevented via the mutex
  
  #if defined(ARMA_USE_OPENMP)
    {
    if(sync_state == 0)
      {
      #pragma omp critical (arma_SpMat_cache)
        {
        sync_cache_simple();
        }
      }
    }
  #elif (!defined(ARMA_DONT_USE_STD_MUTEX))
    {
    if(sync_state == 0)
      {
      cache_mutex.lock();
      
      sync_cache_simple();
      
      cache_mutex.unlock();
      }
    }
  #else
    {
    sync_cache_simple();
    }
  #endif
  }




template<typename eT>
inline
void
SpMat<eT>::sync_cache_simple() const
  {
  arma_extra_debug_sigprint();
  
  if(sync_state == 0)
    {
    cache      = (*this);
    sync_state = 2;
    }
  }




template<typename eT>
inline
void
SpMat<eT>::sync_csc() const
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_OPENMP)
    if(sync_state == 1)
      {
      #pragma omp critical (arma_SpMat_cache)
        {
        sync_csc_simple();
        }
      }
  #elif (!defined(ARMA_DONT_USE_STD_MUTEX))
    if(sync_state == 1)
      {
      cache_mutex.lock();
      
      sync_csc_simple();
      
      cache_mutex.unlock();
      }
  #else
    {
    sync_csc_simple();
    }
  #endif
  }



template<typename eT>
inline
void
SpMat<eT>::sync_csc_simple() const
  {
  arma_extra_debug_sigprint();
  
  // method:
  // 1. construct temporary matrix to prevent the cache from getting zapped
  // 2. steal memory from the temporary matrix
  
  // sync_state is only set to 1 by non-const element access operators,
  // so the shenanigans with const_cast are to satisfy the compiler
  
  // see also the note in sync_cache() above
  
  if(sync_state == 1)
    {
    SpMat<eT>& x = const_cast< SpMat<eT>& >(*this);
    
    SpMat<eT> tmp(cache);
    
    x.steal_mem_simple(tmp);
    
    sync_state = 2;
    }
  }




//
// SpMat_aux



template<typename eT, typename T1>
inline
void
SpMat_aux::set_real(SpMat<eT>& out, const SpBase<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_spmat<T1> tmp(X.get_ref());
  const SpMat<eT>&   A = tmp.M;
  
  arma_debug_assert_same_size( out, A, "SpMat::set_real()" );
  
  out = A;
  }



template<typename eT, typename T1>
inline
void
SpMat_aux::set_imag(SpMat<eT>&, const SpBase<eT,T1>&)
  {
  arma_extra_debug_sigprint();
  }



template<typename T, typename T1>
inline
void
SpMat_aux::set_real(SpMat< std::complex<T> >& out, const SpBase<T,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  const unwrap_spmat<T1> U(X.get_ref());
  const SpMat<T>&    Y = U.M;
  
  arma_debug_assert_same_size(out, Y, "SpMat::set_real()");
  
  SpMat<eT> tmp(Y,arma::imag(out));  // arma:: prefix required due to bugs in GCC 4.4 - 4.6
  
  out.steal_mem(tmp);
  }



template<typename T, typename T1>
inline
void
SpMat_aux::set_imag(SpMat< std::complex<T> >& out, const SpBase<T,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  const unwrap_spmat<T1> U(X.get_ref());
  const SpMat<T>&    Y = U.M;
  
  arma_debug_assert_same_size(out, Y, "SpMat::set_imag()");
  
  SpMat<eT> tmp(arma::real(out),Y);  // arma:: prefix required due to bugs in GCC 4.4 - 4.6
  
  out.steal_mem(tmp);
  }



#if defined(ARMA_EXTRA_SPMAT_MEAT)
  #include ARMA_INCFILE_WRAP(ARMA_EXTRA_SPMAT_MEAT)
#endif



//! @}
