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


//! \addtogroup op_median
//! @{



template<typename T1>
inline
void
op_median::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_median>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> U(expr.m);
  
  const uword dim = expr.aux_uword_a;
  
  arma_debug_check( U.M.has_nan(), "median(): detected NaN"                   );
  arma_debug_check( (dim > 1),     "median(): parameter 'dim' must be 0 or 1" );
  
  if(U.is_alias(out))
    {
    Mat<eT> tmp;
    
    op_median::apply_noalias(out, U.M, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_median::apply_noalias(out, U.M, dim);
    }
  }



template<typename eT>
inline
void
op_median::apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim, const typename arma_not_cx<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  if(dim == 0)  // in each column
    {
    arma_extra_debug_print("op_median::apply(): dim = 0");
    
    out.set_size((X_n_rows > 0) ? 1 : 0, X_n_cols);
    
    if(X_n_rows > 0)
      {
      std::vector<eT> tmp_vec(X_n_rows);
      
      for(uword col=0; col < X_n_cols; ++col)
        {
        arrayops::copy( &(tmp_vec[0]), X.colptr(col), X_n_rows );
        
        out[col] = op_median::direct_median(tmp_vec);
        }
      }
    }
  else
  if(dim == 1)  // in each row
    {
    arma_extra_debug_print("op_median::apply(): dim = 1");
    
    out.set_size(X_n_rows, (X_n_cols > 0) ? 1 : 0);
    
    if(X_n_cols > 0)
      {
      std::vector<eT> tmp_vec(X_n_cols);
        
      for(uword row=0; row < X_n_rows; ++row)
        {
        for(uword col=0; col < X_n_cols; ++col)  { tmp_vec[col] = X.at(row,col); }
        
        out[row] = op_median::direct_median(tmp_vec);
        }
      }
    }
  }



template<typename eT>
inline
void
op_median::apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim, const typename arma_cx_only<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename get_pod_type<eT>::result T;
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  if(dim == 0)  // in each column
    {
    arma_extra_debug_print("op_median::apply(): dim = 0");
    
    out.set_size((X_n_rows > 0) ? 1 : 0, X_n_cols);
    
    if(X_n_rows > 0)
      {
      std::vector< arma_cx_median_packet<T> > tmp_vec(X_n_rows);
      
      for(uword col=0; col<X_n_cols; ++col)
        {
        const eT* colmem = X.colptr(col);
        
        for(uword row=0; row<X_n_rows; ++row)
          {
          tmp_vec[row].val   = std::abs(colmem[row]);
          tmp_vec[row].index = row;
          }
        
        uword index1 = 0;
        uword index2 = 0;
        op_median::direct_cx_median_index(index1, index2, tmp_vec);
          
        out[col] = op_mean::robust_mean(colmem[index1], colmem[index2]);
        }
      }
    }
  else
  if(dim == 1)  // in each row
    {
    arma_extra_debug_print("op_median::apply(): dim = 1");
    
    out.set_size(X_n_rows, (X_n_cols > 0) ? 1 : 0);
    
    if(X_n_cols > 0)
      {
      std::vector< arma_cx_median_packet<T> > tmp_vec(X_n_cols);
      
      for(uword row=0; row<X_n_rows; ++row)
        {
        for(uword col=0; col<X_n_cols; ++col)
          {
          tmp_vec[col].val   = std::abs(X.at(row,col));
          tmp_vec[col].index = col;
          }
        
        uword index1 = 0;
        uword index2 = 0;
        op_median::direct_cx_median_index(index1, index2, tmp_vec);
        
        out[row] = op_mean::robust_mean( X.at(row,index1), X.at(row,index2) );
        }
      }
    }
  }



template<typename T1>
inline
typename T1::elem_type
op_median::median_vec
  (
  const T1& X,
  const typename arma_not_cx<typename T1::elem_type>::result* junk
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> U(X);
  
  const uword n_elem = U.M.n_elem;
  
  if(n_elem == 0)
    {
    arma_debug_check(true, "median(): object has no elements");
    
    return Datum<eT>::nan;
    }
  
  arma_debug_check( U.M.has_nan(), "median(): detected NaN" );
  
  std::vector<eT> tmp_vec(n_elem);
  
  arrayops::copy( &(tmp_vec[0]), U.M.memptr(), n_elem );
  
  return op_median::direct_median(tmp_vec);
  }



template<typename T1>
inline
typename T1::elem_type
op_median::median_vec
  (
  const T1& X,
  const typename arma_cx_only<typename T1::elem_type>::result* junk
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const quasi_unwrap<T1> U(X);
  
  const uword n_elem = U.M.n_elem;
  
  if(n_elem == 0)
    {
    arma_debug_check(true, "median(): object has no elements");
    
    return Datum<eT>::nan;
    }
  
  arma_debug_check( U.M.has_nan(), "median(): detected NaN" );
  
  std::vector< arma_cx_median_packet<T> > tmp_vec(n_elem);
  
  const eT* A = U.M.memptr();
  
  for(uword i=0; i<n_elem; ++i)
    {
    tmp_vec[i].val   = std::abs( A[i] );
    tmp_vec[i].index = i;
    }
  
  uword index1 = 0;
  uword index2 = 0;
  op_median::direct_cx_median_index(index1, index2, tmp_vec);
  
  return op_mean::robust_mean( A[index1], A[index2] );
  }



template<typename eT>
inline 
eT
op_median::direct_median(std::vector<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const uword n_elem = uword(X.size());
  const uword half   = n_elem/2;
  
  typename std::vector<eT>::iterator first    = X.begin();
  typename std::vector<eT>::iterator nth      = first + half;
  typename std::vector<eT>::iterator pastlast = X.end();
  
  std::nth_element(first, nth, pastlast);
  
  if((n_elem % 2) == 0)  // even number of elements
    {
    typename std::vector<eT>::iterator start   = X.begin();
    typename std::vector<eT>::iterator pastend = start + half;
    
    const eT val1 = (*nth);
    const eT val2 = (*(std::max_element(start, pastend)));
    
    return op_mean::robust_mean(val1, val2);
    }
  else  // odd number of elements
    {
    return (*nth);
    }
  }



template<typename T>
inline 
void
op_median::direct_cx_median_index
  (
  uword& out_index1, 
  uword& out_index2, 
  std::vector< arma_cx_median_packet<T> >& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef arma_cx_median_packet<T> eT;
  
  const uword n_elem = uword(X.size());
  const uword half   = n_elem/2;
  
  typename std::vector<eT>::iterator first    = X.begin();
  typename std::vector<eT>::iterator nth      = first + half;
  typename std::vector<eT>::iterator pastlast = X.end();
  
  std::nth_element(first, nth, pastlast);
  
  out_index1 = (*nth).index;
  
  if((n_elem % 2) == 0)  // even number of elements
    {
    typename std::vector<eT>::iterator start   = X.begin();
    typename std::vector<eT>::iterator pastend = start + half;
    
    out_index2 = (*(std::max_element(start, pastend))).index;
    }
  else  // odd number of elements
    {
    out_index2 = out_index1;
    }
  }



//! @}

