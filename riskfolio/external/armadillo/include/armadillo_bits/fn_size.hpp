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


//! \addtogroup fn_size
//! @{



arma_warn_unused
inline
const SizeMat
size(const uword n_rows, const uword n_cols)
  {
  arma_extra_debug_sigprint();
  
  return SizeMat(n_rows, n_cols);
  }



template<typename T1>
arma_warn_unused
inline
const SizeMat
size(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const Proxy<T1> P(X.get_ref());
  
  return SizeMat( P.get_n_rows(), P.get_n_cols() );
  }



// explicit overload to workround ADL issues with C++17 std::size()
template<typename eT>
arma_warn_unused
inline
const SizeMat
size(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return SizeMat( X.n_rows, X.n_cols );
  }



// explicit overload to workround ADL issues with C++17 std::size()
template<typename eT>
arma_warn_unused
inline
const SizeMat
size(const Row<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return SizeMat( X.n_rows, X.n_cols );
  }



// explicit overload to workround ADL issues with C++17 std::size()
template<typename eT>
arma_warn_unused
inline
const SizeMat
size(const Col<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return SizeMat( X.n_rows, X.n_cols );
  }



arma_warn_unused
inline
const SizeMat
size(const arma::span& row_span, const arma::span& col_span)
  {
  arma_extra_debug_sigprint();
  
  uword n_rows = 0;
  uword n_cols = 0;
  
  if(row_span.whole || col_span.whole)
    {
    arma_debug_check(true, "size(): span::all not supported");
    }
  else
    {
    if((row_span.a > row_span.b) || (col_span.a > col_span.b))
      {
      arma_debug_check_bounds(true, "size(): span indices incorrectly used");
      }
    else
      {
      n_rows = row_span.b - row_span.a + 1;
      n_cols = col_span.b - col_span.a + 1;
      }
    }
  
  return SizeMat(n_rows, n_cols);
  }



template<typename T1>
arma_warn_unused
inline
uword
size(const Base<typename T1::elem_type,T1>& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  const Proxy<T1> P(X.get_ref());
  
  return SizeMat( P.get_n_rows(), P.get_n_cols() )( dim );
  }



arma_warn_unused
inline
const SizeCube
size(const uword n_rows, const uword n_cols, const uword n_slices)
  {
  arma_extra_debug_sigprint();
  
  return SizeCube(n_rows, n_cols, n_slices);
  }



template<typename T1>
arma_warn_unused
inline
const SizeCube
size(const BaseCube<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const ProxyCube<T1> P(X.get_ref());
  
  return SizeCube( P.get_n_rows(), P.get_n_cols(), P.get_n_slices() );
  }



// explicit overload to workround ADL issues with C++17 std::size()
template<typename eT>
arma_warn_unused
inline
const SizeCube
size(const Cube<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return SizeCube( X.n_rows, X.n_cols, X.n_slices );
  }



template<typename T1>
arma_warn_unused
inline
uword
size(const BaseCube<typename T1::elem_type, T1>& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  const ProxyCube<T1> P(X.get_ref());
  
  return SizeCube( P.get_n_rows(), P.get_n_cols(), P.get_n_slices() )( dim );
  }



arma_warn_unused
inline
const SizeCube
size(const arma::span& row_span, const arma::span& col_span, const arma::span& slice_span)
  {
  arma_extra_debug_sigprint();
  
  uword n_rows   = 0;
  uword n_cols   = 0;
  uword n_slices = 0;
  
  if(row_span.whole || col_span.whole || slice_span.whole)
    {
    arma_debug_check(true, "size(): span::all not supported");
    }
  else
    {
    if((row_span.a > row_span.b) || (col_span.a > col_span.b) || (slice_span.a > slice_span.b))
      {
      arma_debug_check_bounds(true, "size(): span indices incorrectly used");
      }
    else
      {
      n_rows   =   row_span.b -   row_span.a + 1;
      n_cols   =   col_span.b -   col_span.a + 1;
      n_slices = slice_span.b - slice_span.a + 1;
      }
    }
  
  return SizeCube(n_rows, n_cols, n_slices);
  }



template<typename T1>
arma_warn_unused
inline
const SizeMat
size(const SpBase<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<T1> P(X.get_ref());
  
  return SizeMat( P.get_n_rows(), P.get_n_cols() );
  }



// explicit overload to workround ADL issues with C++17 std::size()
template<typename eT>
arma_warn_unused
inline
const SizeMat
size(const SpMat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return SizeMat( X.n_rows, X.n_cols );
  }



template<typename T1>
arma_warn_unused
inline
uword
size(const SpBase<typename T1::elem_type,T1>& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<T1> P(X.get_ref());
  
  return SizeMat( P.get_n_rows(), P.get_n_cols() )( dim );
  }




template<typename oT>
arma_warn_unused
inline
const SizeCube
size(const field<oT>& X)
  {
  arma_extra_debug_sigprint();
  
  return SizeCube( X.n_rows, X.n_cols, X.n_slices );
  }



template<typename oT>
arma_warn_unused
inline
uword
size(const field<oT>& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return SizeCube( X.n_rows, X.n_cols, X.n_slices )( dim );
  }



template<typename oT>
arma_warn_unused
inline
const SizeCube
size(const subview_field<oT>& X)
  {
  arma_extra_debug_sigprint();
  
  return SizeCube( X.n_rows, X.n_cols, X.n_slices );
  }



template<typename oT>
arma_warn_unused
inline
uword
size(const subview_field<oT>& X, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  return SizeCube( X.n_rows, X.n_cols, X.n_slices )( dim );
  }



//! @}
