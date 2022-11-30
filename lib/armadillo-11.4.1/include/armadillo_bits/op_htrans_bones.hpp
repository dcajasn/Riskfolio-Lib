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


//! \addtogroup op_htrans
//! @{


//! 'hermitian transpose' operation

class op_htrans
  {
  public:
  
  template<typename T1>
  struct traits
    {
    static constexpr bool is_row  = T1::is_col;  // deliberately swapped
    static constexpr bool is_col  = T1::is_row;
    static constexpr bool is_xvec = T1::is_xvec;
    };
  
  template<typename eT>
  arma_hot inline static void apply_mat_noalias(Mat<eT>& out, const Mat<eT>& A, const typename arma_not_cx<eT>::result* junk = nullptr);
  
  template<typename eT>
  arma_hot inline static void apply_mat_noalias(Mat<eT>& out, const Mat<eT>& A, const typename arma_cx_only<eT>::result* junk = nullptr);
  
  //
  
  template<typename T>
  arma_hot inline static void block_worker(std::complex<T>* Y, const std::complex<T>* X, const uword X_n_rows, const uword Y_n_rows, const uword n_rows, const uword n_cols);
  
  template<typename T>
  arma_hot inline static void apply_mat_noalias_large(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& A);
  
  //
  
  template<typename eT>
  arma_hot inline static void apply_mat_inplace(Mat<eT>& out, const typename arma_not_cx<eT>::result* junk = nullptr);
  
  template<typename eT>
  arma_hot inline static void apply_mat_inplace(Mat<eT>& out, const typename arma_cx_only<eT>::result* junk = nullptr);
  
  //
  
  template<typename eT>
  inline static void apply_mat(Mat<eT>& out, const Mat<eT>& A, const typename arma_not_cx<eT>::result* junk = nullptr);
  
  template<typename eT>
  inline static void apply_mat(Mat<eT>& out, const Mat<eT>& A, const typename arma_cx_only<eT>::result* junk = nullptr);
  
  //
  
  template<typename T1>
  inline static void apply_proxy(Mat<typename T1::elem_type>& out, const Proxy<T1>& P);
  
  //
  
  template<typename T1>
  inline static void apply_direct(Mat<typename T1::elem_type>& out, const T1& X);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans>& in, const typename arma_not_cx<typename T1::elem_type>::result* junk = nullptr);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans>& in, const typename arma_cx_only<typename T1::elem_type>::result* junk = nullptr);
  };



class op_htrans2
  {
  public:
  
  template<typename T1>
  struct traits
    {
    static constexpr bool is_row  = T1::is_col;  // deliberately swapped
    static constexpr bool is_col  = T1::is_row;
    static constexpr bool is_xvec = T1::is_xvec;
    };
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans2>& in, const typename arma_not_cx<typename T1::elem_type>::result* junk = nullptr);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans2>& in, const typename arma_cx_only<typename T1::elem_type>::result* junk = nullptr);
  };



//! @}
