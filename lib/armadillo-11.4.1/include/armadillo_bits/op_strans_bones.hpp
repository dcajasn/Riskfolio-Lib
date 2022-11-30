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


//! \addtogroup op_strans
//! @{


//! 'matrix transpose' operation (simple transpose, ie. without taking the conjugate of the elements)

class op_strans
  {
  public:
  
  template<typename T1>
  struct traits
    {
    static constexpr bool is_row  = T1::is_col;  // deliberately swapped
    static constexpr bool is_col  = T1::is_row;
    static constexpr bool is_xvec = T1::is_xvec;
    };
  
  template<const bool do_flip, const uword row, const uword col>
  struct pos
    {
    static constexpr uword n2 = (do_flip == false) ? (row + col*2) : (col + row*2);
    static constexpr uword n3 = (do_flip == false) ? (row + col*3) : (col + row*3);
    static constexpr uword n4 = (do_flip == false) ? (row + col*4) : (col + row*4);
    };
  
  template<typename eT, typename TA>
  arma_cold inline static void apply_mat_noalias_tinysq(Mat<eT>& out, const TA& A);
  
  template<typename eT>
  arma_hot inline static void block_worker(eT* Y, const eT* X, const uword X_n_rows, const uword Y_n_rows, const uword n_rows, const uword n_cols);
  
  template<typename eT>
  arma_hot inline static void apply_mat_noalias_large(Mat<eT>& out, const Mat<eT>& A);
  
  template<typename eT, typename TA>
  arma_hot inline static void apply_mat_noalias(Mat<eT>& out, const TA& A);
  
  template<typename eT>
  arma_hot inline static void apply_mat_inplace(Mat<eT>& out);
  
  template<typename eT, typename TA>
  inline static void apply_mat(Mat<eT>& out, const TA& A);
  
  template<typename T1>
  inline static void apply_proxy(Mat<typename T1::elem_type>& out, const Proxy<T1>& P);
  
  template<typename T1>
  inline static void apply_direct(Mat<typename T1::elem_type>& out, const T1& X);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_strans>& in);
  };



class op_strans_cube
  {
  public:
  
  template<typename eT>
  inline static void apply_noalias(Cube<eT>& out, const Cube<eT>& X);
  };



//! @}
