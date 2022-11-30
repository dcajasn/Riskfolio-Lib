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


//! \addtogroup op_det
//! @{



template<typename T1>
inline
bool
op_det::apply_direct(typename T1::elem_type& out_val, const Base<typename T1::elem_type,T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  if(strip_diagmat<T1>::do_diagmat)
    {
    const strip_diagmat<T1> strip(expr.get_ref());
    
    out_val = op_det::apply_diagmat(strip.M);
    
    return true;
    }
  
  if(strip_trimat<T1>::do_trimat)
    {
    const strip_trimat<T1> strip(expr.get_ref());
    
    out_val = op_det::apply_trimat(strip.M);
    
    return true;
    }
  
  Mat<eT> A(expr.get_ref());
  
  arma_debug_check( (A.is_square() == false), "det(): given matrix must be square sized" );
  
  const uword N = A.n_rows;
  
  if(N == 0)  { out_val = eT(1); return true; }
  if(N == 1)  { out_val =  A[0]; return true; }
  
  if((is_cx<eT>::no) && (N <= 4))
    {
    constexpr T det_min =        std::numeric_limits<T>::epsilon();
    constexpr T det_max = T(1) / std::numeric_limits<T>::epsilon();
    
    eT det_val = eT(0);
    
    if(N == 2)  { det_val = op_det::apply_tiny_2x2(A); }
    if(N == 3)  { det_val = op_det::apply_tiny_3x3(A); }
    if(N == 4)  { det_val = op_det::apply_tiny_4x4(A); }
    
    const T abs_det_val = std::abs(det_val);
    
    if((abs_det_val > det_min) && (abs_det_val < det_max))  { out_val = det_val; return true; }
    
    // fallthrough if det_val is suspect
    }
  
  if(A.is_diagmat())  { out_val = op_det::apply_diagmat(A); return true; }
  
  const bool is_triu =                   trimat_helper::is_triu(A);
  const bool is_tril = is_triu ? false : trimat_helper::is_tril(A);
  
  if(is_triu || is_tril)  { out_val = op_det::apply_trimat(A); return true; }
  
  return auxlib::det(out_val, A);
  }



template<typename T1>
inline
typename T1::elem_type
op_det::apply_diagmat(const Base<typename T1::elem_type,T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const diagmat_proxy<T1> A(expr.get_ref());
  
  arma_debug_check( (A.n_rows != A.n_cols), "det(): given matrix must be square sized" );
  
  const uword N = (std::min)(A.n_rows, A.n_cols);
  
  eT val = eT(1);
  
  for(uword i=0; i<N; ++i)  { val *= A[i]; }
  
  return val;
  }



template<typename T1>
inline
typename T1::elem_type
op_det::apply_trimat(const Base<typename T1::elem_type,T1>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(expr.get_ref());
  
  const uword N = P.get_n_rows();
  
  arma_debug_check( (N != P.get_n_cols()), "det(): given matrix must be square sized" );
  
  eT val = eT(1);
  
  for(uword i=0; i<N; ++i)  { val *= P.at(i,i); }
  
  return val;
  }



template<typename eT>
arma_cold
inline
eT
op_det::apply_tiny_2x2(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const eT* Xm = X.memptr();
  
  return ( Xm[pos<0,0>::n2]*Xm[pos<1,1>::n2] - Xm[pos<0,1>::n2]*Xm[pos<1,0>::n2] );
  }



template<typename eT>
arma_cold
inline
eT
op_det::apply_tiny_3x3(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const eT* Xm = X.memptr();
  
  // const double tmp1 = X.at(0,0) * X.at(1,1) * X.at(2,2);
  // const double tmp2 = X.at(0,1) * X.at(1,2) * X.at(2,0);
  // const double tmp3 = X.at(0,2) * X.at(1,0) * X.at(2,1);
  // const double tmp4 = X.at(2,0) * X.at(1,1) * X.at(0,2);
  // const double tmp5 = X.at(2,1) * X.at(1,2) * X.at(0,0);
  // const double tmp6 = X.at(2,2) * X.at(1,0) * X.at(0,1);
  // return (tmp1+tmp2+tmp3) - (tmp4+tmp5+tmp6);
  
  const eT val1 = Xm[pos<0,0>::n3]*(Xm[pos<2,2>::n3]*Xm[pos<1,1>::n3] - Xm[pos<2,1>::n3]*Xm[pos<1,2>::n3]);
  const eT val2 = Xm[pos<1,0>::n3]*(Xm[pos<2,2>::n3]*Xm[pos<0,1>::n3] - Xm[pos<2,1>::n3]*Xm[pos<0,2>::n3]);
  const eT val3 = Xm[pos<2,0>::n3]*(Xm[pos<1,2>::n3]*Xm[pos<0,1>::n3] - Xm[pos<1,1>::n3]*Xm[pos<0,2>::n3]);
  
  return ( val1 - val2 + val3 );
  }



template<typename eT>
arma_cold
inline
eT
op_det::apply_tiny_4x4(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const eT* Xm = X.memptr();
  
  const eT val_03_12 = Xm[pos<0,3>::n4] * Xm[pos<1,2>::n4];
  const eT val_02_13 = Xm[pos<0,2>::n4] * Xm[pos<1,3>::n4];
  const eT val_03_11 = Xm[pos<0,3>::n4] * Xm[pos<1,1>::n4];
  
  const eT val_01_13 = Xm[pos<0,1>::n4] * Xm[pos<1,3>::n4];
  const eT val_02_11 = Xm[pos<0,2>::n4] * Xm[pos<1,1>::n4];
  const eT val_01_12 = Xm[pos<0,1>::n4] * Xm[pos<1,2>::n4];
  
  const eT val_03_10 = Xm[pos<0,3>::n4] * Xm[pos<1,0>::n4];
  const eT val_00_13 = Xm[pos<0,0>::n4] * Xm[pos<1,3>::n4];
  const eT val_02_10 = Xm[pos<0,2>::n4] * Xm[pos<1,0>::n4];
  const eT val_00_12 = Xm[pos<0,0>::n4] * Xm[pos<1,2>::n4];
  
  const eT val_01_10 = Xm[pos<0,1>::n4] * Xm[pos<1,0>::n4];
  const eT val_00_11 = Xm[pos<0,0>::n4] * Xm[pos<1,1>::n4];
  
  const eT val_21_30 = Xm[pos<2,1>::n4] * Xm[pos<3,0>::n4];
  const eT val_22_30 = Xm[pos<2,2>::n4] * Xm[pos<3,0>::n4];
  const eT val_23_30 = Xm[pos<2,3>::n4] * Xm[pos<3,0>::n4];
  
  const eT val_20_31 = Xm[pos<2,0>::n4] * Xm[pos<3,1>::n4];
  const eT val_22_31 = Xm[pos<2,2>::n4] * Xm[pos<3,1>::n4];
  const eT val_23_31 = Xm[pos<2,3>::n4] * Xm[pos<3,1>::n4];
  
  const eT val_20_32 = Xm[pos<2,0>::n4] * Xm[pos<3,2>::n4];
  const eT val_21_32 = Xm[pos<2,1>::n4] * Xm[pos<3,2>::n4];
  const eT val_23_32 = Xm[pos<2,3>::n4] * Xm[pos<3,2>::n4];
  
  const eT val_20_33 = Xm[pos<2,0>::n4] * Xm[pos<3,3>::n4];
  const eT val_21_33 = Xm[pos<2,1>::n4] * Xm[pos<3,3>::n4];
  const eT val_22_33 = Xm[pos<2,2>::n4] * Xm[pos<3,3>::n4];
  
  const eT val = \
      val_03_12 * val_21_30 \
    - val_02_13 * val_21_30 \
    - val_03_11 * val_22_30 \
    + val_01_13 * val_22_30 \
    + val_02_11 * val_23_30 \
    - val_01_12 * val_23_30 \
    - val_03_12 * val_20_31 \
    + val_02_13 * val_20_31 \
    + val_03_10 * val_22_31 \
    - val_00_13 * val_22_31 \
    - val_02_10 * val_23_31 \
    + val_00_12 * val_23_31 \
    + val_03_11 * val_20_32 \
    - val_01_13 * val_20_32 \
    - val_03_10 * val_21_32 \
    + val_00_13 * val_21_32 \
    + val_01_10 * val_23_32 \
    - val_00_11 * val_23_32 \
    - val_02_11 * val_20_33 \
    + val_01_12 * val_20_33 \
    + val_02_10 * val_21_33 \
    - val_00_12 * val_21_33 \
    - val_01_10 * val_22_33 \
    + val_00_11 * val_22_33 \
    ;
  
  return val;
  }



//! @}
