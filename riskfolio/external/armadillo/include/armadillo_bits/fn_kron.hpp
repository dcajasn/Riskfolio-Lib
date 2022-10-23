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


//! \addtogroup fn_kron
//! @{



template<typename T1, typename T2>
arma_warn_unused
arma_inline
const Glue<T1,T2,glue_kron>
kron(const Base<typename T1::elem_type,T1>& A, const Base<typename T1::elem_type,T2>& B)
  {
  arma_extra_debug_sigprint();

  return Glue<T1, T2, glue_kron>(A.get_ref(), B.get_ref());
  }



template<typename T, typename T1, typename T2>
arma_warn_unused
inline
Mat<typename eT_promoter<T1,T2>::eT>
kron(const Base<std::complex<T>,T1>& X, const Base<T,T2>& Y)
  {
  arma_extra_debug_sigprint();

  typedef typename std::complex<T> eT1;

  promote_type<eT1,T>::check();
  
  const quasi_unwrap<T1> tmp1(X.get_ref());
  const quasi_unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT1>& A = tmp1.M;
  const Mat<T  >& B = tmp2.M;

  Mat<eT1> out;
  
  glue_kron::direct_kron(out, A, B);
  
  return out;
  }



template<typename T, typename T1, typename T2>
arma_warn_unused
inline
Mat<typename eT_promoter<T1,T2>::eT>
kron(const Base<T,T1>& X, const Base<std::complex<T>,T2>& Y)
  {
  arma_extra_debug_sigprint();

  typedef typename std::complex<T> eT2;  

  promote_type<T,eT2>::check();
  
  const quasi_unwrap<T1> tmp1(X.get_ref());
  const quasi_unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<T  >& A = tmp1.M;
  const Mat<eT2>& B = tmp2.M;

  Mat<eT2> out;
  
  glue_kron::direct_kron(out, A, B);
  
  return out;
  }



template<typename T1, typename T2>
arma_warn_unused
arma_inline
const SpGlue<T1, T2, spglue_kron>
kron(const SpBase<typename T1::elem_type,T1>& A, const SpBase<typename T1::elem_type,T2>& B)
  {
  arma_extra_debug_sigprint();
  
  return SpGlue<T1, T2, spglue_kron>(A.get_ref(), B.get_ref());
  }



//! @}
