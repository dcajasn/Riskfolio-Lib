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


//! \addtogroup fn_symmat
//! @{


template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::no, const Op<T1, op_symmatu> >::result
symmatu(const Base<typename T1::elem_type,T1>& X, const bool do_conj = false)
  {
  arma_extra_debug_sigprint();
  arma_ignore(do_conj);
  
  return Op<T1, op_symmatu>(X.get_ref());
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::no, const Op<T1, op_symmatl> >::result
symmatl(const Base<typename T1::elem_type,T1>& X, const bool do_conj = false)
  {
  arma_extra_debug_sigprint();
  arma_ignore(do_conj);
  
  return Op<T1, op_symmatl>(X.get_ref());
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::yes, const Op<T1, op_symmatu_cx> >::result
symmatu(const Base<typename T1::elem_type,T1>& X, const bool do_conj = true)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_symmatu_cx>(X.get_ref(), 0, (do_conj ? 1 : 0));
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::yes, const Op<T1, op_symmatl_cx> >::result
symmatl(const Base<typename T1::elem_type,T1>& X, const bool do_conj = true)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_symmatl_cx>(X.get_ref(), 0, (do_conj ? 1 : 0));
  }



//



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::no, const SpOp<T1, spop_symmat> >::result
symmatu(const SpBase<typename T1::elem_type,T1>& X, const bool do_conj = false)
  {
  arma_extra_debug_sigprint();
  arma_ignore(do_conj);
  
  return SpOp<T1, spop_symmat>(X.get_ref(), 0, 0);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::no, const SpOp<T1, spop_symmat> >::result
symmatl(const SpBase<typename T1::elem_type,T1>& X, const bool do_conj = false)
  {
  arma_extra_debug_sigprint();
  arma_ignore(do_conj);
  
  return SpOp<T1, spop_symmat>(X.get_ref(), 1, 0);
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::yes, const SpOp<T1, spop_symmat_cx> >::result
symmatu(const SpBase<typename T1::elem_type,T1>& X, const bool do_conj = true)
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1, spop_symmat_cx>(X.get_ref(), 0, (do_conj ? 1 : 0));
  }



template<typename T1>
arma_warn_unused
arma_inline
typename enable_if2< is_cx<typename T1::elem_type>::yes, const SpOp<T1, spop_symmat_cx> >::result
symmatl(const SpBase<typename T1::elem_type,T1>& X, const bool do_conj = true)
  {
  arma_extra_debug_sigprint();
  
  return SpOp<T1, spop_symmat_cx>(X.get_ref(), 1, (do_conj ? 1 : 0));
  }



//! @}
