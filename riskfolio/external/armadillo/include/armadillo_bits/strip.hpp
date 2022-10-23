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


//! \addtogroup strip
//! @{



template<typename T1>
struct strip_diagmat
  {
  typedef T1 stored_type;
  
  inline
  strip_diagmat(const T1& X)
    : M(X)
    {
    arma_extra_debug_sigprint();
    }
  
  static constexpr bool do_diagmat = false;
  
  const T1& M;
  };



template<typename T1>
struct strip_diagmat< Op<T1, op_diagmat> >
  {
  typedef T1 stored_type;
  
  inline
  strip_diagmat(const Op<T1, op_diagmat>& X)
    : M(X.m)
    {
    arma_extra_debug_sigprint();
    }
  
  static constexpr bool do_diagmat = true;
  
  const T1& M;
  };



template<typename T1>
struct strip_inv
  {
  typedef T1 stored_type;
  
  inline
  strip_inv(const T1& X)
    : M(X)
    {
    arma_extra_debug_sigprint();
    }
  
  const T1& M;
  
  static constexpr bool do_inv_gen = false;
  static constexpr bool do_inv_spd = false;
  };



template<typename T1>
struct strip_inv< Op<T1, op_inv_gen_default> >
  {
  typedef T1 stored_type;
  
  inline
  strip_inv(const Op<T1, op_inv_gen_default>& X)
    : M(X.m)
    {
    arma_extra_debug_sigprint();
    }
  
  const T1& M;
  
  static constexpr bool do_inv_gen = true;
  static constexpr bool do_inv_spd = false;
  };



template<typename T1>
struct strip_inv< Op<T1, op_inv_spd_default> >
  {
  typedef T1 stored_type;
  
  inline
  strip_inv(const Op<T1, op_inv_spd_default>& X)
    : M(X.m)
    {
    arma_extra_debug_sigprint();
    }
  
  const T1& M;
  
  static constexpr bool do_inv_gen = false;
  static constexpr bool do_inv_spd = true;
  };



template<typename T1>
struct strip_trimat
  {
  typedef T1 stored_type;
  
  const T1& M;
  
  static constexpr bool do_trimat = false;
  static constexpr bool do_triu   = false;
  static constexpr bool do_tril   = false;
  
  inline
  strip_trimat(const T1& X)
    : M(X)
    {
    arma_extra_debug_sigprint();
    }
  };



template<typename T1>
struct strip_trimat< Op<T1, op_trimat> >
  {
  typedef T1 stored_type;
  
  const T1& M;
  
  static constexpr bool do_trimat = true;
  
  const bool do_triu;
  const bool do_tril;
  
  inline
  strip_trimat(const Op<T1, op_trimat>& X)
    : M(X.m)
    , do_triu(X.aux_uword_a == 0)
    , do_tril(X.aux_uword_a == 1)
    {
    arma_extra_debug_sigprint();
    }
  };



//! @}
