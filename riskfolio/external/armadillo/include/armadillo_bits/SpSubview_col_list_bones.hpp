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


//! \addtogroup SpSubview_col_list
//! @{



template<typename eT, typename T1>
class SpSubview_col_list : public SpBase< eT, SpSubview_col_list<eT,T1> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const SpMat<eT>&       m;
  const quasi_unwrap<T1> U_ci;
  
  
  protected:
  
  arma_inline SpSubview_col_list(const SpMat<eT>& in_m, const Base<uword,T1>& in_ci);
  
  
  public:
  
  inline ~SpSubview_col_list();
  inline  SpSubview_col_list() = delete;
  
  template<typename functor> inline void  for_each(functor F);
  template<typename functor> inline void  for_each(functor F) const;
  
  template<typename functor> inline void transform(functor F);
  
  inline void replace(const eT old_val, const eT new_val);
  
  inline void clean(const pod_type threshold);
  
  inline void fill(const eT val);
  inline void zeros();
  inline void ones();
  
  inline void operator+= (const eT val);
  inline void operator-= (const eT val);
  inline void operator*= (const eT val);
  inline void operator/= (const eT val);
  
  template<typename expr> inline void operator= (const Base<eT, expr>& x);
  template<typename expr> inline void operator+=(const Base<eT, expr>& x);
  template<typename expr> inline void operator-=(const Base<eT, expr>& x);
  template<typename expr> inline void operator%=(const Base<eT, expr>& x);
  template<typename expr> inline void operator/=(const Base<eT, expr>& x);
  
                        inline void operator= (const SpSubview_col_list<eT,T1>& x);
  template<typename T2> inline void operator= (const SpSubview_col_list<eT,T2>& x);
  
  template<typename expr> inline void operator=  (const SpBase<eT,expr>& x);
  template<typename expr> inline void operator+= (const SpBase<eT,expr>& x);
  template<typename expr> inline void operator-= (const SpBase<eT,expr>& x);
  template<typename expr> inline void operator%= (const SpBase<eT,expr>& x);
  template<typename expr> inline void operator/= (const SpBase<eT,expr>& x);
  
  inline static void extract(SpMat<eT>& out, const SpSubview_col_list& in);
  
  inline static void  plus_inplace(SpMat<eT>& out, const SpSubview_col_list& in);
  inline static void minus_inplace(SpMat<eT>& out, const SpSubview_col_list& in);
  inline static void schur_inplace(SpMat<eT>& out, const SpSubview_col_list& in);
  inline static void   div_inplace(SpMat<eT>& out, const SpSubview_col_list& in);
  
  
  friend class SpMat<eT>;
  };



//! @}
