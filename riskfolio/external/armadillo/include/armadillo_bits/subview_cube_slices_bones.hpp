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


//! \addtogroup subview_cube_slices
//! @{



template<typename eT, typename T1>
class subview_cube_slices : public BaseCube< eT, subview_cube_slices<eT,T1> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  arma_aligned const Cube<eT>&       m;
  arma_aligned const Base<uword,T1>& base_si;
  
  
  protected:
  
  arma_inline subview_cube_slices(const Cube<eT>& in_m, const Base<uword,T1>& in_si);
  
  
  public:
  
  inline ~subview_cube_slices();
  inline  subview_cube_slices() = delete;
  
  inline void inplace_rand(const uword rand_mode);
  
  template<typename op_type>                inline void inplace_op(const eT val);
  template<typename op_type, typename expr> inline void inplace_op(const BaseCube<eT,expr>& x);
  
  inline void fill(const eT val);
  inline void zeros();
  inline void ones();
  inline void randu();
  inline void randn();
  
  inline void operator+= (const eT val);
  inline void operator-= (const eT val);
  inline void operator*= (const eT val);
  inline void operator/= (const eT val);
  
  
  // deliberately returning void
  template<typename T2> inline void operator_equ(const subview_cube_slices<eT,T2>& x);
  template<typename T2> inline void operator=   (const subview_cube_slices<eT,T2>& x);
                        inline void operator=   (const subview_cube_slices<eT,T1>& x);
  
  template<typename T2> inline void operator+=  (const subview_cube_slices<eT,T2>& x);
  template<typename T2> inline void operator-=  (const subview_cube_slices<eT,T2>& x);
  template<typename T2> inline void operator%=  (const subview_cube_slices<eT,T2>& x);
  template<typename T2> inline void operator/=  (const subview_cube_slices<eT,T2>& x);
  
  template<typename expr> inline void operator=  (const BaseCube<eT,expr>& x);
  template<typename expr> inline void operator+= (const BaseCube<eT,expr>& x);
  template<typename expr> inline void operator-= (const BaseCube<eT,expr>& x);
  template<typename expr> inline void operator%= (const BaseCube<eT,expr>& x);
  template<typename expr> inline void operator/= (const BaseCube<eT,expr>& x);
  
  inline static void extract(Cube<eT>& out, const subview_cube_slices& in);
  
  inline static void  plus_inplace(Cube<eT>& out, const subview_cube_slices& in);
  inline static void minus_inplace(Cube<eT>& out, const subview_cube_slices& in);
  inline static void schur_inplace(Cube<eT>& out, const subview_cube_slices& in);
  inline static void   div_inplace(Cube<eT>& out, const subview_cube_slices& in);
  
  
  friend class Cube<eT>;
  };



//! @}
