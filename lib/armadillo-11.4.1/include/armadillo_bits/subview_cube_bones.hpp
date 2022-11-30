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


//! \addtogroup subview_cube
//! @{


//! Class for storing data required to construct or apply operations to a subcube
//! (ie. where the subcube starts and ends as well as a reference/pointer to the original cube),
template<typename eT>
class subview_cube : public BaseCube< eT, subview_cube<eT> >
  {
  public:    
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  arma_aligned const Cube<eT>& m;
  
  const uword aux_row1;
  const uword aux_col1;
  const uword aux_slice1;
  
  const uword n_rows;
  const uword n_cols;
  const uword n_elem_slice;
  const uword n_slices;
  const uword n_elem;
  
  
  protected:
  
  arma_inline subview_cube(const Cube<eT>& in_m, const uword in_row1, const uword in_col1, const uword in_slice1, const uword in_n_rows, const uword in_n_cols, const uword in_n_slices);
  
  
  public:
  
  inline ~subview_cube();
  inline  subview_cube() = delete;
  
  inline  subview_cube(const subview_cube&  in);
  inline  subview_cube(      subview_cube&& in);
  
  template<typename op_type             > inline void inplace_op(const eT                val                        );
  template<typename op_type, typename T1> inline void inplace_op(const BaseCube<eT,T1>&  x,   const char* identifier);
  template<typename op_type             > inline void inplace_op(const subview_cube<eT>& x,   const char* identifier);
  
  inline void operator=  (const eT val);
  inline void operator+= (const eT val);
  inline void operator-= (const eT val);
  inline void operator*= (const eT val);
  inline void operator/= (const eT val);
  
  // deliberately returning void
  template<typename T1> inline void operator=  (const BaseCube<eT,T1>& x);
  template<typename T1> inline void operator+= (const BaseCube<eT,T1>& x);
  template<typename T1> inline void operator-= (const BaseCube<eT,T1>& x);
  template<typename T1> inline void operator%= (const BaseCube<eT,T1>& x);
  template<typename T1> inline void operator/= (const BaseCube<eT,T1>& x);
  
  inline void operator=  (const subview_cube& x);
  inline void operator+= (const subview_cube& x);
  inline void operator-= (const subview_cube& x);
  inline void operator%= (const subview_cube& x);
  inline void operator/= (const subview_cube& x);
  
  template<typename T1> inline void operator=  (const Base<eT,T1>& x);
  template<typename T1> inline void operator+= (const Base<eT,T1>& x);
  template<typename T1> inline void operator-= (const Base<eT,T1>& x);
  template<typename T1> inline void operator%= (const Base<eT,T1>& x);
  template<typename T1> inline void operator/= (const Base<eT,T1>& x);
  
  template<typename gen_type> inline void operator=(const GenCube<eT,gen_type>& x);
  
  inline static void       extract(Cube<eT>& out, const subview_cube& in);
  inline static void  plus_inplace(Cube<eT>& out, const subview_cube& in);
  inline static void minus_inplace(Cube<eT>& out, const subview_cube& in);
  inline static void schur_inplace(Cube<eT>& out, const subview_cube& in);
  inline static void   div_inplace(Cube<eT>& out, const subview_cube& in);
  
  inline static void       extract(Mat<eT>& out, const subview_cube& in);
  inline static void  plus_inplace(Mat<eT>& out, const subview_cube& in);
  inline static void minus_inplace(Mat<eT>& out, const subview_cube& in);
  inline static void schur_inplace(Mat<eT>& out, const subview_cube& in);
  inline static void   div_inplace(Mat<eT>& out, const subview_cube& in);
  
  template<typename functor> inline void  for_each(functor F);
  template<typename functor> inline void  for_each(functor F) const;
  
  template<typename functor> inline void transform(functor F);
  template<typename functor> inline void     imbue(functor F);
  
  inline void each_slice(const std::function< void(      Mat<eT>&) >& F);
  inline void each_slice(const std::function< void(const Mat<eT>&) >& F) const;
  
  inline void replace(const eT old_val, const eT new_val);
  
  inline void clean(const pod_type threshold);
  
  inline void clamp(const eT min_val, const eT max_val);
  
  inline void fill(const eT val);
  inline void zeros();
  inline void ones();
  inline void randu();
  inline void randn();
  
  inline arma_warn_unused bool is_finite() const;
  inline arma_warn_unused bool is_zero(const pod_type tol = 0) const;
  
  inline arma_warn_unused bool has_inf() const;
  inline arma_warn_unused bool has_nan() const;
  
  inline eT  at_alt    (const uword i) const;
  
  inline eT& operator[](const uword i);
  inline eT  operator[](const uword i) const;
  
  inline eT& operator()(const uword i);
  inline eT  operator()(const uword i) const;
  
  arma_inline eT& operator()(const uword in_row, const uword in_col, const uword in_slice);
  arma_inline eT  operator()(const uword in_row, const uword in_col, const uword in_slice) const;
  
  arma_inline eT&         at(const uword in_row, const uword in_col, const uword in_slice);
  arma_inline eT          at(const uword in_row, const uword in_col, const uword in_slice) const;
  
  arma_inline       eT* slice_colptr(const uword in_slice, const uword in_col);
  arma_inline const eT* slice_colptr(const uword in_slice, const uword in_col) const;
  
  template<typename eT2>
  inline bool check_overlap(const subview_cube<eT2>& x) const;
  
  inline bool check_overlap(const Mat<eT>&           x) const;
  
  
  class const_iterator;
  
  class iterator
    {
    public:
    
    inline iterator();
    inline iterator(const iterator& X);
    inline iterator(subview_cube<eT>& in_sv, const uword in_row, const uword in_col, const uword in_slice);
    
    inline arma_warn_unused eT& operator*();
    
    inline                  iterator& operator++();
    inline arma_warn_unused iterator  operator++(int);
    
    inline arma_warn_unused bool operator==(const       iterator& rhs) const;
    inline arma_warn_unused bool operator!=(const       iterator& rhs) const;
    inline arma_warn_unused bool operator==(const const_iterator& rhs) const;
    inline arma_warn_unused bool operator!=(const const_iterator& rhs) const;
    
    typedef std::forward_iterator_tag iterator_category;
    typedef eT                        value_type;
    typedef std::ptrdiff_t            difference_type;  // TODO: not certain on this one
    typedef eT*                       pointer;
    typedef eT&                       reference;
    
    arma_aligned Cube<eT>* M;
    arma_aligned eT*       current_ptr;
    arma_aligned uword     current_row;
    arma_aligned uword     current_col;
    arma_aligned uword     current_slice;
    
    arma_aligned const uword aux_row1;
    arma_aligned const uword aux_col1;
    
    arma_aligned const uword aux_row2_p1;
    arma_aligned const uword aux_col2_p1;
    };
  
  
  class const_iterator
    {
    public:
    
    inline const_iterator();
    inline const_iterator(const       iterator& X);
    inline const_iterator(const const_iterator& X);
    inline const_iterator(const subview_cube<eT>& in_sv, const uword in_row, const uword in_col, const uword in_slice);
    
    inline arma_warn_unused const eT& operator*();
    
    inline                  const_iterator& operator++();
    inline arma_warn_unused const_iterator  operator++(int);
    
    inline arma_warn_unused bool operator==(const       iterator& rhs) const;
    inline arma_warn_unused bool operator!=(const       iterator& rhs) const;
    inline arma_warn_unused bool operator==(const const_iterator& rhs) const;
    inline arma_warn_unused bool operator!=(const const_iterator& rhs) const;
    
    // So that we satisfy the STL iterator types.
    typedef std::forward_iterator_tag iterator_category;
    typedef eT                        value_type;
    typedef std::ptrdiff_t            difference_type;  // TODO: not certain on this one
    typedef const eT*                 pointer;
    typedef const eT&                 reference;
    
    arma_aligned const Cube<eT>* M;
    arma_aligned const eT*       current_ptr;
    arma_aligned       uword     current_row;
    arma_aligned       uword     current_col;
    arma_aligned       uword     current_slice;
    
    arma_aligned const uword aux_row1;
    arma_aligned const uword aux_col1;
    
    arma_aligned const uword aux_row2_p1;
    arma_aligned const uword aux_col2_p1;
    };
  
  
  inline       iterator  begin();
  inline const_iterator  begin() const;
  inline const_iterator cbegin() const;
  
  inline       iterator  end();
  inline const_iterator  end() const;
  inline const_iterator cend() const;
  
  
  friend class  Mat<eT>;
  friend class Cube<eT>;
  };



//! @}
