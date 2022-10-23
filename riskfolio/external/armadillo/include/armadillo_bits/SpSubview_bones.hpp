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


//! \addtogroup SpSubview
//! @{


template<typename eT>
class SpSubview : public SpBase< eT, SpSubview<eT> >
  {
  public:
  
  const SpMat<eT>& m;
  
  typedef eT elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const uword aux_row1;
  const uword aux_col1;
  const uword n_rows;
  const uword n_cols;
  const uword n_elem;
  const uword n_nonzero;
  
  protected:
  
  inline SpSubview(const SpMat<eT>& in_m, const uword in_row1, const uword in_col1, const uword in_n_rows, const uword in_n_cols);
  
  public:
  
  inline ~SpSubview();
  inline  SpSubview() = delete;
  
  inline  SpSubview(const SpSubview&  in);
  inline  SpSubview(      SpSubview&& in);
  
  inline const SpSubview& operator+= (const eT val);
  inline const SpSubview& operator-= (const eT val);
  inline const SpSubview& operator*= (const eT val);
  inline const SpSubview& operator/= (const eT val);
  
  inline const SpSubview& operator=(const SpSubview& x);
  
  template<typename T1> inline const SpSubview& operator= (const Base<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator+=(const Base<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator-=(const Base<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator*=(const Base<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator%=(const Base<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator/=(const Base<eT, T1>& x);
  
  template<typename T1> inline const SpSubview& operator_equ_common(const SpBase<eT, T1>& x);
  
  template<typename T1> inline const SpSubview& operator= (const SpBase<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator+=(const SpBase<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator-=(const SpBase<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator*=(const SpBase<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator%=(const SpBase<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator/=(const SpBase<eT, T1>& x);
  
  /*
  inline static void extract(SpMat<eT>& out, const SpSubview& in);
  
  inline static void  plus_inplace(Mat<eT>& out, const subview& in);
  inline static void minus_inplace(Mat<eT>& out, const subview& in);
  inline static void schur_inplace(Mat<eT>& out, const subview& in);
  inline static void   div_inplace(Mat<eT>& out, const subview& in);
  */
  
  template<typename functor> inline void  for_each(functor F);
  template<typename functor> inline void  for_each(functor F) const;
  
  template<typename functor> inline void transform(functor F);
  
  inline void replace(const eT old_val, const eT new_val);
  
  inline void clean(const pod_type threshold);
  
  inline void clamp(const eT min_val, const eT max_val);
  
  inline void fill(const eT val);
  inline void zeros();
  inline void ones();
  inline void eye();
  inline void randu();
  inline void randn();
  
  
  arma_hot inline SpSubview_MapMat_val<eT> operator[](const uword i);
  arma_hot inline eT                       operator[](const uword i) const;
  
  arma_hot inline SpSubview_MapMat_val<eT> operator()(const uword i);
  arma_hot inline eT                       operator()(const uword i) const;
  
  arma_hot inline SpSubview_MapMat_val<eT> operator()(const uword in_row, const uword in_col);
  arma_hot inline eT                       operator()(const uword in_row, const uword in_col) const;
  
  arma_hot inline SpSubview_MapMat_val<eT> at(const uword i);
  arma_hot inline eT                       at(const uword i) const;
  
  arma_hot inline SpSubview_MapMat_val<eT> at(const uword in_row, const uword in_col);
  arma_hot inline eT                       at(const uword in_row, const uword in_col) const;
  
  inline bool check_overlap(const SpSubview& x) const;
  
  inline bool is_vec() const;
  
  inline       SpSubview_row<eT> row(const uword row_num);
  inline const SpSubview_row<eT> row(const uword row_num) const;
  
  inline       SpSubview_col<eT> col(const uword col_num);
  inline const SpSubview_col<eT> col(const uword col_num) const;
  
  inline       SpSubview rows(const uword in_row1, const uword in_row2);
  inline const SpSubview rows(const uword in_row1, const uword in_row2) const;
  
  inline       SpSubview cols(const uword in_col1, const uword in_col2);
  inline const SpSubview cols(const uword in_col1, const uword in_col2) const;
  
  inline       SpSubview submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2);
  inline const SpSubview submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2) const;
  
  inline       SpSubview submat(const span& row_span, const span& col_span);
  inline const SpSubview submat(const span& row_span, const span& col_span) const;
  
  inline       SpSubview operator()(const uword row_num, const span& col_span);
  inline const SpSubview operator()(const uword row_num, const span& col_span) const;
  
  inline       SpSubview operator()(const span& row_span, const uword col_num);
  inline const SpSubview operator()(const span& row_span, const uword col_num) const;
  
  inline       SpSubview operator()(const span& row_span, const span& col_span);
  inline const SpSubview operator()(const span& row_span, const span& col_span) const;
  
  
  inline void swap_rows(const uword in_row1, const uword in_row2);
  inline void swap_cols(const uword in_col1, const uword in_col2);
  
  // Forward declarations.
  class iterator_base;
  class const_iterator;
  class iterator;
  class const_row_iterator;
  class row_iterator;
  
  // Similar to SpMat iterators but automatically iterates past and ignores values not in the subview.
  class iterator_base
    {
    public:
    
    inline iterator_base(const SpSubview& in_M);
    inline iterator_base(const SpSubview& in_M, const uword col, const uword pos);
    
    arma_inline uword col() const { return internal_col; }
    arma_inline uword pos() const { return internal_pos; }
    
    arma_aligned const SpSubview* M;
    arma_aligned       uword      internal_col;
    arma_aligned       uword      internal_pos;
    
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef eT                              value_type;
    typedef std::ptrdiff_t                  difference_type;  // TODO: not certain on this one
    typedef const eT*                       pointer;
    typedef const eT&                       reference;
    };
  
  class const_iterator : public iterator_base
    {
    public:
    
    inline const_iterator(const SpSubview& in_M, uword initial_pos = 0);
    inline const_iterator(const SpSubview& in_M, uword in_row, uword in_col);
    inline const_iterator(const SpSubview& in_M, uword in_row, uword in_col, uword in_pos, uword skip_pos);
    inline const_iterator(const const_iterator& other);
    
    arma_inline eT operator*() const;
    
    // Don't hold location internally; call "dummy" methods to get that information.
    arma_inline uword row() const { return iterator_base::M->m.row_indices[iterator_base::internal_pos + skip_pos] - iterator_base::M->aux_row1; }
    
    inline arma_hot         const_iterator& operator++();
    inline arma_warn_unused const_iterator  operator++(int);
    
    inline arma_hot         const_iterator& operator--();
    inline arma_warn_unused const_iterator  operator--(int);
    
    inline arma_hot bool operator!=(const const_iterator& rhs) const;
    inline arma_hot bool operator==(const const_iterator& rhs) const;
    
    inline arma_hot bool operator!=(const typename SpMat<eT>::const_iterator& rhs) const;
    inline arma_hot bool operator==(const typename SpMat<eT>::const_iterator& rhs) const;
    
    inline arma_hot bool operator!=(const const_row_iterator& rhs) const;
    inline arma_hot bool operator==(const const_row_iterator& rhs) const;
    
    inline arma_hot bool operator!=(const typename SpMat<eT>::const_row_iterator& rhs) const;
    inline arma_hot bool operator==(const typename SpMat<eT>::const_row_iterator& rhs) const;
    
    arma_aligned uword skip_pos; // not used in row_iterator or const_row_iterator
    };
  
  class iterator : public const_iterator
    {
    public:
    
    inline iterator(SpSubview& in_M, const uword initial_pos = 0) : const_iterator(in_M, initial_pos) { }
    inline iterator(SpSubview& in_M, const uword in_row, const uword in_col) : const_iterator(in_M, in_row, in_col) { }
    inline iterator(SpSubview& in_M, const uword in_row, const uword in_col, const uword in_pos, const uword in_skip_pos) : const_iterator(in_M, in_row, in_col, in_pos, in_skip_pos) { }
    inline iterator(const iterator& other) : const_iterator(other) { }
    
    inline arma_hot SpValProxy< SpSubview<eT> > operator*();
    
    // overloads needed for return type correctness
    inline arma_hot         iterator& operator++();
    inline arma_warn_unused iterator  operator++(int);
    
    inline arma_hot         iterator& operator--();
    inline arma_warn_unused iterator  operator--(int);
    
    // This has a different value_type than iterator_base.
    typedef       SpValProxy< SpSubview<eT> >  value_type;
    typedef const SpValProxy< SpSubview<eT> >* pointer;
    typedef const SpValProxy< SpSubview<eT> >& reference;
    };
  
  class const_row_iterator : public iterator_base
    {
    public:
    
    inline const_row_iterator();
    inline const_row_iterator(const SpSubview& in_M, uword initial_pos = 0);
    inline const_row_iterator(const SpSubview& in_M, uword in_row, uword in_col);
    inline const_row_iterator(const const_row_iterator& other);
    
    inline arma_hot         const_row_iterator& operator++();
    inline arma_warn_unused const_row_iterator  operator++(int);
    
    inline arma_hot         const_row_iterator& operator--();
    inline arma_warn_unused const_row_iterator  operator--(int);
    
    uword internal_row; // Hold row internally because we use internal_pos differently.
    uword actual_pos; // Actual position in subview's parent matrix.
    
    arma_inline eT operator*() const { return iterator_base::M->m.values[actual_pos]; }
    
    arma_inline uword row() const { return internal_row; }
    
    inline arma_hot bool operator!=(const const_iterator& rhs) const;
    inline arma_hot bool operator==(const const_iterator& rhs) const;
    
    inline arma_hot bool operator!=(const typename SpMat<eT>::const_iterator& rhs) const;
    inline arma_hot bool operator==(const typename SpMat<eT>::const_iterator& rhs) const;
    
    inline arma_hot bool operator!=(const const_row_iterator& rhs) const;
    inline arma_hot bool operator==(const const_row_iterator& rhs) const;
    
    inline arma_hot bool operator!=(const typename SpMat<eT>::const_row_iterator& rhs) const;
    inline arma_hot bool operator==(const typename SpMat<eT>::const_row_iterator& rhs) const;
    };
  
  class row_iterator : public const_row_iterator
    {
    public:
    
    inline row_iterator(SpSubview& in_M, uword initial_pos = 0) : const_row_iterator(in_M, initial_pos) { }
    inline row_iterator(SpSubview& in_M, uword in_row, uword in_col) : const_row_iterator(in_M, in_row, in_col) { }
    inline row_iterator(const row_iterator& other) : const_row_iterator(other) { }
    
    inline arma_hot SpValProxy< SpSubview<eT> > operator*();
    
    // overloads needed for return type correctness
    inline arma_hot         row_iterator& operator++();
    inline arma_warn_unused row_iterator  operator++(int);
    
    inline arma_hot         row_iterator& operator--();
    inline arma_warn_unused row_iterator  operator--(int);
    
    // This has a different value_type than iterator_base.
    typedef       SpValProxy< SpSubview<eT> >  value_type;
    typedef const SpValProxy< SpSubview<eT> >* pointer;
    typedef const SpValProxy< SpSubview<eT> >& reference;
    };
  
  inline       iterator     begin();
  inline const_iterator     begin() const;
  inline const_iterator    cbegin() const;
  
  inline       iterator     begin_col(const uword col_num);
  inline const_iterator     begin_col(const uword col_num) const;
  
  inline       row_iterator begin_row(const uword row_num = 0);
  inline const_row_iterator begin_row(const uword row_num = 0) const;
  
  inline       iterator     end();
  inline const_iterator     end() const;
  inline const_iterator    cend() const;
  
  inline       row_iterator end_row();
  inline const_row_iterator end_row() const;
  
  inline       row_iterator end_row(const uword row_num);
  inline const_row_iterator end_row(const uword row_num) const;
  
  //! don't use this unless you're writing internal Armadillo code
  arma_inline bool is_alias(const SpMat<eT>& X) const;
  
  
  private:
  
  friend class SpMat<eT>;
  friend class SpSubview_col<eT>;
  friend class SpSubview_row<eT>;
  friend class SpValProxy< SpSubview<eT> >;   // allow SpValProxy to call insert_element() and delete_element()
  
  inline arma_warn_unused eT&  insert_element(const uword in_row, const uword in_col, const eT in_val = eT(0));
  inline                  void delete_element(const uword in_row, const uword in_col);
  
  inline void invalidate_cache() const;
  };



template<typename eT>
class SpSubview_col : public SpSubview<eT>
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  inline void operator= (const SpSubview<eT>& x);
  inline void operator= (const SpSubview_col& x);
  
  template<typename T1> inline void operator= (const SpBase<eT,T1>& x);
  template<typename T1> inline void operator= (const   Base<eT,T1>& x);
  
  inline arma_warn_unused const SpOp<SpSubview_col<eT>,spop_htrans>  t() const;
  inline arma_warn_unused const SpOp<SpSubview_col<eT>,spop_htrans> ht() const;
  inline arma_warn_unused const SpOp<SpSubview_col<eT>,spop_strans> st() const;
  
  
  protected:
  
  inline SpSubview_col(const SpMat<eT>& in_m, const uword in_col);
  inline SpSubview_col(const SpMat<eT>& in_m, const uword in_col, const uword in_row1, const uword in_n_rows);
  inline SpSubview_col() = delete;
  
  
  private:
  
  friend class SpMat<eT>;
  friend class SpSubview<eT>;
  };



template<typename eT>
class SpSubview_row : public SpSubview<eT>
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static constexpr bool is_row  = true;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  inline void operator= (const SpSubview<eT>& x);
  inline void operator= (const SpSubview_row& x);
  
  template<typename T1> inline void operator= (const SpBase<eT,T1>& x);
  template<typename T1> inline void operator= (const   Base<eT,T1>& x);
  
  inline arma_warn_unused const SpOp<SpSubview_row<eT>,spop_htrans>  t() const;
  inline arma_warn_unused const SpOp<SpSubview_row<eT>,spop_htrans> ht() const;
  inline arma_warn_unused const SpOp<SpSubview_row<eT>,spop_strans> st() const;
  
  
  protected:
  
  inline SpSubview_row(const SpMat<eT>& in_m, const uword in_row);
  inline SpSubview_row(const SpMat<eT>& in_m, const uword in_row, const uword in_col1, const uword in_n_cols);
  inline SpSubview_row() = delete;
  
  
  private:
  
  friend class SpMat<eT>;
  friend class SpSubview<eT>;
  };



//! @}
