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


//! \addtogroup subview
//! @{


//! Class for storing data required to construct or apply operations to a submatrix
//! (ie. where the submatrix starts and ends as well as a reference/pointer to the original matrix),
template<typename eT>
class subview : public Base< eT, subview<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  arma_aligned const Mat<eT>& m;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const uword aux_row1;
  const uword aux_col1;
  
  const uword n_rows;
  const uword n_cols;
  const uword n_elem;
  
  protected:
  
  arma_inline subview(const Mat<eT>& in_m, const uword in_row1, const uword in_col1, const uword in_n_rows, const uword in_n_cols);
  
  
  public:
  
  inline ~subview();
  inline  subview() = delete;
  
  inline  subview(const subview&  in);
  inline  subview(      subview&& in);
  
  template<typename op_type             > inline void inplace_op(const eT           val                        );
  template<typename op_type, typename T1> inline void inplace_op(const Base<eT,T1>& x,   const char* identifier);
  template<typename op_type             > inline void inplace_op(const subview<eT>& x,   const char* identifier);
  
  // deliberately returning void
  
  inline void operator=  (const eT val);
  inline void operator+= (const eT val);
  inline void operator-= (const eT val);
  inline void operator*= (const eT val);
  inline void operator/= (const eT val);
  
  inline void operator=  (const subview& x);
  inline void operator+= (const subview& x);
  inline void operator-= (const subview& x);
  inline void operator%= (const subview& x);
  inline void operator/= (const subview& x);
  
  template<typename T1> inline void operator=  (const Base<eT,T1>& x);
  template<typename T1> inline void operator+= (const Base<eT,T1>& x);
  template<typename T1> inline void operator-= (const Base<eT,T1>& x);
  template<typename T1> inline void operator%= (const Base<eT,T1>& x);
  template<typename T1> inline void operator/= (const Base<eT,T1>& x);
  
  template<typename T1> inline void operator=  (const SpBase<eT,T1>& x);
  template<typename T1> inline void operator+= (const SpBase<eT,T1>& x);
  template<typename T1> inline void operator-= (const SpBase<eT,T1>& x);
  template<typename T1> inline void operator%= (const SpBase<eT,T1>& x);
  template<typename T1> inline void operator/= (const SpBase<eT,T1>& x);

  template<typename T1, typename gen_type>
  inline typename enable_if2< is_same_type<typename T1::elem_type, eT>::value, void>::result operator=(const Gen<T1,gen_type>& x);
  
  inline void operator=(const std::initializer_list<eT>& list);
  inline void operator=(const std::initializer_list< std::initializer_list<eT> >& list);
  
  
  inline static void extract(Mat<eT>& out, const subview& in);
  
  inline static void  plus_inplace(Mat<eT>& out, const subview& in);
  inline static void minus_inplace(Mat<eT>& out, const subview& in);
  inline static void schur_inplace(Mat<eT>& out, const subview& in);
  inline static void   div_inplace(Mat<eT>& out, const subview& in);
  
  template<typename functor> inline void  for_each(functor F);
  template<typename functor> inline void  for_each(functor F) const;
  
  template<typename functor> inline void transform(functor F);
  template<typename functor> inline void     imbue(functor F);
  
  inline void replace(const eT old_val, const eT new_val);
  
  inline void clean(const pod_type threshold);
  
  inline void clamp(const eT min_val, const eT max_val);
  
  inline void fill(const eT val);
  inline void zeros();
  inline void ones();
  inline void eye();
  inline void randu();
  inline void randn();
  
  inline arma_warn_unused eT  at_alt    (const uword ii) const;
  
  inline arma_warn_unused eT& operator[](const uword ii);
  inline arma_warn_unused eT  operator[](const uword ii) const;
  
  inline arma_warn_unused eT& operator()(const uword ii);
  inline arma_warn_unused eT  operator()(const uword ii) const;
  
  inline arma_warn_unused eT& operator()(const uword in_row, const uword in_col);
  inline arma_warn_unused eT  operator()(const uword in_row, const uword in_col) const;
  
  inline arma_warn_unused eT&         at(const uword in_row, const uword in_col);
  inline arma_warn_unused eT          at(const uword in_row, const uword in_col) const;
  
  inline arma_warn_unused eT& front();
  inline arma_warn_unused eT  front() const;
  
  inline arma_warn_unused eT& back();
  inline arma_warn_unused eT  back() const;
  
  arma_inline       eT* colptr(const uword in_col);
  arma_inline const eT* colptr(const uword in_col) const;
  
  template<typename eT2>
  inline bool check_overlap(const subview<eT2>& x) const;
  
  inline arma_warn_unused bool is_vec()    const;
  inline arma_warn_unused bool is_finite() const;
  inline arma_warn_unused bool is_zero(const pod_type tol = 0) const;
  
  inline arma_warn_unused bool has_inf() const;
  inline arma_warn_unused bool has_nan() const;
  
  inline       subview_row<eT> row(const uword row_num);
  inline const subview_row<eT> row(const uword row_num) const;
  
  inline            subview_row<eT> operator()(const uword row_num, const span& col_span);
  inline      const subview_row<eT> operator()(const uword row_num, const span& col_span) const;
  
  inline       subview_col<eT> col(const uword col_num);
  inline const subview_col<eT> col(const uword col_num) const;
  
  inline            subview_col<eT> operator()(const span& row_span, const uword col_num);
  inline      const subview_col<eT> operator()(const span& row_span, const uword col_num) const;
  
  inline            Col<eT>  unsafe_col(const uword col_num);
  inline      const Col<eT>  unsafe_col(const uword col_num) const;
  
  inline       subview<eT> rows(const uword in_row1, const uword in_row2);
  inline const subview<eT> rows(const uword in_row1, const uword in_row2) const;
  
  inline       subview<eT> cols(const uword in_col1, const uword in_col2);
  inline const subview<eT> cols(const uword in_col1, const uword in_col2) const;
  
  inline       subview<eT> submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2);
  inline const subview<eT> submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2) const;
  
  inline            subview<eT> submat    (const span& row_span, const span& col_span);
  inline      const subview<eT> submat    (const span& row_span, const span& col_span) const;
  
  inline            subview<eT> operator()(const span& row_span, const span& col_span);
  inline      const subview<eT> operator()(const span& row_span, const span& col_span) const;
  
  inline subview_each1< subview<eT>, 0 > each_col();
  inline subview_each1< subview<eT>, 1 > each_row();
  
  template<typename T1> inline subview_each2< subview<eT>, 0, T1 > each_col(const Base<uword, T1>& indices);
  template<typename T1> inline subview_each2< subview<eT>, 1, T1 > each_row(const Base<uword, T1>& indices);
  
  inline void each_col(const std::function< void(      Col<eT>&) >& F);
  inline void each_col(const std::function< void(const Col<eT>&) >& F) const;
  
  inline void each_row(const std::function< void(      Row<eT>&) >& F);
  inline void each_row(const std::function< void(const Row<eT>&) >& F) const;
  
  inline       diagview<eT> diag(const sword in_id = 0);
  inline const diagview<eT> diag(const sword in_id = 0) const;
  
  inline void swap_rows(const uword in_row1, const uword in_row2);
  inline void swap_cols(const uword in_col1, const uword in_col2);
  
  
  class const_iterator;
  
  class iterator
    {
    public:
    
    inline iterator();
    inline iterator(const iterator& X);
    inline iterator(subview<eT>& in_sv, const uword in_row, const uword in_col);
    
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
    
    arma_aligned Mat<eT>* M;
    arma_aligned eT*      current_ptr;
    arma_aligned uword    current_row;
    arma_aligned uword    current_col;
    
    arma_aligned const uword aux_row1;
    arma_aligned const uword aux_row2_p1;
    };
  
  
  class const_iterator
    {
    public:
    
    inline const_iterator();
    inline const_iterator(const       iterator& X);
    inline const_iterator(const const_iterator& X);
    inline const_iterator(const subview<eT>& in_sv, const uword in_row, const uword in_col);
    
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
    
    arma_aligned const Mat<eT>* M;
    arma_aligned const eT*      current_ptr;
    arma_aligned       uword    current_row;
    arma_aligned       uword    current_col;
    
    arma_aligned const uword aux_row1;
    arma_aligned const uword aux_row2_p1;
    };
  
  
  class const_row_iterator;
  
  class row_iterator
    {
    public:
    
    inline row_iterator();
    inline row_iterator(const row_iterator& X);
    inline row_iterator(subview<eT>& in_sv, const uword in_row, const uword in_col);
    
    inline arma_warn_unused eT& operator* ();
    
    inline                  row_iterator& operator++();
    inline arma_warn_unused row_iterator  operator++(int);
    
    inline arma_warn_unused bool operator!=(const       row_iterator& X) const;
    inline arma_warn_unused bool operator==(const       row_iterator& X) const;
    inline arma_warn_unused bool operator!=(const const_row_iterator& X) const;
    inline arma_warn_unused bool operator==(const const_row_iterator& X) const;
    
    typedef std::forward_iterator_tag iterator_category;
    typedef eT                        value_type;
    typedef std::ptrdiff_t            difference_type;  // TODO: not certain on this one
    typedef eT*                       pointer;
    typedef eT&                       reference;
    
    arma_aligned Mat<eT>* M;
    arma_aligned uword    current_row;
    arma_aligned uword    current_col;
    
    arma_aligned const uword aux_col1;
    arma_aligned const uword aux_col2_p1;
    };
  
  
  class const_row_iterator
    {
    public:
    
    inline const_row_iterator();
    inline const_row_iterator(const       row_iterator& X);
    inline const_row_iterator(const const_row_iterator& X);
    inline const_row_iterator(const subview<eT>& in_sv, const uword in_row, const uword in_col);
    
    inline arma_warn_unused const eT& operator*() const;
    
    inline                  const_row_iterator& operator++();
    inline arma_warn_unused const_row_iterator  operator++(int);
    
    inline arma_warn_unused bool operator!=(const       row_iterator& X) const;
    inline arma_warn_unused bool operator==(const       row_iterator& X) const;
    inline arma_warn_unused bool operator!=(const const_row_iterator& X) const;
    inline arma_warn_unused bool operator==(const const_row_iterator& X) const;
    
    typedef std::forward_iterator_tag iterator_category;
    typedef eT                        value_type;
    typedef std::ptrdiff_t            difference_type;  // TODO: not certain on this one
    typedef const eT*                 pointer;
    typedef const eT&                 reference;
    
    arma_aligned const Mat<eT>* M;
    arma_aligned       uword    current_row;
    arma_aligned       uword    current_col;
    
    arma_aligned const uword aux_col1;
    arma_aligned const uword aux_col2_p1;
    };
  
  
  
  inline       iterator  begin();
  inline const_iterator  begin() const;
  inline const_iterator cbegin() const;
  
  inline       iterator  end();
  inline const_iterator  end() const;
  inline const_iterator cend() const;
  
  
  friend class Mat<eT>;
  };



template<typename eT>
class subview_col : public subview<eT>
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const eT* colmem;
  
  inline void operator= (const subview<eT>& x);
  inline void operator= (const subview_col& x);
  inline void operator= (const eT val);
  inline void operator= (const std::initializer_list<eT>& list);
  
  template<typename T1> inline void operator= (const   Base<eT,T1>& x);
  template<typename T1> inline void operator= (const SpBase<eT,T1>& x);
  
  template<typename T1, typename gen_type>
  inline typename enable_if2< is_same_type<typename T1::elem_type, eT>::value, void>::result operator=(const Gen<T1,gen_type>& x);
  
  arma_inline arma_warn_unused const Op<subview_col<eT>,op_htrans>  t() const;
  arma_inline arma_warn_unused const Op<subview_col<eT>,op_htrans> ht() const;
  arma_inline arma_warn_unused const Op<subview_col<eT>,op_strans> st() const;
  
  arma_inline arma_warn_unused const Op<subview_col<eT>,op_strans> as_row() const;
  
  inline void fill(const eT val);
  inline void zeros();
  inline void ones();
  
  arma_inline eT  at_alt    (const uword i) const;
  
  arma_inline eT& operator[](const uword i);
  arma_inline eT  operator[](const uword i) const;
  
  inline eT& operator()(const uword i);
  inline eT  operator()(const uword i) const;
  
  inline eT& operator()(const uword in_row, const uword in_col);
  inline eT  operator()(const uword in_row, const uword in_col) const;
  
  inline eT&         at(const uword in_row, const uword in_col);
  inline eT          at(const uword in_row, const uword in_col) const;
  
  arma_inline       eT* colptr(const uword in_col);
  arma_inline const eT* colptr(const uword in_col) const;
  
  inline       subview_col<eT> rows(const uword in_row1, const uword in_row2);
  inline const subview_col<eT> rows(const uword in_row1, const uword in_row2) const;
  
  inline       subview_col<eT> subvec(const uword in_row1, const uword in_row2);
  inline const subview_col<eT> subvec(const uword in_row1, const uword in_row2) const;
  
  inline       subview_col<eT> subvec(const uword start_row, const SizeMat& s);
  inline const subview_col<eT> subvec(const uword start_row, const SizeMat& s) const;
  
  inline       subview_col<eT> head(const uword N);
  inline const subview_col<eT> head(const uword N) const;
  
  inline       subview_col<eT> tail(const uword N);
  inline const subview_col<eT> tail(const uword N) const;
  
  inline arma_warn_unused eT min() const;
  inline arma_warn_unused eT max() const;
  
  inline eT min(uword& index_of_min_val) const;
  inline eT max(uword& index_of_max_val) const;
  
  inline arma_warn_unused uword index_min() const;
  inline arma_warn_unused uword index_max() const;
  
  inline  subview_col(const subview_col&  in);
  inline  subview_col(      subview_col&& in);
  
  
  protected:
  
  inline subview_col(const Mat<eT>& in_m, const uword in_col);
  inline subview_col(const Mat<eT>& in_m, const uword in_col, const uword in_row1, const uword in_n_rows);
  inline subview_col() = delete;
  
  
  friend class Mat<eT>;
  friend class Col<eT>;
  friend class subview<eT>;
  };



template<typename eT>
class subview_cols : public subview<eT>
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  inline  subview_cols(const subview_cols&  in);
  inline  subview_cols(      subview_cols&& in);
  
  inline void operator= (const subview<eT>&  x);
  inline void operator= (const subview_cols& x);
  inline void operator= (const eT val);
  inline void operator= (const std::initializer_list<eT>& list);
  inline void operator= (const std::initializer_list< std::initializer_list<eT> >& list);
  
  template<typename T1> inline void operator= (const   Base<eT,T1>& x);
  template<typename T1> inline void operator= (const SpBase<eT,T1>& x);
  
  template<typename T1, typename gen_type>
  inline typename enable_if2< is_same_type<typename T1::elem_type, eT>::value, void>::result operator=(const Gen<T1,gen_type>& x);
  
  arma_inline arma_warn_unused const Op<subview_cols<eT>,op_htrans>  t() const;
  arma_inline arma_warn_unused const Op<subview_cols<eT>,op_htrans> ht() const;
  arma_inline arma_warn_unused const Op<subview_cols<eT>,op_strans> st() const;
  
  arma_inline arma_warn_unused const Op<subview_cols<eT>,op_vectorise_col> as_col() const;
  
  inline arma_warn_unused eT  at_alt    (const uword ii) const;
  
  inline arma_warn_unused eT& operator[](const uword ii);
  inline arma_warn_unused eT  operator[](const uword ii) const;
  
  inline arma_warn_unused eT& operator()(const uword ii);
  inline arma_warn_unused eT  operator()(const uword ii) const;
  
  inline arma_warn_unused eT& operator()(const uword in_row, const uword in_col);
  inline arma_warn_unused eT  operator()(const uword in_row, const uword in_col) const;
  
  inline arma_warn_unused eT&         at(const uword in_row, const uword in_col);
  inline arma_warn_unused eT          at(const uword in_row, const uword in_col) const;
  
  arma_inline       eT* colptr(const uword in_col);
  arma_inline const eT* colptr(const uword in_col) const;
  
  protected:
  
  inline subview_cols(const Mat<eT>& in_m, const uword in_col1, const uword in_n_cols);
  inline subview_cols() = delete;
  
  friend class Mat<eT>;
  friend class subview<eT>;
  };



template<typename eT>
class subview_row : public subview<eT>
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static constexpr bool is_row  = true;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  inline void operator= (const subview<eT>& x);
  inline void operator= (const subview_row& x);
  inline void operator= (const eT val);
  inline void operator= (const std::initializer_list<eT>& list);
  
  template<typename T1> inline void operator= (const   Base<eT,T1>& x);
  template<typename T1> inline void operator= (const SpBase<eT,T1>& x);
  
  template<typename T1, typename gen_type>
  inline typename enable_if2< is_same_type<typename T1::elem_type, eT>::value, void>::result operator=(const Gen<T1,gen_type>& x);
  
  arma_inline arma_warn_unused const Op<subview_row<eT>,op_htrans>  t() const;
  arma_inline arma_warn_unused const Op<subview_row<eT>,op_htrans> ht() const;
  arma_inline arma_warn_unused const Op<subview_row<eT>,op_strans> st() const;
  
  arma_inline arma_warn_unused const Op<subview_row<eT>,op_strans> as_col() const;
  
  inline eT  at_alt    (const uword i) const;
  
  inline eT& operator[](const uword i);
  inline eT  operator[](const uword i) const;
  
  inline eT& operator()(const uword i);
  inline eT  operator()(const uword i) const;
  
  inline eT& operator()(const uword in_row, const uword in_col);
  inline eT  operator()(const uword in_row, const uword in_col) const;
  
  inline eT&         at(const uword in_row, const uword in_col);
  inline eT          at(const uword in_row, const uword in_col) const;
  
  inline       subview_row<eT> cols(const uword in_col1, const uword in_col2);
  inline const subview_row<eT> cols(const uword in_col1, const uword in_col2) const;
  
  inline       subview_row<eT> subvec(const uword in_col1, const uword in_col2);
  inline const subview_row<eT> subvec(const uword in_col1, const uword in_col2) const;
  
  inline       subview_row<eT> subvec(const uword start_col, const SizeMat& s);
  inline const subview_row<eT> subvec(const uword start_col, const SizeMat& s) const;
  
  inline       subview_row<eT> head(const uword N);
  inline const subview_row<eT> head(const uword N) const;
  
  inline       subview_row<eT> tail(const uword N);
  inline const subview_row<eT> tail(const uword N) const;
  
  inline arma_warn_unused uword index_min() const;
  inline arma_warn_unused uword index_max() const;
  
  inline typename subview<eT>::row_iterator        begin();
  inline typename subview<eT>::const_row_iterator  begin() const;
  inline typename subview<eT>::const_row_iterator cbegin() const;
  
  inline typename subview<eT>::row_iterator        end();
  inline typename subview<eT>::const_row_iterator  end() const;
  inline typename subview<eT>::const_row_iterator cend() const;
  
  inline  subview_row(const subview_row&  in);
  inline  subview_row(      subview_row&& in);
  
  
  protected:
  
  inline subview_row(const Mat<eT>& in_m, const uword in_row);
  inline subview_row(const Mat<eT>& in_m, const uword in_row, const uword in_col1, const uword in_n_cols);
  inline subview_row() = delete;
  
  
  friend class Mat<eT>;
  friend class Row<eT>;
  friend class subview<eT>;
  };



template<typename eT>
class subview_row_strans : public Base< eT, subview_row_strans<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  arma_aligned const subview_row<eT>& sv_row;
  
  const uword n_rows;     // equal to n_elem
  const uword n_elem;
  
  static constexpr uword n_cols = 1;
  
  
  inline explicit subview_row_strans(const subview_row<eT>& in_sv_row);
  
  inline void extract(Mat<eT>& out) const;
  
  inline eT  at_alt    (const uword i) const;
  
  inline eT  operator[](const uword i) const;
  inline eT  operator()(const uword i) const;
  
  inline eT  operator()(const uword in_row, const uword in_col) const;
  inline eT          at(const uword in_row, const uword in_col) const;
  };



template<typename eT>
class subview_row_htrans : public Base< eT, subview_row_htrans<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  arma_aligned const subview_row<eT>& sv_row;
  
  const uword n_rows;     // equal to n_elem
  const uword n_elem;
  
  static constexpr uword n_cols = 1;
  
  
  inline explicit subview_row_htrans(const subview_row<eT>& in_sv_row);
  
  inline void extract(Mat<eT>& out) const;
  
  inline eT  at_alt    (const uword i) const;
  
  inline eT  operator[](const uword i) const;
  inline eT  operator()(const uword i) const;
  
  inline eT  operator()(const uword in_row, const uword in_col) const;
  inline eT          at(const uword in_row, const uword in_col) const;
  };



//! @}
