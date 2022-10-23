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


//! \addtogroup Cube
//! @{



struct Cube_prealloc
  {
  static constexpr uword mat_ptrs_size = 4;
  static constexpr uword mem_n_elem    = 64;
  };



//! Dense cube class

template<typename eT>
class Cube : public BaseCube< eT, Cube<eT> >
  {
  public:
  
  typedef eT                                elem_type; //!< the type of elements stored in the cube
  typedef typename get_pod_type<eT>::result  pod_type; //!< if eT is std::complex<T>, pod_type is T; otherwise pod_type is eT
  
  const uword n_rows;       //!< number of rows     in each slice (read-only)
  const uword n_cols;       //!< number of columns  in each slice (read-only)
  const uword n_elem_slice; //!< number of elements in each slice (read-only)
  const uword n_slices;     //!< number of slices   in the cube   (read-only)
  const uword n_elem;       //!< number of elements in the cube   (read-only)
  const uword n_alloc;      //!< number of allocated elements     (read-only); NOTE: n_alloc can be 0, even if n_elem > 0
  const uword mem_state;
  
  // mem_state = 0: normal cube which manages its own memory
  // mem_state = 1: use auxiliary memory until a size change
  // mem_state = 2: use auxiliary memory and don't allow the number of elements to be changed
  // mem_state = 3: fixed size (eg. via template based size specification)
  
  arma_aligned const eT* const mem;  //!< pointer to the memory used for storing elements (memory is read-only)
  
  
  protected:
  
  using mat_type = Mat<eT>;
  
  #if defined(ARMA_USE_OPENMP)
    using    raw_mat_ptr_type = mat_type*;
    using atomic_mat_ptr_type = mat_type*;
  #elif (!defined(ARMA_DONT_USE_STD_MUTEX))
    using    raw_mat_ptr_type = mat_type*;
    using atomic_mat_ptr_type = std::atomic<mat_type*>;
  #else
    using    raw_mat_ptr_type = mat_type*;
    using atomic_mat_ptr_type = mat_type*;
  #endif
  
  atomic_mat_ptr_type* mat_ptrs = nullptr;
  
  #if (!defined(ARMA_DONT_USE_STD_MUTEX))
    mutable std::mutex mat_mutex;   // required for slice()
  #endif
  
  arma_aligned   atomic_mat_ptr_type mat_ptrs_local[ Cube_prealloc::mat_ptrs_size ];
  arma_align_mem eT                       mem_local[ Cube_prealloc::mem_n_elem    ];  // local storage, for small cubes
  
  
  public:
  
  inline ~Cube();
  inline  Cube();
  
  inline explicit Cube(const uword in_n_rows, const uword in_n_cols, const uword in_n_slices);
  inline explicit Cube(const SizeCube& s);
  
  template<bool do_zeros> inline explicit Cube(const uword in_n_rows, const uword in_n_cols, const uword in_n_slices, const arma_initmode_indicator<do_zeros>&);
  template<bool do_zeros> inline explicit Cube(const SizeCube& s,                                                     const arma_initmode_indicator<do_zeros>&);
  
  template<typename fill_type> inline Cube(const uword in_n_rows, const uword in_n_cols, const uword in_n_slices, const fill::fill_class<fill_type>& f);
  template<typename fill_type> inline Cube(const SizeCube& s,                                                     const fill::fill_class<fill_type>& f);
  
  inline Cube(const uword in_rows, const uword in_cols, const uword in_slices, const fill::scalar_holder<eT> f);
  inline Cube(const SizeCube& s,                                               const fill::scalar_holder<eT> f);
  
  inline            Cube(Cube&& m);
  inline Cube& operator=(Cube&& m);
  
  inline Cube(      eT* aux_mem, const uword aux_n_rows, const uword aux_n_cols, const uword aux_n_slices, const bool copy_aux_mem = true, const bool strict = false, const bool prealloc_mat = false);
  inline Cube(const eT* aux_mem, const uword aux_n_rows, const uword aux_n_cols, const uword aux_n_slices);
  
  inline Cube& operator= (const eT val);
  inline Cube& operator+=(const eT val);
  inline Cube& operator-=(const eT val);
  inline Cube& operator*=(const eT val);
  inline Cube& operator/=(const eT val);
  
  inline             Cube(const Cube& m);
  inline Cube& operator= (const Cube& m);
  inline Cube& operator+=(const Cube& m);
  inline Cube& operator-=(const Cube& m);
  inline Cube& operator%=(const Cube& m);
  inline Cube& operator/=(const Cube& m);
  
  template<typename T1, typename T2>
  inline explicit Cube(const BaseCube<pod_type,T1>& A, const BaseCube<pod_type,T2>& B);
  
  inline             Cube(const subview_cube<eT>& X);
  inline Cube& operator= (const subview_cube<eT>& X);
  inline Cube& operator+=(const subview_cube<eT>& X);
  inline Cube& operator-=(const subview_cube<eT>& X);
  inline Cube& operator%=(const subview_cube<eT>& X);
  inline Cube& operator/=(const subview_cube<eT>& X);

  template<typename T1> inline             Cube(const subview_cube_slices<eT,T1>& X);
  template<typename T1> inline Cube& operator= (const subview_cube_slices<eT,T1>& X);
  template<typename T1> inline Cube& operator+=(const subview_cube_slices<eT,T1>& X);
  template<typename T1> inline Cube& operator-=(const subview_cube_slices<eT,T1>& X);
  template<typename T1> inline Cube& operator%=(const subview_cube_slices<eT,T1>& X);
  template<typename T1> inline Cube& operator/=(const subview_cube_slices<eT,T1>& X);

  arma_inline       subview_cube<eT> row(const uword in_row);
  arma_inline const subview_cube<eT> row(const uword in_row) const;
  
  arma_inline       subview_cube<eT> col(const uword in_col);
  arma_inline const subview_cube<eT> col(const uword in_col) const;
  
  inline       Mat<eT>& slice(const uword in_slice);
  inline const Mat<eT>& slice(const uword in_slice) const;
  
  arma_inline       subview_cube<eT> rows(const uword in_row1, const uword in_row2);
  arma_inline const subview_cube<eT> rows(const uword in_row1, const uword in_row2) const;
  
  arma_inline       subview_cube<eT> cols(const uword in_col1, const uword in_col2);
  arma_inline const subview_cube<eT> cols(const uword in_col1, const uword in_col2) const;

  arma_inline       subview_cube<eT> slices(const uword in_slice1, const uword in_slice2);
  arma_inline const subview_cube<eT> slices(const uword in_slice1, const uword in_slice2) const;
  
  arma_inline       subview_cube<eT> subcube(const uword in_row1, const uword in_col1, const uword in_slice1, const uword in_row2, const uword in_col2, const uword in_slice2);
  arma_inline const subview_cube<eT> subcube(const uword in_row1, const uword in_col1, const uword in_slice1, const uword in_row2, const uword in_col2, const uword in_slice2) const;
  
  inline            subview_cube<eT> subcube(const uword in_row1, const uword in_col1, const uword in_slice1, const SizeCube& s);
  inline      const subview_cube<eT> subcube(const uword in_row1, const uword in_col1, const uword in_slice1, const SizeCube& s) const;
  
  inline            subview_cube<eT> subcube(const span& row_span, const span& col_span, const span& slice_span);
  inline      const subview_cube<eT> subcube(const span& row_span, const span& col_span, const span& slice_span) const;
  
  inline            subview_cube<eT> operator()(const span& row_span, const span& col_span, const span& slice_span);
  inline      const subview_cube<eT> operator()(const span& row_span, const span& col_span, const span& slice_span) const;
  
  inline            subview_cube<eT> operator()(const uword in_row1, const uword in_col1, const uword in_slice1, const SizeCube& s);
  inline      const subview_cube<eT> operator()(const uword in_row1, const uword in_col1, const uword in_slice1, const SizeCube& s) const;
  
  arma_inline       subview_cube<eT> tube(const uword in_row1, const uword in_col1);
  arma_inline const subview_cube<eT> tube(const uword in_row1, const uword in_col1) const;
  
  arma_inline       subview_cube<eT> tube(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2);
  arma_inline const subview_cube<eT> tube(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2) const;
  
  arma_inline       subview_cube<eT> tube(const uword in_row1, const uword in_col1, const SizeMat& s);
  arma_inline const subview_cube<eT> tube(const uword in_row1, const uword in_col1, const SizeMat& s) const;
  
  inline            subview_cube<eT> tube(const span& row_span, const span& col_span);
  inline      const subview_cube<eT> tube(const span& row_span, const span& col_span) const;
  
  inline            subview_cube<eT> head_slices(const uword N);
  inline      const subview_cube<eT> head_slices(const uword N) const;
  
  inline            subview_cube<eT> tail_slices(const uword N);
  inline      const subview_cube<eT> tail_slices(const uword N) const;
  
  template<typename T1> arma_inline       subview_elem1<eT,T1> elem(const Base<uword,T1>& a);
  template<typename T1> arma_inline const subview_elem1<eT,T1> elem(const Base<uword,T1>& a) const;
  
  template<typename T1> arma_inline       subview_elem1<eT,T1> operator()(const Base<uword,T1>& a);
  template<typename T1> arma_inline const subview_elem1<eT,T1> operator()(const Base<uword,T1>& a) const;
  
  
  arma_inline       subview_cube_each1<eT> each_slice();
  arma_inline const subview_cube_each1<eT> each_slice() const;
  
  template<typename T1> inline       subview_cube_each2<eT, T1> each_slice(const Base<uword, T1>& indices);
  template<typename T1> inline const subview_cube_each2<eT, T1> each_slice(const Base<uword, T1>& indices) const;
  
  inline const Cube& each_slice(const std::function< void(      Mat<eT>&) >& F);
  inline const Cube& each_slice(const std::function< void(const Mat<eT>&) >& F) const;
  
  inline const Cube& each_slice(const std::function< void(      Mat<eT>&) >& F, const bool use_mp);
  inline const Cube& each_slice(const std::function< void(const Mat<eT>&) >& F, const bool use_mp) const;
  
  
  template<typename T1> arma_inline       subview_cube_slices<eT,T1> slices(const Base<uword,T1>& indices);
  template<typename T1> arma_inline const subview_cube_slices<eT,T1> slices(const Base<uword,T1>& indices) const;
  
  
  inline void shed_row(const uword row_num);
  inline void shed_col(const uword col_num);
  inline void shed_slice(const uword slice_num);
  
  inline void shed_rows(const uword in_row1, const uword in_row2);
  inline void shed_cols(const uword in_col1, const uword in_col2);
  inline void shed_slices(const uword in_slice1, const uword in_slice2);
  
  template<typename T1> inline void shed_slices(const Base<uword, T1>& indices);
  
  arma_deprecated inline void insert_rows(const uword row_num, const uword N, const bool set_to_zero);
  arma_deprecated inline void insert_cols(const uword row_num, const uword N, const bool set_to_zero);
  arma_deprecated inline void insert_slices(const uword slice_num, const uword N, const bool set_to_zero);
  
  inline void insert_rows(const uword row_num, const uword N);
  inline void insert_cols(const uword row_num, const uword N);
  inline void insert_slices(const uword slice_num, const uword N);
  
  template<typename T1> inline void insert_rows(const uword row_num, const BaseCube<eT,T1>& X);
  template<typename T1> inline void insert_cols(const uword col_num, const BaseCube<eT,T1>& X);
  template<typename T1> inline void insert_slices(const uword slice_num, const BaseCube<eT,T1>& X);
  template<typename T1> inline void insert_slices(const uword slice_num, const     Base<eT,T1>& X);
  
  
  template<typename gen_type> inline             Cube(const GenCube<eT, gen_type>& X);
  template<typename gen_type> inline Cube& operator= (const GenCube<eT, gen_type>& X);
  template<typename gen_type> inline Cube& operator+=(const GenCube<eT, gen_type>& X);
  template<typename gen_type> inline Cube& operator-=(const GenCube<eT, gen_type>& X);
  template<typename gen_type> inline Cube& operator%=(const GenCube<eT, gen_type>& X);
  template<typename gen_type> inline Cube& operator/=(const GenCube<eT, gen_type>& X);
  
  template<typename T1, typename op_type> inline             Cube(const OpCube<T1, op_type>& X);
  template<typename T1, typename op_type> inline Cube& operator= (const OpCube<T1, op_type>& X);
  template<typename T1, typename op_type> inline Cube& operator+=(const OpCube<T1, op_type>& X);
  template<typename T1, typename op_type> inline Cube& operator-=(const OpCube<T1, op_type>& X);
  template<typename T1, typename op_type> inline Cube& operator%=(const OpCube<T1, op_type>& X);
  template<typename T1, typename op_type> inline Cube& operator/=(const OpCube<T1, op_type>& X);
  
  template<typename T1, typename eop_type> inline             Cube(const eOpCube<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline Cube& operator= (const eOpCube<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline Cube& operator+=(const eOpCube<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline Cube& operator-=(const eOpCube<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline Cube& operator%=(const eOpCube<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline Cube& operator/=(const eOpCube<T1, eop_type>& X);
  
  template<typename T1, typename op_type> inline             Cube(const mtOpCube<eT, T1, op_type>& X);
  template<typename T1, typename op_type> inline Cube& operator= (const mtOpCube<eT, T1, op_type>& X);
  template<typename T1, typename op_type> inline Cube& operator+=(const mtOpCube<eT, T1, op_type>& X);
  template<typename T1, typename op_type> inline Cube& operator-=(const mtOpCube<eT, T1, op_type>& X);
  template<typename T1, typename op_type> inline Cube& operator%=(const mtOpCube<eT, T1, op_type>& X);
  template<typename T1, typename op_type> inline Cube& operator/=(const mtOpCube<eT, T1, op_type>& X);
  
  template<typename T1, typename T2, typename glue_type> inline             Cube(const GlueCube<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline Cube& operator= (const GlueCube<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline Cube& operator+=(const GlueCube<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline Cube& operator-=(const GlueCube<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline Cube& operator%=(const GlueCube<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline Cube& operator/=(const GlueCube<T1, T2, glue_type>& X);
  
  template<typename T1, typename T2, typename eglue_type> inline             Cube(const eGlueCube<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline Cube& operator= (const eGlueCube<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline Cube& operator+=(const eGlueCube<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline Cube& operator-=(const eGlueCube<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline Cube& operator%=(const eGlueCube<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline Cube& operator/=(const eGlueCube<T1, T2, eglue_type>& X);
  
  template<typename T1, typename T2, typename glue_type> inline             Cube(const mtGlueCube<eT, T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline Cube& operator= (const mtGlueCube<eT, T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline Cube& operator+=(const mtGlueCube<eT, T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline Cube& operator-=(const mtGlueCube<eT, T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline Cube& operator%=(const mtGlueCube<eT, T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline Cube& operator/=(const mtGlueCube<eT, T1, T2, glue_type>& X);
  
  
  arma_inline arma_warn_unused const eT& at_alt     (const uword i) const;
  
  arma_inline arma_warn_unused       eT& operator[] (const uword i);
  arma_inline arma_warn_unused const eT& operator[] (const uword i) const;
  
  arma_inline arma_warn_unused       eT& at(const uword i);
  arma_inline arma_warn_unused const eT& at(const uword i) const;
  
  arma_inline arma_warn_unused       eT& operator() (const uword i);
  arma_inline arma_warn_unused const eT& operator() (const uword i) const;
  
  #if defined(__cpp_multidimensional_subscript)
  arma_inline arma_warn_unused       eT& operator[] (const uword in_row, const uword in_col, const uword in_slice);
  arma_inline arma_warn_unused const eT& operator[] (const uword in_row, const uword in_col, const uword in_slice) const;
  #endif
  
  arma_inline arma_warn_unused       eT& at         (const uword in_row, const uword in_col, const uword in_slice);
  arma_inline arma_warn_unused const eT& at         (const uword in_row, const uword in_col, const uword in_slice) const;
  
  arma_inline arma_warn_unused       eT& operator() (const uword in_row, const uword in_col, const uword in_slice);
  arma_inline arma_warn_unused const eT& operator() (const uword in_row, const uword in_col, const uword in_slice) const;
  
  arma_inline const Cube& operator++();
  arma_inline void        operator++(int);
  
  arma_inline const Cube& operator--();
  arma_inline void        operator--(int);
  
       inline arma_warn_unused bool is_finite() const;
  arma_inline arma_warn_unused bool is_empty()  const;
  
  inline arma_warn_unused bool has_inf()       const;
  inline arma_warn_unused bool has_nan()       const;
  inline arma_warn_unused bool has_nonfinite() const;
  
  arma_inline arma_warn_unused bool in_range(const uword i) const;
  arma_inline arma_warn_unused bool in_range(const span& x) const;
  
  arma_inline arma_warn_unused bool in_range(const uword   in_row, const uword   in_col, const uword   in_slice) const;
       inline arma_warn_unused bool in_range(const span& row_span, const span& col_span, const span& slice_span) const;
  
       inline arma_warn_unused bool in_range(const uword   in_row, const uword   in_col, const uword   in_slice, const SizeCube& s) const;
  
  arma_inline arma_warn_unused       eT* memptr();
  arma_inline arma_warn_unused const eT* memptr() const;
  
  arma_inline arma_warn_unused       eT* slice_memptr(const uword slice);
  arma_inline arma_warn_unused const eT* slice_memptr(const uword slice) const;
  
  arma_inline arma_warn_unused       eT* slice_colptr(const uword in_slice, const uword in_col);
  arma_inline arma_warn_unused const eT* slice_colptr(const uword in_slice, const uword in_col) const;
  
  inline void set_size(const uword new_n_rows, const uword new_n_cols, const uword new_n_slices);
  inline void set_size(const SizeCube& s);
  
  inline void reshape(const uword new_n_rows, const uword new_n_cols, const uword new_n_slices);
  inline void reshape(const SizeCube& s);
                  
  inline void resize(const uword new_n_rows, const uword new_n_cols, const uword new_n_slices);
  inline void resize(const SizeCube& s);
  
  
  template<typename eT2> inline void copy_size(const Cube<eT2>& m);
  
  template<typename functor> inline const Cube&  for_each(functor F);
  template<typename functor> inline const Cube&  for_each(functor F) const;
  
  template<typename functor> inline const Cube& transform(functor F);
  template<typename functor> inline const Cube&     imbue(functor F);
  
  inline const Cube& replace(const eT old_val, const eT new_val);
  
  inline const Cube& clean(const pod_type threshold);
  
  inline const Cube& clamp(const eT min_val, const eT max_val);
  
  inline const Cube& fill(const eT val);
  
  inline const Cube& zeros();
  inline const Cube& zeros(const uword new_n_rows, const uword new_n_cols, const uword new_n_slices);
  inline const Cube& zeros(const SizeCube& s);
  
  inline const Cube& ones();
  inline const Cube& ones(const uword new_n_rows, const uword new_n_cols, const uword new_n_slices);
  inline const Cube& ones(const SizeCube& s);
  
  inline const Cube& randu();
  inline const Cube& randu(const uword new_n_rows, const uword new_n_cols, const uword new_n_slices);
  inline const Cube& randu(const SizeCube& s);
  
  inline const Cube& randn();
  inline const Cube& randn(const uword new_n_rows, const uword new_n_cols, const uword new_n_slices);
  inline const Cube& randn(const SizeCube& s);
  
  inline void      reset();
  inline void soft_reset();
  
  
  template<typename T1> inline void set_real(const BaseCube<pod_type,T1>& X);
  template<typename T1> inline void set_imag(const BaseCube<pod_type,T1>& X);
  
  
  inline arma_warn_unused eT min() const;
  inline arma_warn_unused eT max() const;
  
  inline eT min(uword& index_of_min_val) const;
  inline eT max(uword& index_of_max_val) const;
  
  inline eT min(uword& row_of_min_val, uword& col_of_min_val, uword& slice_of_min_val) const;
  inline eT max(uword& row_of_max_val, uword& col_of_max_val, uword& slice_of_max_val) const;
  
  
  inline arma_cold bool save(const std::string   name, const file_type type = arma_binary) const;
  inline arma_cold bool save(const hdf5_name&    spec, const file_type type = hdf5_binary) const;
  inline arma_cold bool save(      std::ostream& os,   const file_type type = arma_binary) const;
  
  inline arma_cold bool load(const std::string   name, const file_type type = auto_detect);
  inline arma_cold bool load(const hdf5_name&    spec, const file_type type = hdf5_binary);
  inline arma_cold bool load(      std::istream& is,   const file_type type = auto_detect);
  
  inline arma_cold bool quiet_save(const std::string   name, const file_type type = arma_binary) const;
  inline arma_cold bool quiet_save(const hdf5_name&    spec, const file_type type = hdf5_binary) const;
  inline arma_cold bool quiet_save(      std::ostream& os,   const file_type type = arma_binary) const;
  
  inline arma_cold bool quiet_load(const std::string   name, const file_type type = auto_detect);
  inline arma_cold bool quiet_load(const hdf5_name&    spec, const file_type type = hdf5_binary);
  inline arma_cold bool quiet_load(      std::istream& is,   const file_type type = auto_detect);
  
  
  // iterators
  
  typedef       eT*       iterator;
  typedef const eT* const_iterator;
  
  typedef       eT*       slice_iterator;
  typedef const eT* const_slice_iterator;
  
  inline       iterator  begin();
  inline const_iterator  begin() const;
  inline const_iterator cbegin() const;
  
  inline       iterator  end();
  inline const_iterator  end() const;
  inline const_iterator cend() const;
  
  inline       slice_iterator begin_slice(const uword slice_num);
  inline const_slice_iterator begin_slice(const uword slice_num) const;
  
  inline       slice_iterator end_slice(const uword slice_num);
  inline const_slice_iterator end_slice(const uword slice_num)   const;
  
  inline void  clear();
  inline bool  empty() const;
  inline uword size()  const;
  
  inline arma_warn_unused       eT& front();
  inline arma_warn_unused const eT& front() const;
  
  inline arma_warn_unused       eT& back();
  inline arma_warn_unused const eT& back() const;
  
  inline void swap(Cube& B);
  
  inline void steal_mem(Cube& X);                      //!< don't use this unless you're writing code internal to Armadillo
  inline void steal_mem(Cube& X, const bool is_move);  //!< don't use this unless you're writing code internal to Armadillo
  
  template<uword fixed_n_rows, uword fixed_n_cols, uword fixed_n_slices> class fixed;
  
  
  protected:
  
  inline void init_cold();
  inline void init_warm(const uword in_n_rows, const uword in_n_cols, const uword in_n_slices);
  
  template<typename T1, typename T2>
  inline void init(const BaseCube<pod_type,T1>& A, const BaseCube<pod_type,T2>& B);
  
  inline void delete_mat();
  inline void create_mat();
  
  inline void     create_mat_ptr(const uword in_slice) const;
  inline Mat<eT>*    get_mat_ptr(const uword in_slice) const;
  
  friend class glue_join;
  friend class op_reshape;
  friend class op_resize;
  friend class subview_cube<eT>;
  
  
  public:
  
  #if defined(ARMA_EXTRA_CUBE_PROTO)
    #include ARMA_INCFILE_WRAP(ARMA_EXTRA_CUBE_PROTO)
  #endif
  };



template<typename eT>
template<uword fixed_n_rows, uword fixed_n_cols, uword fixed_n_slices>
class Cube<eT>::fixed : public Cube<eT>
  {
  private:
  
  static constexpr uword fixed_n_elem       = fixed_n_rows * fixed_n_cols * fixed_n_slices;
  static constexpr uword fixed_n_elem_slice = fixed_n_rows * fixed_n_cols;
  
  static constexpr bool use_extra = (fixed_n_elem > Cube_prealloc::mem_n_elem);
  
  arma_aligned   atomic_mat_ptr_type mat_ptrs_local_extra[ (fixed_n_slices > Cube_prealloc::mat_ptrs_size) ? fixed_n_slices : 1 ];
  arma_align_mem eT                       mem_local_extra[ use_extra                                       ? fixed_n_elem   : 1 ];
  
  arma_inline void mem_setup();
  
  
  public:
  
  inline fixed();
  inline fixed(const fixed<fixed_n_rows, fixed_n_cols, fixed_n_slices>& X);
  
                                     inline fixed(const fill::scalar_holder<eT> f);
  template<typename fill_type>       inline fixed(const fill::fill_class<fill_type>& f);
  template<typename T1>              inline fixed(const BaseCube<eT,T1>& A);
  template<typename T1, typename T2> inline fixed(const BaseCube<pod_type,T1>& A, const BaseCube<pod_type,T2>& B);
  
  using Cube<eT>::operator=;
  using Cube<eT>::operator();
  
  inline Cube& operator=(const fixed<fixed_n_rows, fixed_n_cols, fixed_n_slices>& X);
  
  
  arma_inline arma_warn_unused       eT& operator[] (const uword i);
  arma_inline arma_warn_unused const eT& operator[] (const uword i) const;
  
  arma_inline arma_warn_unused       eT& at         (const uword i);
  arma_inline arma_warn_unused const eT& at         (const uword i) const;
  
  arma_inline arma_warn_unused       eT& operator() (const uword i);
  arma_inline arma_warn_unused const eT& operator() (const uword i) const;
  
  #if defined(__cpp_multidimensional_subscript)
  arma_inline arma_warn_unused       eT& operator[] (const uword in_row, const uword in_col, const uword in_slice);
  arma_inline arma_warn_unused const eT& operator[] (const uword in_row, const uword in_col, const uword in_slice) const;
  #endif
  
  arma_inline arma_warn_unused       eT& at         (const uword in_row, const uword in_col, const uword in_slice);
  arma_inline arma_warn_unused const eT& at         (const uword in_row, const uword in_col, const uword in_slice) const;
  
  arma_inline arma_warn_unused       eT& operator() (const uword in_row, const uword in_col, const uword in_slice);
  arma_inline arma_warn_unused const eT& operator() (const uword in_row, const uword in_col, const uword in_slice) const;
  };



class Cube_aux
  {
  public:
  
  template<typename eT> arma_inline static void prefix_pp(Cube<eT>& x);
  template<typename T>  arma_inline static void prefix_pp(Cube< std::complex<T> >& x);
  
  template<typename eT> arma_inline static void postfix_pp(Cube<eT>& x);
  template<typename T>  arma_inline static void postfix_pp(Cube< std::complex<T> >& x);
  
  template<typename eT> arma_inline static void prefix_mm(Cube<eT>& x);
  template<typename T>  arma_inline static void prefix_mm(Cube< std::complex<T> >& x);
  
  template<typename eT> arma_inline static void postfix_mm(Cube<eT>& x);
  template<typename T>  arma_inline static void postfix_mm(Cube< std::complex<T> >& x);
  
  template<typename eT, typename T1> inline static void set_real(Cube<eT>&                out, const BaseCube<eT,T1>& X);
  template<typename eT, typename T1> inline static void set_imag(Cube<eT>&                out, const BaseCube<eT,T1>& X);
  
  template<typename T,  typename T1> inline static void set_real(Cube< std::complex<T> >& out, const BaseCube< T,T1>& X);
  template<typename T,  typename T1> inline static void set_imag(Cube< std::complex<T> >& out, const BaseCube< T,T1>& X);
  };



//! @}
