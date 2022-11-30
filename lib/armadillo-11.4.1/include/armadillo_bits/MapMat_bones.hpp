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


//! \addtogroup MapMat
//! @{



// this class is for internal use only; subject to change and/or removal without notice
template<typename eT>
class MapMat
  {
  public:
  
  typedef eT                                elem_type;  //!< the type of elements stored in the matrix
  typedef typename get_pod_type<eT>::result  pod_type;  //!< if eT is std::complex<T>, pod_type is T; otherwise pod_type is eT
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const uword n_rows;    //!< number of rows     (read-only)
  const uword n_cols;    //!< number of columns  (read-only)
  const uword n_elem;    //!< number of elements (read-only)
  
  
  private:
  
  typedef typename std::map<uword, eT> map_type;
  
  arma_aligned map_type* map_ptr;
  
  
  public:
  
  inline ~MapMat();
  inline  MapMat();
  
  inline explicit MapMat(const uword in_n_rows, const uword in_n_cols);
  inline explicit MapMat(const SizeMat& s);
  
  inline          MapMat(const MapMat<eT>& x);
  inline void  operator=(const MapMat<eT>& x);
  
  inline explicit MapMat(const SpMat<eT>& x);
  inline void  operator=(const SpMat<eT>& x);
  
  inline          MapMat(MapMat<eT>&& x);
  inline void  operator=(MapMat<eT>&& x);
  
  inline void reset();
  inline void set_size(const uword in_n_rows);
  inline void set_size(const uword in_n_rows, const uword in_n_cols);
  inline void set_size(const SizeMat& s);
  
  inline void zeros();
  inline void zeros(const uword in_n_rows);
  inline void zeros(const uword in_n_rows, const uword in_n_cols);
  inline void zeros(const SizeMat& s);
  
  inline void eye();
  inline void eye(const uword in_n_rows, const uword in_n_cols);
  inline void eye(const SizeMat& s);
  
  inline void speye();
  inline void speye(const uword in_n_rows, const uword in_n_cols);
  inline void speye(const SizeMat& s);
  
  arma_inline arma_warn_unused MapMat_val<eT> operator[](const uword index);
       inline arma_warn_unused            eT  operator[](const uword index) const;
  
  arma_inline arma_warn_unused MapMat_val<eT> operator()(const uword index);
       inline arma_warn_unused            eT  operator()(const uword index) const;
  
  arma_inline arma_warn_unused MapMat_val<eT>         at(const uword in_row, const uword in_col);
       inline arma_warn_unused            eT          at(const uword in_row, const uword in_col) const;
  
  arma_inline arma_warn_unused MapMat_val<eT> operator()(const uword in_row, const uword in_col);
       inline arma_warn_unused            eT  operator()(const uword in_row, const uword in_col) const;
  
  inline arma_warn_unused bool is_empty()  const;
  inline arma_warn_unused bool is_vec()    const;
  inline arma_warn_unused bool is_rowvec() const;
  inline arma_warn_unused bool is_colvec() const;
  inline arma_warn_unused bool is_square() const;
  
  
  inline void sprandu(const uword in_n_rows, const uword in_n_cols, const double density);
  
  inline void print(const std::string& extra_text) const;
  
  inline uword get_n_nonzero() const;
  inline void  get_locval_format(umat& locs, Col<eT>& vals) const;
  
  
  private:
  
  inline void init_cold();
  inline void init_warm(const uword in_n_rows, const uword in_n_cols);
  
  arma_inline void   set_val(const uword index, const eT& in_val);
       inline void erase_val(const uword index);
  
  
  friend class                SpMat<eT>;
  friend class           MapMat_val<eT>;
  friend class     SpMat_MapMat_val<eT>;
  friend class SpSubview_MapMat_val<eT>;
  };



template<typename eT>
class MapMat_val
  {
  private:
  
  arma_aligned MapMat<eT>& parent;
  
  arma_aligned const uword index;
  
  inline MapMat_val(MapMat<eT>& in_parent, const uword in_index);
  
  friend class MapMat<eT>;
  
  
  public:
  
  arma_inline operator eT() const;
  
  arma_inline typename get_pod_type<eT>::result real() const;
  arma_inline typename get_pod_type<eT>::result imag() const;
  
  arma_inline void operator= (const MapMat_val<eT>& x);
  arma_inline void operator= (const eT in_val);
  arma_inline void operator+=(const eT in_val);
  arma_inline void operator-=(const eT in_val);
  arma_inline void operator*=(const eT in_val);
  arma_inline void operator/=(const eT in_val);
  
  arma_inline void operator++();
  arma_inline void operator++(int);
  
  arma_inline void operator--();
  arma_inline void operator--(int);
  };



template<typename eT>
class SpMat_MapMat_val
  {
  private:
  
  arma_aligned  SpMat<eT>& s_parent;
  arma_aligned MapMat<eT>& m_parent;
  
  arma_aligned const uword row;
  arma_aligned const uword col;
  
  inline SpMat_MapMat_val(SpMat<eT>& in_s_parent, MapMat<eT>& in_m_parent, const uword in_row, const uword in_col);
  
  friend class                SpMat<eT>;
  friend class               MapMat<eT>;
  friend class SpSubview_MapMat_val<eT>;
  
  
  public:
  
  inline operator eT() const;
  
  inline typename get_pod_type<eT>::result real() const;
  inline typename get_pod_type<eT>::result imag() const;
  
  inline SpMat_MapMat_val<eT>& operator= (const SpMat_MapMat_val<eT>& x);
  
  inline SpMat_MapMat_val<eT>& operator= (const eT in_val);
  inline SpMat_MapMat_val<eT>& operator+=(const eT in_val);
  inline SpMat_MapMat_val<eT>& operator-=(const eT in_val);
  inline SpMat_MapMat_val<eT>& operator*=(const eT in_val);
  inline SpMat_MapMat_val<eT>& operator/=(const eT in_val);
  
  inline                  SpMat_MapMat_val<eT>& operator++();
  inline arma_warn_unused eT                    operator++(int);
  
  inline                  SpMat_MapMat_val<eT>& operator--();
  inline arma_warn_unused eT                    operator--(int);
  
  inline void set(const eT in_val);
  inline void add(const eT in_val);
  inline void sub(const eT in_val);
  inline void mul(const eT in_val);
  inline void div(const eT in_val);
  };



template<typename eT>
class SpSubview_MapMat_val : public SpMat_MapMat_val<eT>
  {
  private:
  
  arma_inline SpSubview_MapMat_val(SpSubview<eT>& in_sv_parent, MapMat<eT>& in_m_parent, const uword in_row, const uword in_col);
  
  arma_aligned SpSubview<eT>& sv_parent;
  
  friend class            SpMat<eT>;
  friend class           MapMat<eT>;
  friend class        SpSubview<eT>;
  friend class SpMat_MapMat_val<eT>;
  
  
  public:
  
  inline SpSubview_MapMat_val<eT>& operator= (const SpSubview_MapMat_val<eT>& x);
  
  inline SpSubview_MapMat_val<eT>& operator= (const eT in_val);
  inline SpSubview_MapMat_val<eT>& operator+=(const eT in_val);
  inline SpSubview_MapMat_val<eT>& operator-=(const eT in_val);
  inline SpSubview_MapMat_val<eT>& operator*=(const eT in_val);
  inline SpSubview_MapMat_val<eT>& operator/=(const eT in_val);
  
  inline                  SpSubview_MapMat_val<eT>& operator++();
  inline arma_warn_unused eT                        operator++(int);
  
  inline                  SpSubview_MapMat_val<eT>& operator--();
  inline arma_warn_unused eT                        operator--(int);
  };



//! @}
