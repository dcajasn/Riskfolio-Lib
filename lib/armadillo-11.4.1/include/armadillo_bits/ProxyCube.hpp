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


//! \addtogroup ProxyCube
//! @{



template<typename T1>
struct ProxyCube
  {
  inline ProxyCube(const T1&)
    {
    arma_type_check(( is_arma_cube_type<T1>::value == false ));
    }
  };



// ea_type is the "element accessor" type,
// which can provide access to elements via operator[]

template<typename eT>
struct ProxyCube< Cube<eT> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Cube<eT>                                 stored_type;
  typedef const eT*                                ea_type;
  typedef const Cube<eT>&                          aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  arma_aligned const Cube<eT>& Q;
  
  inline explicit ProxyCube(const Cube<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows()       const { return Q.n_rows;       }
  arma_inline uword get_n_cols()       const { return Q.n_cols;       }
  arma_inline uword get_n_elem_slice() const { return Q.n_elem_slice; }
  arma_inline uword get_n_slices()     const { return Q.n_slices;     }
  arma_inline uword get_n_elem()       const { return Q.n_elem;       }
  
  arma_inline elem_type operator[] (const uword i)                               const { return Q[i];          }
  arma_inline elem_type at         (const uword r, const uword c, const uword s) const { return Q.at(r, c, s); }
  arma_inline elem_type at_alt     (const uword i)                               const { return Q.at_alt(i);   }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  arma_inline bool is_alias(const Cube<eT2>& X) const { return (void_ptr(&Q) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview_cube<eT2>& X) const { return is_alias(X.m); }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename eT, typename gen_type>
struct ProxyCube< GenCube<eT, gen_type> >
  {
  typedef          eT                              elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef GenCube<eT, gen_type>                    stored_type;
  typedef const GenCube<eT, gen_type>&             ea_type;
  typedef const GenCube<eT, gen_type>&             aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  arma_aligned const GenCube<eT, gen_type>& Q;
  
  inline explicit ProxyCube(const GenCube<eT, gen_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows()       const { return Q.n_rows;                     }
  arma_inline uword get_n_cols()       const { return Q.n_cols;                     }
  arma_inline uword get_n_elem_slice() const { return Q.n_rows*Q.n_cols;            }
  arma_inline uword get_n_slices()     const { return Q.n_slices;                   }
  arma_inline uword get_n_elem()       const { return Q.n_rows*Q.n_cols*Q.n_slices; }
  
  arma_inline elem_type operator[] (const uword i)                               const { return Q[i];          }
  arma_inline elem_type at         (const uword r, const uword c, const uword s) const { return Q.at(r, c, s); }
  arma_inline elem_type at_alt     (const uword i)                               const { return Q[i];          }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  constexpr bool is_alias(const Cube<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview_cube<eT2>&) const { return false; }
  
  constexpr bool is_aligned() const { return GenCube<eT, gen_type>::is_simple; }
  };



template<typename T1, typename op_type>
struct ProxyCube< OpCube<T1, op_type> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Cube<elem_type>                          stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Cube<elem_type>&                   aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  arma_aligned const Cube<elem_type> Q;
  
  inline explicit ProxyCube(const OpCube<T1, op_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows()       const { return Q.n_rows;       }
  arma_inline uword get_n_cols()       const { return Q.n_cols;       }
  arma_inline uword get_n_elem_slice() const { return Q.n_elem_slice; }
  arma_inline uword get_n_slices()     const { return Q.n_slices;     }
  arma_inline uword get_n_elem()       const { return Q.n_elem;       }
  
  arma_inline elem_type operator[] (const uword i)                               const { return Q[i];          }
  arma_inline elem_type at         (const uword r, const uword c, const uword s) const { return Q.at(r, c, s); }
  arma_inline elem_type at_alt     (const uword i)                               const { return Q.at_alt(i);   }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Cube<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview_cube<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1, typename T2, typename glue_type>
struct ProxyCube< GlueCube<T1, T2, glue_type> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Cube<elem_type>                          stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Cube<elem_type>&                   aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  arma_aligned const Cube<elem_type> Q;
  
  inline explicit ProxyCube(const GlueCube<T1, T2, glue_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }

  arma_inline uword get_n_rows()       const { return Q.n_rows;       }
  arma_inline uword get_n_cols()       const { return Q.n_cols;       }
  arma_inline uword get_n_elem_slice() const { return Q.n_elem_slice; }
  arma_inline uword get_n_slices()     const { return Q.n_slices;     }
  arma_inline uword get_n_elem()       const { return Q.n_elem;       }
  
  arma_inline elem_type operator[] (const uword i)                               const { return Q[i];          }
  arma_inline elem_type at         (const uword r, const uword c, const uword s) const { return Q.at(r, c, s); }
  arma_inline elem_type at_alt     (const uword i)                               const { return Q.at_alt(i);   }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Cube<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview_cube<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename eT>
struct ProxyCube< subview_cube<eT> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview_cube<eT>                         stored_type;
  typedef const subview_cube<eT>&                  ea_type;
  typedef const subview_cube<eT>&                  aligned_ea_type;
  
  static constexpr bool use_at      = true;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  arma_aligned const subview_cube<eT>& Q;
  
  inline explicit ProxyCube(const subview_cube<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows()       const { return Q.n_rows;       }
  arma_inline uword get_n_cols()       const { return Q.n_cols;       }
  arma_inline uword get_n_elem_slice() const { return Q.n_elem_slice; }
  arma_inline uword get_n_slices()     const { return Q.n_slices;     }
  arma_inline uword get_n_elem()       const { return Q.n_elem;       }
  
  arma_inline elem_type operator[] (const uword i)                               const { return Q[i];          }
  arma_inline elem_type at         (const uword r, const uword c, const uword s) const { return Q.at(r, c, s); }
  arma_inline elem_type at_alt     (const uword i)                               const { return Q.at_alt(i);   }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Cube<eT2>& X) const { return (void_ptr(&(Q.m)) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview_cube<eT2>& X) const { return Q.check_overlap(X); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename eT, typename T1>
struct ProxyCube< subview_cube_slices<eT,T1> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Cube<eT>                                 stored_type;
  typedef const eT*                                ea_type;
  typedef const Cube<eT>&                          aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  arma_aligned const Cube<eT> Q;
  
  inline explicit ProxyCube(const subview_cube_slices<eT,T1>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows()       const { return Q.n_rows;       }
  arma_inline uword get_n_cols()       const { return Q.n_cols;       }
  arma_inline uword get_n_elem_slice() const { return Q.n_elem_slice; }
  arma_inline uword get_n_slices()     const { return Q.n_slices;     }
  arma_inline uword get_n_elem()       const { return Q.n_elem;       }
  
  arma_inline elem_type operator[] (const uword i)                               const { return Q[i];          }
  arma_inline elem_type at         (const uword r, const uword c, const uword s) const { return Q.at(r, c, s); }
  arma_inline elem_type at_alt     (const uword i)                               const { return Q.at_alt(i);   }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Cube<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview_cube<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1, typename eop_type>
struct ProxyCube< eOpCube<T1, eop_type > >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eOpCube<T1, eop_type>                    stored_type;
  typedef const eOpCube<T1, eop_type>&             ea_type;
  typedef const eOpCube<T1, eop_type>&             aligned_ea_type;
  
  static constexpr bool use_at      = eOpCube<T1, eop_type>::use_at;
  static constexpr bool use_mp      = eOpCube<T1, eop_type>::use_mp;
  static constexpr bool has_subview = eOpCube<T1, eop_type>::has_subview;
  
  arma_aligned const eOpCube<T1, eop_type>& Q;
  
  inline explicit ProxyCube(const eOpCube<T1, eop_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows()       const { return Q.get_n_rows();       }
  arma_inline uword get_n_cols()       const { return Q.get_n_cols();       }
  arma_inline uword get_n_elem_slice() const { return Q.get_n_elem_slice(); }
  arma_inline uword get_n_slices()     const { return Q.get_n_slices();     }
  arma_inline uword get_n_elem()       const { return Q.get_n_elem();       }
  
  arma_inline elem_type operator[] (const uword i)                               const { return Q[i];          }
  arma_inline elem_type at         (const uword r, const uword c, const uword s) const { return Q.at(r, c, s); }
  arma_inline elem_type at_alt     (const uword i)                               const { return Q.at_alt(i);   }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Cube<eT2>& X) const { return Q.P.is_alias(X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview_cube<eT2>& X) const { return Q.P.has_overlap(X); }
  
  arma_inline bool is_aligned() const { return Q.P.is_aligned(); }
  };



template<typename T1, typename T2, typename eglue_type>
struct ProxyCube< eGlueCube<T1, T2, eglue_type > >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eGlueCube<T1, T2, eglue_type>            stored_type;
  typedef const eGlueCube<T1, T2, eglue_type>&     ea_type;
  typedef const eGlueCube<T1, T2, eglue_type>&     aligned_ea_type;
  
  static constexpr bool use_at      = eGlueCube<T1, T2, eglue_type>::use_at;
  static constexpr bool use_mp      = eGlueCube<T1, T2, eglue_type>::use_mp;
  static constexpr bool has_subview = eGlueCube<T1, T2, eglue_type>::has_subview;
  
  arma_aligned const eGlueCube<T1, T2, eglue_type>& Q;
  
  inline explicit ProxyCube(const eGlueCube<T1, T2, eglue_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows()       const { return Q.get_n_rows();       }
  arma_inline uword get_n_cols()       const { return Q.get_n_cols();       }
  arma_inline uword get_n_elem_slice() const { return Q.get_n_elem_slice(); }
  arma_inline uword get_n_slices()     const { return Q.get_n_slices();     }
  arma_inline uword get_n_elem()       const { return Q.get_n_elem();       }
  
  arma_inline elem_type operator[] (const uword i)                               const { return Q[i];          }
  arma_inline elem_type at         (const uword r, const uword c, const uword s) const { return Q.at(r, c, s); }
  arma_inline elem_type at_alt     (const uword i)                               const { return Q.at_alt(i);   }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Cube<eT2>& X) const { return (Q.P1.is_alias(X) || Q.P2.is_alias(X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview_cube<eT2>& X) const { return (Q.P1.has_overlap(X) || Q.P2.has_overlap(X)); }
  
  arma_inline bool is_aligned() const { return Q.P1.is_aligned() && Q.P2.is_aligned(); }
  };



template<typename out_eT, typename T1, typename op_type>
struct ProxyCube< mtOpCube<out_eT, T1, op_type> >
  {
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  typedef          Cube<out_eT>                 stored_type;
  typedef          const elem_type*             ea_type;
  typedef          const Cube<out_eT>&          aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  arma_aligned const Cube<out_eT> Q;
  
  inline explicit ProxyCube(const mtOpCube<out_eT, T1, op_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows()       const { return Q.n_rows;       }
  arma_inline uword get_n_cols()       const { return Q.n_cols;       }
  arma_inline uword get_n_elem_slice() const { return Q.n_elem_slice; }
  arma_inline uword get_n_slices()     const { return Q.n_slices;     }
  arma_inline uword get_n_elem()       const { return Q.n_elem;       }
  
  arma_inline elem_type operator[] (const uword i)                               const { return Q[i];          }
  arma_inline elem_type at         (const uword r, const uword c, const uword s) const { return Q.at(r, c, s); }
  arma_inline elem_type at_alt     (const uword i)                               const { return Q.at_alt(i);   }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Cube<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview_cube<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename out_eT, typename T1, typename T2, typename glue_type>
struct ProxyCube< mtGlueCube<out_eT, T1, T2, glue_type > >
  {
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  typedef          Cube<out_eT>                 stored_type;
  typedef          const elem_type*             ea_type;
  typedef          const Cube<out_eT>&          aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  arma_aligned const Cube<out_eT> Q;
  
  inline explicit ProxyCube(const mtGlueCube<out_eT, T1, T2, glue_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows()       const { return Q.n_rows;       }
  arma_inline uword get_n_cols()       const { return Q.n_cols;       }
  arma_inline uword get_n_elem_slice() const { return Q.n_elem_slice; }
  arma_inline uword get_n_slices()     const { return Q.n_slices;     }
  arma_inline uword get_n_elem()       const { return Q.n_elem;       }
  
  arma_inline elem_type operator[] (const uword i)                               const { return Q[i];          }
  arma_inline elem_type at         (const uword r, const uword c, const uword s) const { return Q.at(r, c, s); }
  arma_inline elem_type at_alt     (const uword i)                               const { return Q.at_alt(i);   }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Cube<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview_cube<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



//! @}
