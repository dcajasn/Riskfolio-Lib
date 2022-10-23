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


//! \addtogroup op_inv_gen
//! @{



class op_inv_gen_default
  : public traits_op_default
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_gen_default>& in);
  
  template<typename T1>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, const char* caller_sig);
  };



class op_inv_gen_full
  : public traits_op_default
  {
  public:
  
  template<const uword row, const uword col>
  struct pos
    {
    static constexpr uword n2 = row + col*2;
    static constexpr uword n3 = row + col*3;
    static constexpr uword n4 = row + col*4;
    };
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_gen_full>& in);
  
  template<typename T1, const bool has_user_flags = true>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, const char* caller_sig, const uword flags);
  
  template<typename eT>
  arma_cold inline static bool apply_tiny_2x2(Mat<eT>& X);
  
  template<typename eT>
  arma_cold inline static bool apply_tiny_3x3(Mat<eT>& X);
  
  template<typename eT>
  arma_cold inline static bool apply_tiny_4x4(Mat<eT>& X);
  };



template<typename T>
struct op_inv_gen_state
  {
  T    rcond   = T(0);
  bool is_diag = false;
  bool is_sym  = false;
  };



class op_inv_gen_rcond
  : public traits_op_default
  {
  public:
  
  template<typename T1>
  inline static bool apply_direct(Mat<typename T1::elem_type>& out_inv, op_inv_gen_state<typename T1::pod_type>& out_state, const Base<typename T1::elem_type,T1>& expr);
  };



namespace inv_opts
  {
  struct opts
    {
    const uword flags;
    
    inline explicit opts(const uword in_flags);
    
    inline const opts operator+(const opts& rhs) const;
    };
  
  inline
  opts::opts(const uword in_flags)
    : flags(in_flags)
    {}
  
  inline
  const opts
  opts::operator+(const opts& rhs) const
    {
    const opts result( flags | rhs.flags );
    
    return result;
    }
  
  // The values below (eg. 1u << 1) are for internal Armadillo use only.
  // The values can change without notice.
  
  static constexpr uword flag_none         = uword(0       );
  static constexpr uword flag_tiny         = uword(1u <<  0);
  static constexpr uword flag_allow_approx = uword(1u <<  1);
  static constexpr uword flag_likely_sympd = uword(1u <<  2);
  static constexpr uword flag_no_sympd     = uword(1u <<  3);
  static constexpr uword flag_no_ugly      = uword(1u <<  4);
  
  struct opts_none         : public opts { inline opts_none()         : opts(flag_none        ) {} };
  struct opts_tiny         : public opts { inline opts_tiny()         : opts(flag_tiny        ) {} };
  struct opts_allow_approx : public opts { inline opts_allow_approx() : opts(flag_allow_approx) {} };
  struct opts_likely_sympd : public opts { inline opts_likely_sympd() : opts(flag_likely_sympd) {} };
  struct opts_no_sympd     : public opts { inline opts_no_sympd()     : opts(flag_no_sympd    ) {} };
  struct opts_no_ugly      : public opts { inline opts_no_ugly()      : opts(flag_no_ugly     ) {} };
  
  static const opts_none         none;
  static const opts_tiny         tiny;
  static const opts_allow_approx allow_approx;
  static const opts_likely_sympd likely_sympd;
  static const opts_no_sympd     no_sympd;
  static const opts_no_ugly      no_ugly;
  }



//! @}
