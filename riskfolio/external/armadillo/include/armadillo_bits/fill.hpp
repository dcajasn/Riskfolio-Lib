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


//! \addtogroup fill
//! @{


namespace fill
  {
  struct fill_none  {};
  struct fill_zeros {};
  struct fill_ones  {};
  struct fill_eye   {};
  struct fill_randu {};
  struct fill_randn {};
  
  template<typename fill_type> 
  struct fill_class { inline fill_class() {} };
  
  static const fill_class<fill_none > none;
  static const fill_class<fill_zeros> zeros;
  static const fill_class<fill_ones > ones;
  static const fill_class<fill_eye  > eye;
  static const fill_class<fill_randu> randu;
  static const fill_class<fill_randn> randn;
  
  //
  
  template<typename from_type, typename to_type>
  struct allow_conversion
    {
    static constexpr bool value = true;
    };

  template<> struct allow_conversion<std::complex<double>, double> { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<double>, float > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<double>, u64   > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<double>, s64   > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<double>, u32   > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<double>, s32   > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<double>, u16   > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<double>, s16   > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<double>, u8    > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<double>, s8    > { static constexpr bool value = false; };

  template<> struct allow_conversion<std::complex<float>, double> { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<float>, float > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<float>, u64   > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<float>, s64   > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<float>, u32   > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<float>, s32   > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<float>, u16   > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<float>, s16   > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<float>, u8    > { static constexpr bool value = false; };
  template<> struct allow_conversion<std::complex<float>, s8    > { static constexpr bool value = false; };
  
  //
  
  template<typename eT> inline bool isfinite_wrapper(eT                )  { return true;                                               }
  template<>            inline bool isfinite_wrapper(float            x)  { return std::isfinite(x);                                   }
  template<>            inline bool isfinite_wrapper(double           x)  { return std::isfinite(x);                                   }
  template<typename  T> inline bool isfinite_wrapper(std::complex<T>& x)  { return std::isfinite(x.real()) && std::isfinite(x.imag()); }
  
  //
  
  template<typename scalar_type1>
  struct scalar_holder
    {
    const scalar_type1 scalar;
    
    inline explicit scalar_holder(const scalar_type1& in_scalar) : scalar(in_scalar) {}
    
    inline scalar_holder() = delete;
    
    template
      <
      typename scalar_type2,
      typename arma::enable_if2<allow_conversion<scalar_type1, scalar_type2>::value, int>::result = 0
      >
    inline
    operator scalar_holder<scalar_type2>() const
      {
      const bool ok_conversion = (std::is_integral<scalar_type2>::value && std::is_floating_point<scalar_type1>::value) ? isfinite_wrapper(scalar) : true;
      
      return scalar_holder<scalar_type2>( ok_conversion ? scalar_type2(scalar) : scalar_type2(0) );
      }
    };
  
  //
  
  template<typename scalar_type>
  inline
  typename enable_if2< is_supported_elem_type<scalar_type>::value, scalar_holder<scalar_type> >::result
  value(const scalar_type& in_scalar)
    {
    return scalar_holder<scalar_type>(in_scalar);
    }
  }


//! @}
