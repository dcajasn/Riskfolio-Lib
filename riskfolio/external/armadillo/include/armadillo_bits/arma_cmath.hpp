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



//! \addtogroup arma_cmath
//! @{



//
// wrappers for isfinite


template<typename eT>
inline
bool
arma_isfinite(eT)
  {
  return true;
  }



template<>
inline
bool
arma_isfinite(float x)
  {
  return std::isfinite(x);
  }



template<>
inline
bool
arma_isfinite(double x)
  {
  return std::isfinite(x);
  }



template<typename T>
inline
bool
arma_isfinite(const std::complex<T>& x)
  {
  return ( arma_isfinite(x.real()) && arma_isfinite(x.imag()) );
  }



//
// wrappers for isinf


template<typename eT>
inline
bool
arma_isinf(eT)
  {
  return false;
  }



template<>
inline
bool
arma_isinf(float x)
  {
  return std::isinf(x);
  }



template<>
inline
bool
arma_isinf(double x)
  {
  return std::isinf(x);
  }



template<typename T>
inline
bool
arma_isinf(const std::complex<T>& x)
  {
  return ( arma_isinf(x.real()) || arma_isinf(x.imag()) );
  }



//
// wrappers for isnan


template<typename eT>
inline
bool
arma_isnan(eT val)
  {
  arma_ignore(val);
    
  return false;
  }



template<>
inline
bool
arma_isnan(float x)
  {
  return std::isnan(x);
  }



template<>
inline
bool
arma_isnan(double x)
  {
  return std::isnan(x);
  }



template<typename T>
inline
bool
arma_isnan(const std::complex<T>& x)
  {
  return ( arma_isnan(x.real()) || arma_isnan(x.imag()) );
  }



//
// implementation of arma_sign()


template<typename eT>
constexpr
typename arma_unsigned_integral_only<eT>::result
arma_sign(const eT x)
  {
  return (x > eT(0)) ? eT(+1) : eT(0);
  }



template<typename eT>
constexpr
typename arma_signed_integral_only<eT>::result
arma_sign(const eT x)
  {
  return (x > eT(0)) ? eT(+1) : ( (x < eT(0)) ? eT(-1) : eT(0) );
  }



template<typename eT>
constexpr
typename arma_real_only<eT>::result
arma_sign(const eT x)
  {
  return (x > eT(0)) ? eT(+1) : ( (x < eT(0)) ? eT(-1) : ((x == eT(0)) ? eT(0) : x) );
  }



template<typename eT>
inline
typename arma_cx_only<eT>::result
arma_sign(const eT& x)
  {
  typedef typename eT::value_type T;
  
  const T abs_x = std::abs(x);
  
  return (abs_x != T(0)) ? (x / abs_x) : x;
  }



//
// wrappers for hypot(x, y) = sqrt(x^2 + y^2)


template<typename eT>
inline
eT
arma_hypot(const eT x, const eT y)
  {
  arma_ignore(x);
  arma_ignore(y);
  
  arma_stop_runtime_error("arma_hypot(): not implemented for integer or complex element types");
  
  return eT(0);
  }



template<>
inline
float
arma_hypot(const float x, const float y)
  {
  return std::hypot(x, y);
  }



template<>
inline
double
arma_hypot(const double x, const double y)
  {
  return std::hypot(x, y);
  }



//
// implementation of arma_sinc()


template<typename eT>
inline
eT
arma_sinc_generic(const eT x)
  {
  typedef typename get_pod_type<eT>::result T;
  
  const eT tmp = Datum<T>::pi * x;
  
  return (tmp == eT(0)) ? eT(1) : eT( std::sin(tmp) / tmp );
  }



template<typename eT>
inline
eT
arma_sinc(const eT x)
  {
  return eT( arma_sinc_generic( double(x) ) );
  }



template<>
inline
float
arma_sinc(const float x)
  {
  return arma_sinc_generic(x);
  }



template<>
inline
double
arma_sinc(const double x)
  {
  return arma_sinc_generic(x);
  }



template<typename T>
inline
std::complex<T>
arma_sinc(const std::complex<T>& x)
  {
  return arma_sinc_generic(x);
  }



//
// wrappers for arg()


template<typename eT>
struct arma_arg
  {
  static
  inline
  eT
  eval(const eT x)
    {
    return eT( std::arg(x) );
    }
  };



template<>
struct arma_arg<float>
  {
  static
  inline
  float
  eval(const float x)
    {
    return std::arg(x);
    }
  };



template<>
struct arma_arg<double>
  {
  static
  inline
  double
  eval(const double x)
    {
    return std::arg(x);
    }
  };



template<>
struct arma_arg< std::complex<float> >
  {
  static
  inline
  float
  eval(const std::complex<float>& x)
    {
    return std::arg(x);
    }
  };



template<>
struct arma_arg< std::complex<double> >
  {
  static
  inline
  double
  eval(const std::complex<double>& x)
    {
    return std::arg(x);
    }
  };



//! @}
