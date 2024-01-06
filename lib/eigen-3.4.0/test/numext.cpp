// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2017 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "main.h"

template<typename T, typename U>
bool check_if_equal_or_nans(const T& actual, const U& expected) {
  return ((actual == expected) || ((numext::isnan)(actual) && (numext::isnan)(expected)));
}

template<typename T, typename U>
bool check_if_equal_or_nans(const std::complex<T>& actual, const std::complex<U>& expected) {
  return check_if_equal_or_nans(numext::real(actual), numext::real(expected))
         && check_if_equal_or_nans(numext::imag(actual), numext::imag(expected));
}

template<typename T, typename U>
bool test_is_equal_or_nans(const T& actual, const U& expected)
{
    if (check_if_equal_or_nans(actual, expected)) {
      return true;
    }

    // false:
    std::cerr
        << "\n    actual   = " << actual
        << "\n    expected = " << expected << "\n\n";
    return false;
}

#define VERIFY_IS_EQUAL_OR_NANS(a, b) VERIFY(test_is_equal_or_nans(a, b))

template<typename T>
void check_abs() {
  typedef typename NumTraits<T>::Real Real;
  Real zero(0);

  if(NumTraits<T>::IsSigned)
    VERIFY_IS_EQUAL(numext::abs(-T(1)), T(1));
  VERIFY_IS_EQUAL(numext::abs(T(0)), T(0));
  VERIFY_IS_EQUAL(numext::abs(T(1)), T(1));

  for(int k=0; k<100; ++k)
  {
    T x = internal::random<T>();
    if(!internal::is_same<T,bool>::value)
      x = x/Real(2);
    if(NumTraits<T>::IsSigned)
    {
      VERIFY_IS_EQUAL(numext::abs(x), numext::abs(-x));
      VERIFY( numext::abs(-x) >= zero );
    }
    VERIFY( numext::abs(x) >= zero );
    VERIFY_IS_APPROX( numext::abs2(x), numext::abs2(numext::abs(x)) );
  }
}

template<typename T>
void check_arg() {
  typedef typename NumTraits<T>::Real Real;
  VERIFY_IS_EQUAL(numext::abs(T(0)), T(0));
  VERIFY_IS_EQUAL(numext::abs(T(1)), T(1));

  for(int k=0; k<100; ++k)
  {
    T x = internal::random<T>();
    Real y = numext::arg(x);
    VERIFY_IS_APPROX( y, std::arg(x) );
  }
}

template<typename T>
struct check_sqrt_impl {
  static void run() {
    for (int i=0; i<1000; ++i) {
      const T x = numext::abs(internal::random<T>());
      const T sqrtx = numext::sqrt(x);
      VERIFY_IS_APPROX(sqrtx*sqrtx, x);
    }

    // Corner cases.
    const T zero = T(0);
    const T one = T(1);
    const T inf = std::numeric_limits<T>::infinity();
    const T nan = std::numeric_limits<T>::quiet_NaN();
    VERIFY_IS_EQUAL(numext::sqrt(zero), zero);
    VERIFY_IS_EQUAL(numext::sqrt(inf), inf);
    VERIFY((numext::isnan)(numext::sqrt(nan)));
    VERIFY((numext::isnan)(numext::sqrt(-one)));
  }
};

template<typename T>
struct check_sqrt_impl<std::complex<T>  > {
  static void run() {
    typedef typename std::complex<T> ComplexT;

    for (int i=0; i<1000; ++i) {
      const ComplexT x = internal::random<ComplexT>();
      const ComplexT sqrtx = numext::sqrt(x);
      VERIFY_IS_APPROX(sqrtx*sqrtx, x);
    }

    // Corner cases.
    const T zero = T(0);
    const T one = T(1);
    const T inf = std::numeric_limits<T>::infinity();
    const T nan = std::numeric_limits<T>::quiet_NaN();

    // Set of corner cases from https://en.cppreference.com/w/cpp/numeric/complex/sqrt
    const int kNumCorners = 20;
    const ComplexT corners[kNumCorners][2] = {
      {ComplexT(zero, zero), ComplexT(zero, zero)},
      {ComplexT(-zero, zero), ComplexT(zero, zero)},
      {ComplexT(zero, -zero), ComplexT(zero, zero)},
      {ComplexT(-zero, -zero), ComplexT(zero, zero)},
      {ComplexT(one, inf), ComplexT(inf, inf)},
      {ComplexT(nan, inf), ComplexT(inf, inf)},
      {ComplexT(one, -inf), ComplexT(inf, -inf)},
      {ComplexT(nan, -inf), ComplexT(inf, -inf)},
      {ComplexT(-inf, one), ComplexT(zero, inf)},
      {ComplexT(inf, one), ComplexT(inf, zero)},
      {ComplexT(-inf, -one), ComplexT(zero, -inf)},
      {ComplexT(inf, -one), ComplexT(inf, -zero)},
      {ComplexT(-inf, nan), ComplexT(nan, inf)},
      {ComplexT(inf, nan), ComplexT(inf, nan)},
      {ComplexT(zero, nan), ComplexT(nan, nan)},
      {ComplexT(one, nan), ComplexT(nan, nan)},
      {ComplexT(nan, zero), ComplexT(nan, nan)},
      {ComplexT(nan, one), ComplexT(nan, nan)},
      {ComplexT(nan, -one), ComplexT(nan, nan)},
      {ComplexT(nan, nan), ComplexT(nan, nan)},
    };

    for (int i=0; i<kNumCorners; ++i) {
      const ComplexT& x = corners[i][0];
      const ComplexT sqrtx = corners[i][1];
      VERIFY_IS_EQUAL_OR_NANS(numext::sqrt(x), sqrtx);
    }
  }
};

template<typename T>
void check_sqrt() {
  check_sqrt_impl<T>::run();
}

template<typename T>
struct check_rsqrt_impl {
  static void run() {
    const T zero = T(0);
    const T one = T(1);
    const T inf = std::numeric_limits<T>::infinity();
    const T nan = std::numeric_limits<T>::quiet_NaN();

    for (int i=0; i<1000; ++i) {
      const T x = numext::abs(internal::random<T>());
      const T rsqrtx = numext::rsqrt(x);
      const T invx = one / x;
      VERIFY_IS_APPROX(rsqrtx*rsqrtx, invx);
    }

    // Corner cases.
    VERIFY_IS_EQUAL(numext::rsqrt(zero), inf);
    VERIFY_IS_EQUAL(numext::rsqrt(inf), zero);
    VERIFY((numext::isnan)(numext::rsqrt(nan)));
    VERIFY((numext::isnan)(numext::rsqrt(-one)));
  }
};

template<typename T>
struct check_rsqrt_impl<std::complex<T> > {
  static void run() {
    typedef typename std::complex<T> ComplexT;
    const T zero = T(0);
    const T one = T(1);
    const T inf = std::numeric_limits<T>::infinity();
    const T nan = std::numeric_limits<T>::quiet_NaN();

    for (int i=0; i<1000; ++i) {
      const ComplexT x = internal::random<ComplexT>();
      const ComplexT invx = ComplexT(one, zero) / x;
      const ComplexT rsqrtx = numext::rsqrt(x);
      VERIFY_IS_APPROX(rsqrtx*rsqrtx, invx);
    }

    // GCC and MSVC differ in their treatment of 1/(0 + 0i)
    //   GCC/clang = (inf, nan)
    //   MSVC = (nan, nan)
    // and 1 / (x + inf i)
    //   GCC/clang = (0, 0)
    //   MSVC = (nan, nan)
    #if (EIGEN_COMP_GNUC)
    {
      const int kNumCorners = 20;
      const ComplexT corners[kNumCorners][2] = {
        // Only consistent across GCC, clang
        {ComplexT(zero, zero), ComplexT(zero, zero)},
        {ComplexT(-zero, zero), ComplexT(zero, zero)},
        {ComplexT(zero, -zero), ComplexT(zero, zero)},
        {ComplexT(-zero, -zero), ComplexT(zero, zero)},
        {ComplexT(one, inf), ComplexT(inf, inf)},
        {ComplexT(nan, inf), ComplexT(inf, inf)},
        {ComplexT(one, -inf), ComplexT(inf, -inf)},
        {ComplexT(nan, -inf), ComplexT(inf, -inf)},
        // Consistent across GCC, clang, MSVC
        {ComplexT(-inf, one), ComplexT(zero, inf)},
        {ComplexT(inf, one), ComplexT(inf, zero)},
        {ComplexT(-inf, -one), ComplexT(zero, -inf)},
        {ComplexT(inf, -one), ComplexT(inf, -zero)},
        {ComplexT(-inf, nan), ComplexT(nan, inf)},
        {ComplexT(inf, nan), ComplexT(inf, nan)},
        {ComplexT(zero, nan), ComplexT(nan, nan)},
        {ComplexT(one, nan), ComplexT(nan, nan)},
        {ComplexT(nan, zero), ComplexT(nan, nan)},
        {ComplexT(nan, one), ComplexT(nan, nan)},
        {ComplexT(nan, -one), ComplexT(nan, nan)},
        {ComplexT(nan, nan), ComplexT(nan, nan)},
      };

      for (int i=0; i<kNumCorners; ++i) {
        const ComplexT& x = corners[i][0];
        const ComplexT rsqrtx = ComplexT(one, zero) / corners[i][1];
        VERIFY_IS_EQUAL_OR_NANS(numext::rsqrt(x), rsqrtx);
      }
    }
    #endif
  }
};

template<typename T>
void check_rsqrt() {
  check_rsqrt_impl<T>::run();
}

EIGEN_DECLARE_TEST(numext) {
  for(int k=0; k<g_repeat; ++k)
  {
    CALL_SUBTEST( check_abs<bool>() );
    CALL_SUBTEST( check_abs<signed char>() );
    CALL_SUBTEST( check_abs<unsigned char>() );
    CALL_SUBTEST( check_abs<short>() );
    CALL_SUBTEST( check_abs<unsigned short>() );
    CALL_SUBTEST( check_abs<int>() );
    CALL_SUBTEST( check_abs<unsigned int>() );
    CALL_SUBTEST( check_abs<long>() );
    CALL_SUBTEST( check_abs<unsigned long>() );
    CALL_SUBTEST( check_abs<half>() );
    CALL_SUBTEST( check_abs<bfloat16>() );
    CALL_SUBTEST( check_abs<float>() );
    CALL_SUBTEST( check_abs<double>() );
    CALL_SUBTEST( check_abs<long double>() );
    CALL_SUBTEST( check_abs<std::complex<float> >() );
    CALL_SUBTEST( check_abs<std::complex<double> >() );

    CALL_SUBTEST( check_arg<std::complex<float> >() );
    CALL_SUBTEST( check_arg<std::complex<double> >() );

    CALL_SUBTEST( check_sqrt<float>() );
    CALL_SUBTEST( check_sqrt<double>() );
    CALL_SUBTEST( check_sqrt<std::complex<float> >() );
    CALL_SUBTEST( check_sqrt<std::complex<double> >() );
    
    CALL_SUBTEST( check_rsqrt<float>() );
    CALL_SUBTEST( check_rsqrt<double>() );
    CALL_SUBTEST( check_rsqrt<std::complex<float> >() );
    CALL_SUBTEST( check_rsqrt<std::complex<double> >() );
  }
}
