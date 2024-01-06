// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <sstream>

#include "main.h"

#include <Eigen/src/Core/arch/Default/Half.h>

#define VERIFY_HALF_BITS_EQUAL(h, bits) \
  VERIFY_IS_EQUAL((numext::bit_cast<numext::uint16_t>(h)), (static_cast<numext::uint16_t>(bits)))

// Make sure it's possible to forward declare Eigen::half
namespace Eigen {
struct half;
}

using Eigen::half;

void test_conversion()
{
  using Eigen::half_impl::__half_raw;

  // Round-trip bit-cast with uint16.
  VERIFY_IS_EQUAL(
    numext::bit_cast<half>(numext::bit_cast<numext::uint16_t>(half(1.0f))),
    half(1.0f));
  VERIFY_IS_EQUAL(
    numext::bit_cast<half>(numext::bit_cast<numext::uint16_t>(half(0.5f))),
    half(0.5f));
  VERIFY_IS_EQUAL(
    numext::bit_cast<half>(numext::bit_cast<numext::uint16_t>(half(-0.33333f))),
    half(-0.33333f));
   VERIFY_IS_EQUAL(
    numext::bit_cast<half>(numext::bit_cast<numext::uint16_t>(half(0.0f))),
    half(0.0f));

  // Conversion from float.
  VERIFY_HALF_BITS_EQUAL(half(1.0f), 0x3c00);
  VERIFY_HALF_BITS_EQUAL(half(0.5f), 0x3800);
  VERIFY_HALF_BITS_EQUAL(half(0.33333f), 0x3555);
  VERIFY_HALF_BITS_EQUAL(half(0.0f), 0x0000);
  VERIFY_HALF_BITS_EQUAL(half(-0.0f), 0x8000);
  VERIFY_HALF_BITS_EQUAL(half(65504.0f), 0x7bff);
  VERIFY_HALF_BITS_EQUAL(half(65536.0f), 0x7c00);  // Becomes infinity.

  // Denormals.
  VERIFY_HALF_BITS_EQUAL(half(-5.96046e-08f), 0x8001);
  VERIFY_HALF_BITS_EQUAL(half(5.96046e-08f), 0x0001);
  VERIFY_HALF_BITS_EQUAL(half(1.19209e-07f), 0x0002);

  // Verify round-to-nearest-even behavior.
  float val1 = float(half(__half_raw(0x3c00)));
  float val2 = float(half(__half_raw(0x3c01)));
  float val3 = float(half(__half_raw(0x3c02)));
  VERIFY_HALF_BITS_EQUAL(half(0.5f * (val1 + val2)), 0x3c00);
  VERIFY_HALF_BITS_EQUAL(half(0.5f * (val2 + val3)), 0x3c02);

  // Conversion from int.
  VERIFY_HALF_BITS_EQUAL(half(-1), 0xbc00);
  VERIFY_HALF_BITS_EQUAL(half(0), 0x0000);
  VERIFY_HALF_BITS_EQUAL(half(1), 0x3c00);
  VERIFY_HALF_BITS_EQUAL(half(2), 0x4000);
  VERIFY_HALF_BITS_EQUAL(half(3), 0x4200);

  // Conversion from bool.
  VERIFY_HALF_BITS_EQUAL(half(false), 0x0000);
  VERIFY_HALF_BITS_EQUAL(half(true), 0x3c00);

  // Conversion to float.
  VERIFY_IS_EQUAL(float(half(__half_raw(0x0000))), 0.0f);
  VERIFY_IS_EQUAL(float(half(__half_raw(0x3c00))), 1.0f);

  // Denormals.
  VERIFY_IS_APPROX(float(half(__half_raw(0x8001))), -5.96046e-08f);
  VERIFY_IS_APPROX(float(half(__half_raw(0x0001))), 5.96046e-08f);
  VERIFY_IS_APPROX(float(half(__half_raw(0x0002))), 1.19209e-07f);

  // NaNs and infinities.
  VERIFY(!(numext::isinf)(float(half(65504.0f))));  // Largest finite number.
  VERIFY(!(numext::isnan)(float(half(0.0f))));
  VERIFY((numext::isinf)(float(half(__half_raw(0xfc00)))));
  VERIFY((numext::isnan)(float(half(__half_raw(0xfc01)))));
  VERIFY((numext::isinf)(float(half(__half_raw(0x7c00)))));
  VERIFY((numext::isnan)(float(half(__half_raw(0x7c01)))));

#if !EIGEN_COMP_MSVC
  // Visual Studio errors out on divisions by 0
  VERIFY((numext::isnan)(float(half(0.0 / 0.0))));
  VERIFY((numext::isinf)(float(half(1.0 / 0.0))));
  VERIFY((numext::isinf)(float(half(-1.0 / 0.0))));
#endif

  // Exactly same checks as above, just directly on the half representation.
  VERIFY(!(numext::isinf)(half(__half_raw(0x7bff))));
  VERIFY(!(numext::isnan)(half(__half_raw(0x0000))));
  VERIFY((numext::isinf)(half(__half_raw(0xfc00))));
  VERIFY((numext::isnan)(half(__half_raw(0xfc01))));
  VERIFY((numext::isinf)(half(__half_raw(0x7c00))));
  VERIFY((numext::isnan)(half(__half_raw(0x7c01))));

#if !EIGEN_COMP_MSVC
  // Visual Studio errors out on divisions by 0
  VERIFY((numext::isnan)(half(0.0 / 0.0)));
  VERIFY((numext::isinf)(half(1.0 / 0.0)));
  VERIFY((numext::isinf)(half(-1.0 / 0.0)));
#endif

  // Conversion to bool
  VERIFY(!static_cast<bool>(half(0.0)));
  VERIFY(!static_cast<bool>(half(-0.0)));
  VERIFY(static_cast<bool>(half(__half_raw(0x7bff))));
  VERIFY(static_cast<bool>(half(-0.33333)));
  VERIFY(static_cast<bool>(half(1.0)));
  VERIFY(static_cast<bool>(half(-1.0)));
  VERIFY(static_cast<bool>(half(-5.96046e-08f)));
}

void test_numtraits()
{
  std::cout << "epsilon       = " << NumTraits<half>::epsilon() << "  (0x" << std::hex << numext::bit_cast<numext::uint16_t>(NumTraits<half>::epsilon()) << ")" << std::endl;
  std::cout << "highest       = " << NumTraits<half>::highest() << "  (0x" << std::hex << numext::bit_cast<numext::uint16_t>(NumTraits<half>::highest()) << ")" << std::endl;
  std::cout << "lowest        = " << NumTraits<half>::lowest() << "  (0x" << std::hex << numext::bit_cast<numext::uint16_t>(NumTraits<half>::lowest()) << ")" << std::endl;
  std::cout << "min           = " << (std::numeric_limits<half>::min)() << "  (0x" << std::hex << numext::bit_cast<numext::uint16_t>(half((std::numeric_limits<half>::min)())) << ")" << std::endl;
  std::cout << "denorm min    = " << (std::numeric_limits<half>::denorm_min)() << "  (0x" << std::hex << numext::bit_cast<numext::uint16_t>(half((std::numeric_limits<half>::denorm_min)())) << ")" << std::endl;
  std::cout << "infinity      = " << NumTraits<half>::infinity() << "  (0x" << std::hex << numext::bit_cast<numext::uint16_t>(NumTraits<half>::infinity()) << ")" << std::endl;
  std::cout << "quiet nan     = " << NumTraits<half>::quiet_NaN() << "  (0x" << std::hex << numext::bit_cast<numext::uint16_t>(NumTraits<half>::quiet_NaN()) << ")" << std::endl;
  std::cout << "signaling nan = " << std::numeric_limits<half>::signaling_NaN() << "  (0x" << std::hex << numext::bit_cast<numext::uint16_t>(std::numeric_limits<half>::signaling_NaN()) << ")" << std::endl;

  VERIFY(NumTraits<half>::IsSigned);

  VERIFY_IS_EQUAL(
    numext::bit_cast<numext::uint16_t>(std::numeric_limits<half>::infinity()),
    numext::bit_cast<numext::uint16_t>(half(std::numeric_limits<float>::infinity())) );
  // There is no guarantee that casting a 32-bit NaN to 16-bit has a precise
  // bit pattern.  We test that it is in fact a NaN, then test the signaling
  // bit (msb of significand is 1 for quiet, 0 for signaling).
  const numext::uint16_t HALF_QUIET_BIT = 0x0200;
  VERIFY(
    (numext::isnan)(std::numeric_limits<half>::quiet_NaN())
    && (numext::isnan)(half(std::numeric_limits<float>::quiet_NaN()))
    && ((numext::bit_cast<numext::uint16_t>(std::numeric_limits<half>::quiet_NaN()) & HALF_QUIET_BIT) > 0)
    && ((numext::bit_cast<numext::uint16_t>(half(std::numeric_limits<float>::quiet_NaN())) & HALF_QUIET_BIT) > 0) );
  // After a cast to half, a signaling NaN may become non-signaling
  // (e.g. in the case of casting float to native __fp16). Thus, we check that
  // both are NaN, and that only the `numeric_limits` version is signaling.
  VERIFY(
    (numext::isnan)(std::numeric_limits<half>::signaling_NaN())
    && (numext::isnan)(half(std::numeric_limits<float>::signaling_NaN()))
    && ((numext::bit_cast<numext::uint16_t>(std::numeric_limits<half>::signaling_NaN()) & HALF_QUIET_BIT) == 0) );

  VERIFY( (std::numeric_limits<half>::min)() > half(0.f) );
  VERIFY( (std::numeric_limits<half>::denorm_min)() > half(0.f) );
  VERIFY( (std::numeric_limits<half>::min)()/half(2) > half(0.f) );
  VERIFY_IS_EQUAL( (std::numeric_limits<half>::denorm_min)()/half(2), half(0.f) );
}

void test_arithmetic()
{
  VERIFY_IS_EQUAL(float(half(2) + half(2)), 4);
  VERIFY_IS_EQUAL(float(half(2) + half(-2)), 0);
  VERIFY_IS_APPROX(float(half(0.33333f) + half(0.66667f)), 1.0f);
  VERIFY_IS_EQUAL(float(half(2.0f) * half(-5.5f)), -11.0f);
  VERIFY_IS_APPROX(float(half(1.0f) / half(3.0f)), 0.33333f);
  VERIFY_IS_EQUAL(float(-half(4096.0f)), -4096.0f);
  VERIFY_IS_EQUAL(float(-half(-4096.0f)), 4096.0f);
  
  half x(3);
  half y = ++x;
  VERIFY_IS_EQUAL(x, half(4));
  VERIFY_IS_EQUAL(y, half(4));
  y = --x;
  VERIFY_IS_EQUAL(x, half(3));
  VERIFY_IS_EQUAL(y, half(3));
  y = x++;
  VERIFY_IS_EQUAL(x, half(4));
  VERIFY_IS_EQUAL(y, half(3));
  y = x--;
  VERIFY_IS_EQUAL(x, half(3));
  VERIFY_IS_EQUAL(y, half(4));
}

void test_comparison()
{
  VERIFY(half(1.0f) > half(0.5f));
  VERIFY(half(0.5f) < half(1.0f));
  VERIFY(!(half(1.0f) < half(0.5f)));
  VERIFY(!(half(0.5f) > half(1.0f)));

  VERIFY(!(half(4.0f) > half(4.0f)));
  VERIFY(!(half(4.0f) < half(4.0f)));

  VERIFY(!(half(0.0f) < half(-0.0f)));
  VERIFY(!(half(-0.0f) < half(0.0f)));
  VERIFY(!(half(0.0f) > half(-0.0f)));
  VERIFY(!(half(-0.0f) > half(0.0f)));

  VERIFY(half(0.2f) > half(-1.0f));
  VERIFY(half(-1.0f) < half(0.2f));
  VERIFY(half(-16.0f) < half(-15.0f));

  VERIFY(half(1.0f) == half(1.0f));
  VERIFY(half(1.0f) != half(2.0f));

  // Comparisons with NaNs and infinities.
#if !EIGEN_COMP_MSVC
  // Visual Studio errors out on divisions by 0
  VERIFY(!(half(0.0 / 0.0) == half(0.0 / 0.0)));
  VERIFY(half(0.0 / 0.0) != half(0.0 / 0.0));

  VERIFY(!(half(1.0) == half(0.0 / 0.0)));
  VERIFY(!(half(1.0) < half(0.0 / 0.0)));
  VERIFY(!(half(1.0) > half(0.0 / 0.0)));
  VERIFY(half(1.0) != half(0.0 / 0.0));

  VERIFY(half(1.0) < half(1.0 / 0.0));
  VERIFY(half(1.0) > half(-1.0 / 0.0));
#endif
}

void test_basic_functions()
{
  VERIFY_IS_EQUAL(float(numext::abs(half(3.5f))), 3.5f);
  VERIFY_IS_EQUAL(float(abs(half(3.5f))), 3.5f);
  VERIFY_IS_EQUAL(float(numext::abs(half(-3.5f))), 3.5f);
  VERIFY_IS_EQUAL(float(abs(half(-3.5f))), 3.5f);

  VERIFY_IS_EQUAL(float(numext::floor(half(3.5f))), 3.0f);
  VERIFY_IS_EQUAL(float(floor(half(3.5f))), 3.0f);
  VERIFY_IS_EQUAL(float(numext::floor(half(-3.5f))), -4.0f);
  VERIFY_IS_EQUAL(float(floor(half(-3.5f))), -4.0f);

  VERIFY_IS_EQUAL(float(numext::ceil(half(3.5f))), 4.0f);
  VERIFY_IS_EQUAL(float(ceil(half(3.5f))), 4.0f);
  VERIFY_IS_EQUAL(float(numext::ceil(half(-3.5f))), -3.0f);
  VERIFY_IS_EQUAL(float(ceil(half(-3.5f))), -3.0f);

  VERIFY_IS_APPROX(float(numext::sqrt(half(0.0f))), 0.0f);
  VERIFY_IS_APPROX(float(sqrt(half(0.0f))), 0.0f);
  VERIFY_IS_APPROX(float(numext::sqrt(half(4.0f))), 2.0f);
  VERIFY_IS_APPROX(float(sqrt(half(4.0f))), 2.0f);

  VERIFY_IS_APPROX(float(numext::pow(half(0.0f), half(1.0f))), 0.0f);
  VERIFY_IS_APPROX(float(pow(half(0.0f), half(1.0f))), 0.0f);
  VERIFY_IS_APPROX(float(numext::pow(half(2.0f), half(2.0f))), 4.0f);
  VERIFY_IS_APPROX(float(pow(half(2.0f), half(2.0f))), 4.0f);

  VERIFY_IS_EQUAL(float(numext::exp(half(0.0f))), 1.0f);
  VERIFY_IS_EQUAL(float(exp(half(0.0f))), 1.0f);
  VERIFY_IS_APPROX(float(numext::exp(half(EIGEN_PI))), 20.f + float(EIGEN_PI));
  VERIFY_IS_APPROX(float(exp(half(EIGEN_PI))), 20.f + float(EIGEN_PI));

  VERIFY_IS_EQUAL(float(numext::expm1(half(0.0f))), 0.0f);
  VERIFY_IS_EQUAL(float(expm1(half(0.0f))), 0.0f);
  VERIFY_IS_APPROX(float(numext::expm1(half(2.0f))), 6.3890561f);
  VERIFY_IS_APPROX(float(expm1(half(2.0f))), 6.3890561f);

  VERIFY_IS_EQUAL(float(numext::log(half(1.0f))), 0.0f);
  VERIFY_IS_EQUAL(float(log(half(1.0f))), 0.0f);
  VERIFY_IS_APPROX(float(numext::log(half(10.0f))), 2.30273f);
  VERIFY_IS_APPROX(float(log(half(10.0f))), 2.30273f);

  VERIFY_IS_EQUAL(float(numext::log1p(half(0.0f))), 0.0f);
  VERIFY_IS_EQUAL(float(log1p(half(0.0f))), 0.0f);
  VERIFY_IS_APPROX(float(numext::log1p(half(10.0f))), 2.3978953f);
  VERIFY_IS_APPROX(float(log1p(half(10.0f))), 2.3978953f);
  
  VERIFY_IS_APPROX(numext::fmod(half(5.3f), half(2.0f)), half(1.3f));
  VERIFY_IS_APPROX(fmod(half(5.3f), half(2.0f)), half(1.3f));
  VERIFY_IS_APPROX(numext::fmod(half(-18.5f), half(-4.2f)), half(-1.7f));
  VERIFY_IS_APPROX(fmod(half(-18.5f), half(-4.2f)), half(-1.7f));
}

void test_trigonometric_functions()
{
  VERIFY_IS_APPROX(numext::cos(half(0.0f)), half(cosf(0.0f)));
  VERIFY_IS_APPROX(cos(half(0.0f)), half(cosf(0.0f)));
  VERIFY_IS_APPROX(numext::cos(half(EIGEN_PI)), half(cosf(EIGEN_PI)));
  // VERIFY_IS_APPROX(numext::cos(half(EIGEN_PI/2)), half(cosf(EIGEN_PI/2)));
  // VERIFY_IS_APPROX(numext::cos(half(3*EIGEN_PI/2)), half(cosf(3*EIGEN_PI/2)));
  VERIFY_IS_APPROX(numext::cos(half(3.5f)), half(cosf(3.5f)));

  VERIFY_IS_APPROX(numext::sin(half(0.0f)), half(sinf(0.0f)));
  VERIFY_IS_APPROX(sin(half(0.0f)), half(sinf(0.0f)));
  //  VERIFY_IS_APPROX(numext::sin(half(EIGEN_PI)), half(sinf(EIGEN_PI)));
  VERIFY_IS_APPROX(numext::sin(half(EIGEN_PI/2)), half(sinf(EIGEN_PI/2)));
  VERIFY_IS_APPROX(numext::sin(half(3*EIGEN_PI/2)), half(sinf(3*EIGEN_PI/2)));
  VERIFY_IS_APPROX(numext::sin(half(3.5f)), half(sinf(3.5f)));

  VERIFY_IS_APPROX(numext::tan(half(0.0f)), half(tanf(0.0f)));
  VERIFY_IS_APPROX(tan(half(0.0f)), half(tanf(0.0f)));
  //  VERIFY_IS_APPROX(numext::tan(half(EIGEN_PI)), half(tanf(EIGEN_PI)));
  //  VERIFY_IS_APPROX(numext::tan(half(EIGEN_PI/2)), half(tanf(EIGEN_PI/2)));
  //VERIFY_IS_APPROX(numext::tan(half(3*EIGEN_PI/2)), half(tanf(3*EIGEN_PI/2)));
  VERIFY_IS_APPROX(numext::tan(half(3.5f)), half(tanf(3.5f)));
}

void test_array()
{
  typedef Array<half,1,Dynamic> ArrayXh;
  Index size = internal::random<Index>(1,10);
  Index i = internal::random<Index>(0,size-1);
  ArrayXh a1 = ArrayXh::Random(size), a2 = ArrayXh::Random(size);
  VERIFY_IS_APPROX( a1+a1, half(2)*a1 );
  VERIFY( (a1.abs() >= half(0)).all() );
  VERIFY_IS_APPROX( (a1*a1).sqrt(), a1.abs() );

  VERIFY( ((a1.min)(a2) <= (a1.max)(a2)).all() );
  a1(i) = half(-10.);
  VERIFY_IS_EQUAL( a1.minCoeff(), half(-10.) );
  a1(i) = half(10.);
  VERIFY_IS_EQUAL( a1.maxCoeff(), half(10.) );

  std::stringstream ss;
  ss << a1;
}

void test_product()
{
  typedef Matrix<half,Dynamic,Dynamic> MatrixXh;
  Index rows  = internal::random<Index>(1,EIGEN_TEST_MAX_SIZE);
  Index cols  = internal::random<Index>(1,EIGEN_TEST_MAX_SIZE);
  Index depth = internal::random<Index>(1,EIGEN_TEST_MAX_SIZE);
  MatrixXh Ah = MatrixXh::Random(rows,depth);
  MatrixXh Bh = MatrixXh::Random(depth,cols);
  MatrixXh Ch = MatrixXh::Random(rows,cols);
  MatrixXf Af = Ah.cast<float>();
  MatrixXf Bf = Bh.cast<float>();
  MatrixXf Cf = Ch.cast<float>();
  VERIFY_IS_APPROX(Ch.noalias()+=Ah*Bh, (Cf.noalias()+=Af*Bf).cast<half>());
}

EIGEN_DECLARE_TEST(half_float)
{
  CALL_SUBTEST(test_numtraits());
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST(test_conversion());
    CALL_SUBTEST(test_arithmetic());
    CALL_SUBTEST(test_comparison());
    CALL_SUBTEST(test_basic_functions());
    CALL_SUBTEST(test_trigonometric_functions());
    CALL_SUBTEST(test_array());
    CALL_SUBTEST(test_product());
  }
}
