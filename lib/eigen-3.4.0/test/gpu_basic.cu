// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2015-2016 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

// workaround issue between gcc >= 4.7 and cuda 5.5
#if (defined __GNUC__) && (__GNUC__>4 || __GNUC_MINOR__>=7)
  #undef _GLIBCXX_ATOMIC_BUILTINS
  #undef _GLIBCXX_USE_INT128
#endif

#define EIGEN_TEST_NO_LONGDOUBLE
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int

#include "main.h"
#include "gpu_common.h"

// Check that dense modules can be properly parsed by nvcc
#include <Eigen/Dense>

// struct Foo{
//   EIGEN_DEVICE_FUNC
//   void operator()(int i, const float* mats, float* vecs) const {
//     using namespace Eigen;
//   //   Matrix3f M(data);
//   //   Vector3f x(data+9);
//   //   Map<Vector3f>(data+9) = M.inverse() * x;
//     Matrix3f M(mats+i/16);
//     Vector3f x(vecs+i*3);
//   //   using std::min;
//   //   using std::sqrt;
//     Map<Vector3f>(vecs+i*3) << x.minCoeff(), 1, 2;// / x.dot(x);//(M.inverse() *  x) / x.x();
//     //x = x*2 + x.y() * x + x * x.maxCoeff() - x / x.sum();
//   }
// };

template<typename T>
struct coeff_wise {
  EIGEN_DEVICE_FUNC
  void operator()(int i, const typename T::Scalar* in, typename T::Scalar* out) const
  {
    using namespace Eigen;
    T x1(in+i);
    T x2(in+i+1);
    T x3(in+i+2);
    Map<T> res(out+i*T::MaxSizeAtCompileTime);
    
    res.array() += (in[0] * x1 + x2).array() * x3.array();
  }
};

template<typename T>
struct complex_sqrt {
  EIGEN_DEVICE_FUNC
  void operator()(int i, const typename T::Scalar* in, typename T::Scalar* out) const
  {
    using namespace Eigen;
    typedef typename T::Scalar ComplexType;
    typedef typename T::Scalar::value_type ValueType;
    const int num_special_inputs = 18;
    
    if (i == 0) {
      const ValueType nan = std::numeric_limits<ValueType>::quiet_NaN();
      typedef Eigen::Vector<ComplexType, num_special_inputs> SpecialInputs;
      SpecialInputs special_in;
      special_in.setZero();
      int idx = 0;
      special_in[idx++] = ComplexType(0, 0);
      special_in[idx++] = ComplexType(-0, 0);
      special_in[idx++] = ComplexType(0, -0);
      special_in[idx++] = ComplexType(-0, -0);
      // GCC's fallback sqrt implementation fails for inf inputs.
      // It is called when _GLIBCXX_USE_C99_COMPLEX is false or if
      // clang includes the GCC header (which temporarily disables
      // _GLIBCXX_USE_C99_COMPLEX)
      #if !defined(_GLIBCXX_COMPLEX) || \
        (_GLIBCXX_USE_C99_COMPLEX && !defined(__CLANG_CUDA_WRAPPERS_COMPLEX))
      const ValueType inf = std::numeric_limits<ValueType>::infinity();
      special_in[idx++] = ComplexType(1.0, inf);
      special_in[idx++] = ComplexType(nan, inf);
      special_in[idx++] = ComplexType(1.0, -inf);
      special_in[idx++] = ComplexType(nan, -inf);
      special_in[idx++] = ComplexType(-inf, 1.0);
      special_in[idx++] = ComplexType(inf, 1.0);
      special_in[idx++] = ComplexType(-inf, -1.0);
      special_in[idx++] = ComplexType(inf, -1.0);
      special_in[idx++] = ComplexType(-inf, nan);
      special_in[idx++] = ComplexType(inf, nan);
      #endif
      special_in[idx++] = ComplexType(1.0, nan);
      special_in[idx++] = ComplexType(nan, 1.0);
      special_in[idx++] = ComplexType(nan, -1.0);
      special_in[idx++] = ComplexType(nan, nan);
      
      Map<SpecialInputs> special_out(out);
      special_out = special_in.cwiseSqrt();
    }
    
    T x1(in + i);
    Map<T> res(out + num_special_inputs + i*T::MaxSizeAtCompileTime);
    res = x1.cwiseSqrt();
  }
};

template<typename T>
struct complex_operators {
  EIGEN_DEVICE_FUNC
  void operator()(int i, const typename T::Scalar* in, typename T::Scalar* out) const
  {
    using namespace Eigen;
    typedef typename T::Scalar ComplexType;
    typedef typename T::Scalar::value_type ValueType;
    const int num_scalar_operators = 24;
    const int num_vector_operators = 23;  // no unary + operator.
    int out_idx = i * (num_scalar_operators + num_vector_operators * T::MaxSizeAtCompileTime);
    
    // Scalar operators.
    const ComplexType a = in[i];
    const ComplexType b = in[i + 1];
    
    out[out_idx++] = +a;
    out[out_idx++] = -a;
    
    out[out_idx++] = a + b;
    out[out_idx++] = a + numext::real(b);
    out[out_idx++] = numext::real(a) + b;
    out[out_idx++] = a - b;
    out[out_idx++] = a - numext::real(b);
    out[out_idx++] = numext::real(a) - b;
    out[out_idx++] = a * b;
    out[out_idx++] = a * numext::real(b);
    out[out_idx++] = numext::real(a) * b;
    out[out_idx++] = a / b;
    out[out_idx++] = a / numext::real(b);
    out[out_idx++] = numext::real(a) / b;
    
    out[out_idx] = a; out[out_idx++] += b;
    out[out_idx] = a; out[out_idx++] -= b;
    out[out_idx] = a; out[out_idx++] *= b;
    out[out_idx] = a; out[out_idx++] /= b;
    
    const ComplexType true_value = ComplexType(ValueType(1), ValueType(0));
    const ComplexType false_value = ComplexType(ValueType(0), ValueType(0));
    out[out_idx++] = (a == b ? true_value : false_value);
    out[out_idx++] = (a == numext::real(b) ? true_value : false_value);
    out[out_idx++] = (numext::real(a) == b ? true_value : false_value);
    out[out_idx++] = (a != b ? true_value : false_value);
    out[out_idx++] = (a != numext::real(b) ? true_value : false_value);
    out[out_idx++] = (numext::real(a) != b ? true_value : false_value);
    
    // Vector versions.
    T x1(in + i);
    T x2(in + i + 1);
    const int res_size = T::MaxSizeAtCompileTime * num_scalar_operators;
    const int size = T::MaxSizeAtCompileTime;
    int block_idx = 0;
    
    Map<VectorX<ComplexType>> res(out + out_idx, res_size);
    res.segment(block_idx, size) = -x1;
    block_idx += size;
    
    res.segment(block_idx, size) = x1 + x2;
    block_idx += size;
    res.segment(block_idx, size) = x1 + x2.real();
    block_idx += size;
    res.segment(block_idx, size) = x1.real() + x2;
    block_idx += size;
    res.segment(block_idx, size) = x1 - x2;
    block_idx += size;
    res.segment(block_idx, size) = x1 - x2.real();
    block_idx += size;
    res.segment(block_idx, size) = x1.real() - x2;
    block_idx += size;
    res.segment(block_idx, size) = x1.array() * x2.array();
    block_idx += size;
    res.segment(block_idx, size) = x1.array() * x2.real().array();
    block_idx += size;
    res.segment(block_idx, size) = x1.real().array() * x2.array();
    block_idx += size;
    res.segment(block_idx, size) = x1.array() / x2.array();
    block_idx += size;
    res.segment(block_idx, size) = x1.array() / x2.real().array();
    block_idx += size;
    res.segment(block_idx, size) = x1.real().array() / x2.array();
    block_idx += size;
    
    res.segment(block_idx, size) = x1; res.segment(block_idx, size) += x2;
    block_idx += size;
    res.segment(block_idx, size) = x1; res.segment(block_idx, size) -= x2;
    block_idx += size;
    res.segment(block_idx, size) = x1; res.segment(block_idx, size).array() *= x2.array();
    block_idx += size;
    res.segment(block_idx, size) = x1; res.segment(block_idx, size).array() /= x2.array();
    block_idx += size;

    const T true_vector = T::Constant(true_value);
    const T false_vector = T::Constant(false_value);
    res.segment(block_idx, size) = (x1 == x2 ? true_vector : false_vector);
    block_idx += size;
    // Mixing types in equality comparison does not work.
    // res.segment(block_idx, size) = (x1 == x2.real() ? true_vector : false_vector);
    // block_idx += size;
    // res.segment(block_idx, size) = (x1.real() == x2 ? true_vector : false_vector);
    // block_idx += size;
    res.segment(block_idx, size) = (x1 != x2 ? true_vector : false_vector);
    block_idx += size;
    // res.segment(block_idx, size) = (x1 != x2.real() ? true_vector : false_vector);
    // block_idx += size;
    // res.segment(block_idx, size) = (x1.real() != x2 ? true_vector : false_vector);
    // block_idx += size;
  }
};

template<typename T>
struct replicate {
  EIGEN_DEVICE_FUNC
  void operator()(int i, const typename T::Scalar* in, typename T::Scalar* out) const
  {
    using namespace Eigen;
    T x1(in+i);
    int step   = x1.size() * 4;
    int stride = 3 * step;
    
    typedef Map<Array<typename T::Scalar,Dynamic,Dynamic> > MapType;
    MapType(out+i*stride+0*step, x1.rows()*2, x1.cols()*2) = x1.replicate(2,2);
    MapType(out+i*stride+1*step, x1.rows()*3, x1.cols()) = in[i] * x1.colwise().replicate(3);
    MapType(out+i*stride+2*step, x1.rows(), x1.cols()*3) = in[i] * x1.rowwise().replicate(3);
  }
};

template<typename T>
struct alloc_new_delete {
  EIGEN_DEVICE_FUNC
  void operator()(int i, const typename T::Scalar* in, typename T::Scalar* out) const
  {
    int offset = 2*i*T::MaxSizeAtCompileTime;
    T* x = new T(in + offset);
    Eigen::Map<T> u(out + offset);
    u = *x;
    delete x;
    
    offset += T::MaxSizeAtCompileTime;
    T* y = new T[1];
    y[0] = T(in + offset);
    Eigen::Map<T> v(out + offset);
    v = y[0];    
    delete[] y;
  }
};

template<typename T>
struct redux {
  EIGEN_DEVICE_FUNC
  void operator()(int i, const typename T::Scalar* in, typename T::Scalar* out) const
  {
    using namespace Eigen;
    int N = 10;
    T x1(in+i);
    out[i*N+0] = x1.minCoeff();
    out[i*N+1] = x1.maxCoeff();
    out[i*N+2] = x1.sum();
    out[i*N+3] = x1.prod();
    out[i*N+4] = x1.matrix().squaredNorm();
    out[i*N+5] = x1.matrix().norm();
    out[i*N+6] = x1.colwise().sum().maxCoeff();
    out[i*N+7] = x1.rowwise().maxCoeff().sum();
    out[i*N+8] = x1.matrix().colwise().squaredNorm().sum();
  }
};

template<typename T1, typename T2>
struct prod_test {
  EIGEN_DEVICE_FUNC
  void operator()(int i, const typename T1::Scalar* in, typename T1::Scalar* out) const
  {
    using namespace Eigen;
    typedef Matrix<typename T1::Scalar, T1::RowsAtCompileTime, T2::ColsAtCompileTime> T3;
    T1 x1(in+i);
    T2 x2(in+i+1);
    Map<T3> res(out+i*T3::MaxSizeAtCompileTime);
    res += in[i] * x1 * x2;
  }
};

template<typename T1, typename T2>
struct diagonal {
  EIGEN_DEVICE_FUNC
  void operator()(int i, const typename T1::Scalar* in, typename T1::Scalar* out) const
  {
    using namespace Eigen;
    T1 x1(in+i);
    Map<T2> res(out+i*T2::MaxSizeAtCompileTime);
    res += x1.diagonal();
  }
};

template<typename T>
struct eigenvalues_direct {
  EIGEN_DEVICE_FUNC
  void operator()(int i, const typename T::Scalar* in, typename T::Scalar* out) const
  {
    using namespace Eigen;
    typedef Matrix<typename T::Scalar, T::RowsAtCompileTime, 1> Vec;
    T M(in+i);
    Map<Vec> res(out+i*Vec::MaxSizeAtCompileTime);
    T A = M*M.adjoint();
    SelfAdjointEigenSolver<T> eig;
    eig.computeDirect(A);
    res = eig.eigenvalues();
  }
};

template<typename T>
struct eigenvalues {
  EIGEN_DEVICE_FUNC
  void operator()(int i, const typename T::Scalar* in, typename T::Scalar* out) const
  {
    using namespace Eigen;
    typedef Matrix<typename T::Scalar, T::RowsAtCompileTime, 1> Vec;
    T M(in+i);
    Map<Vec> res(out+i*Vec::MaxSizeAtCompileTime);
    T A = M*M.adjoint();
    SelfAdjointEigenSolver<T> eig;
    eig.compute(A);
    res = eig.eigenvalues();
  }
};

template<typename T>
struct matrix_inverse {
  EIGEN_DEVICE_FUNC
  void operator()(int i, const typename T::Scalar* in, typename T::Scalar* out) const
  {
    using namespace Eigen;
    T M(in+i);
    Map<T> res(out+i*T::MaxSizeAtCompileTime);
    res = M.inverse();
  }
};

template<typename T>
struct numeric_limits_test {
  EIGEN_DEVICE_FUNC
  void operator()(int i, const typename T::Scalar* in, typename T::Scalar* out) const
  {
    EIGEN_UNUSED_VARIABLE(in)
    int out_idx = i * 5;
    out[out_idx++] = numext::numeric_limits<float>::epsilon();
    out[out_idx++] = (numext::numeric_limits<float>::max)();
    out[out_idx++] = (numext::numeric_limits<float>::min)();
    out[out_idx++] = numext::numeric_limits<float>::infinity();
    out[out_idx++] = numext::numeric_limits<float>::quiet_NaN();
  }
};

template<typename Type1, typename Type2>
bool verifyIsApproxWithInfsNans(const Type1& a, const Type2& b, typename Type1::Scalar* = 0) // Enabled for Eigen's type only
{
  if (a.rows() != b.rows()) {
    return false;
  }
  if (a.cols() != b.cols()) {
    return false;
  }
  for (Index r = 0; r < a.rows(); ++r) {
    for (Index c = 0; c < a.cols(); ++c) {
      if (a(r, c) != b(r, c)
          && !((numext::isnan)(a(r, c)) && (numext::isnan)(b(r, c))) 
          && !test_isApprox(a(r, c), b(r, c))) {
        return false;
      }
    }
  }
  return true;
}

template<typename Kernel, typename Input, typename Output>
void test_with_infs_nans(const Kernel& ker, int n, const Input& in, Output& out)
{
  Output out_ref, out_gpu;
  #if !defined(EIGEN_GPU_COMPILE_PHASE)
  out_ref = out_gpu = out;
  #else
  EIGEN_UNUSED_VARIABLE(in);
  EIGEN_UNUSED_VARIABLE(out);
  #endif
  run_on_cpu (ker, n, in,  out_ref);
  run_on_gpu(ker, n, in, out_gpu);
  #if !defined(EIGEN_GPU_COMPILE_PHASE)
  verifyIsApproxWithInfsNans(out_ref, out_gpu);
  #endif
}

EIGEN_DECLARE_TEST(gpu_basic)
{
  ei_test_init_gpu();
  
  int nthreads = 100;
  Eigen::VectorXf in, out;
  Eigen::VectorXcf cfin, cfout;
  
  #if !defined(EIGEN_GPU_COMPILE_PHASE)
  int data_size = nthreads * 512;
  in.setRandom(data_size);
  out.setConstant(data_size, -1);
  cfin.setRandom(data_size);
  cfout.setConstant(data_size, -1);
  #endif
  
  CALL_SUBTEST( run_and_compare_to_gpu(coeff_wise<Vector3f>(), nthreads, in, out) );
  CALL_SUBTEST( run_and_compare_to_gpu(coeff_wise<Array44f>(), nthreads, in, out) );

#if !defined(EIGEN_USE_HIP)
  // FIXME
  // These subtests result in a compile failure on the HIP platform
  //
  //  eigen-upstream/Eigen/src/Core/Replicate.h:61:65: error:
  //           base class 'internal::dense_xpr_base<Replicate<Array<float, 4, 1, 0, 4, 1>, -1, -1> >::type'
  //           (aka 'ArrayBase<Eigen::Replicate<Eigen::Array<float, 4, 1, 0, 4, 1>, -1, -1> >') has protected default constructor
  CALL_SUBTEST( run_and_compare_to_gpu(replicate<Array4f>(), nthreads, in, out) );
  CALL_SUBTEST( run_and_compare_to_gpu(replicate<Array33f>(), nthreads, in, out) );

  // HIP does not support new/delete on device.
  CALL_SUBTEST( run_and_compare_to_gpu(alloc_new_delete<Vector3f>(), nthreads, in, out) );
#endif
  
  CALL_SUBTEST( run_and_compare_to_gpu(redux<Array4f>(), nthreads, in, out) );
  CALL_SUBTEST( run_and_compare_to_gpu(redux<Matrix3f>(), nthreads, in, out) );
  
  CALL_SUBTEST( run_and_compare_to_gpu(prod_test<Matrix3f,Matrix3f>(), nthreads, in, out) );
  CALL_SUBTEST( run_and_compare_to_gpu(prod_test<Matrix4f,Vector4f>(), nthreads, in, out) );
  
  CALL_SUBTEST( run_and_compare_to_gpu(diagonal<Matrix3f,Vector3f>(), nthreads, in, out) );
  CALL_SUBTEST( run_and_compare_to_gpu(diagonal<Matrix4f,Vector4f>(), nthreads, in, out) );

  CALL_SUBTEST( run_and_compare_to_gpu(matrix_inverse<Matrix2f>(), nthreads, in, out) );
  CALL_SUBTEST( run_and_compare_to_gpu(matrix_inverse<Matrix3f>(), nthreads, in, out) );
  CALL_SUBTEST( run_and_compare_to_gpu(matrix_inverse<Matrix4f>(), nthreads, in, out) );
  
  CALL_SUBTEST( run_and_compare_to_gpu(eigenvalues_direct<Matrix3f>(), nthreads, in, out) );
  CALL_SUBTEST( run_and_compare_to_gpu(eigenvalues_direct<Matrix2f>(), nthreads, in, out) );

  // Test std::complex.
  CALL_SUBTEST( run_and_compare_to_gpu(complex_operators<Vector3cf>(), nthreads, cfin, cfout) );
  CALL_SUBTEST( test_with_infs_nans(complex_sqrt<Vector3cf>(), nthreads, cfin, cfout) );

  // numeric_limits
  CALL_SUBTEST( test_with_infs_nans(numeric_limits_test<Vector3f>(), 1, in, out) );

#if defined(__NVCC__)
  // FIXME
  // These subtests compiles only with nvcc and fail with HIPCC and clang-cuda
  CALL_SUBTEST( run_and_compare_to_gpu(eigenvalues<Matrix4f>(), nthreads, in, out) );
  typedef Matrix<float,6,6> Matrix6f;
  CALL_SUBTEST( run_and_compare_to_gpu(eigenvalues<Matrix6f>(), nthreads, in, out) );
#endif
}
