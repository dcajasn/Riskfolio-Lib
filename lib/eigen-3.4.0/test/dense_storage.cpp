// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2013 Hauke Heibel <hauke.heibel@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "main.h"
#include "AnnoyingScalar.h"
#include "SafeScalar.h"

#include <Eigen/Core>

#if EIGEN_HAS_TYPE_TRAITS && EIGEN_HAS_CXX11
using DenseStorageD3x3 = Eigen::DenseStorage<double, 3, 3, 3, 3>;
static_assert(std::is_trivially_move_constructible<DenseStorageD3x3>::value, "DenseStorage not trivially_move_constructible");
static_assert(std::is_trivially_move_assignable<DenseStorageD3x3>::value, "DenseStorage not trivially_move_assignable");
#if !defined(EIGEN_DENSE_STORAGE_CTOR_PLUGIN)
static_assert(std::is_trivially_copy_constructible<DenseStorageD3x3>::value, "DenseStorage not trivially_copy_constructible");
static_assert(std::is_trivially_copy_assignable<DenseStorageD3x3>::value, "DenseStorage not trivially_copy_assignable");
static_assert(std::is_trivially_copyable<DenseStorageD3x3>::value, "DenseStorage not trivially_copyable");
#endif
#endif

template <typename T, int Size, int Rows, int Cols>
void dense_storage_copy(int rows, int cols)
{
  typedef DenseStorage<T, Size, Rows, Cols, 0> DenseStorageType;
  
  const int size = rows*cols;
  DenseStorageType reference(size, rows, cols);
  T* raw_reference = reference.data();
  for (int i=0; i<size; ++i)
    raw_reference[i] = static_cast<T>(i);
    
  DenseStorageType copied_reference(reference);
  const T* raw_copied_reference = copied_reference.data();
  for (int i=0; i<size; ++i)
    VERIFY_IS_EQUAL(raw_reference[i], raw_copied_reference[i]);
}

template <typename T, int Size, int Rows, int Cols>
void dense_storage_assignment(int rows, int cols)
{
  typedef DenseStorage<T, Size, Rows, Cols, 0> DenseStorageType;
  
  const int size = rows*cols;
  DenseStorageType reference(size, rows, cols);
  T* raw_reference = reference.data();
  for (int i=0; i<size; ++i)
    raw_reference[i] = static_cast<T>(i);
    
  DenseStorageType copied_reference;
  copied_reference = reference;
  const T* raw_copied_reference = copied_reference.data();
  for (int i=0; i<size; ++i)
    VERIFY_IS_EQUAL(raw_reference[i], raw_copied_reference[i]);
}

template <typename T, int Size, int Rows, int Cols>
void dense_storage_swap(int rows0, int cols0, int rows1, int cols1)
{
  typedef DenseStorage<T, Size, Rows, Cols, 0> DenseStorageType;
  
  const int size0 = rows0*cols0;
  DenseStorageType a(size0, rows0, cols0);
  for (int i=0; i<size0; ++i) {
    a.data()[i] = static_cast<T>(i);
  }
  
  const int size1 = rows1*cols1;
  DenseStorageType b(size1, rows1, cols1);
  for (int i=0; i<size1; ++i) {
    b.data()[i] = static_cast<T>(-i);
  }
  
  a.swap(b);
  
  for (int i=0; i<size0; ++i) {
    VERIFY_IS_EQUAL(b.data()[i], static_cast<T>(i));
  }
  
  for (int i=0; i<size1; ++i) {
    VERIFY_IS_EQUAL(a.data()[i], static_cast<T>(-i));
  }
}

template<typename T, int Size, std::size_t Alignment>
void dense_storage_alignment()
{
  #if EIGEN_HAS_ALIGNAS
  
  struct alignas(Alignment) Empty1 {};
  VERIFY_IS_EQUAL(std::alignment_of<Empty1>::value, Alignment);

  struct EIGEN_ALIGN_TO_BOUNDARY(Alignment) Empty2 {};
  VERIFY_IS_EQUAL(std::alignment_of<Empty2>::value, Alignment);

  struct Nested1 { EIGEN_ALIGN_TO_BOUNDARY(Alignment) T data[Size]; };
  VERIFY_IS_EQUAL(std::alignment_of<Nested1>::value, Alignment);

  VERIFY_IS_EQUAL( (std::alignment_of<internal::plain_array<T,Size,AutoAlign,Alignment> >::value), Alignment);

  const std::size_t default_alignment = internal::compute_default_alignment<T,Size>::value;

  VERIFY_IS_EQUAL( (std::alignment_of<DenseStorage<T,Size,1,1,AutoAlign> >::value), default_alignment);
  VERIFY_IS_EQUAL( (std::alignment_of<Matrix<T,Size,1,AutoAlign> >::value), default_alignment);
  struct Nested2 { Matrix<T,Size,1,AutoAlign> mat; };
  VERIFY_IS_EQUAL(std::alignment_of<Nested2>::value, default_alignment);

  #endif
}

template<typename T>
void dense_storage_tests() {
  // Dynamic Storage.
  dense_storage_copy<T,Dynamic,Dynamic,Dynamic>(4, 3);  
  dense_storage_copy<T,Dynamic,Dynamic,3>(4, 3);
  dense_storage_copy<T,Dynamic,4,Dynamic>(4, 3);
  // Fixed Storage.
  dense_storage_copy<T,12,4,3>(4, 3);
  dense_storage_copy<T,12,Dynamic,Dynamic>(4, 3);
  dense_storage_copy<T,12,4,Dynamic>(4, 3);
  dense_storage_copy<T,12,Dynamic,3>(4, 3);
  // Fixed Storage with Uninitialized Elements.
  dense_storage_copy<T,18,Dynamic,Dynamic>(4, 3);
  dense_storage_copy<T,18,4,Dynamic>(4, 3);
  dense_storage_copy<T,18,Dynamic,3>(4, 3);
  
  // Dynamic Storage.
  dense_storage_assignment<T,Dynamic,Dynamic,Dynamic>(4, 3);  
  dense_storage_assignment<T,Dynamic,Dynamic,3>(4, 3);
  dense_storage_assignment<T,Dynamic,4,Dynamic>(4, 3);
  // Fixed Storage.
  dense_storage_assignment<T,12,4,3>(4, 3);
  dense_storage_assignment<T,12,Dynamic,Dynamic>(4, 3);
  dense_storage_assignment<T,12,4,Dynamic>(4, 3);
  dense_storage_assignment<T,12,Dynamic,3>(4, 3);
  // Fixed Storage with Uninitialized Elements.
  dense_storage_assignment<T,18,Dynamic,Dynamic>(4, 3);
  dense_storage_assignment<T,18,4,Dynamic>(4, 3);
  dense_storage_assignment<T,18,Dynamic,3>(4, 3);
  
  // Dynamic Storage.
  dense_storage_swap<T,Dynamic,Dynamic,Dynamic>(4, 3, 4, 3); 
  dense_storage_swap<T,Dynamic,Dynamic,Dynamic>(4, 3, 2, 1);  
  dense_storage_swap<T,Dynamic,Dynamic,Dynamic>(2, 1, 4, 3);
  dense_storage_swap<T,Dynamic,Dynamic,3>(4, 3, 4, 3);
  dense_storage_swap<T,Dynamic,Dynamic,3>(4, 3, 2, 3);
  dense_storage_swap<T,Dynamic,Dynamic,3>(2, 3, 4, 3);
  dense_storage_swap<T,Dynamic,4,Dynamic>(4, 3, 4, 3);
  dense_storage_swap<T,Dynamic,4,Dynamic>(4, 3, 4, 1);
  dense_storage_swap<T,Dynamic,4,Dynamic>(4, 1, 4, 3);
  // Fixed Storage.
  dense_storage_swap<T,12,4,3>(4, 3, 4, 3);
  dense_storage_swap<T,12,Dynamic,Dynamic>(4, 3, 4, 3);
  dense_storage_swap<T,12,Dynamic,Dynamic>(4, 3, 2, 1);
  dense_storage_swap<T,12,Dynamic,Dynamic>(2, 1, 4, 3);
  dense_storage_swap<T,12,4,Dynamic>(4, 3, 4, 3);
  dense_storage_swap<T,12,4,Dynamic>(4, 3, 4, 1);
  dense_storage_swap<T,12,4,Dynamic>(4, 1, 4, 3);
  dense_storage_swap<T,12,Dynamic,3>(4, 3, 4, 3);
  dense_storage_swap<T,12,Dynamic,3>(4, 3, 2, 3);
  dense_storage_swap<T,12,Dynamic,3>(2, 3, 4, 3);
  // Fixed Storage with Uninitialized Elements.
  dense_storage_swap<T,18,Dynamic,Dynamic>(4, 3, 4, 3);
  dense_storage_swap<T,18,Dynamic,Dynamic>(4, 3, 2, 1);
  dense_storage_swap<T,18,Dynamic,Dynamic>(2, 1, 4, 3);
  dense_storage_swap<T,18,4,Dynamic>(4, 3, 4, 3);
  dense_storage_swap<T,18,4,Dynamic>(4, 3, 4, 1);
  dense_storage_swap<T,18,4,Dynamic>(4, 1, 4, 3);
  dense_storage_swap<T,18,Dynamic,3>(4, 3, 4, 3);
  dense_storage_swap<T,18,Dynamic,3>(4, 3, 2, 3);
  dense_storage_swap<T,18,Dynamic,3>(2, 3, 4, 3);
  
  dense_storage_alignment<T,16,8>();
  dense_storage_alignment<T,16,16>();
  dense_storage_alignment<T,16,32>();
  dense_storage_alignment<T,16,64>();
}

EIGEN_DECLARE_TEST(dense_storage)
{
  dense_storage_tests<int>();
  dense_storage_tests<float>();
  dense_storage_tests<SafeScalar<float> >();
  dense_storage_tests<AnnoyingScalar>();
}
