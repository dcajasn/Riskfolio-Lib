// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008-2009 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "main.h"
#include <Eigen/Geometry>

using namespace std;

// NOTE the following workaround was needed on some 32 bits builds to kill extra precision of x87 registers.
// It seems that it is not needed anymore, but let's keep it here, just in case...

template<typename T> EIGEN_DONT_INLINE
void kill_extra_precision(T& /* x */) {
  // This one worked but triggered a warning:
  /* eigen_assert((void*)(&x) != (void*)0); */
  // An alternative could be:
  /* volatile T tmp = x; */
  /* x = tmp; */
}


template<typename BoxType> void alignedbox(const BoxType& box)
{
  /* this test covers the following files:
     AlignedBox.h
  */
  typedef typename BoxType::Scalar Scalar;
  typedef NumTraits<Scalar> ScalarTraits;
  typedef typename ScalarTraits::Real RealScalar;
  typedef Matrix<Scalar, BoxType::AmbientDimAtCompileTime, 1> VectorType;

  const Index dim = box.dim();

  VectorType p0 = VectorType::Random(dim);
  VectorType p1 = VectorType::Random(dim);
  while( p1 == p0 ){
      p1 =  VectorType::Random(dim); }
  RealScalar s1 = internal::random<RealScalar>(0,1);

  BoxType b0(dim);
  BoxType b1(VectorType::Random(dim),VectorType::Random(dim));
  BoxType b2;

  kill_extra_precision(b1);
  kill_extra_precision(p0);
  kill_extra_precision(p1);

  b0.extend(p0);
  b0.extend(p1);
  VERIFY(b0.contains(p0*s1+(Scalar(1)-s1)*p1));
  VERIFY(b0.contains(b0.center()));
  VERIFY_IS_APPROX(b0.center(),(p0+p1)/Scalar(2));

  (b2 = b0).extend(b1);
  VERIFY(b2.contains(b0));
  VERIFY(b2.contains(b1));
  VERIFY_IS_APPROX(b2.clamp(b0), b0);

  // intersection
  BoxType box1(VectorType::Random(dim));
  box1.extend(VectorType::Random(dim));
  BoxType box2(VectorType::Random(dim));
  box2.extend(VectorType::Random(dim));

  VERIFY(box1.intersects(box2) == !box1.intersection(box2).isEmpty());

  // alignment -- make sure there is no memory alignment assertion
  BoxType *bp0 = new BoxType(dim);
  BoxType *bp1 = new BoxType(dim);
  bp0->extend(*bp1);
  delete bp0;
  delete bp1;

  // sampling
  for( int i=0; i<10; ++i )
  {
      VectorType r = b0.sample();
      VERIFY(b0.contains(r));
  }

}

template<typename BoxType> void alignedboxTranslatable(const BoxType& box)
{
  typedef typename BoxType::Scalar Scalar;
  typedef Matrix<Scalar, BoxType::AmbientDimAtCompileTime, 1> VectorType;
  typedef Transform<Scalar, BoxType::AmbientDimAtCompileTime, Isometry> IsometryTransform;
  typedef Transform<Scalar, BoxType::AmbientDimAtCompileTime, Affine> AffineTransform;

  alignedbox(box);

  const VectorType Ones = VectorType::Ones();
  const VectorType UnitX = VectorType::UnitX();
  const Index dim = box.dim();

  // box((-1, -1, -1), (1, 1, 1))
  BoxType a(-Ones, Ones);

  VERIFY_IS_APPROX(a.sizes(), Ones * Scalar(2));

  BoxType b = a;
  VectorType translate = Ones;
  translate[0] = Scalar(2);
  b.translate(translate);
  // translate by (2, 1, 1) -> box((1, 0, 0), (3, 2, 2))

  VERIFY_IS_APPROX(b.sizes(), Ones * Scalar(2));
  VERIFY_IS_APPROX((b.min)(), UnitX);
  VERIFY_IS_APPROX((b.max)(), Ones * Scalar(2) + UnitX);

  // Test transform

  IsometryTransform tf = IsometryTransform::Identity();
  tf.translation() = -translate;

  BoxType c = b.transformed(tf);
  // translate by (-2, -1, -1) -> box((-1, -1, -1), (1, 1, 1))
  VERIFY_IS_APPROX(c.sizes(), a.sizes());
  VERIFY_IS_APPROX((c.min)(), (a.min)());
  VERIFY_IS_APPROX((c.max)(), (a.max)());

  c.transform(tf);
  // translate by (-2, -1, -1) -> box((-3, -2, -2), (-1, 0, 0))
  VERIFY_IS_APPROX(c.sizes(), a.sizes());
  VERIFY_IS_APPROX((c.min)(), Ones * Scalar(-2) - UnitX);
  VERIFY_IS_APPROX((c.max)(), -UnitX);

  // Scaling

  AffineTransform atf = AffineTransform::Identity();
  atf.scale(Scalar(3));
  c.transform(atf);
  // scale by 3 -> box((-9, -6, -6), (-3, 0, 0))
  VERIFY_IS_APPROX(c.sizes(), Scalar(3) * a.sizes());
  VERIFY_IS_APPROX((c.min)(), Ones * Scalar(-6) - UnitX * Scalar(3));
  VERIFY_IS_APPROX((c.max)(), UnitX * Scalar(-3));

  atf = AffineTransform::Identity();
  atf.scale(Scalar(-3));
  c.transform(atf);
  // scale by -3 -> box((27, 18, 18), (9, 0, 0))
  VERIFY_IS_APPROX(c.sizes(), Scalar(9) * a.sizes());
  VERIFY_IS_APPROX((c.min)(), UnitX * Scalar(9));
  VERIFY_IS_APPROX((c.max)(), Ones * Scalar(18) + UnitX * Scalar(9));

  // Check identity transform within numerical precision.
  BoxType transformedC = c.transformed(IsometryTransform::Identity());
  VERIFY_IS_APPROX(transformedC, c);

  for (size_t i = 0; i < 10; ++i)
  {
    VectorType minCorner;
    VectorType maxCorner;
    for (Index d = 0; d < dim; ++d)
    {
      minCorner[d] = internal::random<Scalar>(-10,10);
      maxCorner[d] = minCorner[d] + internal::random<Scalar>(0, 10);
    }

    c = BoxType(minCorner, maxCorner);

    translate = VectorType::Random();
    c.translate(translate);

    VERIFY_IS_APPROX((c.min)(), minCorner + translate);
    VERIFY_IS_APPROX((c.max)(), maxCorner + translate);
  }
}

template<typename Scalar, typename Rotation>
Rotation rotate2D(Scalar angle) {
  return Rotation2D<Scalar>(angle);
}

template<typename Scalar, typename Rotation>
Rotation rotate2DIntegral(typename NumTraits<Scalar>::NonInteger angle) {
  typedef typename NumTraits<Scalar>::NonInteger NonInteger;
  return Rotation2D<NonInteger>(angle).toRotationMatrix().
      template cast<Scalar>();
}

template<typename Scalar, typename Rotation>
Rotation rotate3DZAxis(Scalar angle) {
  return AngleAxis<Scalar>(angle, Matrix<Scalar, 3, 1>(0, 0, 1));
}

template<typename Scalar, typename Rotation>
Rotation rotate3DZAxisIntegral(typename NumTraits<Scalar>::NonInteger angle) {
  typedef typename NumTraits<Scalar>::NonInteger NonInteger;
  return AngleAxis<NonInteger>(angle, Matrix<NonInteger, 3, 1>(0, 0, 1)).
      toRotationMatrix().template cast<Scalar>();
}

template<typename Scalar, typename Rotation>
Rotation rotate4DZWAxis(Scalar angle) {
  Rotation result = Matrix<Scalar, 4, 4>::Identity();
  result.block(0, 0, 3, 3) = rotate3DZAxis<Scalar, AngleAxisd>(angle).toRotationMatrix();
  return result;
}

template <typename MatrixType>
MatrixType randomRotationMatrix()
{
  // algorithm from
  // https://www.isprs-ann-photogramm-remote-sens-spatial-inf-sci.net/III-7/103/2016/isprs-annals-III-7-103-2016.pdf
  const MatrixType rand = MatrixType::Random();
  const MatrixType q = rand.householderQr().householderQ();
  const JacobiSVD<MatrixType> svd = q.jacobiSvd(ComputeFullU | ComputeFullV);
  const typename MatrixType::Scalar det = (svd.matrixU() * svd.matrixV().transpose()).determinant();
  MatrixType diag = rand.Identity();
  diag(MatrixType::RowsAtCompileTime - 1, MatrixType::ColsAtCompileTime - 1) = det;
  const MatrixType rotation = svd.matrixU() * diag * svd.matrixV().transpose();
  return rotation;
}

template <typename Scalar, int Dim>
Matrix<Scalar, Dim, (1<<Dim)> boxGetCorners(const Matrix<Scalar, Dim, 1>& min_, const Matrix<Scalar, Dim, 1>& max_)
{
  Matrix<Scalar, Dim, (1<<Dim) > result;
  for(Index i=0; i<(1<<Dim); ++i)
  {
    for(Index j=0; j<Dim; ++j)
      result(j,i) = (i & (1<<j)) ? min_(j) : max_(j);
  }
  return result;
}

template<typename BoxType, typename Rotation> void alignedboxRotatable(
    const BoxType& box,
    Rotation (*rotate)(typename NumTraits<typename BoxType::Scalar>::NonInteger /*_angle*/))
{
  alignedboxTranslatable(box);

  typedef typename BoxType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::NonInteger NonInteger;
  typedef Matrix<Scalar, BoxType::AmbientDimAtCompileTime, 1> VectorType;
  typedef Transform<Scalar, BoxType::AmbientDimAtCompileTime, Isometry> IsometryTransform;
  typedef Transform<Scalar, BoxType::AmbientDimAtCompileTime, Affine> AffineTransform;

  const VectorType Zero = VectorType::Zero();
  const VectorType Ones = VectorType::Ones();
  const VectorType UnitX = VectorType::UnitX();
  const VectorType UnitY = VectorType::UnitY();
  // this is vector (0, 0, -1, -1, -1, ...), i.e. with zeros at first and second dimensions
  const VectorType UnitZ = Ones - UnitX - UnitY;

  // in this kind of comments the 3D case values will be illustrated
  // box((-1, -1, -1), (1, 1, 1))
  BoxType a(-Ones, Ones);

  // to allow templating this test for both 2D and 3D cases, we always set all
  // but the first coordinate to the same value; so basically 3D case works as
  // if you were looking at the scene from top

  VectorType minPoint = -2 * Ones;
  minPoint[0] = -3;
  VectorType maxPoint = Zero;
  maxPoint[0] = -1;
  BoxType c(minPoint, maxPoint);
  // box((-3, -2, -2), (-1, 0, 0))

  IsometryTransform tf2 = IsometryTransform::Identity();
  // for some weird reason the following statement has to be put separate from
  // the following rotate call, otherwise precision problems arise...
  Rotation rot = rotate(NonInteger(EIGEN_PI));
  tf2.rotate(rot);

  c.transform(tf2);
  // rotate by 180 deg around origin -> box((1, 0, -2), (3, 2, 0))

  VERIFY_IS_APPROX(c.sizes(), a.sizes());
  VERIFY_IS_APPROX((c.min)(), UnitX - UnitZ * Scalar(2));
  VERIFY_IS_APPROX((c.max)(), UnitX * Scalar(3) + UnitY * Scalar(2));

  rot = rotate(NonInteger(EIGEN_PI / 2));
  tf2.setIdentity();
  tf2.rotate(rot);

  c.transform(tf2);
  // rotate by 90 deg around origin ->  box((-2, 1, -2), (0, 3, 0))

  VERIFY_IS_APPROX(c.sizes(), a.sizes());
  VERIFY_IS_APPROX((c.min)(), Ones * Scalar(-2) + UnitY * Scalar(3));
  VERIFY_IS_APPROX((c.max)(), UnitY * Scalar(3));

  // box((-1, -1, -1), (1, 1, 1))
  AffineTransform atf = AffineTransform::Identity();
  atf.linearExt()(0, 1) = Scalar(1);
  c = BoxType(-Ones, Ones);
  c.transform(atf);
  // 45 deg shear in x direction -> box((-2, -1, -1), (2, 1, 1))

  VERIFY_IS_APPROX(c.sizes(), Ones * Scalar(2) + UnitX * Scalar(2));
  VERIFY_IS_APPROX((c.min)(), -Ones - UnitX);
  VERIFY_IS_APPROX((c.max)(), Ones + UnitX);
}

template<typename BoxType, typename Rotation> void alignedboxNonIntegralRotatable(
    const BoxType& box,
    Rotation (*rotate)(typename NumTraits<typename BoxType::Scalar>::NonInteger /*_angle*/))
{
  alignedboxRotatable(box, rotate);

  typedef typename BoxType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::NonInteger NonInteger;
  enum { Dim = BoxType::AmbientDimAtCompileTime };
  typedef Matrix<Scalar, Dim, 1> VectorType;
  typedef Matrix<Scalar, Dim, (1 << Dim)> CornersType;
  typedef Transform<Scalar, Dim, Isometry> IsometryTransform;
  typedef Transform<Scalar, Dim, Affine> AffineTransform;

  const Index dim = box.dim();
  const VectorType Zero = VectorType::Zero();
  const VectorType Ones = VectorType::Ones();

  VectorType minPoint = -2 * Ones;
  minPoint[1] = 1;
  VectorType maxPoint = Zero;
  maxPoint[1] = 3;
  BoxType c(minPoint, maxPoint);
  // ((-2, 1, -2), (0, 3, 0))

  VectorType cornerBL = (c.min)();
  VectorType cornerTR = (c.max)();
  VectorType cornerBR = (c.min)(); cornerBR[0] = cornerTR[0];
  VectorType cornerTL = (c.max)(); cornerTL[0] = cornerBL[0];

  NonInteger angle = NonInteger(EIGEN_PI/3);
  Rotation rot = rotate(angle);
  IsometryTransform tf2;
  tf2.setIdentity();
  tf2.rotate(rot);

  c.transform(tf2);
  // rotate by 60 deg ->  box((-3.59, -1.23, -2), (-0.86, 1.5, 0))

  cornerBL = tf2 * cornerBL;
  cornerBR = tf2 * cornerBR;
  cornerTL = tf2 * cornerTL;
  cornerTR = tf2 * cornerTR;

  VectorType minCorner = Ones * Scalar(-2);
  VectorType maxCorner = Zero;
  minCorner[0] = (min)((min)(cornerBL[0], cornerBR[0]), (min)(cornerTL[0], cornerTR[0]));
  maxCorner[0] = (max)((max)(cornerBL[0], cornerBR[0]), (max)(cornerTL[0], cornerTR[0]));
  minCorner[1] = (min)((min)(cornerBL[1], cornerBR[1]), (min)(cornerTL[1], cornerTR[1]));
  maxCorner[1] = (max)((max)(cornerBL[1], cornerBR[1]), (max)(cornerTL[1], cornerTR[1]));

  for (Index d = 2; d < dim; ++d)
    VERIFY_IS_APPROX(c.sizes()[d], Scalar(2));

  VERIFY_IS_APPROX((c.min)(), minCorner);
  VERIFY_IS_APPROX((c.max)(), maxCorner);

  VectorType minCornerValue = Ones * Scalar(-2);
  VectorType maxCornerValue = Zero;
  minCornerValue[0] = Scalar(Scalar(-sqrt(2*2 + 3*3)) * Scalar(cos(Scalar(atan(2.0/3.0)) - angle/2)));
  minCornerValue[1] = Scalar(Scalar(-sqrt(1*1 + 2*2)) * Scalar(sin(Scalar(atan(2.0/1.0)) - angle/2)));
  maxCornerValue[0] = Scalar(-sin(angle));
  maxCornerValue[1] = Scalar(3 * cos(angle));
  VERIFY_IS_APPROX((c.min)(), minCornerValue);
  VERIFY_IS_APPROX((c.max)(), maxCornerValue);

  // randomized test - translate and rotate the box and compare to a box made of transformed vertices
  for (size_t i = 0; i < 10; ++i)
  {
    for (Index d = 0; d < dim; ++d)
    {
      minCorner[d] = internal::random<Scalar>(-10,10);
      maxCorner[d] = minCorner[d] + internal::random<Scalar>(0, 10);
    }

    c = BoxType(minCorner, maxCorner);

    CornersType corners = boxGetCorners(minCorner, maxCorner);

    typename AffineTransform::LinearMatrixType rotation =
        randomRotationMatrix<typename AffineTransform::LinearMatrixType>();

    tf2.setIdentity();
    tf2.rotate(rotation);
    tf2.translate(VectorType::Random());

    c.transform(tf2);
    corners = tf2 * corners;

    minCorner = corners.rowwise().minCoeff();
    maxCorner = corners.rowwise().maxCoeff();

    VERIFY_IS_APPROX((c.min)(), minCorner);
    VERIFY_IS_APPROX((c.max)(), maxCorner);
  }

  // randomized test - transform the box with a random affine matrix and compare to a box made of transformed vertices
  for (size_t i = 0; i < 10; ++i)
  {
    for (Index d = 0; d < dim; ++d)
    {
      minCorner[d] = internal::random<Scalar>(-10,10);
      maxCorner[d] = minCorner[d] + internal::random<Scalar>(0, 10);
    }

    c = BoxType(minCorner, maxCorner);

    CornersType corners = boxGetCorners(minCorner, maxCorner);

    AffineTransform atf = AffineTransform::Identity();
    atf.linearExt() = AffineTransform::LinearPart::Random();
    atf.translate(VectorType::Random());

    c.transform(atf);
    corners = atf * corners;

    minCorner = corners.rowwise().minCoeff();
    maxCorner = corners.rowwise().maxCoeff();

    VERIFY_IS_APPROX((c.min)(), minCorner);
    VERIFY_IS_APPROX((c.max)(), maxCorner);
  }
}

template<typename BoxType>
void alignedboxCastTests(const BoxType& box)
{
  // casting
  typedef typename BoxType::Scalar Scalar;
  typedef Matrix<Scalar, BoxType::AmbientDimAtCompileTime, 1> VectorType;

  const Index dim = box.dim();

  VectorType p0 = VectorType::Random(dim);
  VectorType p1 = VectorType::Random(dim);

  BoxType b0(dim);

  b0.extend(p0);
  b0.extend(p1);

  const int Dim = BoxType::AmbientDimAtCompileTime;
  typedef typename GetDifferentType<Scalar>::type OtherScalar;
  AlignedBox<OtherScalar,Dim> hp1f = b0.template cast<OtherScalar>();
  VERIFY_IS_APPROX(hp1f.template cast<Scalar>(),b0);
  AlignedBox<Scalar,Dim> hp1d = b0.template cast<Scalar>();
  VERIFY_IS_APPROX(hp1d.template cast<Scalar>(),b0);
}


void specificTest1()
{
    Vector2f m; m << -1.0f, -2.0f;
    Vector2f M; M <<  1.0f,  5.0f;

    typedef AlignedBox2f  BoxType;
    BoxType box( m, M );

    Vector2f sides = M-m;
    VERIFY_IS_APPROX(sides, box.sizes() );
    VERIFY_IS_APPROX(sides[1], box.sizes()[1] );
    VERIFY_IS_APPROX(sides[1], box.sizes().maxCoeff() );
    VERIFY_IS_APPROX(sides[0], box.sizes().minCoeff() );

    VERIFY_IS_APPROX( 14.0f, box.volume() );
    VERIFY_IS_APPROX( 53.0f, box.diagonal().squaredNorm() );
    VERIFY_IS_APPROX( std::sqrt( 53.0f ), box.diagonal().norm() );

    VERIFY_IS_APPROX( m, box.corner( BoxType::BottomLeft ) );
    VERIFY_IS_APPROX( M, box.corner( BoxType::TopRight ) );
    Vector2f bottomRight; bottomRight << M[0], m[1];
    Vector2f topLeft; topLeft << m[0], M[1];
    VERIFY_IS_APPROX( bottomRight, box.corner( BoxType::BottomRight ) );
    VERIFY_IS_APPROX( topLeft, box.corner( BoxType::TopLeft ) );
}


void specificTest2()
{
    Vector3i m; m << -1, -2, 0;
    Vector3i M; M <<  1,  5, 3;

    typedef AlignedBox3i  BoxType;
    BoxType box( m, M );

    Vector3i sides = M-m;
    VERIFY_IS_APPROX(sides, box.sizes() );
    VERIFY_IS_APPROX(sides[1], box.sizes()[1] );
    VERIFY_IS_APPROX(sides[1], box.sizes().maxCoeff() );
    VERIFY_IS_APPROX(sides[0], box.sizes().minCoeff() );

    VERIFY_IS_APPROX( 42, box.volume() );
    VERIFY_IS_APPROX( 62, box.diagonal().squaredNorm() );

    VERIFY_IS_APPROX( m, box.corner( BoxType::BottomLeftFloor ) );
    VERIFY_IS_APPROX( M, box.corner( BoxType::TopRightCeil ) );
    Vector3i bottomRightFloor; bottomRightFloor << M[0], m[1], m[2];
    Vector3i topLeftFloor; topLeftFloor << m[0], M[1], m[2];
    VERIFY_IS_APPROX( bottomRightFloor, box.corner( BoxType::BottomRightFloor ) );
    VERIFY_IS_APPROX( topLeftFloor, box.corner( BoxType::TopLeftFloor ) );
}


EIGEN_DECLARE_TEST(geo_alignedbox)
{
  for(int i = 0; i < g_repeat; i++)
  {
    CALL_SUBTEST_1( (alignedboxNonIntegralRotatable<AlignedBox2f, Rotation2Df>(AlignedBox2f(), &rotate2D)) );
    CALL_SUBTEST_2( alignedboxCastTests(AlignedBox2f()) );

    CALL_SUBTEST_3( (alignedboxNonIntegralRotatable<AlignedBox3f, AngleAxisf>(AlignedBox3f(), &rotate3DZAxis)) );
    CALL_SUBTEST_4( alignedboxCastTests(AlignedBox3f()) );

    CALL_SUBTEST_5( (alignedboxNonIntegralRotatable<AlignedBox4d, Matrix4d>(AlignedBox4d(), &rotate4DZWAxis)) );
    CALL_SUBTEST_6( alignedboxCastTests(AlignedBox4d()) );

    CALL_SUBTEST_7( alignedboxTranslatable(AlignedBox1d()) );
    CALL_SUBTEST_8( alignedboxCastTests(AlignedBox1d()) );

    CALL_SUBTEST_9( alignedboxTranslatable(AlignedBox1i()) );
    CALL_SUBTEST_10( (alignedboxRotatable<AlignedBox2i, Matrix2i>(AlignedBox2i(), &rotate2DIntegral<int, Matrix2i>)) );
    CALL_SUBTEST_11( (alignedboxRotatable<AlignedBox3i, Matrix3i>(AlignedBox3i(), &rotate3DZAxisIntegral<int, Matrix3i>)) );

    CALL_SUBTEST_14( alignedbox(AlignedBox<double,Dynamic>(4)) );
  }
  CALL_SUBTEST_12( specificTest1() );
  CALL_SUBTEST_13( specificTest2() );
}
