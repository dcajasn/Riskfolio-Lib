// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2011 Gael Guennebaud <g.gael@free.fr>
// Copyright (C) 2012 Kolja Brix <brix@igpm.rwth-aaachen.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "../../test/sparse_solver.h"
#include <Eigen/IterativeSolvers>

template<typename T> void test_idrs_T()
{
  IDRS<SparseMatrix<T>, DiagonalPreconditioner<T> > idrs_colmajor_diag;
  IDRS<SparseMatrix<T>, IncompleteLUT<T> >           idrs_colmajor_ilut;

  CALL_SUBTEST( check_sparse_square_solving(idrs_colmajor_diag)  );
  CALL_SUBTEST( check_sparse_square_solving(idrs_colmajor_ilut)     );
}

EIGEN_DECLARE_TEST(idrs)
{
  CALL_SUBTEST_1(test_idrs_T<double>());
  CALL_SUBTEST_2(test_idrs_T<std::complex<double> >());
}
