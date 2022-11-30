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


namespace newarp
{


template<typename eT>
inline
UpperHessenbergEigen<eT>::UpperHessenbergEigen()
  : n_rows(0)
  , computed(false)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
UpperHessenbergEigen<eT>::UpperHessenbergEigen(const Mat<eT>& mat_obj)
  : n_rows(mat_obj.n_rows)
  , computed(false)
  {
  arma_extra_debug_sigprint();
  
  compute(mat_obj);
  }



template<typename eT>
inline
void
UpperHessenbergEigen<eT>::compute(const Mat<eT>& mat_obj)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (mat_obj.is_square() == false), "newarp::UpperHessenbergEigen::compute(): matrix must be square" );
  
  n_rows = mat_obj.n_rows;
  
  mat_Z.set_size(n_rows, n_rows);
  mat_T.set_size(n_rows, n_rows);
  evals.set_size(n_rows);
  
  mat_Z.eye();
  mat_T = mat_obj;
  
  blas_int want_T = blas_int(1);
  blas_int want_Z = blas_int(1);
  
  blas_int n    = blas_int(n_rows);
  blas_int ilo  = blas_int(1);
  blas_int ihi  = blas_int(n_rows);
  blas_int iloz = blas_int(1);
  blas_int ihiz = blas_int(n_rows);
  
  blas_int info = blas_int(0);
  
  podarray<eT> wr(n_rows);
  podarray<eT> wi(n_rows);
  
  arma_extra_debug_print("lapack::lahqr()");
  lapack::lahqr(&want_T, &want_Z, &n, &ilo, &ihi, mat_T.memptr(), &n, wr.memptr(), wi.memptr(), &iloz, &ihiz, mat_Z.memptr(), &n, &info);
  
  if(info != 0)  { arma_stop_runtime_error("lapack::lahqr(): failed to compute all eigenvalues"); return; }
  
  for(uword i=0; i < n_rows; i++)  { evals(i) = std::complex<eT>(wr[i], wi[i]); }
  
  char     side   = 'R';
  char     howmny = 'B';
  blas_int m      = blas_int(0);
  
  podarray<eT> work(3*n);
  
  arma_extra_debug_print("lapack::trevc()");
  lapack::trevc(&side, &howmny, (blas_int*) NULL, &n, mat_T.memptr(), &n, (eT*) NULL, &n, mat_Z.memptr(), &n, &n, &m, work.memptr(), &info);
  
  if(info != 0)  { arma_stop_runtime_error("lapack::trevc(): illegal value"); return; }
  
  computed = true;
  }



template<typename eT>
inline
Col< std::complex<eT> >
UpperHessenbergEigen<eT>::eigenvalues()
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (computed == false), "newarp::UpperHessenbergEigen::eigenvalues(): need to call compute() first" );
  
  return evals;
  }



template<typename eT>
inline
Mat< std::complex<eT> >
UpperHessenbergEigen<eT>::eigenvectors()
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (computed == false), "newarp::UpperHessenbergEigen::eigenvectors(): need to call compute() first" );
  
  // Lapack will set the imaginary parts of real eigenvalues to be exact zero
  Mat< std::complex<eT> > evecs(n_rows, n_rows, arma_zeros_indicator());
  
  std::complex<eT>* col_ptr = evecs.memptr();
  
  for(uword i=0; i < n_rows; i++)
    {
    if(cx_attrib::is_real(evals(i), eT(0)))
      {
      // for real eigenvector, normalise and copy
      const eT z_norm = norm(mat_Z.col(i));
      
      for(uword j=0; j < n_rows; j++)
        {
        col_ptr[j] = std::complex<eT>(mat_Z(j, i) / z_norm, eT(0));
        }
      
      col_ptr += n_rows;
      }
    else
      {
      // complex eigenvectors are stored in consecutive columns
      const eT r2 = dot(mat_Z.col(i  ), mat_Z.col(i  ));
      const eT i2 = dot(mat_Z.col(i+1), mat_Z.col(i+1));
      
      const eT  z_norm = std::sqrt(r2 + i2);
      const eT* z_ptr  = mat_Z.colptr(i);
      
      for(uword j=0; j < n_rows; j++)
        {
        col_ptr[j         ] = std::complex<eT>(z_ptr[j] / z_norm, z_ptr[j + n_rows] / z_norm);
        col_ptr[j + n_rows] = std::conj(col_ptr[j]);
        }
      
      i++;
      col_ptr += 2 * n_rows;
      }
    }
  
  return evecs;
  }


}  // namespace newarp
