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
TridiagEigen<eT>::TridiagEigen()
  : n(0)
  , computed(false)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
TridiagEigen<eT>::TridiagEigen(const Mat<eT>& mat_obj)
  : n(mat_obj.n_rows)
  , computed(false)
  {
  arma_extra_debug_sigprint();
  
  compute(mat_obj);
  }



template<typename eT>
inline
void
TridiagEigen<eT>::compute(const Mat<eT>& mat_obj)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (mat_obj.is_square() == false), "newarp::TridiagEigen::compute(): matrix must be square" );
  
  n = blas_int(mat_obj.n_rows);
  
  main_diag = mat_obj.diag();
  sub_diag  = mat_obj.diag(-1);
  
  evecs.set_size(n, n);
  
  char     compz      = 'I';
  blas_int lwork_min  = 1 + 4*n + n*n;
  blas_int liwork_min = 3 + 5*n;
  blas_int info       = blas_int(0);
  
  blas_int  lwork_proposed = 0;
  blas_int liwork_proposed = 0;
  
  if(n >= 32)
    {
    eT        work_query[2] = {};
    blas_int lwork_query    = blas_int(-1);
    
    blas_int  iwork_query[2] = {};
    blas_int liwork_query    = blas_int(-1);
    
    arma_extra_debug_print("lapack::stedc()");
    lapack::stedc(&compz, &n, main_diag.memptr(), sub_diag.memptr(), evecs.memptr(), &n, &work_query[0], &lwork_query, &iwork_query[0], &liwork_query, &info);
    
    if(info != 0)  { arma_stop_runtime_error("lapack::stedc(): couldn't get size of work arrays"); return; }
    
     lwork_proposed = static_cast<blas_int>( work_query[0] );
    liwork_proposed = iwork_query[0];
    }
  
  blas_int  lwork = (std::max)( lwork_min,  lwork_proposed);
  blas_int liwork = (std::max)(liwork_min, liwork_proposed);
  
  podarray<eT>        work( static_cast<uword>( lwork) );
  podarray<blas_int> iwork( static_cast<uword>(liwork) );
  
  arma_extra_debug_print("lapack::stedc()");
  lapack::stedc(&compz, &n, main_diag.memptr(), sub_diag.memptr(), evecs.memptr(), &n, work.memptr(), &lwork, iwork.memptr(), &liwork, &info);
  
  if(info != 0)  { arma_stop_runtime_error("lapack::stedc(): failed to compute all eigenvalues"); return; }
  
  computed = true;
  }



template<typename eT>
inline
Col<eT>
TridiagEigen<eT>::eigenvalues()
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (computed == false), "newarp::TridiagEigen::eigenvalues(): need to call compute() first" );

  // After calling compute(), main_diag will contain the eigenvalues.
  return main_diag;
  }



template<typename eT>
inline
Mat<eT>
TridiagEigen<eT>::eigenvectors()
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (computed == false), "newarp::TridiagEigen::eigenvectors(): need to call compute() first" );

  return evecs;
  }


}  // namespace newarp
