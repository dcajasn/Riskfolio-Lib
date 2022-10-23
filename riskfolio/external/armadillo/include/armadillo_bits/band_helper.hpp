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


//! \addtogroup band_helper
//! @{


namespace band_helper
{



template<typename eT>
inline
bool
is_band(uword& out_KL, uword& out_KU, const Mat<eT>& A, const uword N_min)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: assuming that A has a square size
  // NOTE: assuming that N_min is >= 4
  
  const uword N = A.n_rows;
  
  if(N < N_min)  { return false; }
  
  // first, quickly check bottom-left and top-right corners
  
  const eT eT_zero = eT(0);
  
  const eT* A_col0 = A.memptr();
  const eT* A_col1 = A_col0 + N;
  
  if( (A_col0[N-2] != eT_zero) || (A_col0[N-1] != eT_zero) || (A_col1[N-2] != eT_zero) || (A_col1[N-1] != eT_zero) )  { return false; }
  
  const eT* A_colNm2 = A.colptr(N-2);
  const eT* A_colNm1 = A_colNm2 + N;
  
  if( (A_colNm2[0] != eT_zero) || (A_colNm2[1] != eT_zero) || (A_colNm1[0] != eT_zero) || (A_colNm1[1] != eT_zero) )  { return false; }
  
  // if we reached this point, go through the entire matrix to work out number of subdiagonals and superdiagonals
  
  const uword n_nonzero_threshold = (N*N)/4;  // empirically determined
  
  uword KL = 0;  // number of   subdiagonals
  uword KU = 0;  // number of superdiagonals
  
  const eT* A_colptr = A.memptr();
  
  for(uword col=0; col < N; ++col)
    {
    uword first_nonzero_row = col;
    uword  last_nonzero_row = col;
    
    for(uword row=0; row < col; ++row)
      {
      if( A_colptr[row] != eT_zero )  { first_nonzero_row = row; break; }
      }
    
    for(uword row=(col+1); row < N; ++row)
      {
      last_nonzero_row = (A_colptr[row] != eT_zero) ? row : last_nonzero_row;
      }
    
    const uword L_count = last_nonzero_row - col;
    const uword U_count = col - first_nonzero_row;
    
    if( (L_count > KL) || (U_count > KU) )
      {
      KL = (std::max)(KL, L_count);
      KU = (std::max)(KU, U_count);
      
      const uword n_nonzero = N*(KL+KU+1) - (KL*(KL+1) + KU*(KU+1))/2;
      
      // return as soon as we know that it's not worth analysing the matrix any further
      
      if(n_nonzero > n_nonzero_threshold)  { return false; }
      }
    
    A_colptr += N;
    }
  
  out_KL = KL;
  out_KU = KU;
  
  return true;
  }



template<typename eT>
inline
bool
is_band_lower(uword& out_KD, const Mat<eT>& A, const uword N_min)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: assuming that A has a square size
  // NOTE: assuming that N_min is >= 4
  
  const uword N = A.n_rows;
  
  if(N < N_min)  { return false; }
  
  // first, quickly check bottom-left corner
  
  const eT eT_zero = eT(0);
  
  const eT* A_col0 = A.memptr();
  const eT* A_col1 = A_col0 + N;
  
  if( (A_col0[N-2] != eT_zero) || (A_col0[N-1] != eT_zero) || (A_col1[N-2] != eT_zero) || (A_col1[N-1] != eT_zero) )  { return false; }
  
  // if we reached this point, go through the bottom triangle to work out number of subdiagonals
  
  const uword n_nonzero_threshold = ( N*N - (N*(N-1))/2 ) / 4;  // empirically determined
  
  uword KL = 0;  // number of subdiagonals
  
  const eT* A_colptr = A.memptr();
  
  for(uword col=0; col < N; ++col)
    {
    uword last_nonzero_row = col;
    
    for(uword row=(col+1); row < N; ++row)
      {
      last_nonzero_row = (A_colptr[row] != eT_zero) ? row : last_nonzero_row;
      }
    
    const uword L_count = last_nonzero_row - col;
    
    if(L_count > KL)
      {
      KL = L_count;
      
      const uword n_nonzero = N*(KL+1) - (KL*(KL+1))/2;
      
      // return as soon as we know that it's not worth analysing the matrix any further
      
      if(n_nonzero > n_nonzero_threshold)  { return false; }
      }
    
    A_colptr += N;
    }
  
  out_KD = KL;
  
  return true;
  }



template<typename eT>
inline
bool
is_band_upper(uword& out_KD, const Mat<eT>& A, const uword N_min)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: assuming that A has a square size
  // NOTE: assuming that N_min is >= 4
  
  const uword N = A.n_rows;
  
  if(N < N_min)  { return false; }
  
  // first, quickly check top-right corner
  
  const eT eT_zero = eT(0);
  
  const eT* A_colNm2 = A.colptr(N-2);
  const eT* A_colNm1 = A_colNm2 + N;
  
  if( (A_colNm2[0] != eT_zero) || (A_colNm2[1] != eT_zero) || (A_colNm1[0] != eT_zero) || (A_colNm1[1] != eT_zero) )  { return false; }
  
  // if we reached this point, go through the entire matrix to work out number of superdiagonals
  
  const uword n_nonzero_threshold = ( N*N - (N*(N-1))/2 ) / 4;  // empirically determined
  
  uword KU = 0;  // number of superdiagonals
  
  const eT* A_colptr = A.memptr();
  
  for(uword col=0; col < N; ++col)
    {
    uword first_nonzero_row = col;
    
    for(uword row=0; row < col; ++row)
      {
      if( A_colptr[row] != eT_zero )  { first_nonzero_row = row; break; }
      }
    
    const uword U_count = col - first_nonzero_row;
    
    if(U_count > KU)
      {
      KU = U_count;
      
      const uword n_nonzero = N*(KU+1) - (KU*(KU+1))/2;
      
      // return as soon as we know that it's not worth analysing the matrix any further
      
      if(n_nonzero > n_nonzero_threshold)  { return false; }
      }
    
    A_colptr += N;
    }
  
  out_KD = KU;
  
  return true;
  }



template<typename eT>
inline
void
compress(Mat<eT>& AB, const Mat<eT>& A, const uword KL, const uword KU, const bool use_offset)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: assuming that A has a square size
  
  // band matrix storage format
  // http://www.netlib.org/lapack/lug/node124.html  
  
  // for ?gbsv,  matrix AB size: 2*KL+KU+1 x N; band representation of A stored in rows KL+1 to 2*KL+KU+1  (note: fortran counts from 1)
  // for ?gbsvx, matrix AB size:   KL+KU+1 x N; band representaiton of A stored in rows    1 to   KL+KU+1  (note: fortran counts from 1)
  //
  // the +1 in the above formulas is to take into account the main diagonal
  
  const uword AB_n_rows = (use_offset) ? uword(2*KL + KU + 1) : uword(KL + KU + 1);
  const uword N         = A.n_rows;
  
  AB.set_size(AB_n_rows, N);
  
  if(A.is_empty())  { AB.zeros(); return; }
  
  if(AB_n_rows == uword(1))
    {
    eT* AB_mem = AB.memptr();
    
    for(uword i=0; i<N; ++i)  { AB_mem[i] = A.at(i,i); }
    }
  else
    {
    AB.zeros();  // paranoia
    
    for(uword j=0; j < N; ++j)
      {
      const uword A_row_start = (j > KU) ? uword(j - KU) : uword(0);
      const uword A_row_endp1 = (std::min)(N, j+KL+1);
      
      const uword length = A_row_endp1 - A_row_start;
      
      const uword AB_row_start = (KU > j) ? (KU - j) : uword(0);
      
      const eT*  A_colptr =  A.colptr(j) +  A_row_start;
            eT* AB_colptr = AB.colptr(j) + AB_row_start + ( (use_offset) ? KL : uword(0) );
      
      arrayops::copy( AB_colptr, A_colptr, length );
      }
    }
  }



template<typename eT>
inline
void
uncompress(Mat<eT>& A, const Mat<eT>& AB, const uword KL, const uword KU, const bool use_offset)
  {
  arma_extra_debug_sigprint();
  
  const uword AB_n_rows = AB.n_rows;
  const uword N         = AB.n_cols;
  
  arma_debug_check( (AB_n_rows != ((use_offset) ? uword(2*KL + KU + 1) : uword(KL + KU + 1))), "band_helper::uncompress(): detected inconsistency" );
  
  A.zeros(N,N);  // assuming there is no aliasing between A and AB
  
  if(AB_n_rows == uword(1))
    {
    const eT* AB_mem = AB.memptr();
    
    for(uword i=0; i<N; ++i)  { A.at(i,i) = AB_mem[i]; }
    }
  else
    {
    for(uword j=0; j < N; ++j)
      {
      const uword A_row_start = (j > KU) ? uword(j - KU) : uword(0);
      const uword A_row_endp1 = (std::min)(N, j+KL+1);
      
      const uword length = A_row_endp1 - A_row_start;
      
      const uword AB_row_start = (KU > j) ? (KU - j) : uword(0);
      
      const eT* AB_colptr = AB.colptr(j) + AB_row_start + ( (use_offset) ? KL : uword(0) );
            eT*  A_colptr =  A.colptr(j) +  A_row_start;
      
      arrayops::copy( A_colptr, AB_colptr, length );
      }
    }
  }



template<typename eT>
inline
void
extract_tridiag(Mat<eT>& out, const Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: assuming that A has a square size and is at least 2x2
  
  const uword N = A.n_rows;
  
  out.set_size(N, 3);  // assuming there is no aliasing between 'out' and 'A'
  
  if(N < 2)  { return; }
  
  eT* DL = out.colptr(0);
  eT* DD = out.colptr(1);
  eT* DU = out.colptr(2);
  
  DD[0] = A[0];
  DL[0] = A[1];
  
  const uword Nm1 = N-1;
  const uword Nm2 = N-2;
  
  for(uword i=0; i < Nm2; ++i)
    {
    const uword ip1 = i+1;
    
    const eT* data = &(A.at(i, ip1));
    
    const eT tmp0 = data[0];
    const eT tmp1 = data[1];
    const eT tmp2 = data[2];
    
    DL[ip1] = tmp2;
    DD[ip1] = tmp1;
    DU[i  ] = tmp0;
    }
  
  const eT* data = &(A.at(Nm2, Nm1));
  
  DL[Nm1] = 0;
  DU[Nm2] = data[0];
  DU[Nm1] = 0;
  DD[Nm1] = data[1]; 
  }



}  // end of namespace band_helper


//! @}
