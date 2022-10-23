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


//! \addtogroup op_symmat
//! @{



template<typename T1>
inline
void
op_symmatu::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_symmatu>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( (A.is_square() == false), "symmatu(): given matrix must be square sized" );
  
  const uword N = A.n_rows;
  
  if(&out != &A)
    {
    out.copy_size(A);
    
    // upper triangular: copy the diagonal and the elements above the diagonal
    
    for(uword i=0; i<N; ++i)
      {
      const eT* A_data   = A.colptr(i);
            eT* out_data = out.colptr(i);
      
      arrayops::copy( out_data, A_data, i+1 );
      }
    }
  
  
  // reflect elements across the diagonal from upper triangle to lower triangle
  
  for(uword col=1; col < N; ++col)
    {
    const eT* coldata = out.colptr(col);
          eT* row_ptr = &(out.at(col,0));
    
    for(uword row=0; row < col; ++row)
      {
      // out.at(col,row) = coldata[row];
      
      (*row_ptr) = coldata[row];
      row_ptr += N;
      }
    }
  }



template<typename T1>
inline
void
op_symmatl::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_symmatl>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( (A.is_square() == false), "symmatl(): given matrix must be square sized" );
  
  const uword N = A.n_rows;
  
  if(&out != &A)
    {
    out.copy_size(A);
    
    // lower triangular: copy the diagonal and the elements below the diagonal
    
    for(uword i=0; i<N; ++i)
      {
      const eT* A_data   = A.colptr(i);
            eT* out_data = out.colptr(i);
      
      arrayops::copy( &out_data[i], &A_data[i], N-i );
      }
    }
  
  
  // reflect elements across the diagonal from lower triangle to upper triangle
  
  for(uword col=0; col < N; ++col)
    {
    const eT* coldata = out.colptr(col);
          eT* row_ptr = &(out.at(col,col+1));
    
    for(uword row=(col+1); row < N; ++row)
      {
      //out.at(col,row) = coldata[row];
      
      (*row_ptr) = coldata[row];
      row_ptr += N;
      }
    }
  }



//



template<typename T1>
inline
void
op_symmatu_cx::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_symmatu_cx>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( (A.is_square() == false), "symmatu(): given matrix must be square sized" );
  
  const uword N  = A.n_rows;
  
  const bool do_conj = (in.aux_uword_b == 1);
  
  if(&out != &A)
    {
    out.copy_size(A);
    
    // upper triangular: copy the diagonal and the elements above the diagonal
    
    for(uword i=0; i<N; ++i)
      {
      const eT* A_data   = A.colptr(i);
            eT* out_data = out.colptr(i);
      
      arrayops::copy( out_data, A_data, i+1 );
      }
    }
  
  
  if(do_conj)
    {
    // reflect elements across the diagonal from upper triangle to lower triangle
    
    for(uword col=1; col < N; ++col)
      {
      const eT* coldata = out.colptr(col);
            eT* row_ptr = &(out.at(col,0));
      
      for(uword row=0; row < col; ++row)
        {
        // out.at(col,row) = std::conj(coldata[row]);
        
        (*row_ptr) = std::conj(coldata[row]);
        row_ptr += N;
        }
      }
    }
  else  // don't do complex conjugation
    {
    // reflect elements across the diagonal from upper triangle to lower triangle
    
    for(uword col=1; col < N; ++col)
      {
      const eT* coldata = out.colptr(col);
            eT* row_ptr = &(out.at(col,0));
      
      for(uword row=0; row < col; ++row)
        {
        // out.at(col,row) = coldata[row];
        
        (*row_ptr) = coldata[row];
        row_ptr += N;
        }
      }
    }
  }



template<typename T1>
inline
void
op_symmatl_cx::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_symmatl_cx>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( (A.is_square() == false), "symmatl(): given matrix must be square sized" );
  
  const uword N  = A.n_rows;
  
  const bool do_conj = (in.aux_uword_b == 1);
  
  if(&out != &A)
    {
    out.copy_size(A);
    
    // lower triangular: copy the diagonal and the elements below the diagonal
    
    for(uword i=0; i<N; ++i)
      {
      const eT* A_data   = A.colptr(i);
            eT* out_data = out.colptr(i);
      
      arrayops::copy( &out_data[i], &A_data[i], N-i );
      }
    }
  
  
  if(do_conj)
    {
    // reflect elements across the diagonal from lower triangle to upper triangle
    
    for(uword col=0; col < N; ++col)
      {
      const eT* coldata = out.colptr(col);
            eT* row_ptr = &(out.at(col,col+1));
      
      for(uword row=(col+1); row < N; ++row)
        {
        // out.at(col,row) = std::conj(coldata[row]);
        
        (*row_ptr) = std::conj(coldata[row]);
        row_ptr += N;
        }
      }
    }
  else  // don't do complex conjugation
    {
    // reflect elements across the diagonal from lower triangle to upper triangle
    
    for(uword col=0; col < N; ++col)
      {
      const eT* coldata = out.colptr(col);
            eT* row_ptr = &(out.at(col,col+1));
      
      for(uword row=(col+1); row < N; ++row)
        {
        // out.at(col,row) = coldata[row];
        
        (*row_ptr) = coldata[row];
        row_ptr += N;
        }
      }
    }
  }



//! @}
