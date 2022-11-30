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



//! \addtogroup op_cor
//! @{



template<typename T1>
inline
void
op_cor::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cor>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword norm_type = in.aux_uword_a;
  
  const unwrap<T1>   U(in.m);
  const Mat<eT>& A = U.M;
  
  if(A.n_elem == 0)
    {
    out.reset();
    return;
    }
  
  if(A.n_elem == 1)
    {
    out.set_size(1,1);
    out[0] = eT(1);
    return;
    }
  
  const Mat<eT>& AA = (A.n_rows == 1) 
                      ? Mat<eT>(const_cast<eT*>(A.memptr()), A.n_cols, A.n_rows, false, false)
                      : Mat<eT>(const_cast<eT*>(A.memptr()), A.n_rows, A.n_cols, false, false);
  
  const uword N        = AA.n_rows;
  const eT    norm_val = (norm_type == 0) ? ( (N > 1) ? eT(N-1) : eT(1) ) : eT(N);
  
  const Mat<eT> tmp = AA.each_row() - mean(AA,0);
  
  out = tmp.t() * tmp;
  out /= norm_val;
  
  const Col<eT> s = sqrt(out.diag());
  
  out /= (s * s.t());  // TODO: check for zeros in s?
  }



template<typename T1>
inline
void
op_cor::apply(Mat<typename T1::elem_type>& out, const Op< Op<T1,op_htrans>, op_cor>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword norm_type = in.aux_uword_a;
  
  if(is_cx<eT>::yes)
    {
    const Mat<eT> tmp = in.m;  // force the evaluation of Op<T1,op_htrans>
    
    out = cor(tmp, norm_type);
    }
  else
    {
    const unwrap<T1>   U(in.m.m);
    const Mat<eT>& A = U.M;
    
    if(A.n_elem == 0)
      {
      out.reset();
      return;
      }
    
    if(A.n_elem == 1)
      {
      out.set_size(1,1);
      out[0] = eT(1);
      return;
      }
    
    const Mat<eT>& AA = (A.n_cols == 1)
                        ? Mat<eT>(const_cast<eT*>(A.memptr()), A.n_cols, A.n_rows, false, false)
                        : Mat<eT>(const_cast<eT*>(A.memptr()), A.n_rows, A.n_cols, false, false);
    
    const uword N        = AA.n_cols;
    const eT    norm_val = (norm_type == 0) ? ( (N > 1) ? eT(N-1) : eT(1) ) : eT(N);
    
    const Mat<eT> tmp = AA.each_col() - mean(AA,1);
    
    out = tmp * tmp.t();
    out /= norm_val;
    
    const Col<eT> s = sqrt(out.diag());
    
    out /= (s * s.t());  // TODO: check for zeros in s?
    }
  }



//! @}
