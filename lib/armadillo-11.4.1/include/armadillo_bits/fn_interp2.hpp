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


//! \addtogroup fn_interp2
//! @{



template<typename eT>
inline
void
interp2_helper_nearest(const Mat<eT>& XG, const Mat<eT>& ZG, const Mat<eT>& XI, Mat<eT>& ZI, const eT extrap_val, const uword mode)
  {
  arma_extra_debug_sigprint();
  
  const eT XG_min = XG.min();
  const eT XG_max = XG.max();
  
  // mode = 0: interpolate across rows     (eg. expand in vertical   direction)
  // mode = 1: interpolate across columns  (eg. expand in horizontal direction)
  
  if(mode == 0)  { ZI.set_size(XI.n_elem, ZG.n_cols); }
  if(mode == 1)  { ZI.set_size(ZG.n_rows, XI.n_elem); }
  
  const eT* XG_mem = XG.memptr();
  const eT* XI_mem = XI.memptr();
  
  const uword NG = XG.n_elem;
  const uword NI = XI.n_elem;
  
  uword best_j = 0;
  
  for(uword i=0; i<NI; ++i)
    {
    eT best_err = Datum<eT>::inf;
    
    const eT XI_val = XI_mem[i];
    
    if((XI_val < XG_min) || (XI_val > XG_max))
      {
      if(mode == 0)  { ZI.row(i).fill(extrap_val); }
      if(mode == 1)  { ZI.col(i).fill(extrap_val); }
      }
    else
      {
      // XG and XI are guaranteed to be sorted in ascending manner,
      // so start searching XG from last known optimum position 
      
      for(uword j=best_j; j<NG; ++j)
        {
        const eT tmp = XG_mem[j] - XI_val;
        const eT err = (tmp >= eT(0)) ? tmp : -tmp;
        
        if(err >= best_err)
          {
          // error is going up, so we have found the optimum position
          break;
          }
        else
          {
          best_err = err;
          best_j   = j;   // remember the optimum position
          }
        }
      
      if(mode == 0)  { ZI.row(i) = ZG.row(best_j); }
      if(mode == 1)  { ZI.col(i) = ZG.col(best_j); }
      }
    }
  }



template<typename eT>
inline
void
interp2_helper_linear(const Mat<eT>& XG, const Mat<eT>& ZG, const Mat<eT>& XI, Mat<eT>& ZI, const eT extrap_val, const uword mode)
  {
  arma_extra_debug_sigprint();
  
  const eT XG_min = XG.min();
  const eT XG_max = XG.max();
  
  // mode = 0: interpolate across rows     (eg. expand in vertical   direction)
  // mode = 1: interpolate across columns  (eg. expand in horizontal direction)
  
  if(mode == 0)  { ZI.set_size(XI.n_elem, ZG.n_cols); }
  if(mode == 1)  { ZI.set_size(ZG.n_rows, XI.n_elem); }
  
  const eT* XG_mem = XG.memptr();
  const eT* XI_mem = XI.memptr();
  
  const uword NG = XG.n_elem;
  const uword NI = XI.n_elem;
  
  uword a_best_j = 0;
  uword b_best_j = 0;
  
  for(uword i=0; i<NI; ++i)
    {
    const eT XI_val = XI_mem[i];
    
    if((XI_val < XG_min) || (XI_val > XG_max))
      {
      if(mode == 0)  { ZI.row(i).fill(extrap_val); }
      if(mode == 1)  { ZI.col(i).fill(extrap_val); }
      }
    else
      {
      // XG and XI are guaranteed to be sorted in ascending manner,
      // so start searching XG from last known optimum position 
      
      eT a_best_err = Datum<eT>::inf;
      eT b_best_err = Datum<eT>::inf;
      
      for(uword j=a_best_j; j<NG; ++j)
        {
        const eT tmp = XG_mem[j] - XI_val;
        const eT err = (tmp >= eT(0)) ? tmp : -tmp;
        
        if(err >= a_best_err)
          {
          break;
          }
        else
          {
          a_best_err = err;
          a_best_j   = j;
          }
        }
      
      if( (XG_mem[a_best_j] - XI_val) <= eT(0) )
        {
        // a_best_j is to the left of the interpolated position
        
        b_best_j = ( (a_best_j+1) < NG) ? (a_best_j+1) : a_best_j; 
        }
      else
        {
        // a_best_j is to the right of the interpolated position
        
        b_best_j = (a_best_j >= 1) ? (a_best_j-1) : a_best_j; 
        }
      
      b_best_err = std::abs( XG_mem[b_best_j] - XI_val );
      
      if(a_best_j > b_best_j)
        {
        std::swap(a_best_j,   b_best_j  );
        std::swap(a_best_err, b_best_err);
        }
      
      const eT weight = (a_best_err > eT(0)) ? (a_best_err / (a_best_err + b_best_err)) : eT(0);
      
      if(mode == 0)  { ZI.row(i) = (eT(1) - weight)*ZG.row(a_best_j) + (weight)*ZG.row(b_best_j); }
      if(mode == 1)  { ZI.col(i) = (eT(1) - weight)*ZG.col(a_best_j) + (weight)*ZG.col(b_best_j); }
      }
    }
  }



template<typename T1, typename T2, typename T3, typename T4, typename T5>
inline
typename
enable_if2< is_real<typename T1::elem_type>::value, void >::result
interp2
  (
  const Base<typename T1::elem_type, T1>& X,
  const Base<typename T1::elem_type, T2>& Y,
  const Base<typename T1::elem_type, T3>& Z,
  const Base<typename T1::elem_type, T4>& XI,
  const Base<typename T1::elem_type, T5>& YI,
         Mat<typename T1::elem_type>&     ZI,
  const char*                             method     = "linear",
  const typename T1::elem_type            extrap_val = Datum<typename T1::elem_type>::nan
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const char sig = (method != nullptr) ? method[0] : char(0);
  
  arma_debug_check( ((sig != 'n') && (sig != 'l')), "interp2(): unsupported interpolation type" ); 
  
  const quasi_unwrap<T1> UXG(  X.get_ref() );
  const quasi_unwrap<T2> UYG(  Y.get_ref() );
  const quasi_unwrap<T3> UZG(  Z.get_ref() );
  const quasi_unwrap<T4> UXI( XI.get_ref() );
  const quasi_unwrap<T5> UYI( YI.get_ref() );
  
  arma_debug_check( (UXG.M.is_vec() == false), "interp2(): X must resolve to a vector" );
  arma_debug_check( (UYG.M.is_vec() == false), "interp2(): Y must resolve to a vector" );
  
  arma_debug_check( (UXI.M.is_vec() == false), "interp2(): XI must resolve to a vector" );
  arma_debug_check( (UYI.M.is_vec() == false), "interp2(): YI must resolve to a vector" );
  
  arma_debug_check( (UXG.M.n_elem < 2), "interp2(): X must have at least two unique elements" );
  arma_debug_check( (UYG.M.n_elem < 2), "interp2(): Y must have at least two unique elements" );
  
  arma_debug_check( (UXG.M.n_elem != UZG.M.n_cols), "interp2(): number of elements in X must equal the number of columns in Z" );
  arma_debug_check( (UYG.M.n_elem != UZG.M.n_rows), "interp2(): number of elements in Y must equal the number of rows in Z"    );
  
  arma_debug_check( (UXG.M.is_sorted("strictascend") == false), "interp2(): X must be monotonically increasing" );
  arma_debug_check( (UYG.M.is_sorted("strictascend") == false), "interp2(): Y must be monotonically increasing" );
  
  arma_debug_check( (UXI.M.is_sorted("strictascend") == false), "interp2(): XI must be monotonically increasing" );
  arma_debug_check( (UYI.M.is_sorted("strictascend") == false), "interp2(): YI must be monotonically increasing" );
  
  Mat<eT> tmp;
  
  if( UXG.is_alias(ZI) || UXI.is_alias(ZI) )
    {
    Mat<eT> out;
    
    if(sig == 'n')
      {
      interp2_helper_nearest(UYG.M, UZG.M, UYI.M, tmp, extrap_val, 0);
      interp2_helper_nearest(UXG.M, tmp,   UXI.M, out, extrap_val, 1);
      }
    else
    if(sig == 'l')
      {
      interp2_helper_linear(UYG.M, UZG.M, UYI.M, tmp, extrap_val, 0);
      interp2_helper_linear(UXG.M, tmp,   UXI.M, out, extrap_val, 1);
      }
    
    ZI.steal_mem(out);
    }
  else
    {
    if(sig == 'n')
      {
      interp2_helper_nearest(UYG.M, UZG.M, UYI.M, tmp, extrap_val, 0);
      interp2_helper_nearest(UXG.M, tmp,   UXI.M,  ZI, extrap_val, 1);
      }
    else
    if(sig == 'l')
      {
      interp2_helper_linear(UYG.M, UZG.M, UYI.M, tmp, extrap_val, 0);
      interp2_helper_linear(UXG.M, tmp,   UXI.M,  ZI, extrap_val, 1);
      }
    }
  }



//! @}
