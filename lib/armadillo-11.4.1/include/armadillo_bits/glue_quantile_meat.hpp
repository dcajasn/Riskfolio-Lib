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


//! \addtogroup glue_quantile
//! @{


template<typename eTa, typename eTb>
inline
void
glue_quantile::worker(eTb* out_mem, Col<eTa>& Y, const Mat<eTb>& P)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: assuming out_mem is an array with P.n_elem elements
  
  // TODO: ignore non-finite values ?
  
  // algorithm based on "Definition 5" in:
  // Rob J. Hyndman and Yanan Fan.
  // Sample Quantiles in Statistical Packages.
  // The American Statistician, Vol. 50, No. 4, pp. 361-365, 1996.
  // http://doi.org/10.2307/2684934
  
  const eTb*  P_mem    = P.memptr();
  const uword P_n_elem = P.n_elem;
  
  const eTb alpha = 0.5;
  const eTb N     = eTb(Y.n_elem);
  const eTb P_min = (eTb(1) - alpha) / N;
  const eTb P_max = (N      - alpha) / N;
  
  for(uword i=0; i < P_n_elem; ++i)
    {
    const eTb P_i = P_mem[i];
    
    eTb out_val = eTb(0);
    
    if(P_i < P_min)
      {
      out_val = (P_i < eTb(0)) ? eTb(-std::numeric_limits<eTb>::infinity()) : eTb(Y.min());
      }
    else
    if(P_i > P_max)
      {
      out_val = (P_i > eTb(1)) ? eTb( std::numeric_limits<eTb>::infinity()) : eTb(Y.max());
      }
    else
      {
      const uword   k = uword(std::floor(N * P_i + alpha));
      const eTb   P_k = (eTb(k) - alpha) / N;
      
      const eTb w = (P_i - P_k) * N;
      
      eTa* Y_k_ptr = Y.begin() + uword(k);
      std::nth_element( Y.begin(), Y_k_ptr, Y.end() );
      const eTa Y_k_val = (*Y_k_ptr);
      
      eTa* Y_km1_ptr = Y.begin() + uword(k-1);
      // std::nth_element( Y.begin(), Y_km1_ptr, Y.end() );
      std::nth_element( Y.begin(), Y_km1_ptr, Y_k_ptr );
      const eTa Y_km1_val = (*Y_km1_ptr);
      
      out_val = ((eTb(1) - w) * Y_km1_val) + (w * Y_k_val);
      }
    
    out_mem[i] = out_val;
    }
  }



template<typename eTa, typename eTb>
inline
void
glue_quantile::apply_noalias(Mat<eTb>& out, const Mat<eTa>& X, const Mat<eTb>& P, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ((P.is_vec() == false) && (P.is_empty() == false)), "quantile(): parameter 'P' must be a vector" );
  
  if(X.is_empty())  { out.reset(); return; }
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword P_n_elem = P.n_elem;
  
  if(dim == 0)
    {
    out.set_size(P_n_elem, X_n_cols);
    
    if(out.is_empty())  { return; }
    
    Col<eTa> Y(X_n_rows, arma_nozeros_indicator());
    
    if(X_n_cols == 1)
      {
      arrayops::copy(Y.memptr(), X.memptr(), X_n_rows);
      
      glue_quantile::worker(out.memptr(), Y, P);
      }
    else
      {
      for(uword col=0; col < X_n_cols; ++col)
        {
        arrayops::copy(Y.memptr(), X.colptr(col), X_n_rows);
        
        glue_quantile::worker(out.colptr(col), Y, P);
        }
      }
    }
  else
  if(dim == 1)
    {
    out.set_size(X_n_rows, P_n_elem);
    
    if(out.is_empty())  { return; }
    
    Col<eTa> Y(X_n_cols, arma_nozeros_indicator());
    
    if(X_n_rows == 1)
      {
      arrayops::copy(Y.memptr(), X.memptr(), X_n_cols);
      
      glue_quantile::worker(out.memptr(), Y, P);
      }
    else
      {
      Col<eTb> tmp(P_n_elem, arma_nozeros_indicator());
      
      eTb* tmp_mem = tmp.memptr();
      
      for(uword row=0; row < X_n_rows; ++row)
        {
        eTa* Y_mem = Y.memptr();
        
        for(uword col=0; col < X_n_cols; ++col)  { Y_mem[col] = X.at(row,col); }
        
        glue_quantile::worker(tmp_mem, Y, P);
        
        for(uword i=0; i < P_n_elem; ++i)  { out.at(row,i) = tmp_mem[i]; }
        }
      }
    }
  }



template<typename T1, typename T2>
inline
void
glue_quantile::apply(Mat<typename T2::elem_type>& out, const mtGlue<typename T2::elem_type,T1,T2,glue_quantile>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T2::elem_type eTb;
  
  const uword dim = expr.aux_uword;
  
  arma_debug_check( (dim > 1), "quantile(): parameter 'dim' must be 0 or 1" );
  
  const quasi_unwrap<T1> UA(expr.A);
  const quasi_unwrap<T2> UB(expr.B);
  
  arma_debug_check((UA.M.has_nan() || UB.M.has_nan()), "quantile(): detected NaN");
  
  if(UA.is_alias(out) || UB.is_alias(out))
    {
    Mat<eTb> tmp;
    
    glue_quantile::apply_noalias(tmp, UA.M, UB.M, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    glue_quantile::apply_noalias(out, UA.M, UB.M, dim);
    }
  }



template<typename T1, typename T2>
inline
void
glue_quantile_default::apply(Mat<typename T2::elem_type>& out, const mtGlue<typename T2::elem_type,T1,T2,glue_quantile_default>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T2::elem_type eTb;
  
  const quasi_unwrap<T1> UA(expr.A);
  const quasi_unwrap<T2> UB(expr.B);
  
  const uword dim = (T1::is_xvec) ? uword(UA.M.is_rowvec() ? 1 : 0) : uword((T1::is_row) ? 1 : 0);
  
  arma_debug_check((UA.M.has_nan() || UB.M.has_nan()), "quantile(): detected NaN");
  
  if(UA.is_alias(out) || UB.is_alias(out))
    {
    Mat<eTb> tmp;
    
    glue_quantile::apply_noalias(tmp, UA.M, UB.M, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    glue_quantile::apply_noalias(out, UA.M, UB.M, dim);
    }
  }


//! @}
