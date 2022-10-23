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


//! \addtogroup op_chol
//! @{



template<typename T1>
inline
void
op_chol::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_chol>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = op_chol::apply_direct(out, X.m, X.aux_uword_a);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("chol(): decomposition failed");
    }
  }



template<typename T1>
inline
bool
op_chol::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& A_expr, const uword layout)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  out = A_expr.get_ref();
  
  arma_debug_check( (out.is_square() == false), "chol(): given matrix must be square sized", [&](){ out.soft_reset(); } );
  
  if(out.is_empty())  { return true; }
  
  if((arma_config::debug) && (auxlib::rudimentary_sym_check(out) == false))
    {
    if(is_cx<eT>::no )  { arma_debug_warn_level(1, "chol(): given matrix is not symmetric"); }
    if(is_cx<eT>::yes)  { arma_debug_warn_level(1, "chol(): given matrix is not hermitian"); }
    }
  
  uword KD = 0;
  
  const bool is_band = arma_config::optimise_band && ((auxlib::crippled_lapack(out)) ? false : ((layout == 0) ? band_helper::is_band_upper(KD, out, uword(32)) : band_helper::is_band_lower(KD, out, uword(32))));
  
  const bool status = (is_band) ? auxlib::chol_band(out, KD, layout) : auxlib::chol(out, layout);
  
  return status;
  }



//! @}
