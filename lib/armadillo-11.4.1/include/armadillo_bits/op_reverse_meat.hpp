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


//! \addtogroup op_reverse
//! @{



template<typename T1>
inline
void
op_reverse::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_reverse>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  
  arma_debug_check( (dim > 1), "reverse(): parameter 'dim' must be 0 or 1" );
  
  if(is_Mat<T1>::value)
    {
    // allow detection of in-place operation
    
    const unwrap<T1> U(in.m);
    
    if(dim == 0)  { op_flipud::apply_direct(out, U.M); }
    if(dim == 1)  { op_fliplr::apply_direct(out, U.M); }
    }
  else
    {
    const Proxy<T1> P(in.m);
    
    if(P.is_alias(out))
      {
      Mat<eT> tmp;
      
      if(dim == 0)  { op_flipud::apply_proxy_noalias(tmp, P); }
      if(dim == 1)  { op_fliplr::apply_proxy_noalias(tmp, P); }
      
      out.steal_mem(tmp);
      }
    else
      {
      if(dim == 0)  { op_flipud::apply_proxy_noalias(out, P); }
      if(dim == 1)  { op_fliplr::apply_proxy_noalias(out, P); }
      }
    }
  }



template<typename T1>
inline
void
op_reverse_vec::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_reverse_vec>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_Mat<T1>::value)
    {
    // allow detection of in-place operation
    
    const unwrap<T1> U(in.m);
    
    if((T1::is_xvec) ? bool(U.M.is_rowvec()) : bool(T1::is_row))
      {
      op_fliplr::apply_direct(out, U.M);
      }
    else
      {
      op_flipud::apply_direct(out, U.M);
      }
    }
  else
    {
    const Proxy<T1> P(in.m);
    
    if(P.is_alias(out))
      {
      Mat<eT> tmp;
      
      if((T1::is_xvec) ? bool(P.get_n_rows() == 1) : bool(T1::is_row))
        {
        op_fliplr::apply_proxy_noalias(tmp, P);
        }
      else
        {
        op_flipud::apply_proxy_noalias(tmp, P);
        }
      
      out.steal_mem(tmp);
      }
    else
      {
      if((T1::is_xvec) ? bool(P.get_n_rows() == 1) : bool(T1::is_row))
        {
        op_fliplr::apply_proxy_noalias(out, P);
        }
      else
        {
        op_flipud::apply_proxy_noalias(out, P);
        }
      }
    }
  }



//! @}
