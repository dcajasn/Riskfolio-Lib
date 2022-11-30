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


//! \addtogroup op_roots
//! @{



template<typename T1>
inline
void
op_roots::apply(Mat< std::complex<typename T1::pod_type> >& out, const mtOp<std::complex<typename T1::pod_type>, T1, op_roots>& expr)
  {
  arma_extra_debug_sigprint();
  
  const bool status = op_roots::apply_direct(out, expr.m);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("roots(): eigen decomposition failed");
    }
  }



template<typename T1>
inline
bool
op_roots::apply_direct(Mat< std::complex<typename T1::pod_type> >& out, const Base<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<typename T1::pod_type> out_eT;
  
  const quasi_unwrap<T1> U(X.get_ref());
  
  bool status = false;
  
  if(U.is_alias(out))
    {
    Mat<out_eT> tmp;
    
    status = op_roots::apply_noalias(tmp, U.M);
    
    out.steal_mem(tmp);
    }
  else
    {
    status = op_roots::apply_noalias(out, U.M);
    }
  
  return status;
  }



template<typename eT>
inline
bool
op_roots::apply_noalias(Mat< std::complex<typename get_pod_type<eT>::result> >& out, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
  typedef std::complex<typename get_pod_type<eT>::result> out_eT;
  
  arma_debug_check( (X.is_vec() == false), "roots(): given object must be a vector" );
  
  if(X.is_finite() == false)  { return false; }
  
  // treat X as a column vector
  
  const Col<eT> Y( const_cast<eT*>(X.memptr()), X.n_elem, false, false);
  
  const T Y_max = (Y.is_empty() == false) ? T(max(abs(Y))) : T(0);
  
  if(Y_max == T(0))  { out.set_size(1,0); return true; }
  
  const uvec indices = find( Y / Y_max );
  
  const uword n_tail_zeros = (indices.n_elem > 0) ? uword( (Y.n_elem-1) - indices[indices.n_elem-1] ) : uword(0);
  
  const Col<eT> Z = Y.subvec( indices[0], indices[indices.n_elem-1] );
  
  if(Z.n_elem >= uword(2))
    {
    Mat<eT> tmp;
    
    if(Z.n_elem == uword(2))
      {
      tmp.set_size(1,1);
      
      tmp[0] = -Z[1] / Z[0];
      }
    else
      {
      tmp = diagmat(ones< Col<eT> >(Z.n_elem - 2), -1);
      
      tmp.row(0) = strans(-Z.subvec(1, Z.n_elem-1) / Z[0]);
      }
    
    Mat<out_eT> junk;
    
    const bool status = auxlib::eig_gen(out, junk, false, tmp);
    
    if(status == false)  { return false; }
    
    if(n_tail_zeros > 0)
      {
      out.resize(out.n_rows + n_tail_zeros, 1);
      }
    }
  else
    {
    out.zeros(n_tail_zeros,1);
    }
  
  return true;
  }



//! @}
