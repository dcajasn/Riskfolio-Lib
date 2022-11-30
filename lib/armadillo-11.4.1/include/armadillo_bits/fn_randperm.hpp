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


//! \addtogroup fn_randperm
//! @{



template<typename obj_type>
inline
void
internal_randperm_helper(obj_type& x, const uword N, const uword N_keep)
  {
  arma_extra_debug_sigprint();
  
  typedef typename obj_type::elem_type eT;
  
  // see op_sort_index_bones.hpp for the definition of arma_sort_index_packet
  // and the associated comparison functor
  
  typedef arma_sort_index_packet<int> packet;
  
  std::vector<packet> packet_vec(N);
  
  for(uword i=0; i < N; ++i)
    {
    packet_vec[i].val   = int(arma_rng::randi<int>());
    packet_vec[i].index = i;
    }
  
  arma_sort_index_helper_ascend<int> comparator;
  
  if(N >= 2)
    {
    if(N_keep < N)
      {
      typename std::vector<packet>::iterator first    = packet_vec.begin();
      typename std::vector<packet>::iterator nth      = first + N_keep;
      typename std::vector<packet>::iterator pastlast = packet_vec.end();
      
      std::partial_sort(first, nth, pastlast, comparator);
      }
    else
      {
      std::sort( packet_vec.begin(), packet_vec.end(), comparator );
      }
    }
  
  if(is_Row<obj_type>::value)
    {
    x.set_size(1,N_keep);
    }
  else
    {
    x.set_size(N_keep,1);
    }
  
  eT* x_mem = x.memptr();
  
  for(uword i=0; i < N_keep; ++i)
    {
    x_mem[i] = eT( packet_vec[i].index );
    }
  }



template<typename obj_type>
arma_warn_unused
inline
typename enable_if2< is_Mat<obj_type>::value, obj_type >::result
randperm(const uword N)
  {
  arma_extra_debug_sigprint();
  
  obj_type x;
  
  if(N > 0)  { internal_randperm_helper(x, N, N); }
  
  return x;
  }



arma_warn_unused
inline
uvec
randperm(const uword N)
  {
  arma_extra_debug_sigprint();
  
  uvec x;
  
  if(N > 0)  { internal_randperm_helper(x, N, N); }
  
  return x;
  }



template<typename obj_type>
arma_warn_unused
inline
typename enable_if2< is_Mat<obj_type>::value, obj_type >::result
randperm(const uword N, const uword M)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (M > N), "randperm(): 'M' must be less than or equal to 'N'" );
  
  obj_type x;
  
  if( (N > 0) && (M > 0) )  { internal_randperm_helper(x, N, M); }
  
  return x;
  }



arma_warn_unused
inline
uvec
randperm(const uword N, const uword M)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (M > N), "randperm(): 'M' must be less than or equal to 'N'" );
  
  uvec x;
  
  if( (N > 0) && (M > 0) )  { internal_randperm_helper(x, N, M); }
  
  return x;
  }



//! @}
