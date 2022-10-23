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


//! \addtogroup wall_clock
//! @{


inline
wall_clock::wall_clock()
  {
  arma_extra_debug_sigprint();
  }



inline
wall_clock::~wall_clock()
  {
  arma_extra_debug_sigprint();
  }



inline
void
wall_clock::tic()
  {
  arma_extra_debug_sigprint();
  
  chrono_time1 = std::chrono::steady_clock::now();
  valid = true;
  }



inline
double
wall_clock::toc()
  {
  arma_extra_debug_sigprint();
  
  if(valid)
    {
    const std::chrono::steady_clock::time_point chrono_time2 = std::chrono::steady_clock::now();
    
    typedef std::chrono::duration<double> duration_type;  // TODO: check this
    
    const duration_type chrono_span = std::chrono::duration_cast< duration_type >(chrono_time2 - chrono_time1);
    
    return chrono_span.count();
    }
  
  return 0.0;
  }


//! @}
