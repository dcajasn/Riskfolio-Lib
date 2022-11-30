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


//! \addtogroup typedef_elem
//! @{


namespace junk
  {
  struct arma_elem_size_test
    {
    arma_static_check( (sizeof(u8) != 1), "error: type 'u8' has unsupported size" );
    arma_static_check( (sizeof(s8) != 1), "error: type 's8' has unsupported size" );
    
    arma_static_check( (sizeof(u16) != 2), "error: type 'u16' has unsupported size" );
    arma_static_check( (sizeof(s16) != 2), "error: type 's16' has unsupported size" );
    
    arma_static_check( (sizeof(u32) != 4), "error: type 'u32' has unsupported size" );
    arma_static_check( (sizeof(s32) != 4), "error: type 's32' has unsupported size" );
    
    arma_static_check( (sizeof(u64) != 8), "error: type 'u64' has unsupported size" );
    arma_static_check( (sizeof(s64) != 8), "error: type 's64' has unsupported size" );
    
    arma_static_check( (sizeof(float)  != 4), "error: type 'float' has unsupported size" );
    arma_static_check( (sizeof(double) != 8), "error: type 'double' has unsupported size" );
    
    arma_static_check( (sizeof(std::complex<float>)  != 8),  "type 'std::complex<float>' has unsupported size" );
    arma_static_check( (sizeof(std::complex<double>) != 16), "type 'std::complex<double>' has unsupported size" );
    };
  }


//! @}
