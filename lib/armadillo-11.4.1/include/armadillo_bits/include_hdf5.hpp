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


#if defined(ARMA_USE_HDF5)
  
  #undef  H5_USE_110_API
  #define H5_USE_110_API
  
  #if !defined(ARMA_HDF5_INCLUDE_DIR)
    #if defined(__has_include)
      #if __has_include(<hdf5.h>)
        #include <hdf5.h>
      #else
        #undef ARMA_USE_HDF5
        #undef ARMA_USE_HDF5_CMAKE
        #pragma message ("WARNING: use of HDF5 disabled; hdf5.h header not found")
      #endif
    #else
      #include <hdf5.h>
    #endif
  #else
    #undef ARMA_STR1
    #undef ARMA_STR2
    #undef ARMA_HDF5_HEADER
    
    #define ARMA_STR1(x) x
    #define ARMA_STR2(x) ARMA_STR1(x)
    
    #define ARMA_HDF5_HEADER ARMA_STR2(ARMA_HDF5_INCLUDE_DIR)ARMA_STR2(hdf5.h)
    
    #if defined(__has_include)
      #if __has_include(ARMA_INCFILE_WRAP(ARMA_HDF5_HEADER))
        #include ARMA_INCFILE_WRAP(ARMA_HDF5_HEADER)
      #else
        #undef ARMA_USE_HDF5
        #undef ARMA_USE_HDF5_CMAKE
        #pragma message ("WARNING: use of HDF5 disabled; hdf5.h header not found")
      #endif
    #else
      #include ARMA_INCFILE_WRAP(ARMA_HDF5_HEADER)
    #endif
    
    #undef ARMA_STR1
    #undef ARMA_STR2
    #undef ARMA_HDF5_HEADER
  #endif
  
  #if defined(H5_USE_16_API) || defined(H5_USE_16_API_DEFAULT)
    #pragma message ("WARNING: use of HDF5 disabled; incompatible configuration: H5_USE_16_API or H5_USE_16_API_DEFAULT")
    #undef ARMA_USE_HDF5
    #undef ARMA_USE_HDF5_CMAKE
  #endif
  
#endif
