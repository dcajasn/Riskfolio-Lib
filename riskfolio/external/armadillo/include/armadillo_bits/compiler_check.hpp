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


#undef ARMA_HAVE_CXX11
#undef ARMA_HAVE_CXX14
#undef ARMA_HAVE_CXX17
#undef ARMA_HAVE_CXX20

#if (__cplusplus >= 201103L)
  #define ARMA_HAVE_CXX11
#endif

#if (__cplusplus >= 201402L)
  #define ARMA_HAVE_CXX14
#endif

#if (__cplusplus >= 201703L)
  #define ARMA_HAVE_CXX17
#endif

#if (__cplusplus >= 202002L)
  #define ARMA_HAVE_CXX20
#endif


// MS really can't get its proverbial shit together
#if defined(_MSVC_LANG)
  
  #if (_MSVC_LANG >= 201402L)
  #undef  ARMA_HAVE_CXX11
  #undef  ARMA_HAVE_CXX14
  
  #define ARMA_HAVE_CXX11
  #define ARMA_HAVE_CXX14
  #endif
  
  #if (_MSVC_LANG >= 201703L)
    #undef  ARMA_HAVE_CXX17
    #define ARMA_HAVE_CXX17
  #endif
  
  #if (_MSVC_LANG >= 202002L)
    #undef  ARMA_HAVE_CXX20
    #define ARMA_HAVE_CXX20
  #endif
  
#endif


// warn about ignored option used in old versions of Armadillo
#if defined(ARMA_DONT_USE_CXX11)
  #pragma message ("WARNING: option ARMA_DONT_USE_CXX11 ignored")
#endif


#if !defined(ARMA_HAVE_CXX11)
  #error "*** C++11 compiler required; enable C++11 mode in your compiler, or use an earlier version of Armadillo"
#endif


// for compatibility with earlier versions of Armadillo
#undef  ARMA_USE_CXX11
#define ARMA_USE_CXX11
