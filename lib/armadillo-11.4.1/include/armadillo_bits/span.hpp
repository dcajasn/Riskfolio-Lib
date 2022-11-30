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



//! \addtogroup span
//! @{


struct span_alt {};


template<typename junk = int>
class span_base
  {
  public:
  static const span_alt all;
  };


template<typename junk>
const span_alt span_base<junk>::all = span_alt();


class span : public span_base<>
  {
  public:

  uword a;
  uword b;
  bool  whole;
  
  inline
  span()
    : a(0)
    , b(0)
    , whole(true)
    {
    }
  
  
  inline
  span(const span_alt&)
    : a(0)
    , b(0)
    , whole(true)
    {
    }
  
  
  inline
  explicit
  span(const uword in_a)
    : a(in_a)
    , b(in_a)
    , whole(false)
    {
    }
  
  
  // the "explicit" keyword is required here to prevent automatic conversion of {a,b}
  // into an instance of span() when submatrices are specified
  inline
  explicit
  span(const uword in_a, const uword in_b)
    : a(in_a)
    , b(in_b)
    , whole(false)
    {
    }

  };



//! @}
