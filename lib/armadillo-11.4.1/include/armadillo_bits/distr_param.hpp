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



//! \addtogroup distr_param
//! @{



class distr_param
  {
  public:
  
  const uword state;
  
  private:
  
  int a_int;
  int b_int;
    
  double a_double;
  double b_double;
  
  public: 
  
  inline distr_param()
    : state   (0)
    , a_int   (0)
    , b_int   (0)
    , a_double(0)
    , b_double(0)
    {
    }
  
  
  inline explicit distr_param(const int a, const int b)
    : state   (1)
    , a_int   (a)
    , b_int   (b)
    , a_double(double(a))
    , b_double(double(b))
    {
    }
  
  
  inline explicit distr_param(const double a, const double b)
    : state   (2)
    , a_int   (int(a))
    , b_int   (int(b))
    , a_double(a)
    , b_double(b)
    {
    }
  
  
  inline void get_int_vals(int& out_a, int& out_b) const
    {
    if(state == 0)  { return; }
    
    out_a = a_int;
    out_b = b_int;
    }
  
  
  inline void get_double_vals(double& out_a, double& out_b) const
    {
    if(state == 0)  { return; }
    
    out_a = a_double;
    out_b = b_double;
    }
  };



//! @}
