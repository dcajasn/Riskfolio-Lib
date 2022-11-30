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


//! \addtogroup diskio
//! @{


namespace hdf5_opts
  {
  typedef unsigned int flag_type;
  
  struct opts
    {
    const flag_type flags;
    
    inline explicit opts(const flag_type in_flags);
    
    inline const opts operator+(const opts& rhs) const;
    };
  
  inline
  opts::opts(const flag_type in_flags)
    : flags(in_flags)
    {}
  
  inline
  const opts
  opts::operator+(const opts& rhs) const
    {
    const opts result( flags | rhs.flags );
    
    return result;
    }
  
  // The values below (eg. 1u << 0) are for internal Armadillo use only.
  // The values can change without notice.
  
  static const flag_type flag_none    = flag_type(0      );
  static const flag_type flag_trans   = flag_type(1u << 0);
  static const flag_type flag_append  = flag_type(1u << 1);
  static const flag_type flag_replace = flag_type(1u << 2);
  
  struct opts_none    : public opts { inline opts_none()    : opts(flag_none   ) {} };
  struct opts_trans   : public opts { inline opts_trans()   : opts(flag_trans  ) {} };
  struct opts_append  : public opts { inline opts_append()  : opts(flag_append ) {} };
  struct opts_replace : public opts { inline opts_replace() : opts(flag_replace) {} };
  
  static const opts_none    none;
  static const opts_trans   trans;
  static const opts_append  append;
  static const opts_replace replace;
  }


struct hdf5_name
  {
  const std::string     filename;
  const std::string     dsname;
  const hdf5_opts::opts opts;
  
  inline
  hdf5_name(const std::string& in_filename)
    : filename(in_filename    )
    , dsname  (std::string()  )
    , opts    (hdf5_opts::none)
    {}
  
  inline
  hdf5_name(const std::string& in_filename, const std::string& in_dsname, const hdf5_opts::opts& in_opts = hdf5_opts::none)
    : filename(in_filename)
    , dsname  (in_dsname  )
    , opts    (in_opts    )
    {}
  };


//! @}
