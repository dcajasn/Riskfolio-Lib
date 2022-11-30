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


//! Generate the first line of the header used for saving matrices in text format.
//! Format: "ARMA_MAT_TXT_ABXYZ".
//! A is one of: I (for integral types) or F (for floating point types).
//! B is one of: U (for unsigned types), S (for signed types), N (for not applicable) or C (for complex types).
//! XYZ specifies the width of each element in terms of bytes, eg. "008" indicates eight bytes.
template<typename eT>
inline
arma_cold
std::string
diskio::gen_txt_header(const Mat<eT>&)
  {
  arma_type_check(( is_supported_elem_type<eT>::value == false ));
  
  const char* ARMA_MAT_TXT_IU001 = "ARMA_MAT_TXT_IU001";
  const char* ARMA_MAT_TXT_IS001 = "ARMA_MAT_TXT_IS001";
  const char* ARMA_MAT_TXT_IU002 = "ARMA_MAT_TXT_IU002";
  const char* ARMA_MAT_TXT_IS002 = "ARMA_MAT_TXT_IS002";
  const char* ARMA_MAT_TXT_IU004 = "ARMA_MAT_TXT_IU004";
  const char* ARMA_MAT_TXT_IS004 = "ARMA_MAT_TXT_IS004";
  const char* ARMA_MAT_TXT_IU008 = "ARMA_MAT_TXT_IU008";
  const char* ARMA_MAT_TXT_IS008 = "ARMA_MAT_TXT_IS008";
  const char* ARMA_MAT_TXT_FN004 = "ARMA_MAT_TXT_FN004";
  const char* ARMA_MAT_TXT_FN008 = "ARMA_MAT_TXT_FN008";
  const char* ARMA_MAT_TXT_FC008 = "ARMA_MAT_TXT_FC008";
  const char* ARMA_MAT_TXT_FC016 = "ARMA_MAT_TXT_FC016";
  
  char* header = nullptr;
  
       if(       is_u8<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_IU001); }
  else if(       is_s8<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_IS001); }
  else if(      is_u16<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_IU002); }
  else if(      is_s16<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_IS002); }
  else if(      is_u32<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_IU004); }
  else if(      is_s32<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_IS004); }
  else if(      is_u64<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_IU008); }
  else if(      is_s64<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_IS008); }
  else if(is_ulng_t_32<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_IU004); }
  else if(is_slng_t_32<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_IS004); }
  else if(is_ulng_t_64<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_IU008); }
  else if(is_slng_t_64<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_IS008); }
  else if(    is_float<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_FN004); }
  else if(   is_double<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_FN008); }
  else if( is_cx_float<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_FC008); }
  else if(is_cx_double<eT>::value)  { header = const_cast<char*>(ARMA_MAT_TXT_FC016); }
  
  return std::string(header);
  }



//! Generate the first line of the header used for saving matrices in binary format.
//! Format: "ARMA_MAT_BIN_ABXYZ".
//! A is one of: I (for integral types) or F (for floating point types).
//! B is one of: U (for unsigned types), S (for signed types), N (for not applicable) or C (for complex types).
//! XYZ specifies the width of each element in terms of bytes, eg. "008" indicates eight bytes.
template<typename eT>
inline
arma_cold
std::string
diskio::gen_bin_header(const Mat<eT>&)
  {
  arma_type_check(( is_supported_elem_type<eT>::value == false ));
  
  const char* ARMA_MAT_BIN_IU001 = "ARMA_MAT_BIN_IU001";
  const char* ARMA_MAT_BIN_IS001 = "ARMA_MAT_BIN_IS001";
  const char* ARMA_MAT_BIN_IU002 = "ARMA_MAT_BIN_IU002";
  const char* ARMA_MAT_BIN_IS002 = "ARMA_MAT_BIN_IS002";
  const char* ARMA_MAT_BIN_IU004 = "ARMA_MAT_BIN_IU004";
  const char* ARMA_MAT_BIN_IS004 = "ARMA_MAT_BIN_IS004";
  const char* ARMA_MAT_BIN_IU008 = "ARMA_MAT_BIN_IU008";
  const char* ARMA_MAT_BIN_IS008 = "ARMA_MAT_BIN_IS008";
  const char* ARMA_MAT_BIN_FN004 = "ARMA_MAT_BIN_FN004";
  const char* ARMA_MAT_BIN_FN008 = "ARMA_MAT_BIN_FN008";
  const char* ARMA_MAT_BIN_FC008 = "ARMA_MAT_BIN_FC008";
  const char* ARMA_MAT_BIN_FC016 = "ARMA_MAT_BIN_FC016";  
  
  char* header = nullptr;
  
       if(       is_u8<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_IU001); }
  else if(       is_s8<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_IS001); }
  else if(      is_u16<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_IU002); }
  else if(      is_s16<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_IS002); }
  else if(      is_u32<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_IU004); }
  else if(      is_s32<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_IS004); }
  else if(      is_u64<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_IU008); }
  else if(      is_s64<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_IS008); }
  else if(is_ulng_t_32<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_IU004); }
  else if(is_slng_t_32<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_IS004); }
  else if(is_ulng_t_64<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_IU008); }
  else if(is_slng_t_64<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_IS008); }
  else if(    is_float<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_FN004); }
  else if(   is_double<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_FN008); }
  else if( is_cx_float<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_FC008); }
  else if(is_cx_double<eT>::value)  { header = const_cast<char*>(ARMA_MAT_BIN_FC016); }
  
  return std::string(header);
  }



//! Generate the first line of the header used for saving matrices in binary format.
//! Format: "ARMA_SPM_BIN_ABXYZ".
//! A is one of: I (for integral types) or F (for floating point types).
//! B is one of: U (for unsigned types), S (for signed types), N (for not applicable) or C (for complex types).
//! XYZ specifies the width of each element in terms of bytes, eg. "008" indicates eight bytes.
template<typename eT>
inline
arma_cold
std::string
diskio::gen_bin_header(const SpMat<eT>&)
  {
  arma_type_check(( is_supported_elem_type<eT>::value == false ));
  
  const char* ARMA_SPM_BIN_IU001 = "ARMA_SPM_BIN_IU001";
  const char* ARMA_SPM_BIN_IS001 = "ARMA_SPM_BIN_IS001";
  const char* ARMA_SPM_BIN_IU002 = "ARMA_SPM_BIN_IU002";
  const char* ARMA_SPM_BIN_IS002 = "ARMA_SPM_BIN_IS002";
  const char* ARMA_SPM_BIN_IU004 = "ARMA_SPM_BIN_IU004";
  const char* ARMA_SPM_BIN_IS004 = "ARMA_SPM_BIN_IS004";
  const char* ARMA_SPM_BIN_IU008 = "ARMA_SPM_BIN_IU008";
  const char* ARMA_SPM_BIN_IS008 = "ARMA_SPM_BIN_IS008";
  const char* ARMA_SPM_BIN_FN004 = "ARMA_SPM_BIN_FN004";
  const char* ARMA_SPM_BIN_FN008 = "ARMA_SPM_BIN_FN008";
  const char* ARMA_SPM_BIN_FC008 = "ARMA_SPM_BIN_FC008";
  const char* ARMA_SPM_BIN_FC016 = "ARMA_SPM_BIN_FC016";  
  
  char* header = nullptr;
  
       if(       is_u8<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_IU001); }
  else if(       is_s8<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_IS001); }
  else if(      is_u16<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_IU002); }
  else if(      is_s16<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_IS002); }
  else if(      is_u32<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_IU004); }
  else if(      is_s32<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_IS004); }
  else if(      is_u64<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_IU008); }
  else if(      is_s64<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_IS008); }
  else if(is_ulng_t_32<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_IU004); }
  else if(is_slng_t_32<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_IS004); }
  else if(is_ulng_t_64<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_IU008); }
  else if(is_slng_t_64<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_IS008); }
  else if(    is_float<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_FN004); }
  else if(   is_double<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_FN008); }
  else if( is_cx_float<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_FC008); }
  else if(is_cx_double<eT>::value)  { header = const_cast<char*>(ARMA_SPM_BIN_FC016); }
  
  return std::string(header);
  }


//! Generate the first line of the header used for saving cubes in text format.
//! Format: "ARMA_CUB_TXT_ABXYZ".
//! A is one of: I (for integral types) or F (for floating point types).
//! B is one of: U (for unsigned types), S (for signed types), N (for not applicable) or C (for complex types).
//! XYZ specifies the width of each element in terms of bytes, eg. "008" indicates eight bytes.
template<typename eT>
inline
arma_cold
std::string
diskio::gen_txt_header(const Cube<eT>&)
  {
  arma_type_check(( is_supported_elem_type<eT>::value == false ));

  const char* ARMA_CUB_TXT_IU001 = "ARMA_CUB_TXT_IU001";
  const char* ARMA_CUB_TXT_IS001 = "ARMA_CUB_TXT_IS001";
  const char* ARMA_CUB_TXT_IU002 = "ARMA_CUB_TXT_IU002";
  const char* ARMA_CUB_TXT_IS002 = "ARMA_CUB_TXT_IS002";
  const char* ARMA_CUB_TXT_IU004 = "ARMA_CUB_TXT_IU004";
  const char* ARMA_CUB_TXT_IS004 = "ARMA_CUB_TXT_IS004";
  const char* ARMA_CUB_TXT_IU008 = "ARMA_CUB_TXT_IU008";
  const char* ARMA_CUB_TXT_IS008 = "ARMA_CUB_TXT_IS008";
  const char* ARMA_CUB_TXT_FN004 = "ARMA_CUB_TXT_FN004";
  const char* ARMA_CUB_TXT_FN008 = "ARMA_CUB_TXT_FN008";
  const char* ARMA_CUB_TXT_FC008 = "ARMA_CUB_TXT_FC008";
  const char* ARMA_CUB_TXT_FC016 = "ARMA_CUB_TXT_FC016";
  
  char* header = nullptr;
  
       if(       is_u8<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_IU001); }
  else if(       is_s8<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_IS001); }
  else if(      is_u16<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_IU002); }
  else if(      is_s16<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_IS002); }
  else if(      is_u32<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_IU004); }
  else if(      is_s32<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_IS004); }
  else if(      is_u64<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_IU008); }
  else if(      is_s64<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_IS008); }
  else if(is_ulng_t_32<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_IU004); }
  else if(is_slng_t_32<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_IS004); }
  else if(is_ulng_t_64<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_IU008); }
  else if(is_slng_t_64<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_IS008); }
  else if(    is_float<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_FN004); }
  else if(   is_double<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_FN008); }
  else if( is_cx_float<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_FC008); }
  else if(is_cx_double<eT>::value)  { header = const_cast<char*>(ARMA_CUB_TXT_FC016); }
  
  return std::string(header);
  }



//! Generate the first line of the header used for saving cubes in binary format.
//! Format: "ARMA_CUB_BIN_ABXYZ".
//! A is one of: I (for integral types) or F (for floating point types).
//! B is one of: U (for unsigned types), S (for signed types), N (for not applicable) or C (for complex types).
//! XYZ specifies the width of each element in terms of bytes, eg. "008" indicates eight bytes.
template<typename eT>
inline
arma_cold
std::string
diskio::gen_bin_header(const Cube<eT>&)
  {
  arma_type_check(( is_supported_elem_type<eT>::value == false ));
  
  const char* ARMA_CUB_BIN_IU001 = "ARMA_CUB_BIN_IU001";
  const char* ARMA_CUB_BIN_IS001 = "ARMA_CUB_BIN_IS001";
  const char* ARMA_CUB_BIN_IU002 = "ARMA_CUB_BIN_IU002";
  const char* ARMA_CUB_BIN_IS002 = "ARMA_CUB_BIN_IS002";
  const char* ARMA_CUB_BIN_IU004 = "ARMA_CUB_BIN_IU004";
  const char* ARMA_CUB_BIN_IS004 = "ARMA_CUB_BIN_IS004";
  const char* ARMA_CUB_BIN_IU008 = "ARMA_CUB_BIN_IU008";
  const char* ARMA_CUB_BIN_IS008 = "ARMA_CUB_BIN_IS008";
  const char* ARMA_CUB_BIN_FN004 = "ARMA_CUB_BIN_FN004";
  const char* ARMA_CUB_BIN_FN008 = "ARMA_CUB_BIN_FN008";
  const char* ARMA_CUB_BIN_FC008 = "ARMA_CUB_BIN_FC008";
  const char* ARMA_CUB_BIN_FC016 = "ARMA_CUB_BIN_FC016";
  
  char* header = nullptr;
  
       if(       is_u8<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_IU001); }
  else if(       is_s8<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_IS001); }
  else if(      is_u16<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_IU002); }
  else if(      is_s16<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_IS002); }
  else if(      is_u32<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_IU004); }
  else if(      is_s32<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_IS004); }
  else if(      is_u64<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_IU008); }
  else if(      is_s64<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_IS008); }
  else if(is_ulng_t_32<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_IU004); }
  else if(is_slng_t_32<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_IS004); }
  else if(is_ulng_t_64<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_IU008); }
  else if(is_slng_t_64<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_IS008); }
  else if(    is_float<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_FN004); }
  else if(   is_double<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_FN008); }
  else if( is_cx_float<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_FC008); }
  else if(is_cx_double<eT>::value)  { header = const_cast<char*>(ARMA_CUB_BIN_FC016); }
  
  return std::string(header);
  }



inline
arma_deprecated
file_type
diskio::guess_file_type(std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  return diskio::guess_file_type_internal(f);
  }



inline
arma_cold
file_type
diskio::guess_file_type_internal(std::istream& f)
  {
  arma_extra_debug_sigprint();
  
  f.clear();
  const std::fstream::pos_type pos1 = f.tellg();
  
  f.clear();
  f.seekg(0, ios::end);
  
  f.clear();
  const std::fstream::pos_type pos2 = f.tellg();
  
  const uword N_max = ( (pos1 >= 0) && (pos2 >= 0) && (pos2 > pos1) ) ? uword(pos2 - pos1) : uword(0);
  
  f.clear();
  f.seekg(pos1);
  
  if(N_max == 0)  { return file_type_unknown; }
  
  const uword N_use = (std::min)(N_max, uword(4096));
  
  podarray<unsigned char> data(N_use);
  data.zeros();
  
  unsigned char* data_mem = data.memptr();
  
  f.clear();
  f.read( reinterpret_cast<char*>(data_mem), std::streamsize(N_use) );
  
  const bool load_okay = f.good();
  
  f.clear();
  f.seekg(pos1);
  
  if(load_okay == false)  { return file_type_unknown; }
  
  bool has_binary    = false;
  bool has_bracket   = false;
  bool has_comma     = false;
  bool has_semicolon = false;
  
  for(uword i=0; i<N_use; ++i)
    {
    const unsigned char val = data_mem[i];
    
    if( (val <=   8) || (val >= 123) )  { has_binary    = true; break; }  // the range checking can be made more elaborate
    
    if( (val == '(') || (val == ')') )  { has_bracket   = true;        }
    
    if( (val == ';')                 )  { has_semicolon = true;        }
    
    if( (val == ',')                 )  { has_comma     = true;        }
    }
  
  if(has_binary)  { return raw_binary; }
  
  // ssv_ascii has to be before csv_ascii;
  // if the data has semicolons, it suggests a CSV file with semicolon as the separating character;
  // the semicolon may be used to allow the comma character to represent the decimal point (eg. 1,2345 vs 1.2345)
  
  if(has_semicolon && (has_bracket == false))  { return ssv_ascii; }
  
  if(has_comma     && (has_bracket == false))  { return csv_ascii; }
  
  return raw_ascii;
  }



//! Append a quasi-random string to the given filename.
//! Avoiding use of rand() to preserve its state. 
inline
arma_cold
std::string
diskio::gen_tmp_name(const std::string& x)
  {
  union { uword val; void* ptr; } u;
  
  u.val = uword(0);
  u.ptr = const_cast<std::string*>(&x);
  
  const u16 a = u16( (u.val >> 8)   & 0xFFFF );
  const u16 b = u16( (std::clock()) & 0xFFFF );
  
  std::ostringstream ss;
  
  ss << x << ".tmp_";
  
  ss.setf(std::ios_base::hex, std::ios_base::basefield);
  
  ss.width(4);
  ss.fill('0');
  ss << a;
  
  ss.width(4);
  ss.fill('0');
  ss << b;
  
  return ss.str();
  }



//! Safely rename a file.
//! Before renaming, test if we can write to the final file.
//! This should prevent:
//! (i)  overwriting files that are write protected,
//! (ii) overwriting directories.
inline
arma_cold
bool
diskio::safe_rename(const std::string& old_name, const std::string& new_name)
  {
  const char* new_name_c_str = new_name.c_str();
  
  std::fstream f(new_name_c_str, std::fstream::out | std::fstream::app);
  f.put(' ');
  
  if(f.good()) { f.close(); } else { return false; }
  
  if(std::remove(                  new_name_c_str) != 0)  { return false; }
  if(std::rename(old_name.c_str(), new_name_c_str) != 0)  { return false; }
  
  return true;
  }



template<typename eT>
inline
bool
diskio::convert_token(eT& val, const std::string& token)
  {
  const size_t N = size_t(token.length());
  
  if(N == 0)  { val = eT(0); return true; }
  
  const char* str = token.c_str();
  
  if( (N == 3) || (N == 4) )
    {
    const bool neg = (str[0] == '-');
    const bool pos = (str[0] == '+');
    
    const size_t offset = ( (neg || pos) && (N == 4) ) ? 1 : 0;
    
    const char sig_a = str[offset  ];
    const char sig_b = str[offset+1];
    const char sig_c = str[offset+2];
    
    if( ((sig_a == 'i') || (sig_a == 'I')) && ((sig_b == 'n') || (sig_b == 'N')) && ((sig_c == 'f') || (sig_c == 'F')) )
      {
      val = neg ? cond_rel< is_signed<eT>::value >::make_neg(Datum<eT>::inf) : Datum<eT>::inf;
      
      return true;
      }
    else
    if( ((sig_a == 'n') || (sig_a == 'N')) && ((sig_b == 'a') || (sig_b == 'A')) && ((sig_c == 'n') || (sig_c == 'N')) )
      {
      val = Datum<eT>::nan;
      
      return true;
      }
    }
  
  
  char* endptr = nullptr;
  
  if(is_real<eT>::value)
    {
    val = eT( std::strtod(str, &endptr) );
    }
  else
    {
    if(is_signed<eT>::value)
      {
      // signed integer
      
      val = eT( std::strtoll(str, &endptr, 10) );
      }
    else
      {
      // unsigned integer
      
      if((str[0] == '-') && (N >= 2))
        {
        val = eT(0);  
        
        if((str[1] == '-') || (str[1] == '+')) { return false; }
        
        const char* str_offset1 = &(str[1]);
        
        std::strtoull(str_offset1, &endptr, 10);
        
        if(str_offset1 == endptr)  { return false; }
        
        return true;
        }
      
      val = eT( std::strtoull(str, &endptr, 10) );
      }
    }
  
  if(str == endptr)  { return false; }
  
  return true;
  }



template<typename T>
inline
bool
diskio::convert_token(std::complex<T>& val, const std::string& token)
  {
  const size_t N   = size_t(token.length());
  const size_t Nm1 = N-1;
  
  if(N == 0)  { val = std::complex<T>(0); return true; }
  
  const char* str = token.c_str();
  
  // valid complex number formats:
  // (real,imag)
  // (real)
  // ()
  
  if( (token[0] != '(') || (token[Nm1] != ')') )
    {
    // no brackets, so treat the token as a non-complex number
    
    T val_real;
    
    const bool state = diskio::convert_token(val_real, token);  // use the non-complex version of this function
    
    val = std::complex<T>(val_real);
    
    return state;
    }
  
  // does the token contain only the () brackets?
  if(N <= 2)  { val = std::complex<T>(0); return true; }
  
  size_t comma_loc   = 0;
  bool   comma_found = false;
  
  for(size_t i=0; i<N; ++i)  { if(str[i] == ',')  { comma_loc = i; comma_found = true; break; } }
  
  bool state = false;
  
  if(comma_found == false)
    {
    // only the real part is available
    
    const std::string token_real( &(str[1]), (Nm1 - 1) );
    
    T val_real;
    
    state = diskio::convert_token(val_real, token_real);  // use the non-complex version of this function
    
    val = std::complex<T>(val_real);
    }
  else
    {
    const std::string token_real( &(str[1]),           (comma_loc - 1      ) );
    const std::string token_imag( &(str[comma_loc+1]), (Nm1 - 1 - comma_loc) );
    
    T val_real;
    T val_imag;
    
    const bool state_real = diskio::convert_token(val_real, token_real);
    const bool state_imag = diskio::convert_token(val_imag, token_imag);
    
    state = (state_real && state_imag);
    
    val = std::complex<T>(val_real, val_imag);
    }
  
  return state;
  }



template<typename eT>
inline
std::streamsize
diskio::prepare_stream(std::ostream& f)
  {
  std::streamsize cell_width = f.width();
  
  if(is_real<eT>::value)
    {
    f.unsetf(ios::fixed);
    f.setf(ios::scientific);
    f.fill(' ');
    
    f.precision(16);
    cell_width = 24;
    
    // NOTE: for 'float' the optimum settings are f.precision(8) and cell_width = 15
    // NOTE: however, to avoid introducing errors in case single precision data is loaded as double precision,
    // NOTE: the same settings must be used for both 'float' and 'double'
    }
  else
  if(is_cx<eT>::value)
    {
    f.unsetf(ios::fixed);
    f.setf(ios::scientific);
    
    f.precision(16);
    }
  
  return cell_width;
  }
  



//! Save a matrix as raw text (no header, human readable).
//! Matrices can be loaded in Matlab and Octave, as long as they don't have complex elements.
template<typename eT>
inline
bool
diskio::save_raw_ascii(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::fstream f(tmp_name.c_str(), std::fstream::out);
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_raw_ascii(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



//! Save a matrix as raw text (no header, human readable).
//! Matrices can be loaded in Matlab and Octave, as long as they don't have complex elements.
template<typename eT>
inline
bool
diskio::save_raw_ascii(const Mat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(f);
  
  const std::streamsize cell_width = diskio::prepare_stream<eT>(f);
  
  for(uword row=0; row < x.n_rows; ++row)
    {
    for(uword col=0; col < x.n_cols; ++col)
      {
      f.put(' ');
      
      if(is_real<eT>::value)  { f.width(cell_width); }
      
      arma_ostream::raw_print_elem(f, x.at(row,col));
      }
    
    f.put('\n');
    }
  
  const bool save_okay = f.good();
  
  stream_state.restore(f);
  
  return save_okay;
  }



//! Save a matrix as raw binary (no header)
template<typename eT>
inline
bool
diskio::save_raw_binary(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str(), std::fstream::binary);
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_raw_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



template<typename eT>
inline
bool
diskio::save_raw_binary(const Mat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  f.write( reinterpret_cast<const char*>(x.mem), std::streamsize(x.n_elem*sizeof(eT)) );
  
  return f.good();
  }



//! Save a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_ascii(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str());
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_arma_ascii(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



//! Save a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_ascii(const Mat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(f);
  
  f << diskio::gen_txt_header(x) << '\n';
  f << x.n_rows << ' ' << x.n_cols << '\n';
  
  const std::streamsize cell_width = diskio::prepare_stream<eT>(f);
  
  for(uword row=0; row < x.n_rows; ++row)
    {
    for(uword col=0; col < x.n_cols; ++col)
      {
      f.put(' ');
      
      if(is_real<eT>::value)  { f.width(cell_width); }
      
      arma_ostream::raw_print_elem(f, x.at(row,col));
      }
    
    f.put('\n');
    }
  
  const bool save_okay = f.good();
  
  stream_state.restore(f);
  
  return save_okay;
  }



//! Save a matrix in CSV text format (human readable)
template<typename eT>
inline
bool
diskio::save_csv_ascii(const Mat<eT>& x, const std::string& final_name, const field<std::string>& header, const bool with_header, const char separator)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str());
  
  bool save_okay = f.is_open();
  
  if(save_okay == false)  { return false; }
  
  if(with_header)
    {
    arma_extra_debug_print("diskio::save_csv_ascii(): writing header");
    
    for(uword i=0; i < header.n_elem; ++i)
      {
      f << header.at(i);
      
      if(i != (header.n_elem-1))  { f.put(separator); }
      }
    
    f.put('\n');
    
    save_okay = f.good();
    }
  
  if(save_okay)  { save_okay = diskio::save_csv_ascii(x, f, separator); }
  
  f.flush();
  f.close();
  
  if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
  
  return save_okay;
  }



//! Save a matrix in CSV text format (human readable)
template<typename eT>
inline
bool
diskio::save_csv_ascii(const Mat<eT>& x, std::ostream& f, const char separator)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(f);
  
  diskio::prepare_stream<eT>(f);
  
  uword x_n_rows = x.n_rows;
  uword x_n_cols = x.n_cols;
  
  for(uword row=0; row < x_n_rows; ++row)
    {
    for(uword col=0; col < x_n_cols; ++col)
      {
      arma_ostream::raw_print_elem(f, x.at(row,col));
      
      if( col < (x_n_cols-1) )  { f.put(separator); }
      }
    
    f.put('\n');
    }
  
  const bool save_okay = f.good();
  
  stream_state.restore(f);
  
  return save_okay;
  }



//! Save a matrix in CSV text format (human readable); complex numbers stored in "a+bi" format
template<typename T>
inline
bool
diskio::save_csv_ascii(const Mat< std::complex<T> >& x, std::ostream& f, const char separator)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  const arma_ostream_state stream_state(f);
  
  diskio::prepare_stream<eT>(f);
  
  uword x_n_rows = x.n_rows;
  uword x_n_cols = x.n_cols;
  
  for(uword row=0; row < x_n_rows; ++row)
    {
    for(uword col=0; col < x_n_cols; ++col)
      {
      const eT& val = x.at(row,col);
      
      const T    tmp_r     = std::real(val);
      const T    tmp_i     = std::imag(val);
      const T    tmp_i_abs = (tmp_i < T(0)) ? T(-tmp_i) : T(tmp_i);
      const char tmp_sign  = (tmp_i < T(0)) ? char('-') : char('+');
      
      arma_ostream::raw_print_elem(f, tmp_r    );
      f.put(tmp_sign);
      arma_ostream::raw_print_elem(f, tmp_i_abs);
      f.put('i');
      
      if( col < (x_n_cols-1) )  { f.put(separator); }
      }
    
    f.put('\n');
    }
  
  const bool save_okay = f.good();
  
  stream_state.restore(f);
  
  return save_okay;
  }



template<typename eT>
inline
bool
diskio::save_coord_ascii(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str());
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_coord_ascii(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



template<typename eT>
inline
bool
diskio::save_coord_ascii(const Mat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(f);
  
  diskio::prepare_stream<eT>(f);
  
  for(uword col=0; col < x.n_cols; ++col)
  for(uword row=0; row < x.n_rows; ++row)
    {
    const eT val = x.at(row,col);
    
    if(val != eT(0))
      {
      f << row << ' ' << col << ' ' << val << '\n';
      }
    }
  
  // make sure it's possible to figure out the matrix size later
  if( (x.n_rows > 0) && (x.n_cols > 0) )
    {
    const uword max_row = (x.n_rows > 0) ? x.n_rows-1 : 0;
    const uword max_col = (x.n_cols > 0) ? x.n_cols-1 : 0;
    
    if( x.at(max_row, max_col) == eT(0) )
      {
      f << max_row << ' ' << max_col << " 0\n";
      }
    }
  
  const bool save_okay = f.good();
  
  stream_state.restore(f);
  
  return save_okay;
  }



template<typename T>
inline
bool
diskio::save_coord_ascii(const Mat< std::complex<T> >& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  const arma_ostream_state stream_state(f);
  
  diskio::prepare_stream<eT>(f);
  
  const eT eT_zero = eT(0);
  
  for(uword col=0; col < x.n_cols; ++col)
  for(uword row=0; row < x.n_rows; ++row)
    {
    const eT val = x.at(row,col);
    
    if(val != eT_zero)
      {
      f << row << ' ' << col << ' ' << val.real() << ' ' << val.imag() << '\n';
      }
    }
  
  // make sure it's possible to figure out the matrix size later
  if( (x.n_rows > 0) && (x.n_cols > 0) )
    {
    const uword max_row = (x.n_rows > 0) ? x.n_rows-1 : 0;
    const uword max_col = (x.n_cols > 0) ? x.n_cols-1 : 0;
    
    if( x.at(max_row, max_col) == eT_zero )
      {
      f << max_row << ' ' << max_col << " 0 0\n";
      }
    }
  
  const bool save_okay = f.good();
  
  stream_state.restore(f);
  
  return save_okay;
  }



//! Save a matrix in binary format,
//! with a header that stores the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_binary(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str(), std::fstream::binary);
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_arma_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



//! Save a matrix in binary format,
//! with a header that stores the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_binary(const Mat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  f << diskio::gen_bin_header(x) << '\n';
  f << x.n_rows << ' ' << x.n_cols << '\n';
  
  f.write( reinterpret_cast<const char*>(x.mem), std::streamsize(x.n_elem*sizeof(eT)) );
  
  return f.good();
  }



//! Save a matrix as a PGM greyscale image
template<typename eT>
inline
bool
diskio::save_pgm_binary(const Mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::fstream f(tmp_name.c_str(), std::fstream::out | std::fstream::binary);
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_pgm_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



//
// TODO:
// add functionality to save the image in a normalised format,
// ie. scaled so that every value falls in the [0,255] range.

//! Save a matrix as a PGM greyscale image
template<typename eT>
inline
bool
diskio::save_pgm_binary(const Mat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  f << "P5" << '\n';
  f << x.n_cols << ' ' << x.n_rows << '\n';
  f << 255 << '\n';
  
  const uword n_elem = x.n_rows * x.n_cols;
  podarray<u8> tmp(n_elem);
  
  uword i = 0;
  
  for(uword row=0; row < x.n_rows; ++row)
  for(uword col=0; col < x.n_cols; ++col)
    {
    tmp[i] = u8( x.at(row,col) );  // TODO: add round() ?
    ++i;
    }
  
  f.write(reinterpret_cast<const char*>(tmp.mem), std::streamsize(n_elem) );
  
  return f.good();
  }



//! Save a matrix as a PGM greyscale image
template<typename T>
inline
bool
diskio::save_pgm_binary(const Mat< std::complex<T> >& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const uchar_mat tmp = conv_to<uchar_mat>::from(x);
  
  return diskio::save_pgm_binary(tmp, final_name);
  }



//! Save a matrix as a PGM greyscale image
template<typename T>
inline
bool
diskio::save_pgm_binary(const Mat< std::complex<T> >& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const uchar_mat tmp = conv_to<uchar_mat>::from(x);
  
  return diskio::save_pgm_binary(tmp, f);
  }



//! Save a matrix as part of a HDF5 file
template<typename eT>
inline 
bool 
diskio::save_hdf5_binary(const Mat<eT>& x, const hdf5_name& spec, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_HDF5)
    {
    hdf5_misc::hdf5_suspend_printing_errors hdf5_print_suspender;
    
    bool save_okay = false;
    
    const bool append  = bool(spec.opts.flags & hdf5_opts::flag_append);
    const bool replace = bool(spec.opts.flags & hdf5_opts::flag_replace);
    
    const bool use_existing_file = ((append || replace) && (arma_H5Fis_hdf5(spec.filename.c_str()) > 0));
    
    const std::string tmp_name = (use_existing_file) ? std::string() : diskio::gen_tmp_name(spec.filename);
    
    // Set up the file according to HDF5's preferences
    hid_t file = (use_existing_file) ? arma_H5Fopen(spec.filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT) : arma_H5Fcreate(tmp_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    if(file < 0)  { return false; }
    
    // We need to create a dataset, datatype, and dataspace
    hsize_t dims[2];
    dims[1] = x.n_rows;
    dims[0] = x.n_cols;
    
    hid_t dataspace = arma_H5Screate_simple(2, dims, NULL);   // treat the matrix as a 2d array dataspace
    hid_t datatype  = hdf5_misc::get_hdf5_type<eT>();
    
    // If this returned something invalid, well, it's time to crash.
    arma_check(datatype == -1, "Mat::save(): unknown datatype for HDF5");
    
    // MATLAB forces the users to specify a name at save time for HDF5;
    // Octave will use the default of 'dataset' unless otherwise specified.
    // If the user hasn't specified a dataset name, we will use 'dataset'
    // We may have to split out the group name from the dataset name.
    std::vector<hid_t> groups;
    std::string full_name = spec.dsname;
    size_t loc;
    while((loc = full_name.find("/")) != std::string::npos)
      {
      // Create another group...
      if(loc != 0) // Ignore the first /, if there is a leading /.
        {
        hid_t gid = arma_H5Gcreate((groups.size() == 0) ? file : groups[groups.size() - 1], full_name.substr(0, loc).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        if((gid < 0) && use_existing_file)
          {
          gid = arma_H5Gopen((groups.size() == 0) ? file : groups[groups.size() - 1], full_name.substr(0, loc).c_str(), H5P_DEFAULT);
          }
        
        groups.push_back(gid);
        }
      
      full_name = full_name.substr(loc + 1);
      }
    
    const std::string dataset_name = full_name.empty() ? std::string("dataset") : full_name;
    
    const hid_t last_group = (groups.size() == 0) ? file : groups[groups.size() - 1];
    
    if(use_existing_file && replace)
      {
      arma_H5Ldelete(last_group, dataset_name.c_str(), H5P_DEFAULT);
      // NOTE: H5Ldelete() in HDF5 v1.8 doesn't reclaim the deleted space; use h5repack to reclaim space: h5repack oldfile.h5 newfile.h5
      // NOTE: has this behaviour changed in HDF5 1.10 ?
      // NOTE: https://lists.hdfgroup.org/pipermail/hdf-forum_lists.hdfgroup.org/2017-August/010482.html
      // NOTE: https://lists.hdfgroup.org/pipermail/hdf-forum_lists.hdfgroup.org/2017-August/010486.html
      }
    
    hid_t dataset = arma_H5Dcreate(last_group, dataset_name.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    if(dataset < 0)
      {
      save_okay = false;
      
      err_msg = "couldn't create dataset";
      }
    else
      {
      save_okay = (arma_H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, x.mem) >= 0);
      
      arma_H5Dclose(dataset);
      }
    
    arma_H5Tclose(datatype);
    arma_H5Sclose(dataspace);
    for(size_t i = 0; i < groups.size(); ++i)  { arma_H5Gclose(groups[i]); }
    arma_H5Fclose(file);
    
    if((use_existing_file == false) && (save_okay == true))  { save_okay = diskio::safe_rename(tmp_name, spec.filename); }
    
    return save_okay;
    }
  #else
    {
    arma_ignore(x);
    arma_ignore(spec);
    arma_ignore(err_msg);
    
    arma_stop_logic_error("Mat::save(): use of HDF5 must be enabled");
    
    return false;
    }
  #endif
  }



//! Load a matrix as raw text (no header, human readable).
//! Can read matrices saved as text in Matlab and Octave.
//! NOTE: this is much slower than reading a file with a header.
template<typename eT>
inline
bool
diskio::load_raw_ascii(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in);
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_raw_ascii(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



//! Load a matrix as raw text (no header, human readable).
//! Can read matrices saved as text in Matlab and Octave.
//! NOTE: this is much slower than reading a file with a header.
template<typename eT>
inline
bool
diskio::load_raw_ascii(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = f.good();
  
  f.clear();
  const std::fstream::pos_type pos1 = f.tellg();
  
  //
  // work out the size
  
  uword f_n_rows = 0;
  uword f_n_cols = 0;
  
  bool f_n_cols_found = false;
  
  std::string       line_string;
  std::stringstream line_stream;
  
  std::string token;
  
  while( f.good() && load_okay )
    {
    std::getline(f, line_string);
    
    // TODO: does it make sense to stop processing the file if an empty line is found ?
    if(line_string.size() == 0)  { break; }
    
    line_stream.clear();
    line_stream.str(line_string);
    
    uword line_n_cols = 0;
    
    while(line_stream >> token)  { ++line_n_cols; }
    
    if(f_n_cols_found == false)
      {
      f_n_cols       = line_n_cols;
      f_n_cols_found = true;
      }
    else
      {
      if(line_n_cols != f_n_cols)
        {
        load_okay = false;
        err_msg = "inconsistent number of columns";
        }
      }
    
    ++f_n_rows;
    }
  
  
  if(load_okay)
    {
    f.clear();
    f.seekg(pos1);
    
    try { x.set_size(f_n_rows, f_n_cols); } catch(...) { err_msg = "not enough memory"; return false; }
    
    for(uword row=0; ((row < x.n_rows) && load_okay); ++row)
    for(uword col=0; ((col < x.n_cols) && load_okay); ++col)
      {
      f >> token;
      
      if(diskio::convert_token(x.at(row,col), token) == false)
        {
        load_okay = false;
        err_msg = "couldn't interpret data";
        }
      }
    }
  
  
  // an empty file indicates an empty matrix
  if( (f_n_cols_found == false) && (load_okay == true) )  { x.reset(); }
  
  
  return load_okay;
  }



//! Load a matrix in binary format (no header);
//! the matrix is assumed to have one column
template<typename eT>
inline
bool
diskio::load_raw_binary(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f;
  f.open(name.c_str(), std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_raw_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::load_raw_binary(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  f.clear();
  const std::streampos pos1 = f.tellg();
  
  f.clear();
  f.seekg(0, ios::end);

  f.clear();
  const std::streampos pos2 = f.tellg();
  
  const uword N = ( (pos1 >= 0) && (pos2 >= 0) ) ? uword(pos2 - pos1) : 0;
  
  f.clear();
  //f.seekg(0, ios::beg);
  f.seekg(pos1);
  
  try { x.set_size(N / uword(sizeof(eT)), 1); } catch(...) { err_msg = "not enough memory"; return false; }
  
  f.clear();
  f.read( reinterpret_cast<char *>(x.memptr()), std::streamsize(x.n_elem * uword(sizeof(eT))) );
  
  return f.good();
  }



//! Load a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::load_arma_ascii(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f(name.c_str());
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_arma_ascii(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



//! Load a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::load_arma_ascii(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::streampos pos = f.tellg();
  
  bool load_okay = true;
  
  std::string f_header;
  uword       f_n_rows;
  uword       f_n_cols;
  
  f >> f_header;
  f >> f_n_rows;
  f >> f_n_cols;
  
  if(f_header == diskio::gen_txt_header(x))
    {
    try { x.zeros(f_n_rows, f_n_cols); } catch(...) { err_msg = "not enough memory"; return false; }
    
    std::string token;
    
    for(uword row=0; row < x.n_rows; ++row)
    for(uword col=0; col < x.n_cols; ++col)
      {
      f >> token;
      
      diskio::convert_token( x.at(row,col), token );
      }
    
    load_okay = f.good();
    }
  else
    {
    load_okay = false;
    err_msg = "incorrect header";
    }
  
  
  // allow automatic conversion of u32/s32 matrices into u64/s64 matrices
  
  if(load_okay == false)
    {
    if( (sizeof(eT) == 8) && is_same_type<uword,eT>::yes )
      {
      Mat<u32>    tmp;
      std::string junk;
      
      f.clear();
      f.seekg(pos);
      
      load_okay = diskio::load_arma_ascii(tmp, f, junk);
      
      if(load_okay)  { x = conv_to< Mat<eT> >::from(tmp); }
      }
    else
    if( (sizeof(eT) == 8) && is_same_type<sword,eT>::yes )
      {
      Mat<s32>    tmp;
      std::string junk;
      
      f.clear();
      f.seekg(pos);
      
      load_okay = diskio::load_arma_ascii(tmp, f, junk);
      
      if(load_okay)  { x = conv_to< Mat<eT> >::from(tmp); }
      }
    }
  
  return load_okay;
  }



//! Load a matrix in CSV text format (human readable)
template<typename eT>
inline
bool
diskio::load_csv_ascii(Mat<eT>& x, const std::string& name, std::string& err_msg, field<std::string>& header, const bool with_header, const char separator)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in);
  
  bool load_okay = f.is_open();
  
  if(load_okay == false)  { return false; }
  
  if(with_header)
    {
    arma_extra_debug_print("diskio::load_csv_ascii(): reading header");
    
    std::string              header_line;
    std::stringstream        header_stream;
    std::vector<std::string> header_tokens;
    
    std::getline(f, header_line);
    
    load_okay = f.good();
    
    if(load_okay)
      {
      std::string token;
      
      header_stream.clear();
      header_stream.str(header_line);
      
      uword header_n_tokens = 0;
      
      while(header_stream.good())
        {
        std::getline(header_stream, token, separator);
        ++header_n_tokens;
        header_tokens.push_back(token);
        }
      
      if(header_n_tokens == uword(0))
        {
        header.reset();
        }
      else
        {
        header.set_size(1,header_n_tokens);
        
        for(uword i=0; i < header_n_tokens; ++i)  { header.at(i) = header_tokens[i]; }
        }
      }
    }
  
  if(load_okay)
    {
    load_okay = diskio::load_csv_ascii(x, f, err_msg, separator);
    }
  
  f.close();
  
  return load_okay;
  }



//! Load a matrix in CSV text format (human readable)
template<typename eT>
inline
bool
diskio::load_csv_ascii(Mat<eT>& x, std::istream& f, std::string& err_msg, const char separator)
  {
  arma_extra_debug_sigprint();
  
  // TODO: replace with more efficient implementation
  
  if(f.good() == false)  { return false; }
  
  f.clear();
  const std::fstream::pos_type pos1 = f.tellg();
  
  //
  // work out the size
  
  uword f_n_rows = 0;
  uword f_n_cols = 0;
  
  std::string       line_string;
  std::stringstream line_stream;
  
  std::string token;
  
  while(f.good())
    {
    std::getline(f, line_string);
    
    if(line_string.size() == 0)  { break; }
    
    line_stream.clear();
    line_stream.str(line_string);
    
    uword line_n_cols = 0;
    
    while(line_stream.good())
      {
      std::getline(line_stream, token, separator);
      ++line_n_cols;
      }
    
    if(f_n_cols < line_n_cols)  { f_n_cols = line_n_cols; }
    
    ++f_n_rows;
    }
  
  f.clear();
  f.seekg(pos1);
  
  try { x.zeros(f_n_rows, f_n_cols); } catch(...) { err_msg = "not enough memory"; return false; }
  
  const bool use_mp = (arma_config::openmp) && (f_n_rows >= 2) && (f_n_cols >= 64);
  
  field<std::string> token_array;
  
  bool token_array_ok = false;
  
  if(use_mp)
    {
    try
      {
      token_array.set_size(f_n_cols);
      
      for(uword i=0; i < f_n_cols; ++i)  { token_array(i).reserve(32); }
      
      token_array_ok = true;
      }
    catch(...)
      {
      token_array.reset();
      }
    }
  
  if(use_mp && token_array_ok)
    {
    #if defined(ARMA_USE_OPENMP)
      {
      uword row = 0;
      
      while(f.good())
        {
        std::getline(f, line_string);
        
        if(line_string.size() == 0)  { break; }
        
        line_stream.clear();
        line_stream.str(line_string);
        
        for(uword i=0; i < f_n_cols; ++i)  { token_array(i).clear(); }
        
        uword line_stream_col = 0;
        
        while(line_stream.good())
          {
          std::getline(line_stream, token_array(line_stream_col), separator);
          
          ++line_stream_col;
          }
        
        const int n_threads = mp_thread_limit::get();
        
        #pragma omp parallel for schedule(static) num_threads(n_threads)
        for(uword col=0; col < line_stream_col; ++col)
          {
          diskio::convert_token( x.at(row,col), token_array(col) );
          }
        
        ++row;
        }
      }
    #endif
    }
  else  // serial implementation
    {
    uword row = 0;
    
    while(f.good())
      {
      std::getline(f, line_string);
      
      if(line_string.size() == 0)  { break; }
      
      line_stream.clear();
      line_stream.str(line_string);
      
      uword col = 0;
      
      while(line_stream.good())
        {
        std::getline(line_stream, token, separator);
        
        diskio::convert_token( x.at(row,col), token );
        
        ++col;
        }
      
      ++row;
      }
    }
  
  return true;
  }



//! Load a matrix in CSV text format (human readable); complex numbers stored in "a+bi" format
template<typename T>
inline
bool
diskio::load_csv_ascii(Mat< std::complex<T> >& x, std::istream& f, std::string& err_msg, const char separator)
  {
  arma_extra_debug_sigprint();
  
  // TODO: replace with more efficient implementation
  
  if(f.good() == false)  { return false; }
  
  f.clear();
  const std::fstream::pos_type pos1 = f.tellg();
  
  //
  // work out the size
  
  uword f_n_rows = 0;
  uword f_n_cols = 0;
  
  std::string       line_string;
  std::stringstream line_stream;
  
  std::string token;
  
  while(f.good())
    {
    std::getline(f, line_string);
    
    if(line_string.size() == 0)  { break; }
    
    line_stream.clear();
    line_stream.str(line_string);
    
    uword line_n_cols = 0;
    
    while(line_stream.good())
      {
      std::getline(line_stream, token, separator);
      ++line_n_cols;
      }
    
    if(f_n_cols < line_n_cols)  { f_n_cols = line_n_cols; }
    
    ++f_n_rows;
    }
  
  f.clear();
  f.seekg(pos1);
  
  try { x.zeros(f_n_rows, f_n_cols); } catch(...) { err_msg = "not enough memory"; return false; }
  
  uword row = 0;
  
  std::string str_real;
  std::string str_imag;
  
  while(f.good())
    {
    std::getline(f, line_string);
    
    if(line_string.size() == 0)  { break; }
    
    line_stream.clear();
    line_stream.str(line_string);
    
    uword col = 0;
    
    while(line_stream.good())
      {
      std::getline(line_stream, token, separator);
      
      // remove spaces and tabs
      if(token.length() > 0)
        {
        const char c_front = token.front();
        const char c_back  = token.back();
        
        if( (c_front == ' ') || (c_front == '\t') || (c_back == ' ') || (c_back == '\t') )
          {
          token.erase(std::remove_if(token.begin(), token.end(), [](char c) { return ((c == ' ') || (c == '\t')); }), token.end());
          }
        }
      
      const size_t token_len = size_t( token.length() );
      
      if(token_len == 0)  { col++; continue; }
      
      // handle special cases: inf and nan, without the imaginary part
      if( (token_len == 3) || (token_len == 4) )
        {
        const char* str = token.c_str();
        
        const bool neg = (str[0] == '-');
        const bool pos = (str[0] == '+');
        
        const size_t offset = ( (neg || pos) && (token_len == 4) ) ? 1 : 0;
        
        const char sig_a = str[offset  ];
        const char sig_b = str[offset+1];
        const char sig_c = str[offset+2];
        
        bool found_val_real = false;
        T          val_real = T(0);
        
        if( ((sig_a == 'i') || (sig_a == 'I')) && ((sig_b == 'n') || (sig_b == 'N')) && ((sig_c == 'f') || (sig_c == 'F')) )
          {
          val_real = (neg) ? -(Datum<T>::inf) : Datum<T>::inf;
          
          found_val_real = true;
          }
        else
        if( ((sig_a == 'n') || (sig_a == 'N')) && ((sig_b == 'a') || (sig_b == 'A')) && ((sig_c == 'n') || (sig_c == 'N')) )
          {
          val_real = Datum<T>::nan;
          
          found_val_real = true;
          }
        
        if(found_val_real)
          {
          x.at(row,col) = std::complex<T>(val_real, T(0));
          
          col++; continue;  // get next token
          }
        }
      
      bool found_x = false;
      std::string::size_type loc_x = 0;  // location of the separator (+ or -) between the real and imaginary part
      
      std::string::size_type loc_i = token.find_last_of('i');  // location of the imaginary part indicator
      
      if(loc_i == std::string::npos)
        {
        str_real = token;
        str_imag.clear();
        }
      else
        {
        bool found_plus  = false;
        bool found_minus = false;
        
        std::string::size_type loc_plus = token.find_last_of('+');
        
        if(loc_plus != std::string::npos)
          {
          if(loc_plus >= 1)
            {
            const char prev_char = token.at(loc_plus-1);
            
            // make sure we're not looking at the sign of the exponent
            if( (prev_char != 'e') && (prev_char != 'E') )
              {
              found_plus = true;
              }
            else
              {
              // search again, omitting the exponent
              loc_plus = token.find_last_of('+', loc_plus-1);
              
              if(loc_plus != std::string::npos)  { found_plus = true; }
              }
            }
          else
            {
            // loc_plus == 0, meaning we're at the start of the string
            found_plus = true;
            }
          }
        
        std::string::size_type loc_minus = token.find_last_of('-');
        
        if(loc_minus != std::string::npos)
          {
          if(loc_minus >= 1)
            {
            const char prev_char = token.at(loc_minus-1);
            
            // make sure we're not looking at the sign of the exponent
            if( (prev_char != 'e') && (prev_char != 'E') )
              {
              found_minus = true;
              }
            else
              {
              // search again, omitting the exponent
              loc_minus = token.find_last_of('-', loc_minus-1);
              
              if(loc_minus != std::string::npos)  { found_minus = true; }
              }
            }
          else
            {
            // loc_minus == 0, meaning we're at the start of the string
            found_minus = true;
            }
          }
        
        if(found_plus && found_minus)
          {
          if( (loc_i > loc_plus) && (loc_i > loc_minus) )
            {
            // choose the sign closest to the "i" to be the separator between the real and imaginary part
            loc_x = ( (loc_i - loc_plus) < (loc_i - loc_minus) ) ? loc_plus : loc_minus;
            found_x = true;
            }
          }
        else if(found_plus )  { loc_x = loc_plus;  found_x = true; }
        else if(found_minus)  { loc_x = loc_minus; found_x = true; }
        
        if(found_x)
          {
          if( loc_x    > 0           ) { str_real = token.substr(0,loc_x);                     } else { str_real.clear(); }
          if((loc_x+1) < token.size()) { str_imag = token.substr(loc_x, token.size()-loc_x-1); } else { str_imag.clear(); }
          }
        else
          {
          str_real.clear();
          str_imag.clear();
          }
        }
      
      T val_real = T(0);
      T val_imag = T(0);
      
      diskio::convert_token(val_real, str_real);
      diskio::convert_token(val_imag, str_imag);
      
      x.at(row,col) = std::complex<T>(val_real, val_imag);
      
      ++col;
      }
    
    ++row;
    }
  
  return true;
  }



template<typename eT>
inline
bool
diskio::load_coord_ascii(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in);
  
  bool load_okay = f.is_open();
  
  if(load_okay == false)  { return false; }
  
  if(load_okay)
    {
    load_okay = diskio::load_coord_ascii(x, f, err_msg);
    }
  
  f.close();
  
  return load_okay;
  }



//! Load a matrix in CSV text format (human readable)
template<typename eT>
inline
bool
diskio::load_coord_ascii(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  if(f.good() == false)  { return false; }
  
  f.clear();
  const std::fstream::pos_type pos1 = f.tellg();
  
  // work out the size
  
  uword f_n_rows = 0;
  uword f_n_cols = 0;
  
  bool size_found = false;
  
  std::string       line_string;
  std::stringstream line_stream;
  
  std::string token;
  
  while(f.good())
    {
    std::getline(f, line_string);
    
    if(line_string.size() == 0)  { break; }
    
    line_stream.clear();
    line_stream.str(line_string);
    
    uword line_row = 0;
    uword line_col = 0;
    
    // a valid line in co-ord format has at least 2 entries
    
    line_stream >> line_row;
    
    if(line_stream.good() == false)  { err_msg = "incorrect format"; return false; }
    
    line_stream >> line_col;
    
    size_found = true;
    
    if(f_n_rows < line_row)  { f_n_rows = line_row; }
    if(f_n_cols < line_col)  { f_n_cols = line_col; }
    }
  
  // take into account that indices start at 0
  if(size_found)  { ++f_n_rows;  ++f_n_cols; }
  
  f.clear();
  f.seekg(pos1);
  
  try
    {
    Mat<eT> tmp(f_n_rows, f_n_cols, arma_zeros_indicator());
    
    while(f.good())
      {
      std::getline(f, line_string);
      
      if(line_string.size() == 0)  { break; }
      
      line_stream.clear();
      line_stream.str(line_string);
      
      uword line_row = 0;
      uword line_col = 0;
      
      line_stream >> line_row;
      line_stream >> line_col;
      
      eT val = eT(0);
      
      line_stream >> token;
      
      if(line_stream.fail() == false)
        {
        diskio::convert_token( val, token );
        }
      
      if(val != eT(0))  { tmp(line_row,line_col) = val; }
      }
    
    x.steal_mem(tmp);
    }
  catch(...)
    {
    err_msg = "not enough memory";
    return false;
    }
  
  return true;
  }



template<typename T>
inline
bool
diskio::load_coord_ascii(Mat< std::complex<T> >& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  if(f.good() == false)  { return false; }
  
  f.clear();
  const std::fstream::pos_type pos1 = f.tellg();
  
  // work out the size
  
  uword f_n_rows = 0;
  uword f_n_cols = 0;
  
  bool size_found = false;
  
  std::string       line_string;
  std::stringstream line_stream;
  
  std::string token_real;
  std::string token_imag;
  
  while(f.good())
    {
    std::getline(f, line_string);
    
    if(line_string.size() == 0)  { break; }
    
    line_stream.clear();
    line_stream.str(line_string);
    
    uword line_row = 0;
    uword line_col = 0;
    
    // a valid line in co-ord format has at least 2 entries
    
    line_stream >> line_row;
    
    if(line_stream.good() == false)  { err_msg = "incorrect format"; return false; }
    
    line_stream >> line_col;
    
    size_found = true;
    
    if(f_n_rows < line_row)  f_n_rows = line_row;
    if(f_n_cols < line_col)  f_n_cols = line_col;
    }
  
  // take into account that indices start at 0
  if(size_found)  { ++f_n_rows;  ++f_n_cols; }
  
  f.clear();
  f.seekg(pos1);
  
  try
    {
    Mat< std::complex<T> > tmp(f_n_rows, f_n_cols, arma_zeros_indicator());
    
    while(f.good())
      {
      std::getline(f, line_string);
      
      if(line_string.size() == 0)  { break; }
      
      line_stream.clear();
      line_stream.str(line_string);
      
      uword line_row = 0;
      uword line_col = 0;
      
      line_stream >> line_row;
      line_stream >> line_col;
      
      T val_real = T(0);
      T val_imag = T(0);
      
      line_stream >> token_real;
      
      if(line_stream.fail() == false)
        {
        diskio::convert_token( val_real, token_real );
        }
      
      
      line_stream >> token_imag;
      
      if(line_stream.fail() == false)
        {
        diskio::convert_token( val_imag, token_imag );
        }
      
      if( (val_real != T(0)) || (val_imag != T(0)) )
        {
        tmp(line_row,line_col) = std::complex<T>(val_real, val_imag);
        }
      }
    
    x.steal_mem(tmp);
    }
  catch(...)
    {
    err_msg = "not enough memory";
    return false;
    }
  
  return true;
  }



//! Load a matrix in binary format,
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::load_arma_binary(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f;
  f.open(name.c_str(), std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_arma_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::load_arma_binary(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::streampos pos = f.tellg();
  
  bool load_okay = true;
  
  std::string f_header;
  uword       f_n_rows;
  uword       f_n_cols;
  
  f >> f_header;
  f >> f_n_rows;
  f >> f_n_cols;
  
  if(f_header == diskio::gen_bin_header(x))
    {
    //f.seekg(1, ios::cur);  // NOTE: this may not be portable, as on a Windows machine a newline could be two characters
    f.get();
    
    try { x.set_size(f_n_rows,f_n_cols); } catch(...) { err_msg = "not enough memory"; return false; }
    
    f.read( reinterpret_cast<char *>(x.memptr()), std::streamsize(x.n_elem*sizeof(eT)) );
    
    load_okay = f.good();
    }
  else
    {
    load_okay = false;
    err_msg = "incorrect header";
    }
  
  
  // allow automatic conversion of u32/s32 matrices into u64/s64 matrices
  
  if(load_okay == false)
    {
    if( (sizeof(eT) == 8) && is_same_type<uword,eT>::yes )
      {
      Mat<u32>    tmp;
      std::string junk;
      
      f.clear();
      f.seekg(pos);
      
      load_okay = diskio::load_arma_binary(tmp, f, junk);
      
      if(load_okay)  { x = conv_to< Mat<eT> >::from(tmp); }
      }
    else
    if( (sizeof(eT) == 8) && is_same_type<sword,eT>::yes )
      {
      Mat<s32>    tmp;
      std::string junk;
      
      f.clear();
      f.seekg(pos);
      
      load_okay = diskio::load_arma_binary(tmp, f, junk);
      
      if(load_okay)  { x = conv_to< Mat<eT> >::from(tmp); }
      }
    }
  
  return load_okay;
  }



inline
void
diskio::pnm_skip_comments(std::istream& f)
  {
  while( isspace(f.peek()) )
    {
    while( isspace(f.peek()) )  { f.get(); }
    
    if(f.peek() == '#')
      {
      while( (f.peek() != '\r') && (f.peek() != '\n') )  { f.get(); }
      }
    }
  }



//! Load a PGM greyscale image as a matrix
template<typename eT>
inline
bool
diskio::load_pgm_binary(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_pgm_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



//! Load a PGM greyscale image as a matrix
template<typename eT>
inline
bool
diskio::load_pgm_binary(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  bool load_okay = true;
  
  std::string f_header;
  
  f >> f_header;
  
  if(f_header == "P5")
    {
    uword f_n_rows = 0;
    uword f_n_cols = 0;
    int f_maxval = 0;
    
    diskio::pnm_skip_comments(f);
    
    f >> f_n_cols;
    diskio::pnm_skip_comments(f);
    
    f >> f_n_rows;
    diskio::pnm_skip_comments(f);
    
    f >> f_maxval;
    f.get();
    
    if( (f_maxval > 0) && (f_maxval <= 65535) )
      {
      try { x.set_size(f_n_rows,f_n_cols); } catch(...) { err_msg = "not enough memory"; return false; }
      
      if(f_maxval <= 255)
        {
        const uword n_elem = f_n_cols*f_n_rows;
        podarray<u8> tmp(n_elem);
        
        f.read( reinterpret_cast<char*>(tmp.memptr()), std::streamsize(n_elem) );
        
        uword i = 0;
        
        //cout << "f_n_cols = " << f_n_cols << endl;
        //cout << "f_n_rows = " << f_n_rows << endl;
        
        for(uword row=0; row < f_n_rows; ++row)
        for(uword col=0; col < f_n_cols; ++col)
          {
          x.at(row,col) = eT(tmp[i]);
          ++i;
          }
        }
      else
        {
        const uword n_elem = f_n_cols*f_n_rows;
        podarray<u16> tmp(n_elem);
        
        f.read( reinterpret_cast<char *>(tmp.memptr()), std::streamsize(n_elem*2) );
        
        uword i = 0;
        
        for(uword row=0; row < f_n_rows; ++row)
        for(uword col=0; col < f_n_cols; ++col)
          {
          x.at(row,col) = eT(tmp[i]);
          ++i;
          }
        }
      }
    else
      {
      load_okay = false;
      err_msg = "functionality unimplemented";
      }
    
    if(f.good() == false)  { load_okay = false; }
    }
  else
    {
    load_okay = false;
    err_msg = "unsupported header";
    }
  
  return load_okay;
  }



//! Load a PGM greyscale image as a matrix
template<typename T>
inline
bool
diskio::load_pgm_binary(Mat< std::complex<T> >& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  uchar_mat tmp;
  const bool load_okay = diskio::load_pgm_binary(tmp, name, err_msg);
  
  x = conv_to< Mat< std::complex<T> > >::from(tmp);
  
  return load_okay;
  }



//! Load a PGM greyscale image as a matrix
template<typename T>
inline
bool
diskio::load_pgm_binary(Mat< std::complex<T> >& x, std::istream& is, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  uchar_mat tmp;
  const bool load_okay = diskio::load_pgm_binary(tmp, is, err_msg);
  
  x = conv_to< Mat< std::complex<T> > >::from(tmp);
  
  return load_okay;
  }



//! Load a HDF5 file as a matrix
template<typename eT>
inline
bool
diskio::load_hdf5_binary(Mat<eT>& x, const hdf5_name& spec, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_HDF5)
    {
    hdf5_misc::hdf5_suspend_printing_errors hdf5_print_suspender;
    
    bool load_okay = false;
    
    hid_t fid = arma_H5Fopen(spec.filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    
    if(fid >= 0)
      {
      // MATLAB HDF5 dataset names are user-specified;
      // Octave tends to store the datasets in a group, with the actual dataset being referred to as "value".
      // If the user hasn't specified a dataset, we will search for "dataset" and "value",
      // and if those are not found we will take the first dataset we do find.
      
      std::vector<std::string> searchNames;
      
      const bool exact = (spec.dsname.empty() == false);
      
      if(exact)
        {
        searchNames.push_back(spec.dsname);
        }
      else
        {
        searchNames.push_back("dataset");
        searchNames.push_back("value"  );
        }
      
      hid_t dataset = hdf5_misc::search_hdf5_file(searchNames, fid, 2, exact);
      
      if(dataset >= 0)
        {
        hid_t filespace = arma_H5Dget_space(dataset);
        
        // This must be <= 2 due to our search rules.
        const int ndims = arma_H5Sget_simple_extent_ndims(filespace);
        
        hsize_t dims[2];
        const herr_t query_status = arma_H5Sget_simple_extent_dims(filespace, dims, NULL);
        
        // arma_check(query_status < 0, "Mat::load(): cannot get size of HDF5 dataset");
        if(query_status < 0)
          {
          err_msg = "cannot get size of HDF5 dataset";
          
          arma_H5Sclose(filespace);
          arma_H5Dclose(dataset);
          arma_H5Fclose(fid);
          
          return false;
          }
        
        if(ndims == 1) { dims[1] = 1; }  // Vector case; fake second dimension (one column).
        
        try { x.set_size(dims[1], dims[0]); } catch(...) { err_msg = "not enough memory"; return false; }
        
        // Now we have to see what type is stored to figure out how to load it.
        hid_t datatype = arma_H5Dget_type(dataset);
        hid_t mat_type = hdf5_misc::get_hdf5_type<eT>();
        
        // If these are the same type, it is simple.
        if(arma_H5Tequal(datatype, mat_type) > 0)
          {
          // Load directly; H5S_ALL used so that we load the entire dataset.
          hid_t read_status = arma_H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, void_ptr(x.memptr()));
          
          if(read_status >= 0) { load_okay = true; }
          }
        else
          {
          // Load into another array and convert its type accordingly.
          hid_t read_status = hdf5_misc::load_and_convert_hdf5(x.memptr(), dataset, datatype, x.n_elem);
          
          if(read_status >= 0) { load_okay = true; }
          }
        
        // Now clean up.
        arma_H5Tclose(datatype);
        arma_H5Tclose(mat_type);
        arma_H5Sclose(filespace);
        }
      
      arma_H5Dclose(dataset);
      
      arma_H5Fclose(fid);
      
      if(load_okay == false)
        {
        err_msg = "unsupported or missing HDF5 data";
        }
      }
    else
      {
      err_msg = "cannot open";
      }
    
    return load_okay;
    }
  #else
    {
    arma_ignore(x);
    arma_ignore(spec);
    arma_ignore(err_msg);
    
    arma_stop_logic_error("Mat::load(): use of HDF5 must be enabled");
    
    return false;
    }
  #endif
  }



//! Try to load a matrix by automatically determining its type
template<typename eT>
inline
bool
diskio::load_auto_detect(Mat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_HDF5)
    // We're currently using the C bindings for the HDF5 library, which don't support C++ streams
    if( arma_H5Fis_hdf5(name.c_str()) ) { return load_hdf5_binary(x, name, err_msg); }
  #endif
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_auto_detect(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



//! Try to load a matrix by automatically determining its type
template<typename eT>
inline
bool
diskio::load_auto_detect(Mat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  const char* ARMA_MAT_TXT_str = "ARMA_MAT_TXT";
  const char* ARMA_MAT_BIN_str = "ARMA_MAT_BIN";
  const char*           P5_str = "P5";
  
  const uword ARMA_MAT_TXT_len = uword(12);
  const uword ARMA_MAT_BIN_len = uword(12);
  const uword           P5_len = uword(2);
  
  podarray<char> header(ARMA_MAT_TXT_len + 1);
  
  char* header_mem = header.memptr();
  
  std::streampos pos = f.tellg();
    
  f.read( header_mem, std::streamsize(ARMA_MAT_TXT_len) );
  f.clear();
  f.seekg(pos);
  
  header_mem[ARMA_MAT_TXT_len] = '\0';
  
  if( std::strncmp(ARMA_MAT_TXT_str, header_mem, size_t(ARMA_MAT_TXT_len)) == 0 )
    {
    return load_arma_ascii(x, f, err_msg);
    }
  else
  if( std::strncmp(ARMA_MAT_BIN_str, header_mem, size_t(ARMA_MAT_BIN_len)) == 0 )
    {
    return load_arma_binary(x, f, err_msg);
    }
  else
  if( std::strncmp(P5_str, header_mem, size_t(P5_len)) == 0 )
    {
    return load_pgm_binary(x, f, err_msg);
    }
  else
    {
    const file_type ft = guess_file_type_internal(f);
    
    switch(ft)
      {
      case csv_ascii:
        return load_csv_ascii(x, f, err_msg, char(','));
        break;
      
      case ssv_ascii:
        return load_csv_ascii(x, f, err_msg, char(';'));
        break;
      
      case raw_binary:
        return load_raw_binary(x, f, err_msg);
        break;
        
      case raw_ascii:
        return load_raw_ascii(x, f, err_msg);
        break;
      
      default:
        err_msg = "unknown data";
        return false;
      }
    }
  
  return false;
  }



//
// sparse matrices
//



//! Save a sparse matrix in CSV format
template<typename eT>
inline
bool
diskio::save_csv_ascii(const SpMat<eT>& x, const std::string& final_name, const field<std::string>& header, const bool with_header, const char separator)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str());
  
  bool save_okay = f.is_open();
  
  if(save_okay == false)  { return false; }
  
  if(with_header)
    {
    arma_extra_debug_print("diskio::save_csv_ascii(): writing header");
    
    for(uword i=0; i < header.n_elem; ++i)
      {
      f << header(i);
      
      if(i != (header.n_elem-1))  { f.put(separator); }
      }
    
    f.put('\n');
    
    save_okay = f.good();
    }
  
  if(save_okay)  { save_okay = diskio::save_csv_ascii(x, f, separator); }
  
  f.flush();
  f.close();
  
  if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
  
  return save_okay;
  }



//! Save a sparse matrix in CSV format
template<typename eT>
inline
bool
diskio::save_csv_ascii(const SpMat<eT>& x, std::ostream& f, const char separator)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(f);
  
  diskio::prepare_stream<eT>(f);
  
  x.sync();
  
  uword x_n_rows = x.n_rows;
  uword x_n_cols = x.n_cols;
  
  for(uword row=0; row < x_n_rows; ++row)
    {
    for(uword col=0; col < x_n_cols; ++col)
      {
      const eT val = x.at(row,col);
      
      if(val != eT(0))  { arma_ostream::raw_print_elem(f, val); }
      
      if( col < (x_n_cols-1) )  { f.put(separator); }
      }
    
    f.put('\n');
    }
  
  const bool save_okay = f.good();
  
  stream_state.restore(f);
  
  return save_okay;
  }



//! Save a sparse matrix in CSV format (complex numbers)
template<typename T>
inline
bool
diskio::save_csv_ascii(const SpMat< std::complex<T> >& x, std::ostream& f, const char separator)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(x);
  arma_ignore(f);
  arma_ignore(separator);
  
  arma_debug_warn_level(1, "saving complex sparse matrices as csv_ascii not yet implemented");
  
  return false;
  }



//! Save a matrix in ASCII coord format
template<typename eT>
inline
bool
diskio::save_coord_ascii(const SpMat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str());
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_coord_ascii(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



//! Save a matrix in ASCII coord format
template<typename eT>
inline
bool
diskio::save_coord_ascii(const SpMat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(f);
  
  diskio::prepare_stream<eT>(f);
  
  typename SpMat<eT>::const_iterator iter     = x.begin();
  typename SpMat<eT>::const_iterator iter_end = x.end();
  
  for(; iter != iter_end; ++iter)
    {
    const eT val = (*iter);
    
    f << iter.row() << ' ' << iter.col() << ' ' << val << '\n';
    }
  
  
  // make sure it's possible to figure out the matrix size later
  if( (x.n_rows > 0) && (x.n_cols > 0) )
    {
    const uword max_row = (x.n_rows > 0) ? x.n_rows-1 : 0;
    const uword max_col = (x.n_cols > 0) ? x.n_cols-1 : 0;
    
    if( x.at(max_row, max_col) == eT(0) )
      {
      f << max_row << ' ' << max_col << " 0\n";
      }
    }
  
  const bool save_okay = f.good();
  
  stream_state.restore(f);
  
  return save_okay;
  }



//! Save a matrix in ASCII coord format (complex numbers)
template<typename T>
inline
bool
diskio::save_coord_ascii(const SpMat< std::complex<T> >& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  const arma_ostream_state stream_state(f);
  
  diskio::prepare_stream<eT>(f);
  
  typename SpMat<eT>::const_iterator iter     = x.begin();
  typename SpMat<eT>::const_iterator iter_end = x.end();
  
  for(; iter != iter_end; ++iter)
    {
    const eT val = (*iter);
    
    f << iter.row() << ' ' << iter.col() << ' ' << val.real() << ' ' << val.imag() << '\n';
    }
  
  // make sure it's possible to figure out the matrix size later
  if( (x.n_rows > 0) && (x.n_cols > 0) )
    {
    const uword max_row = (x.n_rows > 0) ? x.n_rows-1 : 0;
    const uword max_col = (x.n_cols > 0) ? x.n_cols-1 : 0;
    
    if( x.at(max_row, max_col) == eT(0) )
      {
      f << max_row << ' ' << max_col << " 0 0\n";
      }
    }
  
  const bool save_okay = f.good();
  
  stream_state.restore(f);
  
  return save_okay;
  }



//! Save a matrix in binary format,
//! with a header that stores the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_binary(const SpMat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str(), std::fstream::binary);
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_arma_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



//! Save a matrix in binary format,
//! with a header that stores the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_binary(const SpMat<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  f << diskio::gen_bin_header(x) << '\n';
  f << x.n_rows << ' ' << x.n_cols << ' ' << x.n_nonzero << '\n';
  
  f.write( reinterpret_cast<const char*>(x.values),      std::streamsize(x.n_nonzero*sizeof(eT))     );
  f.write( reinterpret_cast<const char*>(x.row_indices), std::streamsize(x.n_nonzero*sizeof(uword))  );
  f.write( reinterpret_cast<const char*>(x.col_ptrs),    std::streamsize((x.n_cols+1)*sizeof(uword)) );
  
  return f.good();
  }



template<typename eT>
inline
bool
diskio::load_csv_ascii(SpMat<eT>& x, const std::string& name, std::string& err_msg, field<std::string>& header, const bool with_header, const char separator)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in);
  
  bool load_okay = f.is_open();
  
  if(load_okay == false)  { return false; }
  
  if(with_header)
    {
    arma_extra_debug_print("diskio::load_csv_ascii(): reading header");
    
    std::string              header_line;
    std::stringstream        header_stream;
    std::vector<std::string> header_tokens;
    
    std::getline(f, header_line);
    
    load_okay = f.good();
    
    if(load_okay)
      {
      std::string token;
      
      header_stream.clear();
      header_stream.str(header_line);
      
      uword header_n_tokens = 0;
      
      while(header_stream.good())
        {
        std::getline(header_stream, token, separator);
        ++header_n_tokens;
        header_tokens.push_back(token);
        }
      
      if(header_n_tokens == uword(0))
        {
        header.reset();
        }
      else
        {
        header.set_size(1,header_n_tokens);
        
        for(uword i=0; i < header_n_tokens; ++i)  { header.at(i) = header_tokens[i]; }
        }
      }
    }
  
  if(load_okay)
    {
    load_okay = diskio::load_csv_ascii(x, f, err_msg, separator);
    }
  
  f.close();
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::load_csv_ascii(SpMat<eT>& x, std::istream& f, std::string& err_msg, const char separator)
  {
  arma_extra_debug_sigprint();
  
  // TODO: replace with more efficient implementation
  
  if(f.good() == false)  { return false; }
  
  f.clear();
  const std::fstream::pos_type pos1 = f.tellg();
  
  //
  // work out the size
  
  uword f_n_rows = 0;
  uword f_n_cols = 0;
  
  std::string       line_string;
  std::stringstream line_stream;
  
  std::string token;
  
  while(f.good())
    {
    std::getline(f, line_string);
    
    if(line_string.size() == 0)  { break; }
    
    line_stream.clear();
    line_stream.str(line_string);
    
    uword line_n_cols = 0;
    
    while(line_stream.good())
      {
      std::getline(line_stream, token, separator);
      ++line_n_cols;
      }
    
    if(f_n_cols < line_n_cols)  { f_n_cols = line_n_cols; }
    
    ++f_n_rows;
    }
  
  f.clear();
  f.seekg(pos1);
  
  try
    {
    MapMat<eT> tmp(f_n_rows, f_n_cols);
    
    uword row = 0;
    
    while(f.good())
      {
      std::getline(f, line_string);
      
      if(line_string.size() == 0)  { break; }
      
      line_stream.clear();
      line_stream.str(line_string);
      
      uword col = 0;
      
      while(line_stream.good())
        {
        std::getline(line_stream, token, separator);
        
        eT val = eT(0);
        
        diskio::convert_token( val, token );
        
        if(val != eT(0))  { tmp(row,col) = val; }
        
        ++col;
        }
      
      ++row;
      }
    
    x = tmp;
    }
  catch(...)
    {
    err_msg = "not enough memory";
    return false;
    }
  
  return true;
  }



template<typename T>
inline
bool
diskio::load_csv_ascii(SpMat< std::complex<T> >& x, std::istream& f, std::string& err_msg, const char separator)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(x);
  arma_ignore(f);
  arma_ignore(err_msg);
  arma_ignore(separator);
  
  arma_debug_warn_level(1, "loading complex sparse matrices as csv_ascii not yet implemented");
  
  return false;
  }



template<typename eT>
inline
bool
diskio::load_coord_ascii(SpMat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_coord_ascii(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::load_coord_ascii(SpMat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  if(f.good() == false)  { return false; }
  
  f.clear();
  const std::fstream::pos_type pos1 = f.tellg();
  
  // work out the size
  
  uword f_n_rows = 0;
  uword f_n_cols = 0;
  
  bool size_found = false;
  
  std::string       line_string;
  std::stringstream line_stream;
  
  std::string token;
  
  while(f.good())
    {
    std::getline(f, line_string);
    
    if(line_string.size() == 0)  { break; }
    
    line_stream.clear();
    line_stream.str(line_string);
    
    uword line_row = 0;
    uword line_col = 0;
    
    // a valid line in co-ord format has at least 2 entries
    
    line_stream >> line_row;
    
    if(line_stream.good() == false)  { err_msg = "incorrect format"; return false; }
    
    line_stream >> line_col;
    
    size_found = true;
    
    if(f_n_rows < line_row)  { f_n_rows = line_row; }
    if(f_n_cols < line_col)  { f_n_cols = line_col; }
    }
  
  // take into account that indices start at 0
  if(size_found)  { ++f_n_rows;  ++f_n_cols; }
  
  f.clear();
  f.seekg(pos1);
  
  try
    {
    MapMat<eT> tmp(f_n_rows, f_n_cols);
    
    while(f.good())
      {
      std::getline(f, line_string);
      
      if(line_string.size() == 0)  { break; }
      
      line_stream.clear();
      line_stream.str(line_string);
      
      uword line_row = 0;
      uword line_col = 0;
      
      line_stream >> line_row;
      line_stream >> line_col;
      
      eT val = eT(0);
      
      line_stream >> token;
      
      if(line_stream.fail() == false)
        {
        diskio::convert_token( val, token );
        }
      
      if(val != eT(0))  { tmp(line_row,line_col) = val; }
      }
    
    x = tmp;
    }
  catch(...)
    {
    err_msg = "not enough memory";
    return false;
    }
  
  return true;
  }



template<typename T>
inline
bool
diskio::load_coord_ascii(SpMat< std::complex<T> >& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  if(f.good() == false)  { return false; }
  
  f.clear();
  const std::fstream::pos_type pos1 = f.tellg();
  
  // work out the size
  
  uword f_n_rows = 0;
  uword f_n_cols = 0;
  
  bool size_found = false;
  
  std::string       line_string;
  std::stringstream line_stream;
  
  std::string token_real;
  std::string token_imag;
  
  while(f.good())
    {
    std::getline(f, line_string);
    
    if(line_string.size() == 0)  { break; }
    
    line_stream.clear();
    line_stream.str(line_string);
    
    uword line_row = 0;
    uword line_col = 0;
    
    // a valid line in co-ord format has at least 2 entries
    
    line_stream >> line_row;
    
    if(line_stream.good() == false)  { err_msg = "incorrect format"; return false; }
    
    line_stream >> line_col;
    
    size_found = true;
    
    if(f_n_rows < line_row)  f_n_rows = line_row;
    if(f_n_cols < line_col)  f_n_cols = line_col;
    }
  
  // take into account that indices start at 0
  if(size_found)  { ++f_n_rows;  ++f_n_cols; }
  
  f.clear();
  f.seekg(pos1);
  
  try
    {
    MapMat< std::complex<T> > tmp(f_n_rows, f_n_cols);
    
    while(f.good())
      {
      std::getline(f, line_string);
      
      if(line_string.size() == 0)  { break; }
      
      line_stream.clear();
      line_stream.str(line_string);
      
      uword line_row = 0;
      uword line_col = 0;
      
      line_stream >> line_row;
      line_stream >> line_col;
      
      T val_real = T(0);
      T val_imag = T(0);
      
      line_stream >> token_real;
      
      if(line_stream.fail() == false)
        {
        diskio::convert_token( val_real, token_real );
        }
      
      
      line_stream >> token_imag;
      
      if(line_stream.fail() == false)
        {
        diskio::convert_token( val_imag, token_imag );
        }
      
      if( (val_real != T(0)) || (val_imag != T(0)) )
        {
        tmp(line_row,line_col) = std::complex<T>(val_real, val_imag);
        }
      }
    
    x = tmp;
    }
  catch(...)
    {
    err_msg = "not enough memory";
    return false;
    }
  
  return true;
  }



//! Load a matrix in binary format,
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
bool
diskio::load_arma_binary(SpMat<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f;
  f.open(name.c_str(), std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_arma_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::load_arma_binary(SpMat<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  std::string f_header;
  
  f >> f_header;
  
  if(f_header == diskio::gen_bin_header(x))
    {
    uword f_n_rows;
    uword f_n_cols;
    uword f_n_nz;
    
    f >> f_n_rows;
    f >> f_n_cols;
    f >> f_n_nz;
    
    //f.seekg(1, ios::cur);  // NOTE: this may not be portable, as on a Windows machine a newline could be two characters
    f.get();
    
    try { x.reserve(f_n_rows, f_n_cols, f_n_nz); } catch(...) { err_msg = "not enough memory"; return false; }
    
    f.read( reinterpret_cast<char*>(access::rwp(x.values)),      std::streamsize(x.n_nonzero*sizeof(eT))     );
    
    std::streampos pos = f.tellg();
    
    f.read( reinterpret_cast<char*>(access::rwp(x.row_indices)), std::streamsize(x.n_nonzero*sizeof(uword))  );
    f.read( reinterpret_cast<char*>(access::rwp(x.col_ptrs)),    std::streamsize((x.n_cols+1)*sizeof(uword)) );
    
    bool check1 = true;  for(uword i=0; i < x.n_nonzero; ++i)  { if(x.values[i] == eT(0))  { check1 = false; break; } }
    bool check2 = true;  for(uword i=0; i < x.n_cols;    ++i)  { if(x.col_ptrs[i+1] < x.col_ptrs[i])  { check2 = false; break; } }
    bool check3 = (x.col_ptrs[x.n_cols] == x.n_nonzero);
    
    if((check1 == true) && ((check2 == false) || (check3 == false)))
      {
      if(sizeof(uword) == 8)
        {
        arma_extra_debug_print("detected inconsistent data while loading; re-reading integer parts as u32");
        
        // inconstency could be due to a different uword size used during saving,
        // so try loading the row_indices and col_ptrs under the assumption of 32 bit unsigned integers
        
        f.clear();
        f.seekg(pos);
        
        podarray<u32> tmp_a(x.n_nonzero );  tmp_a.zeros();
        podarray<u32> tmp_b(x.n_cols + 1);  tmp_b.zeros();
        
        f.read( reinterpret_cast<char*>(tmp_a.memptr()), std::streamsize( x.n_nonzero   * sizeof(u32)) );
        f.read( reinterpret_cast<char*>(tmp_b.memptr()), std::streamsize((x.n_cols + 1) * sizeof(u32)) );
        
        check2 = true;  for(uword i=0; i < x.n_cols; ++i)  { if(tmp_b[i+1] < tmp_b[i])  { check2 = false; break; } }
        check3 = (tmp_b[x.n_cols] == x.n_nonzero);
        
        load_okay = f.good();
        
        if( load_okay && (check2 == true) && (check3 == true) )
          {
          arma_extra_debug_print("reading integer parts as u32 succeeded");
          
          arrayops::convert(access::rwp(x.row_indices), tmp_a.memptr(), x.n_nonzero );
          arrayops::convert(access::rwp(x.col_ptrs),    tmp_b.memptr(), x.n_cols + 1);
          }
        else
          {
          arma_extra_debug_print("reading integer parts as u32 failed");
          }
        }
      }
    
    if((check1 == false) || (check2 == false) || (check3 == false))
      {
      load_okay = false;
      err_msg = "inconsistent data";
      }
    else
      {
      load_okay = f.good();
      }
    }
  else
    {
    load_okay = false;
    err_msg = "incorrect header";
    }
  
  return load_okay;
  }



// cubes



//! Save a cube as raw text (no header, human readable).
template<typename eT>
inline
bool
diskio::save_raw_ascii(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::fstream f(tmp_name.c_str(), std::fstream::out);
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = save_raw_ascii(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



//! Save a cube as raw text (no header, human readable).
template<typename eT>
inline
bool
diskio::save_raw_ascii(const Cube<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(f);
  
  const std::streamsize cell_width = diskio::prepare_stream<eT>(f);
  
  for(uword slice=0; slice < x.n_slices; ++slice)
    {
    for(uword row=0; row < x.n_rows; ++row)
      {
      for(uword col=0; col < x.n_cols; ++col)
        {
        f.put(' ');
        
        if(is_real<eT>::value)  { f.width(cell_width); }
        
        arma_ostream::raw_print_elem(f, x.at(row,col,slice));
        }
        
      f.put('\n');
      }
    }
  
  const bool save_okay = f.good();
  
  stream_state.restore(f);
  
  return save_okay;
  }



//! Save a cube as raw binary (no header)
template<typename eT>
inline
bool
diskio::save_raw_binary(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str(), std::fstream::binary);
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_raw_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



template<typename eT>
inline
bool
diskio::save_raw_binary(const Cube<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  f.write( reinterpret_cast<const char*>(x.mem), std::streamsize(x.n_elem*sizeof(eT)) );
  
  return f.good();
  }



//! Save a cube in text format (human readable),
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_ascii(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str());
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_arma_ascii(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



//! Save a cube in text format (human readable),
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_ascii(const Cube<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(f);
  
  f << diskio::gen_txt_header(x) << '\n';
  f << x.n_rows << ' ' << x.n_cols << ' ' << x.n_slices << '\n';
  
  const std::streamsize cell_width = diskio::prepare_stream<eT>(f);
  
  for(uword slice=0; slice < x.n_slices; ++slice)
    {
    for(uword row=0; row < x.n_rows; ++row)
      {
      for(uword col=0; col < x.n_cols; ++col)
        {
        f.put(' ');
        
        if(is_real<eT>::value)  { f.width(cell_width); }
        
        arma_ostream::raw_print_elem(f, x.at(row,col,slice));
        }
      
      f.put('\n');
      }
    }
  
  const bool save_okay = f.good();
  
  stream_state.restore(f);
  
  return save_okay;
  }



//! Save a cube in binary format,
//! with a header that stores the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_binary(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str(), std::fstream::binary);
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_arma_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



//! Save a cube in binary format,
//! with a header that stores the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::save_arma_binary(const Cube<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  f << diskio::gen_bin_header(x) << '\n';
  f << x.n_rows << ' ' << x.n_cols << ' ' << x.n_slices << '\n';
  
  f.write( reinterpret_cast<const char*>(x.mem), std::streamsize(x.n_elem*sizeof(eT)) );
  
  return f.good();
  }



//! Save a cube as part of a HDF5 file
template<typename eT>
inline
bool
diskio::save_hdf5_binary(const Cube<eT>& x, const hdf5_name& spec, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_HDF5)
    {
    hdf5_misc::hdf5_suspend_printing_errors hdf5_print_suspender;
    
    bool save_okay = false;
    
    const bool append  = bool(spec.opts.flags & hdf5_opts::flag_append);
    const bool replace = bool(spec.opts.flags & hdf5_opts::flag_replace);
    
    const bool use_existing_file = ((append || replace) && (arma_H5Fis_hdf5(spec.filename.c_str()) > 0));
    
    const std::string tmp_name = (use_existing_file) ? std::string() : diskio::gen_tmp_name(spec.filename);
    
    // Set up the file according to HDF5's preferences
    hid_t file = (use_existing_file) ? arma_H5Fopen(spec.filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT) : arma_H5Fcreate(tmp_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    if(file < 0)  { return false; }
    
    // We need to create a dataset, datatype, and dataspace
    hsize_t dims[3];
    dims[2] = x.n_rows;
    dims[1] = x.n_cols;
    dims[0] = x.n_slices;
    
    hid_t dataspace = arma_H5Screate_simple(3, dims, NULL);   // treat the cube as a 3d array dataspace
    hid_t datatype  = hdf5_misc::get_hdf5_type<eT>();
    
    // If this returned something invalid, well, it's time to crash.
    arma_check(datatype == -1, "Cube::save(): unknown datatype for HDF5");
    
    // MATLAB forces the users to specify a name at save time for HDF5;
    // Octave will use the default of 'dataset' unless otherwise specified.
    // If the user hasn't specified a dataset name, we will use 'dataset'
    // We may have to split out the group name from the dataset name.
    std::vector<hid_t> groups;
    std::string full_name = spec.dsname;
    size_t loc;
    while((loc = full_name.find("/")) != std::string::npos)
      {
      // Create another group...
      if(loc != 0) // Ignore the first /, if there is a leading /.
        {
        hid_t gid = arma_H5Gcreate((groups.size() == 0) ? file : groups[groups.size() - 1], full_name.substr(0, loc).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        if((gid < 0) && use_existing_file)
          {
          gid = arma_H5Gopen((groups.size() == 0) ? file : groups[groups.size() - 1], full_name.substr(0, loc).c_str(), H5P_DEFAULT);
          }
        
        groups.push_back(gid);
        }
      
      full_name = full_name.substr(loc + 1);
      }
    
    const std::string dataset_name = full_name.empty() ? std::string("dataset") : full_name;
    
    const hid_t last_group = (groups.size() == 0) ? file : groups[groups.size() - 1];
    
    if(use_existing_file && replace)
      {
      arma_H5Ldelete(last_group, dataset_name.c_str(), H5P_DEFAULT);
      // NOTE: H5Ldelete() in HDF5 v1.8 doesn't reclaim the deleted space; use h5repack to reclaim space: h5repack oldfile.h5 newfile.h5
      // NOTE: has this behaviour changed in HDF5 1.10 ?
      // NOTE: https://lists.hdfgroup.org/pipermail/hdf-forum_lists.hdfgroup.org/2017-August/010482.html
      // NOTE: https://lists.hdfgroup.org/pipermail/hdf-forum_lists.hdfgroup.org/2017-August/010486.html
      }
    
    hid_t dataset = arma_H5Dcreate(last_group, dataset_name.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    if(dataset < 0)
      {
      save_okay = false;
      
      err_msg = "couldn't create dataset";
      }
    else
      {
      save_okay = (arma_H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, x.mem) >= 0);
      
      arma_H5Dclose(dataset);
      }
    
    arma_H5Tclose(datatype);
    arma_H5Sclose(dataspace);
    for(size_t i = 0; i < groups.size(); ++i)  { arma_H5Gclose(groups[i]); }
    arma_H5Fclose(file);
    
    if((use_existing_file == false) && (save_okay == true))  { save_okay = diskio::safe_rename(tmp_name, spec.filename); }
    
    return save_okay;
    }
  #else
    {
    arma_ignore(x);
    arma_ignore(spec);
    arma_ignore(err_msg);
    
    arma_stop_logic_error("Cube::save(): use of HDF5 must be enabled");
    
    return false;
    }
  #endif
  }



//! Load a cube as raw text (no header, human readable).
//! NOTE: this is much slower than reading a file with a header.
template<typename eT>
inline
bool
diskio::load_raw_ascii(Cube<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp;
  const bool load_okay = diskio::load_raw_ascii(tmp, name, err_msg);
  
  if(load_okay)
    {
    if(tmp.is_empty() == false)
      {
      try { x.set_size(tmp.n_rows, tmp.n_cols, 1); } catch(...) { err_msg = "not enough memory"; return false; }
      
      x.slice(0) = tmp;
      }
    else
      {
      x.reset();
      }
    }
  
  return load_okay;
  }



//! Load a cube as raw text (no header, human readable).
//! NOTE: this is much slower than reading a file with a header.
template<typename eT>
inline
bool
diskio::load_raw_ascii(Cube<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp;
  const bool load_okay = diskio::load_raw_ascii(tmp, f, err_msg);
  
  if(load_okay)
    {
    if(tmp.is_empty() == false)
      {
      try { x.set_size(tmp.n_rows, tmp.n_cols, 1); } catch(...) { err_msg = "not enough memory"; return false; }
      
      x.slice(0) = tmp;
      }
    else
      {
      x.reset();
      }
    }
  
  return load_okay;
  }



//! Load a cube in binary format (no header);
//! the cube is assumed to have one slice with one column
template<typename eT>
inline
bool
diskio::load_raw_binary(Cube<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f;
  f.open(name.c_str(), std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_raw_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::load_raw_binary(Cube<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  f.clear();
  const std::streampos pos1 = f.tellg();
  
  f.clear();
  f.seekg(0, ios::end);
  
  f.clear();
  const std::streampos pos2 = f.tellg();
  
  const uword N = ( (pos1 >= 0) && (pos2 >= 0) ) ? uword(pos2 - pos1) : 0;
  
  f.clear();
  //f.seekg(0, ios::beg);
  f.seekg(pos1);
  
  try { x.set_size(N / uword(sizeof(eT)), 1, 1); } catch(...) { err_msg = "not enough memory"; return false; }
  
  f.clear();
  f.read( reinterpret_cast<char *>(x.memptr()), std::streamsize(x.n_elem * uword(sizeof(eT))) );
  
  return f.good();
  }



//! Load a cube in text format (human readable),
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::load_arma_ascii(Cube<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f(name.c_str());
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_arma_ascii(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }
  


//! Load a cube in text format (human readable),
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::load_arma_ascii(Cube<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::streampos pos = f.tellg();
  
  bool load_okay = true;
  
  std::string f_header;
  uword       f_n_rows;
  uword       f_n_cols;
  uword       f_n_slices;
  
  f >> f_header;
  f >> f_n_rows;
  f >> f_n_cols;
  f >> f_n_slices;
  
  if(f_header == diskio::gen_txt_header(x))
    {
    try { x.set_size(f_n_rows, f_n_cols, f_n_slices); } catch(...) { err_msg = "not enough memory"; return false; }

    for(uword slice = 0; slice < x.n_slices; ++slice)
    for(uword   row = 0;   row < x.n_rows;   ++row  )
    for(uword   col = 0;   col < x.n_cols;   ++col  )
      {
      f >> x.at(row,col,slice);
      }
    
    load_okay = f.good();
    }
  else
    {
    load_okay = false;
    err_msg = "incorrect header";
    }
  
  
  // allow automatic conversion of u32/s32 cubes into u64/s64 cubes
  
  if(load_okay == false)
    {
    if( (sizeof(eT) == 8) && is_same_type<uword,eT>::yes )
      {
      Cube<u32>   tmp;
      std::string junk;
      
      f.clear();
      f.seekg(pos);
      
      load_okay = diskio::load_arma_ascii(tmp, f, junk);
      
      if(load_okay)  { x = conv_to< Cube<eT> >::from(tmp); }
      }
    else
    if( (sizeof(eT) == 8) && is_same_type<sword,eT>::yes )
      {
      Cube<s32>   tmp;
      std::string junk;
      
      f.clear();
      f.seekg(pos);
      
      load_okay = diskio::load_arma_ascii(tmp, f, junk);
      
      if(load_okay)  { x = conv_to< Cube<eT> >::from(tmp); }
      }
    }
  
  return load_okay;
  }



//! Load a cube in binary format,
//! with a header that indicates the cube type as well as its dimensions
template<typename eT>
inline
bool
diskio::load_arma_binary(Cube<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f;
  f.open(name.c_str(), std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_arma_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::load_arma_binary(Cube<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::streampos pos = f.tellg();
    
  bool load_okay = true;
  
  std::string f_header;
  uword       f_n_rows;
  uword       f_n_cols;
  uword       f_n_slices;
  
  f >> f_header;
  f >> f_n_rows;
  f >> f_n_cols;
  f >> f_n_slices;
  
  if(f_header == diskio::gen_bin_header(x))
    {
    //f.seekg(1, ios::cur);  // NOTE: this may not be portable, as on a Windows machine a newline could be two characters
    f.get();
    
    try { x.set_size(f_n_rows, f_n_cols, f_n_slices); } catch(...) { err_msg = "not enough memory"; return false; }
    
    f.read( reinterpret_cast<char *>(x.memptr()), std::streamsize(x.n_elem*sizeof(eT)) );
    
    load_okay = f.good();
    }
  else
    {
    load_okay = false;
    err_msg = "incorrect header";
    }
  
  
  // allow automatic conversion of u32/s32 cubes into u64/s64 cubes
  
  if(load_okay == false)
    {
    if( (sizeof(eT) == 8) && is_same_type<uword,eT>::yes )
      {
      Cube<u32>   tmp;
      std::string junk;
      
      f.clear();
      f.seekg(pos);
      
      load_okay = diskio::load_arma_binary(tmp, f, junk);
      
      if(load_okay)  { x = conv_to< Cube<eT> >::from(tmp); }
      }
    else
    if( (sizeof(eT) == 8) && is_same_type<sword,eT>::yes )
      {
      Cube<s32>   tmp;
      std::string junk;
      
      f.clear();
      f.seekg(pos);
      
      load_okay = diskio::load_arma_binary(tmp, f, junk);
      
      if(load_okay)  { x = conv_to< Cube<eT> >::from(tmp); }
      }
    }
  
  return load_okay;
  }



//! Load a HDF5 file as a cube
template<typename eT>
inline
bool
diskio::load_hdf5_binary(Cube<eT>& x, const hdf5_name& spec, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_HDF5)
    {
    hdf5_misc::hdf5_suspend_printing_errors hdf5_print_suspender;
    
    bool load_okay = false;
    
    hid_t fid = arma_H5Fopen(spec.filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    
    if(fid >= 0)
      {
      // MATLAB HDF5 dataset names are user-specified;
      // Octave tends to store the datasets in a group, with the actual dataset being referred to as "value".
      // If the user hasn't specified a dataset, we will search for "dataset" and "value",
      // and if those are not found we will take the first dataset we do find.
      
      std::vector<std::string> searchNames;
      
      const bool exact = (spec.dsname.empty() == false);
      
      if(exact)
        {
        searchNames.push_back(spec.dsname);
        }
      else
        {
        searchNames.push_back("dataset");
        searchNames.push_back("value"  );
        }
      
      hid_t dataset = hdf5_misc::search_hdf5_file(searchNames, fid, 3, exact);
      
      if(dataset >= 0)
        {
        hid_t filespace = arma_H5Dget_space(dataset);
        
        // This must be <= 3 due to our search rules.
        const int ndims = arma_H5Sget_simple_extent_ndims(filespace);
        
        hsize_t dims[3];
        const herr_t query_status = arma_H5Sget_simple_extent_dims(filespace, dims, NULL);
        
        // arma_check(query_status < 0, "Cube::load(): cannot get size of HDF5 dataset");
        if(query_status < 0)
          {
          err_msg = "cannot get size of HDF5 dataset";
          
          arma_H5Sclose(filespace);
          arma_H5Dclose(dataset);
          arma_H5Fclose(fid);
          
          return false;
          }
        
        if(ndims == 1) { dims[1] = 1; dims[2] = 1; }  // Vector case; one row/colum, several slices
        if(ndims == 2) {              dims[2] = 1; }  // Matrix case; one column, several rows/slices
        
        try { x.set_size(dims[2], dims[1], dims[0]); } catch(...) { err_msg = "not enough memory"; return false; }
        
        // Now we have to see what type is stored to figure out how to load it.
        hid_t datatype = arma_H5Dget_type(dataset);
        hid_t mat_type = hdf5_misc::get_hdf5_type<eT>();
        
        // If these are the same type, it is simple.
        if(arma_H5Tequal(datatype, mat_type) > 0)
          {
          // Load directly; H5S_ALL used so that we load the entire dataset.
          hid_t read_status = arma_H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, void_ptr(x.memptr()));
          
          if(read_status >= 0) { load_okay = true; }
          }
        else
          {
          // Load into another array and convert its type accordingly.
          hid_t read_status = hdf5_misc::load_and_convert_hdf5(x.memptr(), dataset, datatype, x.n_elem);
          
          if(read_status >= 0) { load_okay = true; }
          }
        
        // Now clean up.
        arma_H5Tclose(datatype);
        arma_H5Tclose(mat_type);
        arma_H5Sclose(filespace);
        }
      
      arma_H5Dclose(dataset);
      
      arma_H5Fclose(fid);
      
      if(load_okay == false)
        {
        err_msg = "unsupported or missing HDF5 data";
        }
      }
    else
      {
      err_msg = "cannot open";
      }
    
    return load_okay;
    }
  #else
    {
    arma_ignore(x);
    arma_ignore(spec);
    arma_ignore(err_msg);
    
    arma_stop_logic_error("Cube::load(): use of HDF5 must be enabled");
    
    return false;
    }
  #endif
  }



//! Try to load a cube by automatically determining its type
template<typename eT>
inline
bool
diskio::load_auto_detect(Cube<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_HDF5)
    // We're currently using the C bindings for the HDF5 library, which don't support C++ streams
    if( arma_H5Fis_hdf5(name.c_str()) ) { return load_hdf5_binary(x, name, err_msg); }
  #endif
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_auto_detect(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



//! Try to load a cube by automatically determining its type
template<typename eT>
inline
bool
diskio::load_auto_detect(Cube<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  const char* ARMA_CUB_TXT_str = "ARMA_CUB_TXT";
  const char* ARMA_CUB_BIN_str = "ARMA_CUB_BIN";
  const char*           P6_str = "P6";
  
  const uword ARMA_CUB_TXT_len = uword(12);
  const uword ARMA_CUB_BIN_len = uword(12);
  const uword           P6_len = uword(2);
  
  podarray<char> header(ARMA_CUB_TXT_len + 1);
  
  char* header_mem = header.memptr();
  
  std::streampos pos = f.tellg();
  
  f.read( header_mem, std::streamsize(ARMA_CUB_TXT_len) );
  f.clear();
  f.seekg(pos);
  
  header_mem[ARMA_CUB_TXT_len] = '\0';
  
  if( std::strncmp(ARMA_CUB_TXT_str, header_mem, size_t(ARMA_CUB_TXT_len)) == 0 )
    {
    return load_arma_ascii(x, f, err_msg);
    }
  else
  if( std::strncmp(ARMA_CUB_BIN_str, header_mem, size_t(ARMA_CUB_BIN_len)) == 0 )
    {
    return load_arma_binary(x, f, err_msg);
    }
  else
  if( std::strncmp(P6_str, header_mem, size_t(P6_len)) == 0 )
    {
    return load_ppm_binary(x, f, err_msg);
    }
  else
    {
    const file_type ft = guess_file_type_internal(f);
    
    switch(ft)
      {
      // case csv_ascii:
      //   return load_csv_ascii(x, f, err_msg);
      //   break;
      
      case raw_binary:
        return load_raw_binary(x, f, err_msg);
        break;
        
      case raw_ascii:
        return load_raw_ascii(x, f, err_msg);
        break;
        
      default:
        err_msg = "unknown data";
        return false;
      }
    }
  
  return false;
  }





// fields



template<typename T1>
inline
bool
diskio::save_arma_binary(const field<T1>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f( tmp_name.c_str(), std::fstream::binary );
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_arma_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



template<typename T1>
inline
bool
diskio::save_arma_binary(const field<T1>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( (is_Mat<T1>::value == false) && (is_Cube<T1>::value == false) ));
  
  if(x.n_slices <= 1)
    {
    f << "ARMA_FLD_BIN" << '\n';
    f << x.n_rows       << '\n';
    f << x.n_cols       << '\n';
    }
  else
    {
    f << "ARMA_FL3_BIN" << '\n';
    f << x.n_rows       << '\n';
    f << x.n_cols       << '\n';
    f << x.n_slices     << '\n';
    }
  
  bool save_okay = true;
  
  for(uword i=0; i<x.n_elem; ++i)
    {
    save_okay = diskio::save_arma_binary(x[i], f);
    
    if(save_okay == false)  { break; }
    }
  
  return save_okay;
  }



template<typename T1>
inline
bool
diskio::load_arma_binary(field<T1>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f( name.c_str(), std::fstream::binary );
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_arma_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename T1>
inline
bool
diskio::load_arma_binary(field<T1>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( (is_Mat<T1>::value == false) && (is_Cube<T1>::value == false) ));
  
  bool load_okay = true;
  
  std::string f_type;
  f >> f_type;
  
  if(f_type == "ARMA_FLD_BIN")
    {
    uword f_n_rows;
    uword f_n_cols;
    
    f >> f_n_rows;
    f >> f_n_cols;
    
    try { x.set_size(f_n_rows, f_n_cols); } catch(...) { err_msg = "not enough memory"; return false; }
    
    f.get();
    
    for(uword i=0; i<x.n_elem; ++i)
      {
      load_okay = diskio::load_arma_binary(x[i], f, err_msg);
      
      if(load_okay == false)  { break; }
      }
    }
  else
  if(f_type == "ARMA_FL3_BIN")
    {
    uword f_n_rows;
    uword f_n_cols;
    uword f_n_slices;
    
    f >> f_n_rows;
    f >> f_n_cols;
    f >> f_n_slices;
    
    try { x.set_size(f_n_rows, f_n_cols, f_n_slices); } catch(...) { err_msg = "not enough memory"; return false; }
    
    f.get();
    
    for(uword i=0; i<x.n_elem; ++i)
      {
      load_okay = diskio::load_arma_binary(x[i], f, err_msg);
      
      if(load_okay == false)  { break; }
      }
    }
  else
    {
    load_okay = false;
    err_msg = "unsupported field type";
    }
  
  return load_okay;
  }



inline
bool
diskio::save_std_string(const field<std::string>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f( tmp_name.c_str(), std::fstream::binary );
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_std_string(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



inline
bool
diskio::save_std_string(const field<std::string>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  for(uword row=0; row<x.n_rows; ++row)
  for(uword col=0; col<x.n_cols; ++col)
    {
    f << x.at(row,col);
    
    if(col < x.n_cols-1)
      {
      f << ' ';
      }
    else
      {
      f << '\n';
      }
    }
  
  return f.good();
  }



inline
bool
diskio::load_std_string(field<std::string>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f( name.c_str() );
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_std_string(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



inline
bool
diskio::load_std_string(field<std::string>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  //
  // work out the size
  
  uword f_n_rows = 0;
  uword f_n_cols = 0;
  
  bool f_n_cols_found = false;
  
  std::string line_string;
  std::string token;
  
  while( f.good() && load_okay )
    {
    std::getline(f, line_string);
    
    if(line_string.size() == 0)  { break; }
    
    std::stringstream line_stream(line_string);
    
    uword line_n_cols = 0;
    
    while(line_stream >> token)  { line_n_cols++; }
    
    if(f_n_cols_found == false)
      {
      f_n_cols = line_n_cols;
      f_n_cols_found = true;
      }
    else
      {
      if(line_n_cols != f_n_cols)
        {
        load_okay = false;
        err_msg = "inconsistent number of columns";
        }
      }
    
    ++f_n_rows;
    }
    
  if(load_okay)
    {
    f.clear();
    f.seekg(0, ios::beg);
    //f.seekg(start);
    
    try { x.set_size(f_n_rows, f_n_cols); } catch(...) { err_msg = "not enough memory"; return false; }
    
    for(uword row=0; row < x.n_rows; ++row)
    for(uword col=0; col < x.n_cols; ++col)
      {
      f >> x.at(row,col);
      }
    }
  
  if(f.good() == false)  { load_okay = false; }
  
  return load_okay;
  }



//! Try to load a field by automatically determining its type
template<typename T1>
inline
bool
diskio::load_auto_detect(field<T1>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_auto_detect(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



//! Try to load a field by automatically determining its type
template<typename T1>
inline
bool
diskio::load_auto_detect(field<T1>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_Mat<T1>::value == false ));
  
  static const std::string ARMA_FLD_BIN = "ARMA_FLD_BIN";
  static const std::string ARMA_FL3_BIN = "ARMA_FL3_BIN";
  static const std::string           P6 = "P6";
  
  podarray<char> raw_header(uword(ARMA_FLD_BIN.length()) + 1);
  
  std::streampos pos = f.tellg();
  
  f.read( raw_header.memptr(), std::streamsize(ARMA_FLD_BIN.length()) );
  
  f.clear();
  f.seekg(pos);
  
  raw_header[uword(ARMA_FLD_BIN.length())] = '\0';
  
  const std::string header = raw_header.mem;
  
  if(ARMA_FLD_BIN == header.substr(0, ARMA_FLD_BIN.length()))
    {
    return load_arma_binary(x, f, err_msg);
    }
  else
  if(ARMA_FL3_BIN == header.substr(0, ARMA_FL3_BIN.length()))
    {
    return load_arma_binary(x, f, err_msg);
    }
  else
  if(P6 == header.substr(0, P6.length()))
    {
    return load_ppm_binary(x, f, err_msg);
    }
  else
    {
    err_msg = "unsupported header";
    return false;
    }
  }



//
// handling of PPM images by cubes


template<typename eT>
inline
bool
diskio::load_ppm_binary(Cube<eT>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_ppm_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::load_ppm_binary(Cube<eT>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  std::string f_header;
  
  f >> f_header;
  
  if(f_header == "P6")
    {
    uword f_n_rows = 0;
    uword f_n_cols = 0;
    int f_maxval = 0;
    
    diskio::pnm_skip_comments(f);
    
    f >> f_n_cols;
    diskio::pnm_skip_comments(f);
    
    f >> f_n_rows;
    diskio::pnm_skip_comments(f);
    
    f >> f_maxval;
    f.get();
    
    if( (f_maxval > 0) && (f_maxval <= 65535) )
      {
      try { x.set_size(f_n_rows, f_n_cols, 3); } catch(...) { err_msg = "not enough memory"; return false; }
      
      if(f_maxval <= 255)
        {
        const uword n_elem = 3*f_n_cols*f_n_rows;
        podarray<u8> tmp(n_elem);
        
        f.read( reinterpret_cast<char*>(tmp.memptr()), std::streamsize(n_elem) );
        
        uword i = 0;
        
        //cout << "f_n_cols = " << f_n_cols << endl;
        //cout << "f_n_rows = " << f_n_rows << endl;
        
        for(uword row=0; row < f_n_rows; ++row)
        for(uword col=0; col < f_n_cols; ++col)
          {
          x.at(row,col,0) = eT(tmp[i+0]);
          x.at(row,col,1) = eT(tmp[i+1]);
          x.at(row,col,2) = eT(tmp[i+2]);
          i+=3;
          }
        }
      else
        {
        const uword n_elem = 3*f_n_cols*f_n_rows;
        podarray<u16> tmp(n_elem);
        
        f.read( reinterpret_cast<char *>(tmp.memptr()), std::streamsize(2*n_elem) );
        
        uword i = 0;
        
        for(uword row=0; row < f_n_rows; ++row)
        for(uword col=0; col < f_n_cols; ++col)
          {
          x.at(row,col,0) = eT(tmp[i+0]);
          x.at(row,col,1) = eT(tmp[i+1]);
          x.at(row,col,2) = eT(tmp[i+2]);
          i+=3;
          }
        }
      }
    else
      {
      load_okay = false;
      err_msg = "functionality unimplemented";
      }
    
    if(f.good() == false)  { load_okay = false; }
    }
  else
    {
    load_okay = false;
    err_msg = "unsupported header";
    }
  
  return load_okay;
  }



template<typename eT>
inline
bool
diskio::save_ppm_binary(const Cube<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f( tmp_name.c_str(), std::fstream::binary );
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_ppm_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



template<typename eT>
inline
bool
diskio::save_ppm_binary(const Cube<eT>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (x.n_slices != 3), "diskio::save_ppm_binary(): given cube must have exactly 3 slices" );
  
  const uword n_elem = 3 * x.n_rows * x.n_cols;
  podarray<u8> tmp(n_elem);
  
  uword i = 0;
  for(uword row=0; row < x.n_rows; ++row)
    {
    for(uword col=0; col < x.n_cols; ++col)
      {
      tmp[i+0] = u8( access::tmp_real( x.at(row,col,0) ) );
      tmp[i+1] = u8( access::tmp_real( x.at(row,col,1) ) );
      tmp[i+2] = u8( access::tmp_real( x.at(row,col,2) ) );
      
      i+=3;
      }
    }
  
  f << "P6" << '\n';
  f << x.n_cols << '\n';
  f << x.n_rows << '\n';
  f << 255 << '\n';
  
  f.write( reinterpret_cast<const char*>(tmp.mem), std::streamsize(n_elem) );
  
  return f.good();
  }



//
// handling of PPM images by fields



template<typename T1>
inline
bool
diskio::load_ppm_binary(field<T1>& x, const std::string& name, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = f.is_open();
  
  if(load_okay)
    {
    load_okay = diskio::load_ppm_binary(x, f, err_msg);
    f.close();
    }
  
  return load_okay;
  }



template<typename T1>
inline
bool
diskio::load_ppm_binary(field<T1>& x, std::istream& f, std::string& err_msg)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_Mat<T1>::value == false ));
  typedef typename T1::elem_type eT;
  
  bool load_okay = true;
  
  std::string f_header;
  
  f >> f_header;
  
  if(f_header == "P6")
    {
    uword f_n_rows = 0;
    uword f_n_cols = 0;
    int f_maxval = 0;
    
    diskio::pnm_skip_comments(f);
    
    f >> f_n_cols;
    diskio::pnm_skip_comments(f);
    
    f >> f_n_rows;
    diskio::pnm_skip_comments(f);
    
    f >> f_maxval;
    f.get();
    
    if( (f_maxval > 0) && (f_maxval <= 65535) )
      {
      x.set_size(3);
      Mat<eT>& R = x(0);
      Mat<eT>& G = x(1);
      Mat<eT>& B = x(2);
      
      try { R.set_size(f_n_rows,f_n_cols); } catch(...) { err_msg = "not enough memory"; return false; }
      try { G.set_size(f_n_rows,f_n_cols); } catch(...) { err_msg = "not enough memory"; return false; }
      try { B.set_size(f_n_rows,f_n_cols); } catch(...) { err_msg = "not enough memory"; return false; }
      
      if(f_maxval <= 255)
        {
        const uword n_elem = 3*f_n_cols*f_n_rows;
        podarray<u8> tmp(n_elem);
        
        f.read( reinterpret_cast<char*>(tmp.memptr()), std::streamsize(n_elem) );
        
        uword i = 0;
        
        //cout << "f_n_cols = " << f_n_cols << endl;
        //cout << "f_n_rows = " << f_n_rows << endl;
        
        
        for(uword row=0; row < f_n_rows; ++row)
          {
          for(uword col=0; col < f_n_cols; ++col)
            {
            R.at(row,col) = eT(tmp[i+0]);
            G.at(row,col) = eT(tmp[i+1]);
            B.at(row,col) = eT(tmp[i+2]);
            i+=3;
            }
          
          }
        }
      else
        {
        const uword n_elem = 3*f_n_cols*f_n_rows;
        podarray<u16> tmp(n_elem);
        
        f.read( reinterpret_cast<char *>(tmp.memptr()), std::streamsize(2*n_elem) );
        
        uword i = 0;
        
        for(uword row=0; row < f_n_rows; ++row)
        for(uword col=0; col < f_n_cols; ++col)
          {
          R.at(row,col) = eT(tmp[i+0]);
          G.at(row,col) = eT(tmp[i+1]);
          B.at(row,col) = eT(tmp[i+2]);
          i+=3;
          }
        }
      }
    else
      {
      load_okay = false;
      err_msg = "functionality unimplemented";
      }
    
    if(f.good() == false)  { load_okay = false; }
    }
  else
    {
    load_okay = false;
    err_msg = "unsupported header";
    }
  
  return load_okay;
  }



template<typename T1>
inline
bool
diskio::save_ppm_binary(const field<T1>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  std::ofstream f( tmp_name.c_str(), std::fstream::binary );
  
  bool save_okay = f.is_open();
  
  if(save_okay)
    {
    save_okay = diskio::save_ppm_binary(x, f);
    
    f.flush();
    f.close();
    
    if(save_okay)  { save_okay = diskio::safe_rename(tmp_name, final_name); }
    }
  
  return save_okay;
  }



template<typename T1>
inline
bool
diskio::save_ppm_binary(const field<T1>& x, std::ostream& f)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check(( is_Mat<T1>::value == false ));
  
  typedef typename T1::elem_type eT;
  
  arma_debug_check( (x.n_elem != 3), "diskio::save_ppm_binary(): given field must have exactly 3 matrices of equal size" );
  
  bool same_size = true;
  for(uword i=1; i<3; ++i)
    {
    if( (x(0).n_rows != x(i).n_rows) || (x(0).n_cols != x(i).n_cols) )
      {
      same_size = false;
      break;
      }
    }
  
  arma_debug_check( (same_size != true), "diskio::save_ppm_binary(): given field must have exactly 3 matrices of equal size" );
  
  const Mat<eT>& R = x(0);
  const Mat<eT>& G = x(1);
  const Mat<eT>& B = x(2);
  
  f << "P6" << '\n';
  f << R.n_cols << '\n';
  f << R.n_rows << '\n';
  f << 255 << '\n';

  const uword n_elem = 3 * R.n_rows * R.n_cols;
  podarray<u8> tmp(n_elem);

  uword i = 0;
  for(uword row=0; row < R.n_rows; ++row)
  for(uword col=0; col < R.n_cols; ++col)
    {
    tmp[i+0] = u8( access::tmp_real( R.at(row,col) ) );
    tmp[i+1] = u8( access::tmp_real( G.at(row,col) ) );
    tmp[i+2] = u8( access::tmp_real( B.at(row,col) ) );
    
    i+=3;
    }
  
  f.write( reinterpret_cast<const char*>(tmp.mem), std::streamsize(n_elem) );
  
  return f.good();
  }



//! @}

