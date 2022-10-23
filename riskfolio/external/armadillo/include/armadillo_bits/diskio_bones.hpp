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


//! class for saving and loading matrices and fields - INTERNAL USE ONLY!
class diskio
  {
  public:
  
  inline arma_deprecated static file_type guess_file_type(std::istream& f);
  
  
  private:
  
  template<typename eT> friend class   Mat;
  template<typename eT> friend class  Cube;
  template<typename eT> friend class SpMat;
  template<typename oT> friend class field;
  
  friend class   Mat_aux;
  friend class  Cube_aux;
  friend class SpMat_aux;
  friend class field_aux;
  
  template<typename eT> inline arma_cold static std::string gen_txt_header(const Mat<eT>&);
  template<typename eT> inline arma_cold static std::string gen_bin_header(const Mat<eT>&);
  
  template<typename eT> inline arma_cold static std::string gen_bin_header(const SpMat<eT>&);

  template<typename eT> inline arma_cold static std::string gen_txt_header(const Cube<eT>&);
  template<typename eT> inline arma_cold static std::string gen_bin_header(const Cube<eT>&);
  
  inline arma_cold static file_type guess_file_type_internal(std::istream& f);
  
  inline arma_cold static std::string gen_tmp_name(const std::string& x);
  
  inline arma_cold static bool safe_rename(const std::string& old_name, const std::string& new_name);
  
  template<typename eT> inline static bool convert_token(eT&              val, const std::string& token);
  template<typename  T> inline static bool convert_token(std::complex<T>& val, const std::string& token);
  
  template<typename eT> inline static std::streamsize prepare_stream(std::ostream& f);
  
  
  //
  // matrix saving
  
  template<typename eT> inline static bool save_raw_ascii  (const Mat<eT>&                x, const std::string& final_name);
  template<typename eT> inline static bool save_raw_binary (const Mat<eT>&                x, const std::string& final_name);
  template<typename eT> inline static bool save_arma_ascii (const Mat<eT>&                x, const std::string& final_name);
  template<typename eT> inline static bool save_csv_ascii  (const Mat<eT>&                x, const std::string& final_name, const field<std::string>& header, const bool with_header, const char separator);
  template<typename eT> inline static bool save_coord_ascii(const Mat<eT>&                x, const std::string& final_name);
  template<typename eT> inline static bool save_arma_binary(const Mat<eT>&                x, const std::string& final_name);
  template<typename eT> inline static bool save_pgm_binary (const Mat<eT>&                x, const std::string& final_name);
  template<typename  T> inline static bool save_pgm_binary (const Mat< std::complex<T> >& x, const std::string& final_name);
  template<typename eT> inline static bool save_hdf5_binary(const Mat<eT>&                x, const   hdf5_name& spec, std::string& err_msg);
  
  template<typename eT> inline static bool save_raw_ascii  (const Mat<eT>&                x, std::ostream& f);
  template<typename eT> inline static bool save_raw_binary (const Mat<eT>&                x, std::ostream& f);
  template<typename eT> inline static bool save_arma_ascii (const Mat<eT>&                x, std::ostream& f);
  template<typename eT> inline static bool save_csv_ascii  (const Mat<eT>&                x, std::ostream& f, const char separator);
  template<typename  T> inline static bool save_csv_ascii  (const Mat< std::complex<T> >& x, std::ostream& f, const char separator);
  template<typename eT> inline static bool save_coord_ascii(const Mat<eT>&                x, std::ostream& f);
  template<typename  T> inline static bool save_coord_ascii(const Mat< std::complex<T> >& x, std::ostream& f);
  template<typename eT> inline static bool save_arma_binary(const Mat<eT>&                x, std::ostream& f);
  template<typename eT> inline static bool save_pgm_binary (const Mat<eT>&                x, std::ostream& f);
  template<typename  T> inline static bool save_pgm_binary (const Mat< std::complex<T> >& x, std::ostream& f);
  
  
  //
  // matrix loading
  
  template<typename eT> inline static bool load_raw_ascii  (Mat<eT>&                x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_raw_binary (Mat<eT>&                x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_arma_ascii (Mat<eT>&                x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_csv_ascii  (Mat<eT>&                x, const std::string& name, std::string& err_msg, field<std::string>& header, const bool with_header, const char separator);
  template<typename eT> inline static bool load_coord_ascii(Mat<eT>&                x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_arma_binary(Mat<eT>&                x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_pgm_binary (Mat<eT>&                x, const std::string& name, std::string& err_msg);
  template<typename  T> inline static bool load_pgm_binary (Mat< std::complex<T> >& x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_hdf5_binary(Mat<eT>&                x, const   hdf5_name& spec, std::string& err_msg);
  template<typename eT> inline static bool load_auto_detect(Mat<eT>&                x, const std::string& name, std::string& err_msg);
  
  template<typename eT> inline static bool load_raw_ascii  (Mat<eT>&                x, std::istream& f,  std::string& err_msg);
  template<typename eT> inline static bool load_raw_binary (Mat<eT>&                x, std::istream& f,  std::string& err_msg);
  template<typename eT> inline static bool load_arma_ascii (Mat<eT>&                x, std::istream& f,  std::string& err_msg);
  template<typename eT> inline static bool load_csv_ascii  (Mat<eT>&                x, std::istream& f,  std::string& err_msg, const char separator);
  template<typename  T> inline static bool load_csv_ascii  (Mat< std::complex<T> >& x, std::istream& f,  std::string& err_msg, const char separator);
  template<typename eT> inline static bool load_coord_ascii(Mat<eT>&                x, std::istream& f,  std::string& err_msg);
  template<typename  T> inline static bool load_coord_ascii(Mat< std::complex<T> >& x, std::istream& f,  std::string& err_msg);
  template<typename eT> inline static bool load_arma_binary(Mat<eT>&                x, std::istream& f,  std::string& err_msg);
  template<typename eT> inline static bool load_pgm_binary (Mat<eT>&                x, std::istream& is, std::string& err_msg);
  template<typename  T> inline static bool load_pgm_binary (Mat< std::complex<T> >& x, std::istream& is, std::string& err_msg);
  template<typename eT> inline static bool load_auto_detect(Mat<eT>&                x, std::istream& f,  std::string& err_msg);
  
  inline static void pnm_skip_comments(std::istream& f);
  
  
  //
  // sparse matrix saving
  
  template<typename eT> inline static bool save_csv_ascii  (const SpMat<eT>& x, const std::string& final_name, const field<std::string>& header, const bool with_header, const char separator);
  template<typename eT> inline static bool save_coord_ascii(const SpMat<eT>& x, const std::string& final_name);
  template<typename eT> inline static bool save_arma_binary(const SpMat<eT>& x, const std::string& final_name);
  
  template<typename eT> inline static bool save_csv_ascii  (const SpMat<eT>& x,                std::ostream& f, const char separator);
  template<typename  T> inline static bool save_csv_ascii  (const SpMat< std::complex<T> >& x, std::ostream& f, const char separator);
  template<typename eT> inline static bool save_coord_ascii(const SpMat<eT>& x,                std::ostream& f);
  template<typename  T> inline static bool save_coord_ascii(const SpMat< std::complex<T> >& x, std::ostream& f);
  template<typename eT> inline static bool save_arma_binary(const SpMat<eT>& x,                std::ostream& f);
  
  
  //
  // sparse matrix loading
  
  template<typename eT> inline static bool load_csv_ascii  (SpMat<eT>& x, const std::string& name, std::string& err_msg, field<std::string>& header, const bool with_header, const char separator);
  template<typename eT> inline static bool load_coord_ascii(SpMat<eT>& x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_arma_binary(SpMat<eT>& x, const std::string& name, std::string& err_msg);
  
  template<typename eT> inline static bool load_csv_ascii  (SpMat<eT>& x,                std::istream& f, std::string& err_msg, const char separator);
  template<typename  T> inline static bool load_csv_ascii  (SpMat< std::complex<T> >& x, std::istream& f, std::string& err_msg, const char separator);
  template<typename eT> inline static bool load_coord_ascii(SpMat<eT>& x,                std::istream& f, std::string& err_msg);
  template<typename  T> inline static bool load_coord_ascii(SpMat< std::complex<T> >& x, std::istream& f, std::string& err_msg);
  template<typename eT> inline static bool load_arma_binary(SpMat<eT>& x,                std::istream& f, std::string& err_msg);
  
  
  
  //
  // cube saving
  
  template<typename eT> inline static bool save_raw_ascii  (const Cube<eT>& x, const std::string& name);
  template<typename eT> inline static bool save_raw_binary (const Cube<eT>& x, const std::string& name);
  template<typename eT> inline static bool save_arma_ascii (const Cube<eT>& x, const std::string& name);
  template<typename eT> inline static bool save_arma_binary(const Cube<eT>& x, const std::string& name);
  template<typename eT> inline static bool save_hdf5_binary(const Cube<eT>& x, const   hdf5_name& spec, std::string& err_msg);
  
  template<typename eT> inline static bool save_raw_ascii  (const Cube<eT>& x, std::ostream& f);
  template<typename eT> inline static bool save_raw_binary (const Cube<eT>& x, std::ostream& f);
  template<typename eT> inline static bool save_arma_ascii (const Cube<eT>& x, std::ostream& f);
  template<typename eT> inline static bool save_arma_binary(const Cube<eT>& x, std::ostream& f);
  
  
  //
  // cube loading
  
  template<typename eT> inline static bool load_raw_ascii  (Cube<eT>& x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_raw_binary (Cube<eT>& x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_arma_ascii (Cube<eT>& x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_arma_binary(Cube<eT>& x, const std::string& name, std::string& err_msg);
  template<typename eT> inline static bool load_hdf5_binary(Cube<eT>& x, const   hdf5_name& spec, std::string& err_msg);
  template<typename eT> inline static bool load_auto_detect(Cube<eT>& x, const std::string& name, std::string& err_msg);
  
  template<typename eT> inline static bool load_raw_ascii  (Cube<eT>& x, std::istream& f, std::string& err_msg);
  template<typename eT> inline static bool load_raw_binary (Cube<eT>& x, std::istream& f, std::string& err_msg);
  template<typename eT> inline static bool load_arma_ascii (Cube<eT>& x, std::istream& f, std::string& err_msg);
  template<typename eT> inline static bool load_arma_binary(Cube<eT>& x, std::istream& f, std::string& err_msg);
  template<typename eT> inline static bool load_auto_detect(Cube<eT>& x, std::istream& f, std::string& err_msg);
  
  
  //
  // field saving and loading
  
  template<typename T1> inline static bool save_arma_binary(const field<T1>& x, const std::string&  name);
  template<typename T1> inline static bool save_arma_binary(const field<T1>& x,       std::ostream& f);
  
  template<typename T1> inline static bool load_arma_binary(      field<T1>& x, const std::string&  name, std::string& err_msg);
  template<typename T1> inline static bool load_arma_binary(      field<T1>& x,       std::istream& f,    std::string& err_msg);
  
  template<typename T1> inline static bool load_auto_detect(      field<T1>& x, const std::string&  name, std::string& err_msg);
  template<typename T1> inline static bool load_auto_detect(      field<T1>& x,       std::istream& f,    std::string& err_msg);
  
  inline static bool save_std_string(const field<std::string>& x, const std::string&  name);
  inline static bool save_std_string(const field<std::string>& x,       std::ostream& f);
  
  inline static bool load_std_string(      field<std::string>& x, const std::string&  name, std::string& err_msg);
  inline static bool load_std_string(      field<std::string>& x,       std::istream& f,    std::string& err_msg);
  


  //
  // handling of PPM images by cubes

  template<typename T1> inline static bool save_ppm_binary(const Cube<T1>& x, const std::string&  final_name);
  template<typename T1> inline static bool save_ppm_binary(const Cube<T1>& x,       std::ostream& f);
  
  template<typename T1> inline static bool load_ppm_binary(      Cube<T1>& x, const std::string&  final_name, std::string& err_msg);
  template<typename T1> inline static bool load_ppm_binary(      Cube<T1>& x,       std::istream& f,          std::string& err_msg);


  //
  // handling of PPM images by fields

  template<typename T1> inline static bool save_ppm_binary(const field<T1>& x, const std::string&  final_name);
  template<typename T1> inline static bool save_ppm_binary(const field<T1>& x,       std::ostream& f);
  
  template<typename T1> inline static bool load_ppm_binary(      field<T1>& x, const std::string&  final_name, std::string& err_msg);
  template<typename T1> inline static bool load_ppm_binary(      field<T1>& x,       std::istream& f,          std::string& err_msg);
  


  };



//! @}
