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


//! \addtogroup arma_ostream
//! @{



inline
arma_ostream_state::arma_ostream_state(const std::ostream& o)
  : orig_flags    (o.flags())
  , orig_precision(o.precision())
  , orig_width    (o.width())
  , orig_fill     (o.fill())
  {
  }



inline
void
arma_ostream_state::restore(std::ostream& o) const
  {
  o.flags    (orig_flags);
  o.precision(orig_precision);
  o.width    (orig_width);
  o.fill     (orig_fill);
  }



//
//



template<typename eT>
inline
std::streamsize
arma_ostream::modify_stream(std::ostream& o, const eT* data, const uword n_elem)
  {
  o.unsetf(ios::showbase);
  o.unsetf(ios::uppercase);
  o.unsetf(ios::showpos);
  
  o.fill(' ');
  
  std::streamsize cell_width;
  
  bool use_layout_B = false;
  bool use_layout_C = false;
  bool use_layout_D = false;
  
  for(uword i=0; i<n_elem; ++i)
    {
    const eT val = data[i];
    
    if(arma_isfinite(val) == false)  { continue; }
    
    if(
      ( cond_rel< (sizeof(eT) > 4) && (is_same_type<uword,eT>::yes || is_same_type<sword,eT>::yes) >::geq(val, eT(+10000000000)) )
      ||
      ( cond_rel< (sizeof(eT) > 4) &&  is_same_type<sword,eT>::yes                                 >::leq(val, eT(-10000000000)) )
      )
      {
      use_layout_D = true;
      break;
      }
    
    if(
      ( val >= eT(+100) )
      ||
      //( (is_signed<eT>::value) && (val <= eT(-100)) ) ||
      //( (is_non_integral<eT>::value) && (val > eT(0)) && (val <= eT(+1e-4)) ) ||
      //( (is_non_integral<eT>::value) && (is_signed<eT>::value) && (val < eT(0)) && (val >= eT(-1e-4)) ) 
        (
        cond_rel< is_signed<eT>::value >::leq(val, eT(-100))
        )
      ||
        (
        cond_rel< is_non_integral<eT>::value >::gt(val,  eT(0))
        &&
        cond_rel< is_non_integral<eT>::value >::leq(val, eT(+1e-4))
        )
      ||
        (
        cond_rel< is_non_integral<eT>::value && is_signed<eT>::value >::lt(val, eT(0))
        &&
        cond_rel< is_non_integral<eT>::value && is_signed<eT>::value >::geq(val, eT(-1e-4))
        )
      )
      {
      use_layout_C = true;
      break;
      }
      
    if(
      // (val >= eT(+10)) || ( (is_signed<eT>::value) && (val <= eT(-10)) )
      (val >= eT(+10)) || ( cond_rel< is_signed<eT>::value >::leq(val, eT(-10)) )
      )
      {
      use_layout_B = true;
      }
    }
  
  if(use_layout_D)
    {
    o.setf(ios::scientific);
    o.setf(ios::right);
    o.unsetf(ios::fixed);
    o.precision(4);
    cell_width = 21;
    }
  else
  if(use_layout_C)
    {
    o.setf(ios::scientific);
    o.setf(ios::right);
    o.unsetf(ios::fixed);
    o.precision(4);
    cell_width = 13;
    }
  else
  if(use_layout_B)
    {
    o.unsetf(ios::scientific);
    o.setf(ios::right);
    o.setf(ios::fixed);
    o.precision(4);
    cell_width = 10;
    }
  else
    {
    o.unsetf(ios::scientific);
    o.setf(ios::right);
    o.setf(ios::fixed);
    o.precision(4);
    cell_width = 9;
    }
  
  return cell_width;
  }



//! "better than nothing" settings for complex numbers
template<typename T>
inline
std::streamsize
arma_ostream::modify_stream(std::ostream& o, const std::complex<T>* data, const uword n_elem)
  {
  arma_ignore(data);
  arma_ignore(n_elem);
  
  o.unsetf(ios::showbase);
  o.unsetf(ios::uppercase);
  o.fill(' ');
  
  o.setf(ios::scientific);
  o.setf(ios::showpos);
  o.setf(ios::right);
  o.unsetf(ios::fixed);
  
  std::streamsize cell_width;
  
  o.precision(3);
  cell_width = 2 + 2*(1 + 3 + o.precision() + 5) + 1;
  
  return cell_width;
  }


template<typename eT>
inline
std::streamsize
arma_ostream::modify_stream(std::ostream& o, typename SpMat<eT>::const_iterator begin, const uword n_elem, const typename arma_not_cx<eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  o.unsetf(ios::showbase);
  o.unsetf(ios::uppercase);
  o.unsetf(ios::showpos);
  
  o.fill(' ');
  
  std::streamsize cell_width;
  
  bool use_layout_B  = false;
  bool use_layout_C  = false;
  
  for(typename SpMat<eT>::const_iterator it = begin; it.pos() < n_elem; ++it)
    {
    const eT val = (*it);
    
    if(arma_isfinite(val) == false)  { continue; }
    
    if(
      val >= eT(+100) ||
      ( (is_signed<eT>::value) && (val <= eT(-100)) ) ||
      ( (is_non_integral<eT>::value) && (val > eT(0)) && (val <= eT(+1e-4)) ) ||
      ( (is_non_integral<eT>::value) && (is_signed<eT>::value) && (val < eT(0)) && (val >= eT(-1e-4)) )
      )
      {
      use_layout_C = true;
      break;
      }
    
    if(
      (val >= eT(+10)) || ( (is_signed<eT>::value) && (val <= eT(-10)) )
      )
      {
      use_layout_B = true;
      }
    }
  
  if(use_layout_C)
    {
    o.setf(ios::scientific);
    o.setf(ios::right);
    o.unsetf(ios::fixed);
    o.precision(4);
    cell_width = 13;
    }
  else
  if(use_layout_B)
    {
    o.unsetf(ios::scientific);
    o.setf(ios::right);
    o.setf(ios::fixed);
    o.precision(4);
    cell_width = 10;
    }
  else
    {
    o.unsetf(ios::scientific);
    o.setf(ios::right);
    o.setf(ios::fixed);
    o.precision(4);
    cell_width = 9;
    }
  
  return cell_width;
  }



//! "better than nothing" settings for complex numbers
template<typename eT>
inline
std::streamsize
arma_ostream::modify_stream(std::ostream& o, typename SpMat<eT>::const_iterator begin, const uword n_elem, const typename arma_cx_only<eT>::result* junk)
  {
  arma_ignore(begin);
  arma_ignore(n_elem);
  arma_ignore(junk);
  
  o.unsetf(ios::showbase);
  o.unsetf(ios::uppercase);
  o.fill(' ');
  
  o.setf(ios::scientific);
  o.setf(ios::showpos);
  o.setf(ios::right);
  o.unsetf(ios::fixed);
  
  std::streamsize cell_width;
  
  o.precision(3);
  cell_width = 2 + 2*(1 + 3 + o.precision() + 5) + 1;
  
  return cell_width;
  }



template<typename eT>
inline
void
arma_ostream::print_elem_zero(std::ostream& o, const bool modify)
  {
  typedef typename promote_type<eT, s16>::result promoted_eT;
  
  if(modify)
    {
    const ios::fmtflags   save_flags     = o.flags();
    const std::streamsize save_precision = o.precision();
    
    o.unsetf(ios::scientific);
    o.setf(ios::fixed);
    o.precision(0);
    
    o << promoted_eT(0);
    
    o.flags(save_flags);
    o.precision(save_precision);
    }
  else
    {
    o << promoted_eT(0);
    }
  }



template<typename eT>
inline
void
arma_ostream::print_elem(std::ostream& o, const eT& x, const bool modify)
  {
  if(x == eT(0))
    {
    arma_ostream::print_elem_zero<eT>(o, modify);
    }
  else
    {
    arma_ostream::raw_print_elem(o, x);
    }
  }



template<typename eT>
inline
void
arma_ostream::raw_print_elem(std::ostream& o, const eT& x)
  {
  if(is_signed<eT>::value)
    {
    typedef typename promote_type<eT, s16>::result promoted_eT;
    
    if(arma_isfinite(x))
      {
      o << promoted_eT(x);
      }
    else
      {
      o << ( arma_isinf(x) ? ((x <= eT(0)) ? "-inf" : "inf") : "nan" );
      }
    }
  else
    {
    typedef typename promote_type<eT, u16>::result promoted_eT;
    
    o << promoted_eT(x);
    }
  }



template<typename T>
inline
void
arma_ostream::print_elem(std::ostream& o, const std::complex<T>& x, const bool modify)
  {
  if( (x.real() == T(0)) && (x.imag() == T(0)) && (modify) )
    {
    o << "(0,0)";
    }
  else
    {
    arma_ostream::raw_print_elem(o, x);
    }
  }



template<typename T>
inline
void
arma_ostream::raw_print_elem(std::ostream& o, const std::complex<T>& x)
  {
  std::ostringstream ss;
  ss.flags(o.flags());
  //ss.imbue(o.getloc());
  ss.precision(o.precision());
  
  ss << '(';
  
  const T a = x.real();
  
  if(arma_isfinite(a))
    {
    ss << a;
    }
  else
    {
    ss << ( arma_isinf(a) ? ((a <= T(0)) ? "-inf" : "+inf") : "nan" );
    }
  
  ss << ',';
  
  const T b = x.imag();
  
  if(arma_isfinite(b))
    {
    ss << b;
    }
  else
    {
    ss << ( arma_isinf(b) ? ((b <= T(0)) ? "-inf" : "+inf") : "nan" );
    }
  
  ss << ')';
  
  o << ss.str();
  }



//! Print a matrix to the specified stream
template<typename eT>
arma_cold
inline
void
arma_ostream::print(std::ostream& o, const Mat<eT>& m, const bool modify)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(o);
  
  const std::streamsize cell_width = modify ? arma_ostream::modify_stream(o, m.memptr(), m.n_elem) : o.width();
  
  const uword m_n_rows = m.n_rows;
  const uword m_n_cols = m.n_cols;
  
  if(m.is_empty() == false)
    {
    if(m_n_cols > 0)
      {
      if(cell_width > 0)
        {
        for(uword row=0; row < m_n_rows; ++row)
          {
          for(uword col=0; col < m_n_cols; ++col)
            {
            // the cell width appears to be reset after each element is printed,
            // hence we need to restore it
            o.width(cell_width);
            arma_ostream::print_elem(o, m.at(row,col), modify);
            }
        
          o << '\n';
          }
        }
      else
        {
        for(uword row=0; row < m_n_rows; ++row)
          {
          for(uword col=0; col < m_n_cols-1; ++col)
            {
            arma_ostream::print_elem(o, m.at(row,col), modify);
            o << ' ';
            }
        
          arma_ostream::print_elem(o, m.at(row, m_n_cols-1), modify);
          o << '\n';
          }
        }
      }
    }
  else
    {
    if(modify)
      {
      o.unsetf(ios::showbase);
      o.unsetf(ios::uppercase);
      o.unsetf(ios::showpos);
      o.setf(ios::fixed);
      }
    
    o << "[matrix size: " << m_n_rows << 'x' << m_n_cols << "]\n";
    }
  
  o.flush();
  stream_state.restore(o);
  }



//! Print a cube to the specified stream
template<typename eT>
arma_cold
inline
void
arma_ostream::print(std::ostream& o, const Cube<eT>& x, const bool modify)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(o);
  
  if(x.is_empty() == false)
    {
    for(uword slice=0; slice < x.n_slices; ++slice)
      {
      const Mat<eT> tmp(const_cast<eT*>(x.slice_memptr(slice)), x.n_rows, x.n_cols, false);
      
      o << "[cube slice: " << slice << ']' << '\n';
      arma_ostream::print(o, tmp, modify);
      
      if((slice+1) < x.n_slices)  { o << '\n'; }
      }
    }
  else
    {
    if(modify)
      {
      o.unsetf(ios::showbase);
      o.unsetf(ios::uppercase);
      o.unsetf(ios::showpos);
      o.setf(ios::fixed);
      }
    
    o << "[cube size: " << x.n_rows << 'x' << x.n_cols << 'x' << x.n_slices <<  "]\n";
    }
  
  stream_state.restore(o);
  }




//! Print a field to the specified stream
//! Assumes type oT can be printed, ie. oT has std::ostream& operator<< (std::ostream&, const oT&) 
template<typename oT>
arma_cold
inline
void
arma_ostream::print(std::ostream& o, const field<oT>& x)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(o);
  
  const std::streamsize cell_width = o.width();
  
  const uword x_n_rows   = x.n_rows;
  const uword x_n_cols   = x.n_cols;
  const uword x_n_slices = x.n_slices;
  
  if(x.is_empty() == false)
    {
    if(x_n_slices == 1)
      {
      for(uword col=0; col < x_n_cols; ++col)
        {
        o << "[field column: " << col << ']' << '\n'; 
        
        for(uword row=0; row < x_n_rows; ++row)
          {
          o.width(cell_width);
          o << x.at(row,col) << '\n';
          }
        
        o << '\n';
        }
      }
    else
      {
      for(uword slice=0; slice < x_n_slices; ++slice)
        {
        o << "[field slice: " << slice << ']' << '\n';
        
        for(uword col=0; col < x_n_cols; ++col)
          {
          o << "[field column: " << col << ']' << '\n';
          
          for(uword row=0; row < x_n_rows; ++row)
            {
            o.width(cell_width);
            o << x.at(row,col,slice) << '\n';
            }
          
          o << '\n';
          }
        
        o << '\n';
        }
      }
    }
  else
    {
    o.unsetf(ios::showbase);
    o.unsetf(ios::uppercase);
    o.unsetf(ios::showpos);
    o.setf(ios::fixed);
    
    o << "[field size: " << x_n_rows << 'x' << x_n_cols << 'x' << x_n_slices << "]\n";
    }
  
  o.flush();
  stream_state.restore(o);
  }



//! Print a subfield to the specified stream
//! Assumes type oT can be printed, ie. oT has std::ostream& operator<< (std::ostream&, const oT&) 
template<typename oT>
arma_cold
inline
void
arma_ostream::print(std::ostream& o, const subview_field<oT>& x)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(o);
  
  const std::streamsize cell_width = o.width();
  
  const uword x_n_rows   = x.n_rows;
  const uword x_n_cols   = x.n_cols;
  const uword x_n_slices = x.n_slices;
  
  if(x.is_empty() == false)
    {
    if(x_n_slices == 1)
      {
      for(uword col=0; col < x_n_cols; ++col)
        {
        o << "[field column: " << col << ']' << '\n'; 
        for(uword row=0; row<x_n_rows; ++row)
          {
          o.width(cell_width);
          o << x.at(row,col) << '\n';
          }
        
        o << '\n';
        }
      }
    else
      {
      for(uword slice=0; slice < x_n_slices; ++slice)
        {
        o << "[field slice: " << slice << ']' << '\n';
        
        for(uword col=0; col < x_n_cols; ++col)
          {
          o << "[field column: " << col << ']' << '\n';
          
          for(uword row=0; row < x_n_rows; ++row)
            {
            o.width(cell_width);
            o << x.at(row,col,slice) << '\n';
            }
          
          o << '\n';
          }
        
        o << '\n';
        }
      }
    }
  else
    {
    o.unsetf(ios::showbase);
    o.unsetf(ios::uppercase);
    o.unsetf(ios::showpos);
    o.setf(ios::fixed);
    
    o << "[field size: " << x_n_rows << 'x' << x_n_cols << 'x' << x_n_slices << "]\n";
    }
  
  o.flush();
  stream_state.restore(o);
  }



template<typename eT>
arma_cold
inline
void
arma_ostream::print_dense(std::ostream& o, const SpMat<eT>& m, const bool modify)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(o);
  
  std::streamsize cell_width = o.width();
  
  if(modify)
    {
    if(m.n_nonzero > 0)
      {
      cell_width = arma_ostream::modify_stream<eT>(o, m.begin(), m.n_nonzero);
      }
    else
      {
      eT tmp[1];  tmp[0] = eT(0);
      
      cell_width = arma_ostream::modify_stream(o, &tmp[0], 1);
      }
    }
  
  const uword m_n_rows = m.n_rows;
  const uword m_n_cols = m.n_cols;
  
  if(m.is_empty() == false)
    {
    if(m_n_cols > 0)
      {
      if(cell_width > 0)
        {
        for(uword row=0; row < m_n_rows; ++row)
          {
          for(uword col=0; col < m_n_cols; ++col)
            {
            // the cell width appears to be reset after each element is printed,
            // hence we need to restore it
            o.width(cell_width);
            arma_ostream::print_elem(o, m.at(row,col), modify);
            }
          
          o << '\n';
          }
        }
      else
        {
        for(uword row=0; row < m_n_rows; ++row)
          {
          for(uword col=0; col < m_n_cols-1; ++col)
            {
            arma_ostream::print_elem(o, m.at(row,col), modify);
            o << ' ';
            }
          
          arma_ostream::print_elem(o, m.at(row, m_n_cols-1), modify);
          o << '\n';
          }
        }
      }
    }
  else
    {
    if(modify)
      {
      o.unsetf(ios::showbase);
      o.unsetf(ios::uppercase);
      o.unsetf(ios::showpos);
      o.setf(ios::fixed);
      }
    
    o << "[matrix size: " << m_n_rows << 'x' << m_n_cols << "]\n";
    }
  
  o.flush();
  stream_state.restore(o);
  }



template<typename eT>
arma_cold
inline
void
arma_ostream::print(std::ostream& o, const SpMat<eT>& m, const bool modify)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(o);
  
  o.unsetf(ios::showbase);
  o.unsetf(ios::uppercase);
  o.unsetf(ios::showpos);
  o.unsetf(ios::scientific);
  o.setf(ios::right);
  o.setf(ios::fixed);
  
  const uword  m_n_nonzero = m.n_nonzero;
  const double density     = (m.n_elem > 0) ? (double(m_n_nonzero) / double(m.n_elem) * double(100)) : double(0);
  
  o << "[matrix size: " << m.n_rows << 'x' << m.n_cols << "; n_nonzero: " << m_n_nonzero;
  
  if(density == double(0))
    {
    o.precision(0);
    }
  else
  if(density >= (double(10.0)-std::numeric_limits<double>::epsilon()))
    {
    o.precision(1);
    }
  else
  if(density > (double(0.01)-std::numeric_limits<double>::epsilon()))
    {
    o.precision(2);
    }
  else
  if(density > (double(0.001)-std::numeric_limits<double>::epsilon()))
    {
    o.precision(3);
    }
  else
  if(density > (double(0.0001)-std::numeric_limits<double>::epsilon()))
    {
    o.precision(4);
    }
  else
    {
    o.unsetf(ios::fixed);
    o.setf(ios::scientific);
    o.precision(2);
    }
  
  o << "; density: " << density  << "%]\n\n";
  
  if(modify == false) { stream_state.restore(o); }
  
  if(m_n_nonzero > 0)
    {
    const std::streamsize cell_width = modify ? arma_ostream::modify_stream<eT>(o, m.begin(), m_n_nonzero) : o.width();
    
    typename SpMat<eT>::const_iterator it     = m.begin();
    typename SpMat<eT>::const_iterator it_end = m.end();
    
    while(it != it_end)
      {
      const uword row = it.row();
      const uword col = it.col();
      
      // TODO: change the maximum number of spaces before and after each location to be dependent on n_rows and n_cols
      
           if(row < 10)       { o << "      "; }
      else if(row < 100)      { o << "     ";  }
      else if(row < 1000)     { o << "    ";   }
      else if(row < 10000)    { o << "   ";    }
      else if(row < 100000)   { o << "  ";     }
      else if(row < 1000000)  { o << ' ';      }
      
      o << '(' << row << ", " << col << ") ";
      
           if(col < 10)       { o << "      "; }
      else if(col < 100)      { o << "     ";  }
      else if(col < 1000)     { o << "    ";   }
      else if(col < 10000)    { o << "   ";    }
      else if(col < 100000)   { o << "  ";     }
      else if(col < 1000000)  { o << ' ';      }
      
      if(cell_width > 0) { o.width(cell_width); }
        
      arma_ostream::print_elem(o, eT(*it), modify);
      o << '\n';
      
      ++it;
      }
    
    o << '\n';
    }
  
  o.flush();
  stream_state.restore(o);
  }



arma_cold
inline
void
arma_ostream::print(std::ostream& o, const SizeMat& S)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(o);
  
  o.unsetf(ios::showbase);
  o.unsetf(ios::uppercase);
  o.unsetf(ios::showpos);
  
  o.setf(ios::fixed);
  
  o << S.n_rows << 'x' << S.n_cols;
  
  stream_state.restore(o);
  }



arma_cold
inline
void
arma_ostream::print(std::ostream& o, const SizeCube& S)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(o);
  
  o.unsetf(ios::showbase);
  o.unsetf(ios::uppercase);
  o.unsetf(ios::showpos);
  
  o.setf(ios::fixed);
    
  o << S.n_rows << 'x' << S.n_cols << 'x' << S.n_slices;
  
  stream_state.restore(o);
  }



template<typename eT>
arma_cold
inline
void
arma_ostream::brief_print(std::ostream& o, const Mat<eT>& m, const bool print_size)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(o);
  
  if(print_size)
    {
    o.unsetf(ios::showbase);
    o.unsetf(ios::uppercase);
    o.unsetf(ios::showpos);
    o.setf(ios::fixed);
    
    o << "[matrix size: " << m.n_rows << 'x' << m.n_cols << "]\n";
    }
  
  if(m.n_elem == 0)  { o.flush(); stream_state.restore(o); return; }
  
  if((m.n_rows <= 5) && (m.n_cols <= 5))  { arma_ostream::print(o, m, true); return; }
  
  const bool print_row_ellipsis = (m.n_rows >= 6);
  const bool print_col_ellipsis = (m.n_cols >= 6);
  
  if( (print_row_ellipsis == true) && (print_col_ellipsis == true) )
    {
    Mat<eT> X(4, 4, arma_nozeros_indicator());
    
    X( span(0,2), span(0,2) ) = m( span(0,2),  span(0,2)  );  // top left submatrix
    X( 3,         span(0,2) ) = m( m.n_rows-1, span(0,2)  );  // truncated last row
    X( span(0,2), 3         ) = m( span(0,2),  m.n_cols-1 );  // truncated last column
    X( 3,         3         ) = m( m.n_rows-1, m.n_cols-1 );  // bottom right element
    
    const std::streamsize cell_width = arma_ostream::modify_stream(o, X.memptr(), X.n_elem);
    
    for(uword row=0; row <= 2; ++row)
      {
      for(uword col=0; col <= 2; ++col)
        {
        o.width(cell_width);
        arma_ostream::print_elem(o, X.at(row,col), true);
        }
      
      o.width(6);
      o << "...";
      
      o.width(cell_width);
      arma_ostream::print_elem(o, X.at(row,3), true);
      o << '\n';
      }
    
    for(uword col=0; col <= 2; ++col)
      {
      o.width(cell_width);
      o << ':';
      }
    
    o.width(6);
    o << "...";
    
    o.width(cell_width);
    o << ':' << '\n';
    
    const uword row = 3;
      {
      for(uword col=0; col <= 2; ++col)
        {
        o.width(cell_width);
        arma_ostream::print_elem(o, X.at(row,col), true);
        }
      
      o.width(6);
      o << "...";
      
      o.width(cell_width);
      arma_ostream::print_elem(o, X.at(row,3), true);
      o << '\n';
      }
    }
  
  
  if( (print_row_ellipsis == true) && (print_col_ellipsis == false) )
    {
    Mat<eT> X(4, m.n_cols, arma_nozeros_indicator());
    
    X( span(0,2), span::all ) = m( span(0,2),  span::all );  // top
    X( 3,         span::all ) = m( m.n_rows-1, span::all );  // bottom
    
    const std::streamsize cell_width = arma_ostream::modify_stream(o, X.memptr(), X.n_elem);
    
    for(uword row=0; row <= 2; ++row)  // first 3 rows
      {
      for(uword col=0; col < m.n_cols; ++col)
        {
        o.width(cell_width);
        arma_ostream::print_elem(o, X.at(row,col), true);
        }
      
      o << '\n';
      }
    
    for(uword col=0; col < m.n_cols; ++col)
      {
      o.width(cell_width);
      o << ':';
      }
    
    o.width(cell_width);
    o << '\n';
    
    const uword row = 3;
      {
      for(uword col=0; col < m.n_cols; ++col)
        {
        o.width(cell_width);
        arma_ostream::print_elem(o, X.at(row,col), true);
        }
      }
    
    o << '\n';
    }
  
  
  if( (print_row_ellipsis == false) && (print_col_ellipsis == true) )
    {
    Mat<eT> X(m.n_rows, 4, arma_nozeros_indicator());
    
    X( span::all, span(0,2) ) = m( span::all, span(0,2)  );  // left
    X( span::all, 3         ) = m( span::all, m.n_cols-1 );  // right
    
    const std::streamsize cell_width = arma_ostream::modify_stream(o, X.memptr(), X.n_elem);
    
    for(uword row=0; row < m.n_rows; ++row)
      {
      for(uword col=0; col <= 2; ++col)
        {
        o.width(cell_width);
        arma_ostream::print_elem(o, X.at(row,col), true);
        }
      
      o.width(6);
      o << "...";
      
      o.width(cell_width);
      arma_ostream::print_elem(o, X.at(row,3), true);
      o << '\n';
      }
    }
  
  
  o.flush();
  stream_state.restore(o);
  }



template<typename eT>
arma_cold
inline
void
arma_ostream::brief_print(std::ostream& o, const Cube<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(o);
  
  o.unsetf(ios::showbase);
  o.unsetf(ios::uppercase);
  o.unsetf(ios::showpos);
  o.setf(ios::fixed);
  
  o << "[cube size: " << x.n_rows << 'x' << x.n_cols << 'x' << x.n_slices << "]\n";
  
  if(x.n_elem == 0)  { o.flush(); stream_state.restore(o); return; }
  
  if(x.n_slices <= 3)
    {
    for(uword slice=0; slice < x.n_slices; ++slice)
      {
      const Mat<eT> tmp(const_cast<eT*>(x.slice_memptr(slice)), x.n_rows, x.n_cols, false);
      
      o << "[cube slice: " << slice << ']' << '\n';
      arma_ostream::brief_print(o, tmp, false);
      
      if((slice+1) < x.n_slices)  { o << '\n'; }
      }
    }
  else
    {
    for(uword slice=0; slice <= 1; ++slice)
      {
      const Mat<eT> tmp(const_cast<eT*>(x.slice_memptr(slice)), x.n_rows, x.n_cols, false);
      
      o << "[cube slice: " << slice << ']' << '\n';
      arma_ostream::brief_print(o, tmp, false);
      o << '\n';
      }
      
    o << "[cube slice: ...]\n\n";
    
    const uword slice = x.n_slices-1;
      {
      const Mat<eT> tmp(const_cast<eT*>(x.slice_memptr(slice)), x.n_rows, x.n_cols, false);
      
      o << "[cube slice: " << slice << ']' << '\n';
      arma_ostream::brief_print(o, tmp, false);
      }
    }
  
  stream_state.restore(o);
  }



template<typename eT>
arma_cold
inline
void
arma_ostream::brief_print(std::ostream& o, const SpMat<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  if(m.n_nonzero <= 10)  { arma_ostream::print(o, m, true); return; }
  
  const arma_ostream_state stream_state(o);
  
  o.unsetf(ios::showbase);
  o.unsetf(ios::uppercase);
  o.unsetf(ios::showpos);
  o.unsetf(ios::scientific);
  o.setf(ios::right);
  o.setf(ios::fixed);
  
  const uword  m_n_nonzero = m.n_nonzero;
  const double density     = (m.n_elem > 0) ? (double(m_n_nonzero) / double(m.n_elem) * double(100)) : double(0);
  
  o << "[matrix size: " << m.n_rows << 'x' << m.n_cols << "; n_nonzero: " << m_n_nonzero;
  
  if(density == double(0))
    {
    o.precision(0);
    }
  else
  if(density >= (double(10.0)-std::numeric_limits<double>::epsilon()))
    {
    o.precision(1);
    }
  else
  if(density > (double(0.01)-std::numeric_limits<double>::epsilon()))
    {
    o.precision(2);
    }
  else
  if(density > (double(0.001)-std::numeric_limits<double>::epsilon()))
    {
    o.precision(3);
    }
  else
  if(density > (double(0.0001)-std::numeric_limits<double>::epsilon()))
    {
    o.precision(4);
    }
  else
    {
    o.unsetf(ios::fixed);
    o.setf(ios::scientific);
    o.precision(2);
    }
  
  o << "; density: " << density  << "%]\n\n";
  
  // get the first 9 elements and the last element
  
  typename SpMat<eT>::const_iterator it     = m.begin();
  typename SpMat<eT>::const_iterator it_end = m.end();
  
  uvec    storage_row(10);
  uvec    storage_col(10);
  Col<eT> storage_val(10);
  
  uword count = 0;
  
  while( (it != it_end) && (count < 9) )
    {
    storage_row(count) = it.row();
    storage_col(count) = it.col();
    storage_val(count) = (*it);
    
    ++it;
    ++count;
    }
  
  it = it_end;
  --it;
  
  storage_row(count) = it.row();
  storage_col(count) = it.col();
  storage_val(count) = (*it);
  
  const std::streamsize cell_width = arma_ostream::modify_stream(o, storage_val.memptr(), 10);
  
  for(uword i=0; i < 9; ++i)
    {
    const uword row = storage_row(i);
    const uword col = storage_col(i);
    
         if(row < 10)       { o << "      "; }
    else if(row < 100)      { o << "     ";  }
    else if(row < 1000)     { o << "    ";   }
    else if(row < 10000)    { o << "   ";    }
    else if(row < 100000)   { o << "  ";     }
    else if(row < 1000000)  { o << ' ';      }
    
    o << '(' << row << ", " << col << ") ";
    
         if(col < 10)       { o << "      "; }
    else if(col < 100)      { o << "     ";  }
    else if(col < 1000)     { o << "    ";   }
    else if(col < 10000)    { o << "   ";    }
    else if(col < 100000)   { o << "  ";     }
    else if(col < 1000000)  { o << ' ';      }
    
    if(cell_width > 0)  { o.width(cell_width); }
    
    arma_ostream::print_elem(o, storage_val(i), true);
    o << '\n';
    }
  
  o << "      (:, :)     ";
  if(cell_width > 0)  { o.width(cell_width); }
  o << "...\n";
  
  
  const uword i = 9;
    {
    const uword row = storage_row(i);
    const uword col = storage_col(i);
    
         if(row < 10)       { o << "      "; }
    else if(row < 100)      { o << "     ";  }
    else if(row < 1000)     { o << "    ";   }
    else if(row < 10000)    { o << "   ";    }
    else if(row < 100000)   { o << "  ";     }
    else if(row < 1000000)  { o << ' ';      }
    
    o << '(' << row << ", " << col << ") ";
    
         if(col < 10)       { o << "      "; }
    else if(col < 100)      { o << "     ";  }
    else if(col < 1000)     { o << "    ";   }
    else if(col < 10000)    { o << "   ";    }
    else if(col < 100000)   { o << "  ";     }
    else if(col < 1000000)  { o << ' ';      }
    
    if(cell_width > 0) { o.width(cell_width); }
      
    arma_ostream::print_elem(o, storage_val(i), true);
    o << '\n';
    }
  
  o.flush();
  stream_state.restore(o);
  }



//! @}
