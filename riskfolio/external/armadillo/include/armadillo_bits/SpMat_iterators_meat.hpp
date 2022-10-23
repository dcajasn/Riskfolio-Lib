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


//! \addtogroup SpMat
//! @{


///////////////////////////////////////////////////////////////////////////////
// SpMat::iterator_base implementation                                       //
///////////////////////////////////////////////////////////////////////////////


template<typename eT>
inline
SpMat<eT>::iterator_base::iterator_base()
  : M(nullptr)
  , internal_col(0)
  , internal_pos(0)
  {
  // Technically this iterator is invalid (it does not point to a valid element)
  }



template<typename eT>
inline
SpMat<eT>::iterator_base::iterator_base(const SpMat<eT>& in_M)
  : M(&in_M)
  , internal_col(0)
  , internal_pos(0)
  {
  // Technically this iterator is invalid (it may not point to a valid element)
  }



template<typename eT>
inline
SpMat<eT>::iterator_base::iterator_base(const SpMat<eT>& in_M, const uword in_col, const uword in_pos)
  : M(&in_M)
  , internal_col(in_col)
  , internal_pos(in_pos)
  {
  // Nothing to do.
  }



template<typename eT>
arma_inline
eT
SpMat<eT>::iterator_base::operator*() const
  {
  return M->values[internal_pos];
  }



///////////////////////////////////////////////////////////////////////////////
// SpMat::const_iterator implementation                                      //
///////////////////////////////////////////////////////////////////////////////

template<typename eT>
inline
SpMat<eT>::const_iterator::const_iterator()
  : iterator_base()
  {
  }



template<typename eT>
inline
SpMat<eT>::const_iterator::const_iterator(const SpMat<eT>& in_M, uword initial_pos)
  : iterator_base(in_M, 0, initial_pos)
  {
  // Corner case for empty matrices.
  if(in_M.n_nonzero == 0)
    {
    iterator_base::internal_col = in_M.n_cols;
    return;
    }
  
  // Determine which column we should be in.
  while(iterator_base::M->col_ptrs[iterator_base::internal_col + 1] <= iterator_base::internal_pos)
    {
    iterator_base::internal_col++;
    }
  }



template<typename eT>
inline
SpMat<eT>::const_iterator::const_iterator(const SpMat<eT>& in_M, uword in_row, uword in_col)
  : iterator_base(in_M, in_col, 0)
  {
  // So we have a position we want to be right after.  Skip to the column.
  iterator_base::internal_pos = iterator_base::M->col_ptrs[iterator_base::internal_col];
  
  // Now we have to make sure that is the right column.
  while(iterator_base::M->col_ptrs[iterator_base::internal_col + 1] <= iterator_base::internal_pos)
    {
    iterator_base::internal_col++;
    }
  
  // Now we have to get to the right row.
  while((iterator_base::M->row_indices[iterator_base::internal_pos] < in_row) && (iterator_base::internal_col == in_col))
    {
    ++(*this); // Increment iterator.
    }
  }



template<typename eT>
inline
SpMat<eT>::const_iterator::const_iterator(const SpMat<eT>& in_M, const uword /* in_row */, const uword in_col, const uword in_pos)
  : iterator_base(in_M, in_col, in_pos)
  {
  // Nothing to do.
  }



template<typename eT>
inline
SpMat<eT>::const_iterator::const_iterator(const typename SpMat<eT>::const_iterator& other)
  : iterator_base(*other.M, other.internal_col, other.internal_pos)
  {
  // Nothing to do.
  }



template<typename eT>
inline
arma_hot
typename SpMat<eT>::const_iterator&
SpMat<eT>::const_iterator::operator++()
  {
  ++iterator_base::internal_pos;
  
  if(iterator_base::internal_pos == iterator_base::M->n_nonzero)
    {
    iterator_base::internal_col = iterator_base::M->n_cols;
    return *this;
    }
  
  // Check to see if we moved a column.
  while(iterator_base::M->col_ptrs[iterator_base::internal_col + 1] <= iterator_base::internal_pos)
    {
    ++iterator_base::internal_col;
    }
  
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::const_iterator
SpMat<eT>::const_iterator::operator++(int)
  {
  typename SpMat<eT>::const_iterator tmp(*this);
  
  ++(*this);
  
  return tmp;
  }



template<typename eT>
inline
arma_hot
typename SpMat<eT>::const_iterator&
SpMat<eT>::const_iterator::operator--()
  {
  --iterator_base::internal_pos;
  
  // First, see if we moved back a column.
  while(iterator_base::internal_pos < iterator_base::M->col_ptrs[iterator_base::internal_col])
    {
    --iterator_base::internal_col;
    }
  
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::const_iterator
SpMat<eT>::const_iterator::operator--(int)
  {
  typename SpMat<eT>::const_iterator tmp(*this);
  
  --(*this);
  
  return tmp;
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator==(const const_iterator& rhs) const
  {
  return (rhs.row() == (*this).row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator!=(const const_iterator& rhs) const
  {
  return (rhs.row() != (*this).row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator==(const typename SpSubview<eT>::const_iterator& rhs) const
  {
  return (rhs.row() == (*this).row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator!=(const typename SpSubview<eT>::const_iterator& rhs) const
  {
  return (rhs.row() != (*this).row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator==(const const_row_iterator& rhs) const
  {
  return (rhs.row() == (*this).row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator!=(const const_row_iterator& rhs) const
  {
  return (rhs.row() != (*this).row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator==(const typename SpSubview<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() == (*this).row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_iterator::operator!=(const typename SpSubview<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() != (*this).row()) || (rhs.col() != iterator_base::internal_col);
  }



///////////////////////////////////////////////////////////////////////////////
// SpMat::iterator implementation                                            //
///////////////////////////////////////////////////////////////////////////////

template<typename eT>
inline
arma_hot
SpValProxy< SpMat<eT> >
SpMat<eT>::iterator::operator*()
  {
  return SpValProxy< SpMat<eT> >(
    iterator_base::M->row_indices[iterator_base::internal_pos],
    iterator_base::internal_col,
    access::rw(*iterator_base::M),
    &access::rw(iterator_base::M->values[iterator_base::internal_pos]));
  }



template<typename eT>
inline
arma_hot
typename SpMat<eT>::iterator&
SpMat<eT>::iterator::operator++()
  {
  const_iterator::operator++();
  
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::iterator
SpMat<eT>::iterator::operator++(int)
  {
  typename SpMat<eT>::iterator tmp(*this);
  
  const_iterator::operator++();
  
  return tmp;
  }



template<typename eT>
inline
arma_hot
typename SpMat<eT>::iterator&
SpMat<eT>::iterator::operator--()
  {
  const_iterator::operator--();
  
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::iterator
SpMat<eT>::iterator::operator--(int)
  {
  typename SpMat<eT>::iterator tmp(*this);
  
  const_iterator::operator--();
  
  return tmp;
  }



///////////////////////////////////////////////////////////////////////////////
// SpMat::const_row_iterator implementation                                  //
///////////////////////////////////////////////////////////////////////////////

/**
 * Initialize the const_row_iterator.
 */

template<typename eT>
inline
SpMat<eT>::const_row_iterator::const_row_iterator()
  : iterator_base()
  , internal_row(0)
  , actual_pos(0)
  {
  }



template<typename eT>
inline
SpMat<eT>::const_row_iterator::const_row_iterator(const SpMat<eT>& in_M, uword initial_pos)
  : iterator_base(in_M, 0, initial_pos)
  , internal_row(0)
  , actual_pos(0)
  {
  // Corner case for the end of a matrix.
  if(initial_pos == in_M.n_nonzero)
    {
    iterator_base::internal_col = 0;
    internal_row = in_M.n_rows;
    actual_pos = in_M.n_nonzero;
    iterator_base::internal_pos = in_M.n_nonzero;
    
    return;
    }
  
  // We don't count zeros in our position count, so we have to find the nonzero
  // value corresponding to the given initial position.  We assume initial_pos
  // is valid.
  
  // This is irritating because we don't know where the elements are in each row.
  // What we will do is loop across all columns looking for elements in row 0
  // (and add to our sum), then in row 1, and so forth, until we get to the desired position.
  uword cur_pos = std::numeric_limits<uword>::max(); // Invalid value.
  uword cur_actual_pos = 0;
  
  for(uword row = 0; row < iterator_base::M->n_rows; ++row)
    {
    for(uword col = 0; col < iterator_base::M->n_cols; ++col)
      {
      // Find the first element with row greater than or equal to in_row.
      const uword      col_offset = iterator_base::M->col_ptrs[col    ];
      const uword next_col_offset = iterator_base::M->col_ptrs[col + 1];
      
      const uword* start_ptr = &iterator_base::M->row_indices[     col_offset];
      const uword*   end_ptr = &iterator_base::M->row_indices[next_col_offset];
      
      if(start_ptr != end_ptr)
        {
        const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, row);
        
        // This is the number of elements in the column with row index less than in_row.
        const uword offset = uword(pos_ptr - start_ptr);
        
        if(iterator_base::M->row_indices[col_offset + offset] == row)
          {
          cur_actual_pos = col_offset + offset;
          
          // Increment position portably.
          if(cur_pos == std::numeric_limits<uword>::max())
            { cur_pos = 0; }
          else
            { ++cur_pos; }
          
          // Do we terminate?
          if(cur_pos == initial_pos)
            {
            internal_row = row;
            iterator_base::internal_col = col;
            iterator_base::internal_pos = cur_pos;
            actual_pos = cur_actual_pos;
            
            return;
            }
          }
        }
      }
    }
  
  // If we got to here, then we have gone past the end of the matrix.
  // This shouldn't happen...
  iterator_base::internal_pos = iterator_base::M->n_nonzero;
  iterator_base::internal_col = 0;
  internal_row = iterator_base::M->n_rows;
  actual_pos = iterator_base::M->n_nonzero;
  }



template<typename eT>
inline
SpMat<eT>::const_row_iterator::const_row_iterator(const SpMat<eT>& in_M, uword in_row, uword in_col)
  : iterator_base(in_M, in_col, 0)
  , internal_row(0)
  , actual_pos(0)
  {
  // Start our search in the given row.  We need to find two things:
  //
  //   1. The first nonzero element (iterating by rows) after (in_row, in_col).
  //   2. The number of nonzero elements (iterating by rows) that come before
  //      (in_row, in_col).
  //
  // We'll find these simultaneously, though we will have to loop over all
  // columns.
  
  // This will hold the total number of points with rows less than in_row.
  uword cur_pos = 0;
  uword cur_min_row = iterator_base::M->n_rows;
  uword cur_min_col = 0;
  uword cur_actual_pos = 0;
  
  for(uword col = 0; col < iterator_base::M->n_cols; ++col)
    {
    // Find the first element with row greater than or equal to in_row.
    const uword      col_offset = iterator_base::M->col_ptrs[col    ];
    const uword next_col_offset = iterator_base::M->col_ptrs[col + 1];
    
    const uword* start_ptr = &iterator_base::M->row_indices[     col_offset];
    const uword*   end_ptr = &iterator_base::M->row_indices[next_col_offset];

    if(start_ptr != end_ptr)
      {
      const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, in_row);
      
      // This is the number of elements in the column with row index less than in_row.
      const uword offset = uword(pos_ptr - start_ptr);
      
      cur_pos += offset;
      
      if(pos_ptr != end_ptr)
        {
        // This is the row index of the first element in the column with row index
        // greater than or equal to in_row.
        if((*pos_ptr) < cur_min_row)
          {
          // If we are in the desired row but before the desired column,
          // we can't take this.
          if(col >= in_col)
            {
            cur_min_row = (*pos_ptr);
            cur_min_col = col;
            cur_actual_pos = col_offset + offset;
            }
          }
        }
      }
    }
  
  // Now we know what the minimum row is.
  internal_row = cur_min_row;
  iterator_base::internal_col = cur_min_col;
  iterator_base::internal_pos = cur_pos;
  actual_pos = cur_actual_pos;
  }



/**
 * Initialize the const_row_iterator from another const_row_iterator.
 */
template<typename eT>
inline
SpMat<eT>::const_row_iterator::const_row_iterator(const typename SpMat<eT>::const_row_iterator& other)
  : iterator_base(*other.M, other.internal_col, other.internal_pos)
  , internal_row(other.internal_row)
  , actual_pos(other.actual_pos)
  {
  // Nothing to do.
  }



/**
 * Increment the row_iterator.
 */
template<typename eT>
inline
arma_hot
typename SpMat<eT>::const_row_iterator&
SpMat<eT>::const_row_iterator::operator++()
  {
  // We just need to find the next nonzero element.
  iterator_base::internal_pos++;
  
  if(iterator_base::internal_pos == iterator_base::M->n_nonzero)
    {
    internal_row = iterator_base::M->n_rows;
    iterator_base::internal_col = 0;
    
    return *this;
    }
  
  // Otherwise, we need to search.  We can start in the next column and use
  // lower_bound() to find the next element.
  uword next_min_row = iterator_base::M->n_rows;
  uword next_min_col = iterator_base::M->n_cols;
  uword next_actual_pos = 0;
  
  // Search from the current column to the end of the matrix.
  for(uword col = iterator_base::internal_col + 1; col < iterator_base::M->n_cols; ++col)
    {
    // Find the first element with row greater than or equal to in_row.
    const uword      col_offset = iterator_base::M->col_ptrs[col    ];
    const uword next_col_offset = iterator_base::M->col_ptrs[col + 1];
    
    const uword* start_ptr = &iterator_base::M->row_indices[     col_offset];
    const uword*   end_ptr = &iterator_base::M->row_indices[next_col_offset];
    
    if(start_ptr != end_ptr)
      {
      // Find the first element in the column with row greater than or equal to
      // the current row.
      const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, internal_row);
      
      if(pos_ptr != end_ptr)
        {
        // We found something in the column, but is the row index correct?
        if((*pos_ptr) == internal_row)
          {
          // Exact match---so we are done.
          iterator_base::internal_col = col;
          actual_pos = col_offset + (pos_ptr - start_ptr);
          return *this;
          }
        else if((*pos_ptr) < next_min_row)
          {
          // The first element in this column is in a subsequent row, but it's
          // the minimum row we've seen so far.
          next_min_row = (*pos_ptr);
          next_min_col = col;
          next_actual_pos = col_offset + (pos_ptr - start_ptr);
          }
        else if((*pos_ptr) == next_min_row && col < next_min_col)
          {
          // The first element in this column is in a subsequent row that we
          // already have another element for, but the column index is less so
          // this element will come first.
          next_min_col = col;
          next_actual_pos = col_offset + (pos_ptr - start_ptr);
          }
        }
      }
    }
  
  // Restart the search in the next row.
  for(uword col = 0; col <= iterator_base::internal_col; ++col)
    {
    // Find the first element with row greater than or equal to in_row + 1.
    const uword      col_offset = iterator_base::M->col_ptrs[col    ];
    const uword next_col_offset = iterator_base::M->col_ptrs[col + 1];
    
    const uword* start_ptr = &iterator_base::M->row_indices[     col_offset];
    const uword*   end_ptr = &iterator_base::M->row_indices[next_col_offset];
    
    if(start_ptr != end_ptr)
      {
      const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, internal_row + 1);
      
      if(pos_ptr != end_ptr)
        {
        // We found something in the column, but is the row index correct?
        if((*pos_ptr) == internal_row + 1)
          {
          // Exact match---so we are done.
          iterator_base::internal_col = col;
          internal_row++;
          actual_pos = col_offset + (pos_ptr - start_ptr);
          return *this;
          }
        else if((*pos_ptr) < next_min_row)
          {
          // The first element in this column is in a subsequent row,
          // but it's the minimum row we've seen so far.
          next_min_row = (*pos_ptr);
          next_min_col = col;
          next_actual_pos = col_offset + (pos_ptr - start_ptr);
          }
        else if((*pos_ptr) == next_min_row && col < next_min_col)
          {
          // The first element in this column is in a subsequent row that we
          // already have another element for, but the column index is less so
          // this element will come first.
          next_min_col = col;
          next_actual_pos = col_offset + (pos_ptr - start_ptr);
          }
        }
      }
    }
  
  iterator_base::internal_col = next_min_col;
  internal_row = next_min_row;
  actual_pos = next_actual_pos;
  
  return *this; // Now we are done.
  }



/**
 * Increment the row_iterator (but do not return anything.
 */
template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::const_row_iterator
SpMat<eT>::const_row_iterator::operator++(int)
  {
  typename SpMat<eT>::const_row_iterator tmp(*this);
  
  ++(*this);
  
  return tmp;
  }



/**
 * Decrement the row_iterator.
 */
template<typename eT>
inline
arma_hot
typename SpMat<eT>::const_row_iterator&
SpMat<eT>::const_row_iterator::operator--()
  {
  if(iterator_base::internal_pos == 0)
    {
    // Do nothing; we are already at the beginning.
    return *this;
    }
  
  iterator_base::internal_pos--;
  
  // We have to search backwards.  We'll do this by going backwards over columns
  // and seeing if we find an element in the same row.
  uword max_row = 0;
  uword max_col = 0;
  uword next_actual_pos = 0;
  
  //for(uword col = iterator_base::internal_col; col > 1; --col)
  for(uword col = iterator_base::internal_col; col >= 1; --col)
    {
    // Find the first element with row greater than or equal to in_row + 1.
    const uword      col_offset = iterator_base::M->col_ptrs[col - 1];
    const uword next_col_offset = iterator_base::M->col_ptrs[col    ];
    
    const uword* start_ptr = &iterator_base::M->row_indices[     col_offset];
    const uword*   end_ptr = &iterator_base::M->row_indices[next_col_offset];
    
    if(start_ptr != end_ptr)
      {
      // There are elements in this column.
      const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, internal_row + 1);
      
      if(pos_ptr != start_ptr)
        {
        // The element before pos_ptr is the one we are interested in.
        if(*(pos_ptr - 1) > max_row)
          {
          max_row = *(pos_ptr - 1);
          max_col = col - 1;
          next_actual_pos = col_offset + (pos_ptr - 1 - start_ptr);
          }
        else if(*(pos_ptr - 1) == max_row && (col - 1) > max_col)
          {
          max_col = col - 1;
          next_actual_pos = col_offset + (pos_ptr - 1 - start_ptr);
          }
        }
      }
    }
  
  // Now loop around to the columns at the end of the matrix.
  for(uword col = iterator_base::M->n_cols - 1; col >= iterator_base::internal_col; --col)
    {
    // Find the first element with row greater than or equal to in_row + 1.
    const uword      col_offset = iterator_base::M->col_ptrs[col    ];
    const uword next_col_offset = iterator_base::M->col_ptrs[col + 1];
    
    const uword* start_ptr = &iterator_base::M->row_indices[     col_offset];
    const uword*   end_ptr = &iterator_base::M->row_indices[next_col_offset];
    
    if(start_ptr != end_ptr)
      {
      // There are elements in this column.
      const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, internal_row);
      
      if(pos_ptr != start_ptr)
        {
        // There are elements in this column with row index < internal_row.
        if(*(pos_ptr - 1) > max_row)
          {
          max_row = *(pos_ptr - 1);
          max_col = col;
          next_actual_pos = col_offset + (pos_ptr - 1 - start_ptr);
          }
        else if(*(pos_ptr - 1) == max_row && col > max_col)
          {
          max_col = col;
          next_actual_pos = col_offset + (pos_ptr - 1 - start_ptr);
          }
        }
      }
    
    if(col == 0) // Catch edge case that the loop termination condition won't.
      {
      break;
      }
    }
  
  iterator_base::internal_col = max_col;
  internal_row = max_row;
  actual_pos = next_actual_pos;
  
  return *this;
  }



/**
 * Decrement the row_iterator.
 */
template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::const_row_iterator
SpMat<eT>::const_row_iterator::operator--(int)
  {
  typename SpMat<eT>::const_row_iterator tmp(*this);
  
  --(*this);
  
  return tmp;
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator==(const const_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator!=(const const_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator==(const typename SpSubview<eT>::const_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator!=(const typename SpSubview<eT>::const_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator==(const const_row_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator!=(const const_row_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator==(const typename SpSubview<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpMat<eT>::const_row_iterator::operator!=(const typename SpSubview<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != iterator_base::internal_col);
  }



///////////////////////////////////////////////////////////////////////////////
// SpMat::row_iterator implementation                                        //
///////////////////////////////////////////////////////////////////////////////

template<typename eT>
inline
arma_hot
SpValProxy< SpMat<eT> >
SpMat<eT>::row_iterator::operator*()
  {
  return SpValProxy< SpMat<eT> >(
      const_row_iterator::internal_row,
      iterator_base::internal_col,
      access::rw(*iterator_base::M),
      &access::rw(iterator_base::M->values[const_row_iterator::actual_pos]));
  }



template<typename eT>
inline
arma_hot
typename SpMat<eT>::row_iterator&
SpMat<eT>::row_iterator::operator++()
  {
  const_row_iterator::operator++();
  
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::row_iterator
SpMat<eT>::row_iterator::operator++(int)
  {
  typename SpMat<eT>::row_iterator tmp(*this);
  
  const_row_iterator::operator++();
  
  return tmp;
  }



template<typename eT>
inline
arma_hot
typename SpMat<eT>::row_iterator&
SpMat<eT>::row_iterator::operator--()
  {
  const_row_iterator::operator--();
  
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpMat<eT>::row_iterator
SpMat<eT>::row_iterator::operator--(int)
  {
  typename SpMat<eT>::row_iterator tmp(*this);
  
  const_row_iterator::operator--();
  
  return tmp;
  }


//! @}
