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


//! \addtogroup SpSubview
//! @{


///////////////////////////////////////////////////////////////////////////////
// SpSubview::iterator_base implementation                                   //
///////////////////////////////////////////////////////////////////////////////

template<typename eT>
inline
SpSubview<eT>::iterator_base::iterator_base(const SpSubview<eT>& in_M)
  : M(&in_M)
  , internal_col(0)
  , internal_pos(0)
  {
  // Technically this iterator is invalid (it may not point to a valid element)
  }



template<typename eT>
inline
SpSubview<eT>::iterator_base::iterator_base(const SpSubview<eT>& in_M, const uword in_col, const uword in_pos)
  : M(&in_M)
  , internal_col(in_col)
  , internal_pos(in_pos)
  {
  // Nothing to do.
  }



///////////////////////////////////////////////////////////////////////////////
// SpSubview::const_iterator implementation                                  //
///////////////////////////////////////////////////////////////////////////////

template<typename eT>
inline
SpSubview<eT>::const_iterator::const_iterator(const SpSubview<eT>& in_M, const uword initial_pos)
  : iterator_base(in_M, 0, initial_pos)
  {
  // Corner case for empty subviews.
  if(in_M.n_nonzero == 0)
    {
    iterator_base::internal_col = in_M.n_cols;
    skip_pos                    = in_M.m.n_nonzero;
    return;
    }

  // Figure out the row and column of the position.
  // lskip_pos holds the number of values which aren't part of this subview.
  const uword aux_col = iterator_base::M->aux_col1;
  const uword aux_row = iterator_base::M->aux_row1;
  const uword ln_rows = iterator_base::M->n_rows;
  const uword ln_cols = iterator_base::M->n_cols;

  uword cur_pos   = 0; // off by one because we might be searching for pos 0
  uword lskip_pos = iterator_base::M->m.col_ptrs[aux_col];
  uword cur_col   = 0;

  while(cur_pos < (iterator_base::internal_pos + 1))
    {
    // Have we stepped forward a column (or multiple columns)?
    while(((lskip_pos + cur_pos) >= iterator_base::M->m.col_ptrs[cur_col + aux_col + 1]) && (cur_col < ln_cols))
      {
      ++cur_col;
      }

    // See if the current position is in the subview.
    const uword row_index = iterator_base::M->m.row_indices[cur_pos + lskip_pos];
    if(row_index < aux_row)
      {
      ++lskip_pos; // not valid
      }
    else if(row_index < (aux_row + ln_rows))
      {
      ++cur_pos; // valid, in the subview
      }
    else
      {
      // skip to end of column
      const uword next_colptr = iterator_base::M->m.col_ptrs[cur_col + aux_col + 1];
      lskip_pos += (next_colptr - (cur_pos + lskip_pos));
      }
    }

  iterator_base::internal_col = cur_col;
  skip_pos                    = lskip_pos;
  }



template<typename eT>
inline
SpSubview<eT>::const_iterator::const_iterator(const SpSubview<eT>& in_M, const uword in_row, const uword in_col)
  : iterator_base(in_M, in_col, 0)
  {
  // Corner case for empty subviews.
  if(in_M.n_nonzero == 0)
    {
    // We must be at the last position.
    iterator_base::internal_col = in_M.n_cols;
    skip_pos                    = in_M.m.n_nonzero;
    return;
    }

  // We have a destination we want to be just after, but don't know what position that is.
  // Because we have to count the points in this subview and not in this subview, this becomes a little difficult and slow.
  const uword aux_col = iterator_base::M->aux_col1;
  const uword aux_row = iterator_base::M->aux_row1;
  const uword ln_rows = iterator_base::M->n_rows;
  const uword ln_cols = iterator_base::M->n_cols;

  uword cur_pos = 0;
  skip_pos = iterator_base::M->m.col_ptrs[aux_col];
  uword cur_col = 0;

  // Skip any empty columns.
  while(((skip_pos + cur_pos) >= iterator_base::M->m.col_ptrs[cur_col + aux_col + 1]) && (cur_col < ln_cols))
    {
    ++cur_col;
    }

  while(cur_col < in_col)
    {
    // See if the current position is in the subview.
    const uword row_index = iterator_base::M->m.row_indices[cur_pos + skip_pos];
    if(row_index < aux_row)
      {
      ++skip_pos;
      }
    else if(row_index < (aux_row + ln_rows))
      {
      ++cur_pos;
      }
    else
      {
      // skip to end of column
      const uword next_colptr = iterator_base::M->m.col_ptrs[cur_col + aux_col + 1];
      skip_pos += (next_colptr - (cur_pos + skip_pos));
      }

    // Have we stepped forward a column (or multiple columns)?
    while(((skip_pos + cur_pos) >= iterator_base::M->m.col_ptrs[cur_col + aux_col + 1]) && (cur_col < ln_cols))
      {
      ++cur_col;
      }
    }

  // Now we are either on the right column or ahead of it.
  if(cur_col == in_col)
    {
    // We have to find the right row index.
    uword row_index = iterator_base::M->m.row_indices[cur_pos + skip_pos];
    while((row_index < (in_row + aux_row)))
      {
      if(row_index < aux_row)
        {
        ++skip_pos;
        }
      else
        {
        ++cur_pos;
        }

      // Ensure we didn't step forward a column; if we did, we need to stop.
      while(((skip_pos + cur_pos) >= iterator_base::M->m.col_ptrs[cur_col + aux_col + 1]) && (cur_col < ln_cols))
        {
        ++cur_col;
        }

      if(cur_col != in_col)
        {
        break;
        }

      row_index = iterator_base::M->m.row_indices[cur_pos + skip_pos];
      }
    }

  // Now we need to find the next valid position in the subview.
  uword row_index;
  while(true)
    {
    const uword next_colptr = iterator_base::M->m.col_ptrs[cur_col + aux_col + 1];
    row_index = iterator_base::M->m.row_indices[cur_pos + skip_pos];

    // Are we at the last position?
    if(cur_col >= ln_cols)
      {
      cur_col = ln_cols;
      // Make sure we will be pointing at the last element in the parent matrix.
      skip_pos = iterator_base::M->m.n_nonzero - iterator_base::M->n_nonzero;
      break;
      }

    if(row_index < aux_row)
      {
      ++skip_pos;
      }
    else if(row_index < (aux_row + ln_rows))
      {
      break; // found
      }
    else
      {
      skip_pos += (next_colptr - (cur_pos + skip_pos));
      }

    // Did we move any columns?
    while(((skip_pos + cur_pos) >= iterator_base::M->m.col_ptrs[cur_col + aux_col + 1]) && (cur_col < ln_cols))
      {
      ++cur_col;
      }
    }

  // It is possible we have moved another column.
  while(((skip_pos + cur_pos) >= iterator_base::M->m.col_ptrs[cur_col + aux_col + 1]) && (cur_col < ln_cols))
    {
    ++cur_col;
    }

  iterator_base::internal_pos = cur_pos;
  iterator_base::internal_col = cur_col;
  }



template<typename eT>
inline
SpSubview<eT>::const_iterator::const_iterator(const SpSubview<eT>& in_M, uword in_row, uword in_col, uword in_pos, uword in_skip_pos)
  : iterator_base(in_M, in_col, in_pos)
  , skip_pos(in_skip_pos)
  {
  arma_ignore(in_row);
  
  // Nothing to do.
  }



template<typename eT>
inline
SpSubview<eT>::const_iterator::const_iterator(const const_iterator& other)
  : iterator_base(*other.M, other.internal_col, other.internal_pos)
  , skip_pos(other.skip_pos)
  {
  // Nothing to do.
  }



template<typename eT>
arma_inline
eT
SpSubview<eT>::const_iterator::operator*() const
  {
  return iterator_base::M->m.values[iterator_base::internal_pos + skip_pos];
  }



template<typename eT>
inline
arma_hot
typename SpSubview<eT>::const_iterator&
SpSubview<eT>::const_iterator::operator++()
  {
  const uword aux_col = iterator_base::M->aux_col1;
  const uword aux_row = iterator_base::M->aux_row1;
  const uword ln_rows = iterator_base::M->n_rows;
  const uword ln_cols = iterator_base::M->n_cols;

  uword cur_col   = iterator_base::internal_col;
  uword cur_pos   = iterator_base::internal_pos + 1;
  uword lskip_pos = skip_pos;
  uword row_index;

  while(true)
    {
    const uword next_colptr = iterator_base::M->m.col_ptrs[cur_col + aux_col + 1];
    row_index = iterator_base::M->m.row_indices[cur_pos + lskip_pos];

    // Did we move any columns?
    while((cur_col < ln_cols) && ((lskip_pos + cur_pos) >= iterator_base::M->m.col_ptrs[cur_col + aux_col + 1]))
      {
      ++cur_col;
      }

    // Are we at the last position?
    if(cur_col >= ln_cols)
      {
      cur_col = ln_cols;
      // Make sure we will be pointing at the last element in the parent matrix.
      lskip_pos = iterator_base::M->m.n_nonzero - iterator_base::M->n_nonzero;
      break;
      }

    if(row_index < aux_row)
      {
      ++lskip_pos;
      }
    else if(row_index < (aux_row + ln_rows))
      {
      break; // found
      }
    else
      {
      lskip_pos += (next_colptr - (cur_pos + lskip_pos));
      }
    }

  iterator_base::internal_pos = cur_pos;
  iterator_base::internal_col = cur_col;
  skip_pos                    = lskip_pos;

  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpSubview<eT>::const_iterator
SpSubview<eT>::const_iterator::operator++(int)
  {
  typename SpSubview<eT>::const_iterator tmp(*this);

  ++(*this);

  return tmp;
  }



template<typename eT>
inline
arma_hot
typename SpSubview<eT>::const_iterator&
SpSubview<eT>::const_iterator::operator--()
  {
  const uword aux_col = iterator_base::M->aux_col1;
  const uword aux_row = iterator_base::M->aux_row1;
  const uword ln_rows = iterator_base::M->n_rows;

  uword cur_col  = iterator_base::internal_col;
  uword cur_pos  = iterator_base::internal_pos - 1;

  // Special condition for end of iterator.
  if((skip_pos + cur_pos + 1) == iterator_base::M->m.n_nonzero)
    {
    // We are at the last element.  So we need to set skip_pos back to what it
    // would be if we didn't manually modify it back in operator++().
    skip_pos = iterator_base::M->m.col_ptrs[cur_col + aux_col] - iterator_base::internal_pos;
    }

  uword row_index;

  while(true)
    {
    const uword colptr = iterator_base::M->m.col_ptrs[cur_col + aux_col];
    row_index = iterator_base::M->m.row_indices[cur_pos + skip_pos];

    // Did we move back any columns?
    while((skip_pos + cur_pos) < iterator_base::M->m.col_ptrs[cur_col + aux_col])
      {
      --cur_col;
      }

    if(row_index < aux_row)
      {
      skip_pos -= (colptr - (cur_pos + skip_pos) + 1);
      }
    else if(row_index < (aux_row + ln_rows))
      {
      break; // found
      }
    else
      {
      --skip_pos;
      }
    }

  iterator_base::internal_pos = cur_pos;
  iterator_base::internal_col = cur_col;

  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpSubview<eT>::const_iterator
SpSubview<eT>::const_iterator::operator--(int)
  {
  typename SpSubview<eT>::const_iterator tmp(*this);

  --(*this);

  return tmp;
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_iterator::operator==(const const_iterator& rhs) const
  {
  return (rhs.row() == (*this).row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_iterator::operator!=(const const_iterator& rhs) const
  {
  return (rhs.row() != (*this).row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_iterator::operator==(const typename SpMat<eT>::const_iterator& rhs) const
  {
  return (rhs.row() == (*this).row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_iterator::operator!=(const typename SpMat<eT>::const_iterator& rhs) const
  {
  return (rhs.row() != (*this).row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_iterator::operator==(const const_row_iterator& rhs) const
  {
  return (rhs.row() == (*this).row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_iterator::operator!=(const const_row_iterator& rhs) const
  {
  return (rhs.row() != (*this).row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_iterator::operator==(const typename SpMat<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() == (*this).row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_iterator::operator!=(const typename SpMat<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() != (*this).row()) || (rhs.col() != iterator_base::internal_col);
  }



///////////////////////////////////////////////////////////////////////////////
// SpSubview<eT>::iterator implementation                                    //
///////////////////////////////////////////////////////////////////////////////

template<typename eT>
inline
arma_hot
SpValProxy< SpSubview<eT> >
SpSubview<eT>::iterator::operator*()
  {
  return SpValProxy< SpSubview<eT> >(
    const_iterator::row(),
    iterator_base::col(),
    access::rw(*iterator_base::M),
    &(access::rw(iterator_base::M->m.values[iterator_base::internal_pos + const_iterator::skip_pos])));
  }



template<typename eT>
inline
arma_hot
typename SpSubview<eT>::iterator&
SpSubview<eT>::iterator::operator++()
  {
  const_iterator::operator++();
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpSubview<eT>::iterator
SpSubview<eT>::iterator::operator++(int)
  {
  typename SpSubview<eT>::iterator tmp(*this);

  const_iterator::operator++();

  return tmp;
  }



template<typename eT>
inline
arma_hot
typename SpSubview<eT>::iterator&
SpSubview<eT>::iterator::operator--()
  {
  const_iterator::operator--();
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpSubview<eT>::iterator
SpSubview<eT>::iterator::operator--(int)
  {
  typename SpSubview<eT>::iterator tmp(*this);

  const_iterator::operator--();

  return tmp;
  }



///////////////////////////////////////////////////////////////////////////////
// SpSubview<eT>::const_row_iterator implementation                          //
///////////////////////////////////////////////////////////////////////////////

template<typename eT>
inline
SpSubview<eT>::const_row_iterator::const_row_iterator()
  : iterator_base()
  , internal_row(0)
  , actual_pos(0)
  {
  }



template<typename eT>
inline
SpSubview<eT>::const_row_iterator::const_row_iterator(const SpSubview<eT>& in_M, uword initial_pos)
  : iterator_base(in_M, 0, initial_pos)
  , internal_row(0)
  , actual_pos(0)
  {
  // Corner case for the end of a subview.
  if(initial_pos == in_M.n_nonzero)
    {
    iterator_base::internal_col = 0;
    internal_row = in_M.n_rows;
    return;
    }

  const uword aux_col = iterator_base::M->aux_col1;
  const uword aux_row = iterator_base::M->aux_row1;

  // We don't count zeros in our position count, so we have to find the nonzero
  // value corresponding to the given initial position, and we also have to skip
  // any nonzero elements that aren't a part of the subview.

  uword cur_pos = std::numeric_limits<uword>::max();
  uword cur_actual_pos = 0;

  // Since we don't know where the elements are in each row, we have to loop
  // across all columns looking for elements in row 0 and add to our sum, then
  // in row 1, and so forth, until we get to the desired position.
  for(uword row = 0; row < iterator_base::M->n_rows; ++row)
    {
    for(uword col = 0; col < iterator_base::M->n_cols; ++col)
      {
      // Find the first element with row greater than or equal to row + aux_row.
      const uword      col_offset = iterator_base::M->m.col_ptrs[col + aux_col    ];
      const uword next_col_offset = iterator_base::M->m.col_ptrs[col + aux_col + 1];

      const uword* start_ptr = &iterator_base::M->m.row_indices[     col_offset];
      const uword*   end_ptr = &iterator_base::M->m.row_indices[next_col_offset];

      if(start_ptr != end_ptr)
        {
        const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, row + aux_row);

        const uword offset = uword(pos_ptr - start_ptr);

        if(iterator_base::M->m.row_indices[col_offset + offset] == row + aux_row)
          {
          cur_actual_pos = col_offset + offset;

          // Increment position portably.
          if(cur_pos == std::numeric_limits<uword>::max())
            cur_pos = 0;
          else
            ++cur_pos;

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

  // This shouldn't happen.
  iterator_base::internal_pos = iterator_base::M->n_nonzero;
  iterator_base::internal_col = 0;
  internal_row = iterator_base::M->n_rows;
  actual_pos = iterator_base::M->n_nonzero;
  }



template<typename eT>
inline
SpSubview<eT>::const_row_iterator::const_row_iterator(const SpSubview<eT>& in_M, uword in_row, uword in_col)
  : iterator_base(in_M, in_col, 0)
  , internal_row(0)
  , actual_pos(0)
  {
  // Start our search in the given row.  We need to find two things:
  //
  //  1. The first nonzero element (iterating by rows) after (in_row, in_col).
  //  2. The number of nonzero elements (iterating by rows) that come before
  //     (in_row, in_col).
  //
  // We'll find these simultaneously, though we will have to loop over all
  // columns.

  const uword aux_col = iterator_base::M->aux_col1;
  const uword aux_row = iterator_base::M->aux_row1;

  // This will hold the total number of points in the subview with rows less
  // than in_row.
  uword cur_pos = 0;
  uword cur_min_row = iterator_base::M->n_rows;
  uword cur_min_col = 0;
  uword cur_actual_pos = 0;

  for(uword col = 0; col < iterator_base::M->n_cols; ++col)
    {
    // Find the first element with row greater than or equal to in_row.
    const uword      col_offset = iterator_base::M->m.col_ptrs[col + aux_col    ];
    const uword next_col_offset = iterator_base::M->m.col_ptrs[col + aux_col + 1];

    const uword* start_ptr = &iterator_base::M->m.row_indices[     col_offset];
    const uword*   end_ptr = &iterator_base::M->m.row_indices[next_col_offset];

    if(start_ptr != end_ptr)
      {
      // First let us find the first element that is in the subview.
      const uword* first_subview_ptr = std::lower_bound(start_ptr, end_ptr, aux_row);

      if(first_subview_ptr != end_ptr && (*first_subview_ptr) < aux_row + iterator_base::M->n_rows)
        {
        // There exists at least one element in the subview.
        const uword* pos_ptr = std::lower_bound(first_subview_ptr, end_ptr, aux_row + in_row);

        // This is the number of elements in the subview with row index less
        // than in_row.
        cur_pos += uword(pos_ptr - first_subview_ptr);

        if(pos_ptr != end_ptr && (*pos_ptr) < aux_row + iterator_base::M->n_rows)
          {
          // This is the row index of the first element in the column with row
          // index greater than or equal to in_row + aux_row.
          if((*pos_ptr) - aux_row < cur_min_row)
            {
            // If we are in the desired row but before the desired column, we
            // can't take this.
            if(col >= in_col)
              {
              cur_min_row = (*pos_ptr) - aux_row;
              cur_min_col = col;
              cur_actual_pos = col_offset + (pos_ptr - start_ptr);
              }
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



template<typename eT>
inline
SpSubview<eT>::const_row_iterator::const_row_iterator(const const_row_iterator& other)
  : iterator_base(*other.M, other.internal_col, other.internal_pos)
  , internal_row(other.internal_row)
  , actual_pos(other.actual_pos)
  {
  // Nothing to do.
  }



template<typename eT>
inline
arma_hot
typename SpSubview<eT>::const_row_iterator&
SpSubview<eT>::const_row_iterator::operator++()
  {
  // We just need to find the next nonzero element.
  ++iterator_base::internal_pos;

  // If we have exceeded the bounds, update accordingly.
  if(iterator_base::internal_pos >= iterator_base::M->n_nonzero)
    {
    internal_row = iterator_base::M->n_rows;
    iterator_base::internal_col = 0;
    actual_pos = iterator_base::M->n_nonzero;

    return *this;
    }

  const uword aux_col  = iterator_base::M->aux_col1;
  const uword aux_row  = iterator_base::M->aux_row1;
  const uword M_n_cols = iterator_base::M->n_cols;

  // Otherwise, we need to search.  We have to loop over all of the columns in
  // the subview.
  uword next_min_row = iterator_base::M->n_rows;
  uword next_min_col = 0;
  uword next_actual_pos = 0;

  for(uword col = iterator_base::internal_col + 1; col < M_n_cols; ++col)
    {
    // Find the first element with row greater than or equal to row.
    const uword      col_offset = iterator_base::M->m.col_ptrs[col + aux_col    ];
    const uword next_col_offset = iterator_base::M->m.col_ptrs[col + aux_col + 1];

    const uword* start_ptr = &iterator_base::M->m.row_indices[     col_offset];
    const uword* end_ptr   = &iterator_base::M->m.row_indices[next_col_offset];

    if(start_ptr != end_ptr)
      {
      // Find the first element in the column with row greater than or equal to
      // the current row.  Since this is a subview, it's possible that we may
      // find rows past the end of the subview.
      const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, internal_row + aux_row);

      if(pos_ptr != end_ptr)
        {
        // We found something; is the row index correct?
        if((*pos_ptr) == internal_row + aux_row && (*pos_ptr) < aux_row + iterator_base::M->n_rows)
          {
          // Exact match---so we are done.
          iterator_base::internal_col = col;
          actual_pos = col_offset + (pos_ptr - start_ptr);
          return *this;
          }
        else if((*pos_ptr) < next_min_row + aux_row && (*pos_ptr) < aux_row + iterator_base::M->n_rows)
          {
          // The first element in this column is in a subsequent row, but it's
          // the minimum row we've seen so far.
          next_min_row = (*pos_ptr) - aux_row;
          next_min_col = col;
          next_actual_pos = col_offset + (pos_ptr - start_ptr);
          }
        else if((*pos_ptr) == next_min_row + aux_row && col < next_min_col && (*pos_ptr) < aux_row + iterator_base::M->n_rows)
          {
          // The first element in this column is in a subsequent row that we
          // already have another elemnt for, but the column index is less so
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
    // Find the first element with row greater than or equal to row + 1.
    const uword      col_offset = iterator_base::M->m.col_ptrs[col + aux_col    ];
    const uword next_col_offset = iterator_base::M->m.col_ptrs[col + aux_col + 1];

    const uword* start_ptr = &iterator_base::M->m.row_indices[     col_offset];
    const uword* end_ptr   = &iterator_base::M->m.row_indices[next_col_offset];

    if(start_ptr != end_ptr)
      {
      const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, internal_row + aux_row + 1);

      if(pos_ptr != end_ptr)
        {
        // We found something in the column, but is the row index correct?
        if((*pos_ptr) == internal_row + aux_row + 1 && (*pos_ptr) < aux_row + iterator_base::M->n_rows)
          {
          // Exact match---so we are done.
          iterator_base::internal_col = col;
          internal_row++;
          actual_pos = col_offset + (pos_ptr - start_ptr);
          return *this;
          }
        else if((*pos_ptr) < next_min_row + aux_row && (*pos_ptr) < aux_row + iterator_base::M->n_rows)
          {
          // The first element in this column is in a subsequent row, but it's
          // the minimum row we've seen so far.
          next_min_row = (*pos_ptr) - aux_row;
          next_min_col = col;
          next_actual_pos = col_offset + (pos_ptr - start_ptr);
          }
        else if((*pos_ptr) == next_min_row + aux_row && col < next_min_col && (*pos_ptr) < aux_row + iterator_base::M->n_rows)
          {
          // We've found a better column.
          next_min_col = col;
          next_actual_pos = col_offset + (pos_ptr - start_ptr);
          }
        }
      }
    }

  iterator_base::internal_col = next_min_col;
  internal_row = next_min_row;
  actual_pos = next_actual_pos;

  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpSubview<eT>::const_row_iterator
SpSubview<eT>::const_row_iterator::operator++(int)
  {
  typename SpSubview<eT>::const_row_iterator tmp(*this);

  ++(*this);

  return tmp;
  }



template<typename eT>
inline
arma_hot
typename SpSubview<eT>::const_row_iterator&
SpSubview<eT>::const_row_iterator::operator--()
  {
  if(iterator_base::internal_pos == 0)
    {
    // We are already at the beginning.
    return *this;
    }

  iterator_base::internal_pos--;

  const uword aux_col = iterator_base::M->aux_col1;
  const uword aux_row = iterator_base::M->aux_row1;

  // We have to search backwards.
  uword max_row = 0;
  uword max_col = 0;
  uword next_actual_pos = 0;

  for(uword col = iterator_base::internal_col; col >= 1; --col)
    {
    // Find the first element with row greater than or equal to in_row + 1.
    const uword      col_offset = iterator_base::M->m.col_ptrs[col + aux_col - 1];
    const uword next_col_offset = iterator_base::M->m.col_ptrs[col + aux_col   ];

    const uword* start_ptr = &iterator_base::M->m.row_indices[     col_offset];
    const uword* end_ptr   = &iterator_base::M->m.row_indices[next_col_offset];

    if(start_ptr != end_ptr)
      {
      // There are elements in this column.
      const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, internal_row + aux_row + 1);

      if(pos_ptr != start_ptr)
        {
        if(*(pos_ptr - 1) > max_row + aux_row)
          {
          // There are elements in this column with row index < internal_row.
          max_row = *(pos_ptr - 1) - aux_row;
          max_col = col - 1;
          next_actual_pos = col_offset + (pos_ptr - 1 - start_ptr);
          }
        else if(*(pos_ptr - 1) == max_row + aux_row && (col - 1) >= max_col)
          {
          max_col = col - 1;
          next_actual_pos = col_offset + (pos_ptr - 1 - start_ptr);
          }
        }
      }
    }

  for(uword col = iterator_base::M->n_cols - 1; col >= iterator_base::internal_col; --col)
    {
    // Find the first element with row greater than or equal to row + 1.
    const uword      col_offset = iterator_base::M->m.col_ptrs[col + aux_col    ];
    const uword next_col_offset = iterator_base::M->m.col_ptrs[col + aux_col + 1];

    const uword* start_ptr = &iterator_base::M->m.row_indices[     col_offset];
    const uword*   end_ptr = &iterator_base::M->m.row_indices[next_col_offset];

    if(start_ptr != end_ptr)
      {
      // There are elements in this column.
      const uword* pos_ptr = std::lower_bound(start_ptr, end_ptr, internal_row + aux_row);

      if(pos_ptr != start_ptr)
        {
        // There are elements in this column with row index < internal_row.
        if(*(pos_ptr - 1) > max_row + aux_row)
          {
          max_row = *(pos_ptr - 1) - aux_row;
          max_col = col;
          next_actual_pos = col_offset + (pos_ptr - 1 - start_ptr);
          }
        else if(*(pos_ptr - 1) == max_row + aux_row && col >= max_col)
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



template<typename eT>
inline
arma_warn_unused
typename SpSubview<eT>::const_row_iterator
SpSubview<eT>::const_row_iterator::operator--(int)
  {
  typename SpSubview<eT>::const_row_iterator tmp(*this);

  --(*this);

  return tmp;
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_row_iterator::operator==(const const_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_row_iterator::operator!=(const const_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_row_iterator::operator==(const typename SpMat<eT>::const_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_row_iterator::operator!=(const typename SpMat<eT>::const_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_row_iterator::operator==(const const_row_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_row_iterator::operator!=(const const_row_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_row_iterator::operator==(const typename SpMat<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() == row()) && (rhs.col() == iterator_base::internal_col);
  }



template<typename eT>
inline
arma_hot
bool
SpSubview<eT>::const_row_iterator::operator!=(const typename SpMat<eT>::const_row_iterator& rhs) const
  {
  return (rhs.row() != row()) || (rhs.col() != iterator_base::internal_col);
  }



///////////////////////////////////////////////////////////////////////////////
// SpSubview<eT>::row_iterator implementation                                //
///////////////////////////////////////////////////////////////////////////////

template<typename eT>
inline
arma_hot
SpValProxy< SpSubview<eT> >
SpSubview<eT>::row_iterator::operator*()
  {
  return SpValProxy< SpSubview<eT> >(
    const_row_iterator::internal_row,
    iterator_base::internal_col,
    access::rw(*iterator_base::M),
    &access::rw(iterator_base::M->m.values[const_row_iterator::actual_pos]));
  }



template<typename eT>
inline
arma_hot
typename SpSubview<eT>::row_iterator&
SpSubview<eT>::row_iterator::operator++()
  {
  const_row_iterator::operator++();
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpSubview<eT>::row_iterator
SpSubview<eT>::row_iterator::operator++(int)
  {
  typename SpSubview<eT>::row_iterator tmp(*this);

  ++(*this);

  return tmp;
  }



template<typename eT>
inline
arma_hot
typename SpSubview<eT>::row_iterator&
SpSubview<eT>::row_iterator::operator--()
  {
  const_row_iterator::operator--();
  return *this;
  }



template<typename eT>
inline
arma_warn_unused
typename SpSubview<eT>::row_iterator
SpSubview<eT>::row_iterator::operator--(int)
  {
  typename SpSubview<eT>::row_iterator tmp(*this);

  --(*this);

  return tmp;
  }


//! @}
