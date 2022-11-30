/*  carma/exceptions.h: Exceptions used by Carma
 *  Copyright (c) 2021 Ralph Urlus <rurlus.dev@gmail.com>
 *  All rights reserved. Use of this source code is governed by a
 *  Apache-2.0 license that can be found in the LICENSE file.
 */
#ifndef INCLUDE_CARMA_BITS_EXCEPTIONS_H_
#define INCLUDE_CARMA_BITS_EXCEPTIONS_H_
#include <stdexcept>

namespace carma {

class ConversionError : public std::runtime_error {
 public:
    explicit ConversionError(const char* what) : std::runtime_error(what) {}
};

}  // namespace carma
#endif  // INCLUDE_CARMA_BITS_EXCEPTIONS_H_
