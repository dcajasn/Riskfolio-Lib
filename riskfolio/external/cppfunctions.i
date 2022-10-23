/* File : cppfunctions.i */

/*
 * Copyright (c) 2020-2022, Dany Cajas
 * All rights reserved.
 * This work is licensed under BSD 3-Clause "New" or "Revised" License.
 * License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
 */

%module cppfunctions
%{
/* include C++ header files necessary to compile the interface */
#define SWIG_FILE_WITH_INIT
#include "cppfunctions.h"
%}

/* Now include arma_numpy typemaps */
%include "arma_numpy.i"

/* Parse the header file to generate wrappers */
%include "cppfunctions.h"

