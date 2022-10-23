/* File : cppfunctions.h */

/*
 * Copyright (c) 2020-2022, Dany Cajas
 * All rights reserved.
 * This work is licensed under BSD 3-Clause "New" or "Revised" License.
 * License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
 */

#ifndef _CPPFUNCTIONS_H_
#define _CPPFUNCTIONS_H_
#endif 

#include <iostream>
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>

arma::Mat<double> duplication_matrix(const int &n);

arma::Mat<double> duplication_elimination_matrix(const int &n);

arma::Mat<double> summation_matrix(const int &n);

arma::Mat<double> kurtosis_matrix(const arma::Mat<double> &Y);

arma::Mat<double> semi_kurtosis_matrix(const arma::Mat<double> &Y);