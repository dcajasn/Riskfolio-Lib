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


//! \addtogroup constants_old
//! @{


// DO NOT USE IN NEW CODE !!!
// the Math and Phy classes are kept for compatibility with old code;
// for new code, use the Datum class instead
// eg. instead of math::pi(), use datum::pi

template<typename eT>
class Math
  {
  public:
  
  // the long lengths of the constants are for future support of "long double"
  // and any smart compiler that does high-precision computation at compile-time
  
  //! ratio of any circle's circumference to its diameter
  arma_deprecated static eT pi()        { return eT(Datum<eT>::pi); }  // use datum::pi instead
  
  //! base of the natural logarithm
  arma_deprecated static eT e()         { return eT(Datum<eT>::e); }  // use datum::e instead
  
  //! Euler's constant, aka Euler-Mascheroni constant
  arma_deprecated static eT euler()     { return eT(Datum<eT>::euler); }  // use datum::euler instead
  
  //! golden ratio
  arma_deprecated static eT gratio()    { return eT(Datum<eT>::gratio); }  // use datum::gratio instead
  
  //! square root of 2
  arma_deprecated static eT sqrt2()     { return eT(Datum<eT>::sqrt2); }  // use datum::sqrt2 instead
  
  //! the difference between 1 and the least value greater than 1 that is representable
  arma_deprecated static eT eps()       { return eT(Datum<eT>::eps); }  // use datum::eps instead
  
  //! log of the minimum representable value
  arma_deprecated static eT log_min()   { return eT(Datum<eT>::log_min); }  // use datum::log_min instead
    
  //! log of the maximum representable value
  arma_deprecated static eT log_max()   { return eT(Datum<eT>::log_max); }  // use datum::log_max instead
  
  //! "not a number"
  arma_deprecated static eT nan()       { return eT(Datum<eT>::nan); }  // use datum::nan instead
  
  //! infinity 
  arma_deprecated static eT inf()       { return eT(Datum<eT>::inf); }  // use datum::inf instead
  };



//! Physical constants taken from NIST 2010 CODATA values, and some from WolframAlpha (values provided as of 2009-06-23)
//! http://physics.nist.gov/cuu/Constants
//! http://www.wolframalpha.com
//! See also http://en.wikipedia.org/wiki/Physical_constant
template<typename eT>
class Phy
  {
  public:
  
  //! atomic mass constant (in kg)
  arma_deprecated static eT m_u()       { return eT(Datum<eT>::m_u); }
  
  //! Avogadro constant
  arma_deprecated static eT N_A()       { return eT(Datum<eT>::N_A); }
  
  //! Boltzmann constant (in joules per kelvin)
  arma_deprecated static eT k()         { return eT(Datum<eT>::k); }
  
  //! Boltzmann constant (in eV/K)
  arma_deprecated static eT k_evk()     { return eT(Datum<eT>::k_evk); }
  
  //! Bohr radius (in meters)
  arma_deprecated static eT a_0()       { return eT(Datum<eT>::a_0); }
  
  //! Bohr magneton
  arma_deprecated static eT mu_B()      { return eT(Datum<eT>::mu_B); }
  
  //! characteristic impedance of vacuum (in ohms)
  arma_deprecated static eT Z_0()       { return eT(Datum<eT>::Z_0); }
  
  //! conductance quantum (in siemens)
  arma_deprecated static eT G_0()       { return eT(Datum<eT>::G_0); }
  
  //! Coulomb's constant (in meters per farad)
  arma_deprecated static eT k_e()       { return eT(Datum<eT>::k_e); }
  
  //! electric constant (in farads per meter)
  arma_deprecated static eT eps_0()     { return eT(Datum<eT>::eps_0); }
  
  //! electron mass (in kg)
  arma_deprecated static eT m_e()       { return eT(Datum<eT>::m_e); }
  
  //! electron volt (in joules)
  arma_deprecated static eT eV()        { return eT(Datum<eT>::eV); }
  
  //! elementary charge (in coulombs)
  arma_deprecated static eT e()         { return eT(Datum<eT>::ec); }
  
  //! Faraday constant (in coulombs)
  arma_deprecated static eT F()         { return eT(Datum<eT>::F); }
  
  //! fine-structure constant
  arma_deprecated static eT alpha()     { return eT(Datum<eT>::alpha); }
  
  //! inverse fine-structure constant
  arma_deprecated static eT alpha_inv() { return eT(Datum<eT>::alpha_inv); }
  
  //! Josephson constant
  arma_deprecated static eT K_J()       { return eT(Datum<eT>::K_J); }
  
  //! magnetic constant (in henries per meter)
  arma_deprecated static eT mu_0()      { return eT(Datum<eT>::mu_0); }
  
  //! magnetic flux quantum (in webers)
  arma_deprecated static eT phi_0()     { return eT(Datum<eT>::phi_0); }
  
  //! molar gas constant (in joules per mole kelvin)
  arma_deprecated static eT R()         { return eT(Datum<eT>::R); }
  
  //! Newtonian constant of gravitation (in newton square meters per kilogram squared)
  arma_deprecated static eT G()         { return eT(Datum<eT>::G); }
  
  //! Planck constant (in joule seconds)
  arma_deprecated static eT h()         { return eT(Datum<eT>::h); }
  
  //! Planck constant over 2 pi, aka reduced Planck constant (in joule seconds)
  arma_deprecated static eT h_bar()     { return eT(Datum<eT>::h_bar); }
  
  //! proton mass (in kg)
  arma_deprecated static eT m_p()       { return eT(Datum<eT>::m_p); }
  
  //! Rydberg constant (in reciprocal meters)
  arma_deprecated static eT R_inf()     { return eT(Datum<eT>::R_inf); }
  
  //! speed of light in vacuum (in meters per second)
  arma_deprecated static eT c_0()       { return eT(Datum<eT>::c_0); }
  
  //! Stefan-Boltzmann constant
  arma_deprecated static eT sigma()     { return eT(Datum<eT>::sigma); }
  
  //! von Klitzing constant (in ohms)
  arma_deprecated static eT R_k()       { return eT(Datum<eT>::R_k); }
  
  //! Wien wavelength displacement law constant
  arma_deprecated static eT b()         { return eT(Datum<eT>::b); }
  };



typedef Math<float>  fmath;
typedef Math<double> math;

typedef Phy<float>   fphy;
typedef Phy<double>  phy;



//! @}
