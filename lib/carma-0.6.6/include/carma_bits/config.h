#ifndef INCLUDE_CARMA_BITS_CONFIG_H_
#define INCLUDE_CARMA_BITS_CONFIG_H_
// CARMA Configuration
//
// This is the configuration header for carma, set the option to ON in carma_config.cmake to enable the setting.
// Please see the documentation for details: https://carma.readthedocs.io/en/latest/configuration.html

//-- ENABLE_CARMA_EXTRA_DEBUG --
#if !defined(CARMA_EXTRA_DEBUG)
#cmakedefine ENABLE_CARMA_EXTRA_DEBUG
#cmakedefine CARMA_EXTRA_DEBUG "@CARMA_EXTRA_DEBUG@"
#endif

// -- ENABLE_CARMA_SOFT_STEAL --
// When stealing the data of an array replace it with
// an array containing a single NaN
// This is a safer option compared to HARD_STEAL
// The default approach when staling is to only set the OWNDATA flag to False
#if !defined(CARMA_SOFT_STEAL)
#cmakedefine ENABLE_CARMA_SOFT_STEAL
#cmakedefine CARMA_SOFT_STEAL "@CARMA_SOFT_STEAL@"
#endif

// -- ENABLE_CARMA_HARD_STEAL --
// When stealing the data of an array set nullptr
// NOTE this will cause a segfault when accessing the original array's data
// The default approach when staling is to only set the OWNDATA flag to False
#if !defined(CARMA_HARD_STEAL)
#if !defined(CARMA_SOFT_STEAL)
#cmakedefine ENABLE_CARMA_HARD_STEAL
#cmakedefine CARMA_HARD_STEAL "@CARMA_HARD_STEAL@"
#endif
#endif

// -- REQUIRE_OWNDATA --
// Do NOT copy arrays if the data is not owned by Numpy, default behaviour
// is to copy when OWNDATA is False
#if !defined(CARMA_DONT_REQUIRE_OWNDATA)
#cmakedefine ENABLE_CARMA_DONT_REQUIRE_OWNDATA
#cmakedefine CARMA_DONT_REQUIRE_OWNDATA "@CARMA_DONT_REQUIRE_OWNDATA@"
#endif

// -- REQUIRE_F_CONTIGUOUS --
// Do NOT copy c-style arrays, default behaviour is to copy c-style arrays

#if !defined(CARMA_DONT_REQUIRE_F_CONTIGUOUS)
#cmakedefine ENABLE_CARMA_DONT_REQUIRE_F_CONTIGUOUS
#cmakedefine CARMA_DONT_REQUIRE_F_CONTIGUOUS "@CARMA_DONT_REQUIRE_F_CONTIGUOUS@"
#endif

#endif  // INCLUDE_CARMA_BITS_CONFIG_H_
