//===--- globals.h - Geant-V ------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file globals.h
 * @brief Definition of global variables for Propagator  
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef PROPAGATOR_GLOBALS
#define PROPAGATOR_GLOBALS

#include "GeantPropagator.h"

/** @brief Propagator class */
extern GeantPropagator *gPropagator; /** Propagator class */

#ifdef GEANT_NVCC
__constant__ double gPropagator_fBmag;
#endif

#endif
