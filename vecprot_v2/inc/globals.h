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

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantPropagator;

/** @brief Propagator class */
extern GeantPropagator *gPropagator; /** Propagator class */
#ifdef VECCORE_CUDA
__constant__ double gPropagator_fBmag;
#endif

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
