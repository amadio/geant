//===--- XSecSamplingStage.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file XSecSamplingStage.h
 * @brief This stage deals with proposing the physics step by sampling cross sections
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_XSEC_SAMPLING_STAGE
#define GEANT_XSEC_SAMPLING_STAGE

#include "SimulationStage.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class XSecSamplingStage : public SimulationStage {

private:
  XSecSamplingStage(const XSecSamplingStage &) = delete;
  XSecSamplingStage &operator=(const XSecSamplingStage &) = delete;

protected:

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateHandlers();

public:

  /** @brief Interface to select the handler matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Handler *Select(GeantTrack *track);

public:
  /** @brief Dummy XSecSamplingStage constructor */
  VECCORE_ATT_HOST_DEVICE
  XSecSamplingStage() {}

  /** @brief Standard XSecSamplingStage constructor */
  VECCORE_ATT_HOST_DEVICE
  XSecSamplingStage(GeantPropagator *prop);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  ~XSecSamplingStage();

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
