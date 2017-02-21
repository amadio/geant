//===--- GeomQueryStage.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeomQueryStage.h
 * @brief The geometry query part as simulation stage.
 * @details This simulation stage deals with geometry computation of the
 *          transport length
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_GEOM_QUERY_STAGE
#define GEANT_GEOM_QUERY_STAGE

#include "SimulationStage.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeomQueryStage : public SimulationStage {

private:
  GeomQueryStage(const GeomQueryStage &) = delete;
  GeomQueryStage &operator=(const GeomQueryStage &) = delete;

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
  /** @brief Dummy GeomQueryStage constructor */
  VECCORE_ATT_HOST_DEVICE
  GeomQueryStage() {}

  /** @brief Standard GeomQueryStage constructor */
  VECCORE_ATT_HOST_DEVICE
  GeomQueryStage(GeantPropagator *prop);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  ~GeomQueryStage();

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
