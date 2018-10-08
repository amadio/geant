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

#include "Geant/SimulationStage.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeomQueryStage : public SimulationStage {

protected:
public:
  /** @brief Clone the stage and copy the existing handlers **/
  VECCORE_ATT_HOST_DEVICE
  virtual SimulationStage *Clone() const;
  
  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Handler *Select(Track *track, TaskData *td);

public:
  /** @brief Dummy GeomQueryStage constructor */
  VECCORE_ATT_HOST_DEVICE
  GeomQueryStage() {}

  /** @brief Standard GeomQueryStage constructor */
  VECCORE_ATT_HOST_DEVICE
  GeomQueryStage(Propagator *prop);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~GeomQueryStage() {}

  VECCORE_ATT_HOST_DEVICE
  GeomQueryStage(const GeomQueryStage &);

  VECCORE_ATT_HOST_DEVICE
  GeomQueryStage &operator=(const GeomQueryStage &);

  /** @brief Simulation stage name */
  VECCORE_ATT_HOST_DEVICE
  virtual const char *GetName() const { return "GeometryQuery"; }

  /** @brief Activating basketizing */
  VECCORE_ATT_HOST_DEVICE
  void ActivateBasketizing(bool flag = true);
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
