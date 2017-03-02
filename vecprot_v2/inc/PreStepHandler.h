//===--- PreStepHandler.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file PreStepHandler.h
 * @brief Implementation of the stepping actions handler.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_PRE_STEP_HANDLER
#define GEANT_PRE_STEP_HANDLER

#include "Handler.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Handler invoked once per step per particle. Calling the user actions.
 */
 
class PreStepHandler : public Handler {

protected:  

private:
  PreStepHandler(const PreStepHandler &) = delete;
  PreStepHandler &operator=(const PreStepHandler &) = delete;
    
public:
  /** @brief Default constructor */
  VECCORE_ATT_HOST_DEVICE
  PreStepHandler() : Handler() {}

  /** 
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  VECCORE_ATT_HOST_DEVICE
  PreStepHandler(int threshold, GeantPropagator *propagator);

  /** @brief Geometry filter destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~PreStepHandler();

  /** @brief Scalar DoIt interface */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(GeantTrack *track, Basket& output, GeantTaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(Basket &input, Basket& output, GeantTaskData *td);

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
