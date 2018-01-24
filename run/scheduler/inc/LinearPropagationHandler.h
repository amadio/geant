//===--- LinearPropagationHandler.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file LinearPropagationHandler.h
 * @brief Implementation of the linear propagation as a handler.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_LINEAR_PROPAGATION_HANDLER
#define GEANT_LINEAR_PROPAGATION_HANDLER

#include "Handler.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Handler grouping charged tracks and performing field propagation.
 */
 
class LinearPropagationHandler : public Handler {

protected:  

private:
  LinearPropagationHandler(const LinearPropagationHandler &) = delete;
  LinearPropagationHandler &operator=(const LinearPropagationHandler &) = delete;

  VECCORE_ATT_HOST_DEVICE
  bool IsSameLocation(GeantTrack &track, GeantTaskData *td);
public:
  /** @brief Default constructor */
  VECCORE_ATT_HOST_DEVICE
  LinearPropagationHandler() : Handler() {}

  /** 
   * @brief Default constructor
   * @param threshold Basketizing threshold
   * @param propagator Propagator working with this handler
   */
  VECCORE_ATT_HOST_DEVICE
  LinearPropagationHandler(int threshold, GeantPropagator *propagator);

  /** @brief Geometry filter destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~LinearPropagationHandler();

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
