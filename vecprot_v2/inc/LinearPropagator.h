//===--- LinearPropagator.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file LinearPropagator.h
 * @brief Implementation of the linear propagation as a selector.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_LINEAR_PROPAGATOR
#define GEANT_LINEAR_PROPAGATOR

#include "Selector.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Selector grouping charged tracks and performing field propagation.
 */
 
class LinearPropagator : public Selector {

protected:  

private:
  LinearPropagator(const LinearPropagator &) = delete;
  LinearPropagator &operator=(const LinearPropagator &) = delete;
    
public:
  /** @brief Default constructor */
  VECCORE_ATT_HOST_DEVICE
  LinearPropagator() : Selector() {}

  /** 
   * @brief Default constructor
   * @param threshold Basketizing threshold
   * @param propagator Propagator working with this selector
   */
  VECCORE_ATT_HOST_DEVICE
  LinearPropagator(int threshold, GeantPropagator *propagator);

  /** @brief Geometry filter destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~LinearPropagator();

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
