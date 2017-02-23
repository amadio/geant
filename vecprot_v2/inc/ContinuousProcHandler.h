//===--- ContinuousProcHandler.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file ContinuousProcHandler.h
 * @brief Implementation of the continuous process handler.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_CONTINUOUS_PROC_HANDLER
#define GEANT_CONTINUOUS_PROC_HANDLER

#include "Handler.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Handler grouping charged tracks and performing field propagation.
 */
 
class ContinuousProcHandler : public Handler {

protected:  

private:
  ContinuousProcHandler(const ContinuousProcHandler &) = delete;
  ContinuousProcHandler &operator=(const ContinuousProcHandler &) = delete;
    
public:
  /** @brief Default constructor */
  VECCORE_ATT_HOST_DEVICE
  ContinuousProcHandler() : Handler() {}

  /** 
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  VECCORE_ATT_HOST_DEVICE
  ContinuousProcHandler(int threshold, GeantPropagator *propagator);

  /** @brief Geometry filter destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~ContinuousProcHandler();

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
