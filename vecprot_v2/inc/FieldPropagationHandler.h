//===--- FieldPropagationHandler.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file FieldPropagationHandler.h
 * @brief Implementation of the field propagation as a handler.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_FIELD_PROPAGATION_HANDLER
#define GEANT_FIELD_PROPAGATION_HANDLER

#include "Handler.h"
#include "GeantTaskData.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Handler grouping charged tracks and performing field propagation.
 */
 
class FieldPropagationHandler : public Handler {

protected:
  VECCORE_ATT_HOST_DEVICE
  bool IsSameLocation(GeantTrack &track, GeantTaskData *td);

private:
  FieldPropagationHandler(const FieldPropagationHandler &) = delete;
  FieldPropagationHandler &operator=(const FieldPropagationHandler &) = delete;
  
  
  /** @brief Scalar implementation for magnetic field propagation */
  VECCORE_ATT_HOST_DEVICE
  void PropagateInVolume(GeantTrack &track, double crtstep, GeantTaskData * td);

  /** @brief Vector implementation for magnetic field propagation */
  VECCORE_ATT_HOST_DEVICE
  void PropagateInVolume(TrackVec_t &tracks, const double *crtstep, GeantTaskData *td);

  /** @brief Function that returns safe length */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double SafeLength(const GeantTrack &track, double Bz, double eps = 1.E-4) {
    // Returns the propagation length in field such that the propagated point is
    // shifted less than eps with respect to the linear propagation.
    return 2. * sqrt(eps / track.Curvature(Bz));
  }
  
public:
  /** @brief Default constructor */
  VECCORE_ATT_HOST_DEVICE
  FieldPropagationHandler() : Handler() {}

  /** 
   * @brief Default constructor
   * @param threshold Basketizing threshold
   * @param propagator Propagator working with this handler
   */
  VECCORE_ATT_HOST_DEVICE
  FieldPropagationHandler(int threshold, GeantPropagator *propagator);

  /** @brief Geometry filter destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~FieldPropagationHandler();

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
