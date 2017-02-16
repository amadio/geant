//===--- FieldPropagator.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file FieldPropagator.h
 * @brief Implementation of the field propagation as a selector.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_FIELD_PROPAGATOR
#define GEANT_FIELD_PROPAGATOR

#include "Selector.h"
#include "GeantTaskData.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Selector grouping charged tracks and performing field propagation.
 */
 
class FieldPropagator : public Selector {

protected:  

private:
  FieldPropagator(const FieldPropagator &) = delete;
  FieldPropagator &operator=(const FieldPropagator &) = delete;
  
  
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
  FieldPropagator() : Selector() {}

  /** 
   * @brief Default constructor
   * @param threshold Basketizing threshold
   * @param propagator Propagator working with this selector
   */
  VECCORE_ATT_HOST_DEVICE
  FieldPropagator(int threshold, GeantPropagator *propagator);

  /** @brief Geometry filter destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~FieldPropagator();

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
