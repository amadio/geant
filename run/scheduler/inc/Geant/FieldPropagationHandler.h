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

#include "Geant/Handler.h"
#include "Geant/TaskData.h"

#include "Geant/WorkspaceForFieldPropagation.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Handler grouping charged tracks and performing field propagation.
 */

class FieldPropagationHandler : public Handler {

public:
  /** @brief Scalar DoIt interface */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(Track *track, Basket &output, TaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(Basket &input, Basket &output, TaskData *td);

  /** @brief Default constructor */
  VECCORE_ATT_HOST_DEVICE
  FieldPropagationHandler() : Handler(), fEpsTol(3.0e-5) {}

  /**
   * @brief Default constructor
   * @param threshold Basketizing threshold
   * @param propagator Propagator working with this handler
   */
  VECCORE_ATT_HOST_DEVICE
  FieldPropagationHandler(int threshold, Propagator *propagator, double relativeError = 3.0e-5);

  /** @brief Field Propagation filter destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~FieldPropagationHandler();

  /** @brief maximum acceptable deflection from curved trajectory */
  static const double gEpsDeflection; // = 1.E-2 * geant::units::cm;

protected:
  VECCORE_ATT_HOST_DEVICE
  bool IsSameLocation(Track &track, TaskData *td);

  VECCORE_ATT_HOST_DEVICE
  void CheckTrack(Track &track, const char *msg, double epsilon = 1.0e-5) const;

  VECCORE_ATT_HOST_DEVICE
  static void InitializeStats();

  VECCORE_ATT_HOST_DEVICE
  void PrintStats() const;

  VECCORE_ATT_HOST_DEVICE
  void CheckVsScalar(const vecgeom::Vector3D<double> &startPosition, const vecgeom::Vector3D<double> &startDirection,
                     double charge, double startMomentum, double stepSize, const vecgeom::Vector3D<double> &endPosition,
                     const vecgeom::Vector3D<double> &endDirection, double curvature,
                     int index, // For reporting
                     TaskData *td, bool checkVsHelix = true);

private:
  FieldPropagationHandler(const FieldPropagationHandler &) = delete;
  FieldPropagationHandler &operator=(const FieldPropagationHandler &) = delete;

  /** @brief Relative error acceptable in Runge-Kutta integration */
  const double fEpsTol;

  /** @brief Scalar implementation for magnetic field propagation */
  VECCORE_ATT_HOST_DEVICE
  void PropagateInVolume(Track &track, double crtstep, TaskData *td);

  /** @brief Vector implementation for magnetic field propagation */
  VECCORE_ATT_HOST_DEVICE
  void PropagateInVolume(TrackVec_t &tracks, const double *crtstep, TaskData *td);

  /** @brief Curvature for general field    */
  VECCORE_ATT_HOST_DEVICE
  double Curvature(const Track &track) const;

  /** @brief Function that returns safe length */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double SafeLength(const Track &track, double eps = 1.E-4) const;

  /** @brief Function that return Field Propagator, i.e. the holder of (RK) Integration Driver */
  GEANT_FORCE_INLINE
  GUFieldPropagator *GetFieldPropagator(TaskData *td);

  // - Book keeping methods for task data

  /** @brief Cleanup the thread working buffers */
  VECCORE_ATT_HOST_DEVICE
  void Cleanup(TaskData *td);

  /** @brief Clear the old buffers and create new working buffers */
  VECCORE_ATT_HOST_DEVICE
  void PrepareBuffers(size_t nTracks, TaskData *td);
};

// ---------------------------------------------------------------------------------
// Inline implementation ----

VECCORE_ATT_HOST_DEVICE
GEANT_FORCE_INLINE
double FieldPropagationHandler::SafeLength(const Track &track, double eps) const
{
  // Returns the propagation length in field such that the propagated point is
  // shifted less than eps with respect to the linear propagation.
  // OLD: return 2. * sqrt(eps / track.Curvature(Bz));
  double c   = Curvature(track); //, td);
  double val = 0.0;
  // if (c < 1.E-10) { val= 1.E50; } else
  val = 2. * sqrt(eps / c);
  return val;
}

//______________________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
inline void FieldPropagationHandler::PrepareBuffers(size_t nTracks, TaskData *td)
{
  auto wsp = td->fSpace4FieldProp;
  assert(wsp);
  if (nTracks > wsp->capacity()) {
    //std::cout << "Calling ClearAndResizeBuffers on task/thread " << td->fTid << " with tracks = " << nTracks
    //          << " . Note: capacity = " << wsp->capacity() << std::endl;
    wsp->ClearAndResize(nTracks);
  } else {
    wsp->Resize(0); // Erase the entries, ready for new content!!
  }
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant

#endif
