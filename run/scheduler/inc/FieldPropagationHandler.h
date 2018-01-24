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

#include "WorkspaceForFieldPropagation.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Handler grouping charged tracks and performing field propagation.
 */

class FieldPropagationHandler : public Handler
{

public:

  /** @brief Scalar DoIt interface */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(GeantTrack *track, Basket& output, GeantTaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(Basket &input, Basket& output, GeantTaskData *td);

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

  /** @brief Field Propagation filter destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~FieldPropagationHandler();

  /** @brief maximum acceptable deflection from curved trajectory */
  static const double gEpsDeflection; // = 1.E-2 * geant::cm;

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

  /** @brief Curvature for general field    */
  VECCORE_ATT_HOST_DEVICE
  double Curvature(const GeantTrack &track ) const;

  /** @brief Function that returns safe length */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double SafeLength(const GeantTrack &track, double eps = 1.E-4);

  /** @brief Function that return Field Propagator, i.e. the holder of (RK) Integration Driver */
  GEANT_FORCE_INLINE
  GUFieldPropagator * GetFieldPropagator(GeantTaskData *td);

  // - Book keeping methods for task data

  /** @brief Connect with thread's FieldPropagator & create working buffers */
  VECCORE_ATT_HOST_DEVICE
  GUFieldPropagator * Initialize(GeantTaskData * td);

  /** @brief Cleanup the thread working buffers */
  VECCORE_ATT_HOST_DEVICE
  void Cleanup(GeantTaskData * td);

  /** @brief Clear the old buffers and create new working buffers */
  VECCORE_ATT_HOST_DEVICE
  void PrepareBuffers( size_t nTracks, GeantTaskData *td );
};

// ---------------------------------------------------------------------------------
// Inline implementation ----

VECCORE_ATT_HOST_DEVICE
GEANT_FORCE_INLINE
double
FieldPropagationHandler::
SafeLength(const GeantTrack &track, double eps)
{
   // Returns the propagation length in field such that the propagated point is
   // shifted less than eps with respect to the linear propagation.
   // OLD: return 2. * sqrt(eps / track.Curvature(Bz));
   double c = Curvature(track); //, td);
   double val= 0.0;
   // if (c < 1.E-10) { val= 1.E50; } else
   val = 2. * sqrt(eps / c);
   return val;
}

//______________________________________________________________________________
// VECCORE_ATT_HOST_DEVICE -- not yet
GUFieldPropagator *
FieldPropagationHandler::GetFieldPropagator( GeantTaskData *td )
{
   GUFieldPropagator *fieldPropagator = nullptr;
   bool useRungeKutta = td->fPropagator->fConfig->fUseRungeKutta;

   // static GUFieldPropagatorPool* fieldPropPool= GUFieldPropagatorPool::Instance();
   // if( useRungeKutta && fieldPropPool ){
   //    fieldPropagator = fieldPropPool->GetPropagator(td->fTid);
   //    assert( fieldPropagator );  // To assert, it must be initialised !!
   // }
   if( useRungeKutta ){
      fieldPropagator = td->fFieldPropagator;
      // assert( fieldPropagator );  // To assert, it must be initialised !!
   }
   // GUFieldPropagator *fieldPropagator = useRungeKutta ? td->fFieldPropagator : nullptr;
   return fieldPropagator;
}

//______________________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
inline
void FieldPropagationHandler::PrepareBuffers( size_t nTracks, GeantTaskData *td )
{
   auto wsp = td->fSpace4FieldProp;
   assert(wsp);
   if( nTracks > wsp->capacity() ){
      std::cout << "Calling ClearAndResizeBuffers on task/thread " << td->fTid
                << " with tracks = " << nTracks
                << " . Note: capacity = " << wsp->capacity() << std::endl;
      wsp->ClearAndResize( nTracks );
   } else {
      wsp->Resize(0); // Erase the entries, ready for new content!!
   }
}

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
