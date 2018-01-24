//===--- GeomQueryHandler.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeomQueryHandler.h
 * @brief Implementation of a geometry query for the transport length to next boundary.
 * @details The volume handler performs computation of the allowed geometry step
 *  using both scalar and vector interfaces.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_GEOM_QUERY_HANDLER
#define GEANT_GEOM_QUERY_HANDLER

#include "Handler.h"

#ifdef USE_ROOT
#include "TGeoExtension.h"
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Handler grouping tracks by logical volume and performing geometry actions.
 */

#if !defined(USE_ROOT) || defined(VECCORE_CUDA)
class GeomQueryHandler : public Handler {
#else
class GeomQueryHandler : public Handler,
                         public TGeoExtension {
#endif

protected:  
  Volume_t *fVolume = nullptr;         ///< Associated volume
  int fIndex = -1;                     ///< Handler index in the array of geometry handlers

private:
  GeomQueryHandler(const GeomQueryHandler &) = delete;
  GeomQueryHandler &operator=(const GeomQueryHandler &) = delete;
  
protected:
  VECCORE_ATT_HOST_DEVICE
  void ConnectToVolume();

  VECCORE_ATT_HOST_DEVICE
  void DisconnectVolume();

public:
  /** @brief Default constructor */
  VECCORE_ATT_HOST_DEVICE
  GeomQueryHandler() : Handler() {}

  /** 
   * @brief Volume handler default constructor
   * @param threshold Basketizing threshold
   * @param propagator Propagator working with this handler
   * @param vol Associated volume
   */
  VECCORE_ATT_HOST_DEVICE
  GeomQueryHandler(Volume_t *vol, int threshold, GeantPropagator *propagator, int index = -1);

  /** @brief Geometry filter destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~GeomQueryHandler();

  /** @brief Activate/de-activate the handler */
  VECCORE_ATT_HOST_DEVICE
  virtual void ActivateBasketizing(bool flag);

  /** @brief Scalar DoIt interface */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(GeantTrack *track, Basket& output, GeantTaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(Basket &input, Basket& output, GeantTaskData *td);

  /** @brief Getter for the index in the list of geometry handlers */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetIndex() const { return fIndex; }

  /** @brief Getter for handler number */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  Volume_t *GetVolume() const { return fVolume; }

  /**
   * @brief Grab function for TGeoExtension interface
   * @details Interface of TGeoExtension for getting a reference to this from Volume
   * @return Pointer to the base class
   */
#if defined(USE_ROOT) && !defined(VECCORE_CUDA)
  virtual TGeoExtension *Grab() { return this; }
  /**
   * @brief Release function for TGeoExtension interface
   * @details Interface of TGeoExtension to signal releasing ownership of this from TGeoVolume
   */
  virtual void Release() const {}
#endif  // USE_ROOT.

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
