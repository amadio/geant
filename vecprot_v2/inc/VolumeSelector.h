//===--- VolumeSelector.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file VolumeSelector.h
 * @brief Implementation of a geometry volume selector
 * @details The volume selector performs computation of the allowed geometry step
 *  using both scalar and vector interfaces.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_VOLUME_SELECTOR
#define GEANT_VOLUME_SELECTOR

#include "Selector.h"

#ifdef USE_ROOT
#include "TGeoExtension.h"
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Selector grouping tracks by logical volume and performing geometry actions.
 */
 
#ifndef USE_ROOT
class VolumeSelector : public Selector {
#else
class VolumeSelector : public Selector,
                       public TGeoExtension {
#endif

protected:  
  Volume_t *fVolume = nullptr;         ///< Associated volume
  int fIndex = -1;                     ///< Selector index in the array of geometry selectors

private:
  VolumeSelector(const VolumeSelector &) = delete;
  VolumeSelector &operator=(const VolumeSelector &) = delete;
  
protected:
  VECCORE_ATT_HOST_DEVICE
  void ConnectToVolume();

  VECCORE_ATT_HOST_DEVICE
  void DisconnectVolume();

public:
  /** @brief Default constructor */
  VECCORE_ATT_HOST_DEVICE
  VolumeSelector() : Selector() {}

  /** 
   * @brief NUMA aware volume selector constructor
   *
   * @param threshold Basketizing threshold
   * @param propagator Propagator working with this selector
   * @param vol Associated volume
   * @param node NUMA node where the basket is alocated
   */
  VECCORE_ATT_HOST_DEVICE
  VolumeSelector(Volume_t *vol, int threshold, GeantPropagator *propagator,
         int node = -1, int index = -1);

  /** @brief Geometry filter destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~VolumeSelector();

  /** @brief Activate/de-activate the selector */
  VECCORE_ATT_HOST_DEVICE
  virtual void ActivateBasketizing(bool flag);

  /** @brief Scalar DoIt interface */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(GeantTrack *track, Basket& output, GeantTaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(Basket &input, Basket& output, GeantTaskData *td);

  /** @brief Getter for the index in the list of geometry selectors */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetIndex() const { return fIndex; }

  /** @brief Getter for selector number */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  Volume_t *GetVolume() const { return fVolume; }

  /**
   * @brief Grab function for TGeoExtension interface
   * @details Interface of TGeoExtension for getting a reference to this from Volume
   * @return Pointer to the base class
   */
 #ifdef USE_ROOT
  virtual TGeoExtension *Grab() { return this; }
  /**
   * @brief Release function for TGeoExtension interface
   * @details Interface of TGeoExtension to signal releasing ownership of this from TGeoVolume
   */
  virtual void Release() const {}
 #endif 

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
