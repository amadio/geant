//===--- Filter.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Filter.h
 * @brief Implementation of generic filter using a basketizer
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_FILTER
#define GEANT_FILTER

#ifdef USE_ROOT
#include "TGeoExtension.h"
#endif

#include "Geant/Typedefs.h"
#include "priority_queue.h"
#include "Basketizer.h"
#include "Basket.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTaskData;
class GeantTrack;
class GeantPropagator;
#include "GeantFwd.h"

/**
 * @brief A track filter with a dedicated basketizer.
 * @details Filters can have different locality types.
 */
 
#ifdef USE_ROOT
class Filter : public TGeoExtension {
#else
class Filter {
#endif
  using queue_t = priority_queue<Basket *>;
  using basketizer_t = Basketizer<GeantTrack>;
protected:  
  Basket::ELocality fLocality = Basket::kNone; ///< Locality type
  Volume_t *fVolume = nullptr;         ///< Associated volume if any.
  int fIndex = -1;                     ///< Filter index in the array of geometry filters
  int fNode = -1;                      ///< Numa node for basket allocation
  bool fActive = false;                ///< Activity flag
  int fBcap = 0;                       ///< Minimum capacity for the handled baskets
  std::atomic_int fThreshold;          ///< Basketizing threshold
  basketizer_t *fBasketizer = nullptr; ///< Basketizer for this filter
  queue_t *fFeeder = nullptr;          ///< Queue to which baskets get injected
private:
  Filter(const Filter &) = delete;
  Filter &operator=(const Filter &) = delete;
  
protected:
  VECCORE_ATT_HOST_DEVICE
  void ConnectToVolume();

  VECCORE_ATT_HOST_DEVICE
  void DisconnectVolume();

public:
  /** @brief Default filter constructor */
  VECCORE_ATT_HOST_DEVICE
  Filter() {}

  /** 
   * @brief NUMA aware filter constructor
   *
   * @param threshold Basketizing threshold
   * @param propagator Propagator working with this filter
   * @param vol Associated volume
   * @param node NUMA node where the basket is alocated
   */
  VECCORE_ATT_HOST_DEVICE
  Filter(int threshold, GeantPropagator *propagator,
         Basket::ELocality locality = Basket::kNone, int node = -1,
         int index = -1, Volume_t *vol = nullptr);

  /** @brief Basket destructor */
  VECCORE_ATT_HOST_DEVICE
  ~Filter();

   /** @brief Locality getter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  Basket::ELocality GetLocality() const { return fLocality; }

  /** @brief Threshold getter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetThreshold() const { return fThreshold.load(); }
  
  /** @brief Getter for the index in the list of geometry filters */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetIndex() const { return fIndex; }

  /** @brief NUMA node getter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetNode() const { return fNode; }

  /** @brief Getter for filter number */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  Volume_t *GetVolume() const { return fVolume; }

  /** @brief Getter for the feeder queue */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  queue_t *GetFeeder() const { return fFeeder; }

  /** @brief Setter for the feeder queue */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetFeeder(queue_t *feeder) { fFeeder = feeder; }

  /** @brief Activate/de-activate the filter */
  VECCORE_ATT_HOST_DEVICE
  void SetActive(bool flag = true);

  /** @brief Add a track pointer to filter. */
  VECCORE_ATT_HOST_DEVICE
  bool AddTrack(GeantTrack *track, Basket &collector);

  /** @brief Flush all tracks from the filter into the feeder
   *  @return Number of tracks flushed
   */
  VECCORE_ATT_HOST_DEVICE
  bool Flush(Basket &collector);
  
  /** @brief Get a free basket from the caller thread data storage
   *  @return Fresh basket pointer
   */
  VECCORE_ATT_HOST_DEVICE
  Basket *GetFreeBasket(GeantTaskData *td);

  /**
   * @brief Grab function
   * @details Interface of TGeoExtension for getting a reference to this from Volume
   * @return Pointer to the base class
   */
 #ifdef USE_ROOT
  virtual TGeoExtension *Grab() { return this; }
 #endif 

  /**
   * @brief Release function
   * @details Interface of TGeoExtension to signal releasing ownership of this from TGeoVolume
   */
  virtual void Release() const {}

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
