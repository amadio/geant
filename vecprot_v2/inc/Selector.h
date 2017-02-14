//===--- Selector.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Selector.h
 * @brief Implementation of a generic selector using a basketizer
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_SELECTOR
#define GEANT_SELECTOR

#include "Geant/Typedefs.h"
#include "Basketizer.h"
#include "Basket.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTaskData;
class GeantTrack;
class GeantPropagator;
#include "GeantFwd.h"

/**
 * @brief The base class for track selectors having a dedicated basketizer.
 */
 
class Selector {

public:
  using basketizer_t = Basketizer<GeantTrack>;

protected:  
  int fNode = -1;                      ///< Numa node for basket allocation
  bool fActive = false;                ///< Activity flag
  int fBcap = 0;                       ///< Minimum capacity for the handled baskets
  std::atomic_int fThreshold;          ///< Basketizing threshold
  basketizer_t *fBasketizer = nullptr; ///< Basketizer for this selector
  GeantPropagator *fPropagator = nullptr; ///< Associated propagator
  std::atomic_flag fLock;              ///< Lock for flushing
private:
  Selector(const Selector &) = delete;
  Selector &operator=(const Selector &) = delete;
  
public:
  /** @brief Default selector constructor */
  VECCORE_ATT_HOST_DEVICE
  Selector() {}

  /** 
   * @brief NUMA aware selector constructor
   *
   * @param threshold Basketizing threshold
   * @param propagator Propagator working with this selector
   * @param vol Associated volume
   * @param node NUMA node where the basket is alocated
   */
  VECCORE_ATT_HOST_DEVICE
  Selector(int threshold, GeantPropagator *propagator, int node = -1);

  /** @brief Basket destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~Selector();

  /** @brief Scalar DoIt interface */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(GeantTrack *track, Basket& output, GeantTaskData *td) = 0;

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(Basket &input, Basket& output, GeantTaskData *td);

  /** @brief Threshold getter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetThreshold() const { return fThreshold.load(); }
  
  /** @brief NUMA node getter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetNode() const { return fNode; }

  /** @brief Check if selector is active for basketizing */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsActive() const { return fActive; }

  /** @brief Activate/de-activate the selector */
  VECCORE_ATT_HOST_DEVICE
  virtual void ActivateBasketizing(bool flag = true);

  /** @brief Check if selector is flushing */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsFlushing() {
    if (fLock.test_and_set(std::memory_order_acquire)) return true;
    fLock.clear(std::memory_order_release);
    return false;
  }

  /** @brief Add a track pointer to selector. */
  VECCORE_ATT_HOST_DEVICE
  bool AddTrack(GeantTrack *track, Basket &collector);

  /** @brief Flush all tracks from the selector into a collector basket
   *  @return Number of tracks flushed
   */
  VECCORE_ATT_HOST_DEVICE
  bool Flush(Basket &collector);
  
  /** @brief Get a free basket from the caller thread data storage
   *  @return Fresh basket pointer
   */
  VECCORE_ATT_HOST_DEVICE
  Basket *GetFreeBasket(GeantTaskData *td);

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
