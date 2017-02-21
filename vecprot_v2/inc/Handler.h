//===--- Handler.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Handler.h
 * @brief Implementation of a generic handler using a basketizer
 * @details A handler performs a specific operation in the stepping algorithm
 *  for a single track or for a group of tracks. A handler implementation must
 *  implement at least the scalar DoIt interface and optionally the vector DoIt.
*
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_HANDLER
#define GEANT_HANDLER

#include "Geant/Typedefs.h"
#include "Basketizer.h"
#include "Basket.h"
#include "GeantPropagator.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTaskData;
class GeantTrack;
#include "GeantFwd.h"

/**
 * @brief The base class for track handlers having a dedicated basketizer.
 */
 
class Handler {

public:
  using basketizer_t = Basketizer<GeantTrack>;

protected:  
  bool fActive = false;                ///< Activity flag
  int fBcap = 0;                       ///< Minimum capacity for the handled baskets
  std::atomic_int fThreshold;          ///< Basketizing threshold
  basketizer_t *fBasketizer = nullptr; ///< Basketizer for this handler
  GeantPropagator *fPropagator = nullptr; ///< Associated propagator
  std::atomic_flag fLock;              ///< Lock for flushing
private:
  Handler(const Handler &) = delete;
  Handler &operator=(const Handler &) = delete;
  
public:
  /** @brief Default handler constructor */
  VECCORE_ATT_HOST_DEVICE
  Handler() {}

  /** 
   * @brief Default constructor
   * @param threshold Basketizing threshold
   * @param propagator Propagator working with this handler
   * @param vol Associated volume
   */
  VECCORE_ATT_HOST_DEVICE
  Handler(int threshold, GeantPropagator *propagator);

  /** @brief Basket destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~Handler();

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
  int GetNode() const { return fPropagator->fNuma; }

  /** @brief Check if handler is active for basketizing */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsActive() const { return fActive; }

  /** @brief Activate/de-activate the handler */
  VECCORE_ATT_HOST_DEVICE
  virtual void ActivateBasketizing(bool flag = true);

  /** @brief Check if handler is flushing */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsFlushing() {
    if (fLock.test_and_set(std::memory_order_acquire)) return true;
    fLock.clear(std::memory_order_release);
    return false;
  }

  /** @brief Add a track pointer to handler. */
  VECCORE_ATT_HOST_DEVICE
  bool AddTrack(GeantTrack *track, Basket &collector);

  /** @brief Flush all tracks from the handler into a collector basket
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
