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
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  // On cuda there is one propagator per thread.  So (for now), no need
  // for atomics.
  template <typename T>
  using atomic_t = T;
#else
  template <typename T>
  using atomic_t = std::atomic<T>;
#endif

public:
  using basketizer_t = Basketizer<GeantTrack>;

protected:  
  bool fActive = false;                ///< Activity flag
  int fBcap = 0;                       ///< Minimum capacity for the handled baskets
  atomic_t<int> fThreshold;          ///< Basketizing threshold
  basketizer_t *fBasketizer = nullptr; ///< Basketizer for this handler
  GeantPropagator *fPropagator = nullptr; ///< Associated propagator
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  std::atomic_flag fLock;              ///< Lock for flushing
#endif
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

  /** @brief NUMA node getter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  basketizer_t *GetBasketizer() const { return fBasketizer; }

  /** @brief Threshold getter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetThreshold() const { return fThreshold; }
  
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

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  /** @brief Check if handler is flushing */
  GEANT_FORCE_INLINE
  bool IsFlushing() {
    if (fLock.test_and_set(std::memory_order_acquire)) return true;
    fLock.clear(std::memory_order_release);
    return false;
  }
#endif

  /** @brief Add a track pointer to handler. */
  VECCORE_ATT_HOST_DEVICE
  bool AddTrack(GeantTrack *track, Basket &collector);

  /** @brief Flush all tracks from the handler into a collector basket
   *  @return Number of tracks flushed
   */
  VECCORE_ATT_HOST_DEVICE
  bool Flush(Basket &collector);

#ifndef VECCORE_CUDA
  /** @brief Get a free basket from the caller thread data storage
   *  @return Fresh basket pointer
   */
  Basket *GetFreeBasket(GeantTaskData *td);
#endif

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
