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
#include "Geant/Basketizer.h"
#include "Geant/Basket.h"
#include "Geant/Propagator.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class TaskData;
class Track;
#include "Geant/Fwd.h"

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
  using basketizer_t = Basketizer<Track>;

protected:
  bool fThreadLocal    = false;        ///< Handler is thread local
  bool fActive         = false;        ///< Activity flag
  bool fMayBasketize   = false;        ///< This handler can basketize
  bool fScalarDispatch = false;        ///< Dispatch baskets in scalar mode
  size_t fId           = 0;            ///< Handler id in the stage
  int fBcap            = 0;            ///< Minimum capacity for the handled baskets
  atomic_t<int> fThreshold;            ///< Basketizing threshold
  basketizer_t *fBasketizer = nullptr; ///< Basketizer for this handler
  Propagator *fPropagator   = nullptr; ///< Associated propagator
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  std::atomic_flag fLock; ///< Lock for flushing
  std::string fName;      ///< Handler name
#endif

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
  Handler(int threshold, Propagator *propagator);

  /** @brief Basket destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~Handler();

  /** @brief Copy constructor */
  VECCORE_ATT_HOST_DEVICE
  Handler(const Handler &);

  /** @brief Assignment operator */
  VECCORE_ATT_HOST_DEVICE
  Handler &operator=(const Handler &);

  /** @brief Scalar DoIt interface */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(Track *track, Basket &output, TaskData *td) = 0;

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(Basket &input, Basket &output, TaskData *td);

  /** @brief Scalar emulation of vector DoIt interface. Base class implements it as a loop. */
  VECCORE_ATT_HOST_DEVICE
  void DoItScalar(Basket &input, Basket &output, TaskData *td);

  /** @brief NUMA node getter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  basketizer_t *GetBasketizer() const { return fBasketizer; }

  /** @brief Handler id getter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t GetId() const { return fId; }

  /** @brief Handler id getter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetId(size_t id) { fId = id; }

  /** @brief Threshold getter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetThreshold() const { return fThreshold; }

  /** @brief Name getter */
  const char *GetName() const { return fName.c_str(); }

  /** @brief Name setter */
  void SetName(const char *name) { fName = name; }

  /** @brief Does this handler basketize */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool MayBasketize() const { return fMayBasketize; }

  /** @brief Setter for may basketize flag */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetMayBasketize(bool flag = true) { fMayBasketize = flag; }

  /** @brief Is this handler thread local */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsLocal() const { return fThreadLocal; }

  /** @brief Set this handler thread local */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetLocal(bool flag) { fThreadLocal = flag; }

  /** @brief Is this handler dispatching scalar */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsScalarDispatch() const { return fScalarDispatch; }

  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetScalarDispatch(bool flag = true) { fScalarDispatch = flag; }

  /** @brief NUMA node getter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetNode() const { return fPropagator->fNuma; }

  /** @brief Check if handler has basketized tracks */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  virtual bool HasTracks() const { return fBasketizer->GetNstored() > 0; }

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
  bool IsFlushing()
  {
    if (fLock.test_and_set(std::memory_order_acquire)) return true;
    fLock.clear(std::memory_order_release);
    return false;
  }
#endif

  /** @brief Add a track pointer to handler. */
  VECCORE_ATT_HOST_DEVICE
  virtual bool AddTrack(Track *track, Basket &collector, TaskData *td);

  /** @brief Flush all tracks from the handler into a collector basket
   *  @return Number of tracks flushed
   */
  VECCORE_ATT_HOST_DEVICE
  virtual bool Flush(Basket &collector, TaskData *td);
};

class LocalHandler : public Handler {

private:
  TrackVec_t fLocalBasket; ///< Local basket for thread local handler
  size_t fThr;             ///< Thread-local threshold
  Handler *fHandler;       ///< Handler to which this applies (no state)

  LocalHandler(const LocalHandler &) = delete;
  LocalHandler &operator=(const LocalHandler &) = delete;

public:
  LocalHandler(Handler *handler);
  virtual ~LocalHandler() {}
  /** @brief Scalar DoIt interface */
  VECCORE_ATT_HOST_DEVICE
  Handler *GetHandler() const { return fHandler; }

  /** @brief Activate/de-activate the handler */
  VECCORE_ATT_HOST_DEVICE
  virtual void ActivateBasketizing(bool flag) { fActive = flag; }

  /** @brief Scalar DoIt interface */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(Track *track, Basket &output, TaskData *td) { fHandler->DoIt(track, output, td); }

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  VECCORE_ATT_HOST_DEVICE
  virtual void DoIt(Basket &input, Basket &output, TaskData *td) { fHandler->DoIt(input, output, td); }

  /** @brief Add a track pointer to handler. */
  VECCORE_ATT_HOST_DEVICE
  virtual bool AddTrack(Track *track, Basket &collector, TaskData *td);

  /** @brief Flush all tracks from the handler into a collector basket
   *  @return Number of tracks flushed
   */
  VECCORE_ATT_HOST_DEVICE
  virtual bool Flush(Basket &collector, TaskData *td);

  /** @brief Check if handler has basketized tracks */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  virtual bool HasTracks() const { return fLocalBasket.size() > 0; }
};

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant

#endif
