//===--- Basket.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Basket.h
 * @brief Implementation of baskets of tracks for Geant-V prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_BASKET_NEW
#define GEANT_BASKET_NEW

#ifndef USE_ROOT
#define MIC_BIT(n) (1ULL<<(n))
#endif

#include "Geant/Typedefs.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTaskData;
class GeantTrack;
#include "GeantFwd.h"

/**
 * @brief A basket holding track pointers.
 * @details Basket of tracks in the same volume which are transported by a single thread
 */
 
class Basket {

public:
/** Basket locality types */
enum ELocality {
  kNone,     // No locality criteria
  kVolume,   // Same geometry volume
  kParticle, // Same particle type
  kProcess   // Same physics process
};

/** Basket processing stages */
enum EStage {
  kScheduling,      // Scheduling actions: track fetching, flushing, prioritizing
  kSampleXsec,      // Propose physics step by sampling total Xsec
  kTransportLength, // Compute geometry transport length
  kPropagate,       // Propagation in field
  kAlongStep,       // Along step actions (e.g. continuous energy loss)
  kPostStep         // Post step actions
};

protected:
  ELocality fLocality = Basket::kNone; ///< Locality type
  int fThreshold = 64;                ///< Basket threshold
  int fNode = 0;                      ///< Numa node for basket allocation
  vector_t<GeantTrack *>  fTracks;    ///< Vector of track pointers
  
private:

  Basket(const Basket &) = delete;
  Basket &operator=(const Basket &) = delete;
  
public:
  /** @brief Default basket constructor */
  VECCORE_ATT_HOST_DEVICE
  Basket() {}

  /** 
   * @brief NUMA aware basket constructor
   *
   * @param size Initial basket size
   * @param loc  Initial basket locality type
   * @param node NUMA node where the basket is alocated
   */
  VECCORE_ATT_HOST_DEVICE
  Basket(int size, ELocality loc = Basket::kNone, int threshold = 0, int node = 0);

  /** @brief Basket destructor */
  VECCORE_ATT_HOST_DEVICE
  ~Basket() {}

  /**
   * @brief Add a track pointer to basket. Non thread-safe
   * @return Track index
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void AddTrack(GeantTrack *track) { fTracks.push_back(track); }

  /** @brief Clearing the basket content, no deletion */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void Clear() { fTracks.clear(); }

  /**
   * @brief Check if a basket contains tracks in a given event range
   *
   * @param evstart Start event id.
   * @param nevents Number of events (default 1)
   * @return Truth value
   */
  VECCORE_ATT_HOST_DEVICE
  bool Contains(int evstart, int nevents = 1) const;

  /**
   * @brief Function returning the basket locality type
   * @return Locality type
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  ELocality GetLocality() const { return fLocality; }

  /**
   * @brief Function setting the basket locality type
   * @return Locality type
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetLocality(ELocality loc) { fLocality = loc; }

  /**
   * @brief Function returning the number tracks in the basket
   * @return Number of tracks
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetNtracks() const { return fTracks.size(); }

  /**
   * @brief Function for getting the NUMA locality for the basket
   * @return  Value of NUMA locality
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetNode() const { return fNode; }

  /**
   * @brief Function for getting basket transportability threshold
   * @return  Value of transportability threshold
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetThreshold() const { return fThreshold; }

  /**
   * @brief Function returning the locality type as string.
   * @return Boolean value if the tracks are mixed from several volumes
   */
  VECCORE_ATT_HOST_DEVICE
  const char *GetLocalityString() const;

  /**
   * @brief Function returning the mixed tracks property.
   * @return Boolean value if the tracks are mixed from several volumes
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsMixed() const { return ( fLocality == Basket::kNone ); }

  /**
   * @brief Function returning a reference to the vector of input tracks
   * @return Reference to input vector of tracks
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  vector_t<GeantTrack*> &Tracks() { return fTracks; }

  /**
   * @brief Function returning a const reference to the vector of input tracks
   * @return Reference to input vector of tracks
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  vector_t<GeantTrack*> const &GetTracks() { return fTracks; }

  /**
   * @brief Print the basket content
   */
  void Print(const char *msg = "") const;

  /**
   * @brief Print the parameters for a given track
   *
   * @param itr Track id.
   * @param input Refer to input or output track (default input)
   */
  VECCORE_ATT_HOST_DEVICE
  void PrintTrack(int itr) const;

  /** @brief Recycle this basket */
  VECCORE_ATT_HOST_DEVICE
  void Recycle(GeantTaskData *td);

  /**
   * @brief  Function that providing the size of this basket in bytes
   * @return Size of basket in bytes, not considering the space taken by the pointed tracks
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t SizeOfInstance() const {
    return ( sizeof(this) + fTracks.capacity() * sizeof(GeantTrack*) );
  }

  /**
   * @brief Function to change the transportability threshold for the basket. Not thread safe.
   *
   * @param threshold New threshold value
   */
  void SetThreshold(int threshold);

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
