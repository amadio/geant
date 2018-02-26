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
#include "Track.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class TaskData;
class SimulationStage;
#include "Geant/Fwd.h"

/**
 * @brief A basket holding track pointers.
 * @details Basket of tracks in the same volume which are transported by a single thread
 */
 
class Basket {

protected:
  int fThreshold = 64;                ///< Basket threshold
  int fNode = 0;                      ///< Numa node for basket allocation
  SimulationStage *fStage = nullptr;  ///< Simulation stage to be executed by tracks inside
  TrackVec_t fTracks;                 ///< Vector of track pointers
  
private:

  Basket(const Basket &) = delete;
  Basket &operator=(const Basket &) = delete;
  
public:
  /** @brief Default basket constructor */
  VECCORE_ATT_HOST_DEVICE
  Basket() {}

  /** 
   * @brief Standard basket constructor
   *
   * @param size Initial basket size
   * @param threshold  Initial basket threshold
   */
  VECCORE_ATT_HOST_DEVICE
  Basket(int size, int threshold);

  /** 
   * @brief NUMA aware basket constructor
   *
   * @param size Initial basket size
   * @param threshold  Initial basket threshold
   * @param node NUMA node where the basket is alocated
   */
  VECCORE_ATT_HOST_DEVICE
  Basket(int size, int threshold, int node);

  /** @brief Basket destructor */
  VECCORE_ATT_HOST_DEVICE
  ~Basket() {}

  /**
   * @brief Add a track pointer to basket. Non thread-safe
   * @return Track index
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void AddTrack(Track *track) { fTracks.push_back(track); }

  /**
   * @brief Add several tracks to the basket.
   * @return Track index
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void AddTracks(TrackVec_t &tracks)
  {
#ifndef VECCORE_CUDA
    std::copy(tracks.begin(), tracks.end(), std::back_inserter(Tracks()));
#else
    for (auto track : tracks) fTracks.push_back(track);
#endif
  }

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
   * @brief Function for getting the simulation stage the basket will perform
   * @return  Pointer to simulation stage
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  SimulationStage *GetStage() const { return fStage; }

  /**
   * @brief Function for getting basket transportability threshold
   * @return  Value of transportability threshold
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetThreshold() const { return fThreshold; }

  /** @brief Function checking if a track is already contained */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int HasTrack(Track* const track) const
  {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    return ( std::find(fTracks.begin(), fTracks.end(), track) != fTracks.end() );
#else
    // #!@$f**k%#@! CUDA
    for (size_t i=0; i<fTracks.size(); ++i)
      if (fTracks[i] == track) return true;
    return false;
#endif
  }

  /** @brief Function checking if a track is contained repeatedly */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int HasTrackMany(Track* const track) const
  {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    auto it = std::find(fTracks.begin(), fTracks.end(), track);
    return ( std::find(++it, fTracks.end(), track) != fTracks.end() );
#else
    // #!@$f**k%#@! CUDA
    size_t i;
    for (i=0; i<fTracks.size(); ++i)
      if (fTracks[i] == track) break;
    for (size_t j=i+1; j<fTracks.size(); ++j)
      if (fTracks[j] == track) return true;
    return false;
#endif
  }
  /**
   * @brief Function returning a reference to the vector of input tracks
   * @return Reference to input vector of tracks
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  TrackVec_t &Tracks() { return fTracks; }

  /**
   * @brief Function returning a const reference to the vector of input tracks
   * @return Reference to input vector of tracks
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  TrackVec_t const &GetTracks() { return fTracks; }

  /**
   * @brief Function returning the number of tracks
   * @return Number of tracks
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t size() { return fTracks.size(); }

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
  void Recycle(TaskData *td);

  /**
   * @brief  Function that providing the size of this basket in bytes
   * @return Size of basket in bytes, not considering the space taken by the pointed tracks
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t SizeOfInstance() const {
    return ( sizeof(this) + fTracks.capacity() * sizeof(Track*) );
  }

  /**
   * @brief Function to change the transportability threshold for the basket. Not thread safe.
   *
   * @param threshold New threshold value
   */
  VECCORE_ATT_HOST_DEVICE
  void SetThreshold(int threshold);

  /**
   * @brief Function to change simulation stage
   *
   * @param stage New simulation stage
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetStage(SimulationStage *stage) { fStage = stage; }

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
