//===--- GeantTaskData.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantTaskData.h
 * @brief Implementation of data organized per thread Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TASKDATA
#define GEANT_TASKDATA

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef GEANT_TRACK
#include "GeantTrack.h"
#endif

#include <deque>

class TGeoVolume;
class TRandom;
class GeantBasketMgr;
class GeantBasket;

// SOA3D container from VecGeom
#include "base/SOA3D.h"

/**
 * @brief Class GeantTaskData
 * @details Class descripting data organized per thread
 * 
 */
namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {
class GeantTaskData {
public:
  Int_t fTid;          /** Thread unique id */
  Int_t fNthreads;     /** Number of transport threads */
  Int_t fMaxDepth;     /** Maximum geometry depth */
  Int_t fSizeBool;     /** Size of bool array */
  Int_t fSizeDbl;      /** Size of dbl array */
  Bool_t fToClean;     /** Flag set when the basket queue is to be cleaned */
  TGeoVolume *fVolume; /** Current volume per thread */
  TRandom *fRndm;      /** Random generator for thread */
  Bool_t *fBoolArray;  /** [fSizeBool] Thread array of bools */
  Double_t *fDblArray; /** [fSizeDbl] Thread array of doubles */
  GeantTrack fTrack;   /** Track support for this thread */
  VolumePath_t *fPath; /** Volume path for the thread */
  GeantBasketMgr *fBmgr; /** Basket manager collecting mixed tracks */
#ifdef GEANT_NVCC
   char fPool[sizeof(std::deque<GeantBasket*>)]; // Use the same space ...
#else
   std::deque<GeantBasket*> fPool; /** Pool of empty baskets */
#endif
   vecgeom::SOA3D<double>  *fSOA3Dworkspace1; // Thread SOA3D workspace (to be used for vector navigation)
   vecgeom::SOA3D<double>  *fSOA3Dworkspace2; // SOA3D workspace (to be used for vector navigation)
   int                      fSizeInt;  // current size of IntArray
   int                     *fIntArray; // Thread array of ints (used in vector navigation)

private:
   // a helper function checking internal arrays and allocating more space if necessary
  template <typename T> static void CheckSizeAndAlloc(T *&array, int &currentsize, size_t wantedsize) {
    if (wantedsize < currentsize)
      return;
    T *newarray = new T[wantedsize];
    delete[] array;
    currentsize = wantedsize;
  }

public:

   /** @brief GeantTaskData constructor */
  GeantTaskData();

  /** @brief GeantTaskData constructor */
  GEANT_CUDA_DEVICE_CODE
  GeantTaskData(Int_t nthreads, Int_t maxDepth, Int_t maxPerBasket);

  /** @brief GeantTaskData destructor */
  GEANT_CUDA_DEVICE_CODE
  ~GeantTaskData();

  /**
   * @brief Function that return double array
   * 
   * @param size Size of double array
   */
  Double_t *GetDblArray(Int_t size) {
    CheckSizeAndAlloc<Double_t>(fDblArray, fSizeDbl, size);
    return fDblArray;
  }

  /**
   * @brief Function that return boolean array
   * 
   * @param size Size of boolean array
   */
  Bool_t *GetBoolArray(Int_t size) {
    CheckSizeAndAlloc<Bool_t>(fBoolArray, fSizeBool, size);
    return fBoolArray;
  }

  /**
   * @brief Function that returns int array
   *
   * @param size Size of int array
   */
  int *GetIntArray(Int_t size) {
    CheckSizeAndAlloc<int>(fIntArray, fSizeInt, size);
    return fIntArray;
  }

  /**
   * @brief Function that returns a (per thread/task) preallocated SOA3D workspace
   *
   * @param size Size of container
   */
  vecgeom::SOA3D<double> *GetSOA3DWorkspace1(int size) {
    // TODO: should actually resize the workspace containers together
    // TODO: they should also be close together in memory
    if (size > fSOA3Dworkspace1->size()) {
      delete fSOA3Dworkspace1;
      fSOA3Dworkspace1 = new vecgeom::SOA3D<double>(2 * size);
    }
    return fSOA3Dworkspace1;
  }

  /**
   * @brief Function that returns a (per thread/task) preallocated SOA3D workspace
   *
   * @param size Size of container
   */
  vecgeom::SOA3D<double> *GetSOA3DWorkspace2(int size) {
    if (size > fSOA3Dworkspace2->size()) {
      delete fSOA3Dworkspace2;
      fSOA3Dworkspace2 = new vecgeom::SOA3D<double>(2 * size);
    }
    return fSOA3Dworkspace2;
  }

  /**
   * @brief Function that returns a (per thread/task) preallocated NavigationState object
   *
   */
  VolumePath_t *GetPath() {
    return fPath;
  }

  /**
   * @brief Get the cleared storedtrack
   */
  GeantTrack &GetTrack() { fTrack.Clear(); return fTrack; }

  /**
   * @brief Get next free basket or null if not available
   * @details Get pointer to next free basket
   */
  GeantBasket *GetNextBasket();

  /** @brief Setter for the toclean flag */
  void SetToClean(Bool_t flag) { fToClean = flag; }

  /** @brief Getter for the toclean flag */
  bool NeedsToClean() const { return fToClean; }

  /**
   * @brief Recycles a given basket
   * 
   * @param b Pointer to current GeantBasket for recycling
   */
  void RecycleBasket(GeantBasket *b);

  /**
   * @brief Function cleaning a number of free baskets
   * 
   * @param ntoclean Number of baskets to be cleaned
   * @return Number of baskets actually cleaned
   */
  Int_t CleanBaskets(size_t ntoclean);

private:

  /**
   * @brief Constructor GeantTaskData
   * @todo Still not implemented
   */
  GeantTaskData(const GeantTaskData &);

  /**
   * @brief Operator &operator=
   * @todo Still not implemented
   */
  GeantTaskData &operator=(const GeantTaskData &);

   // ClassDef(GeantTaskData, 1) // Stateful data organized per thread
};
} // GEANT_IMPL_NAMESPACE
} // Geant
#endif
