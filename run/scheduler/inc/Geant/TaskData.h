//===--- TaskData.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file TaskData.h
 * @brief Implementation of data organized per thread Geant-V prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TASKDATA
#define GEANT_TASKDATA

#include <deque>
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
#include <atomic>
#endif
#include "Geant/Typedefs.h"
#include "Geant/Error.h"
#include "Geant/mpmc_bounded_queue.h"
#include "Geant/NumaBlockMgr.h"
#include "Geant/Propagator.h"
#include "Geant/Track.h"

namespace geantphysics {
  class PhysicsData;
}

#include "base/RNG.h"

class GUFieldPropagator;
class VVectorField;

/**
 * @brief Class TaskData
 * @details Class descripting data organized per thread
 *
 */
namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class Basket;
class TrackGeo_v;
class StackLikeBuffer;
class TrackStat;
struct WorkspaceForFieldPropagation;
struct BasketCounters;

class TaskData {
private:
  Track *fTrack = nullptr; /** Blueprint track */

public:

  using NumaTrackBlock_t = NumaBlock<Track>;
  using UserDataVect_t = vector_t<char*>;

  Propagator *fPropagator = nullptr; /** Propagator */
  int fTid = -1;         /** Thread unique id */
  int fNode = -1;        /** Locality node */
  size_t fNthreads = 0;  /** Number of transport threads */
  int fSizeBool = 0;     /** Size of bool array */
  int fSizeInt = 0;      /*  Size of int array */
  int fSizeDbl = 0;      /** Size of dbl array */
  bool fToClean = false; /** Flag set when the basket queue is to be cleaned */
  bool fTaskCompleted;   /** Flag set when a task is completed */
  Volume_t *fVolume = nullptr; /** Current volume per thread */
  vecgeom::RNG *fRndm = nullptr;           /** Random generator for thread */
  bool *fBoolArray = nullptr;              /** [fSizeBool] Thread array of bools */
  double *fDblArray = nullptr;             /** [fSizeDbl] Thread array of doubles */
  int *fIntArray = nullptr;                /** [fSizeInt] Thread array of ints */
  VolumePath_t *fPath = nullptr;           /** Volume path for the thread */
  VolumePath_t **fPathV = nullptr;         /** Volume path for the thread */
  VolumePath_t **fNextpathV = nullptr;     /** Volume path for the thread */
  TrackGeo_v *fGeoTrack = nullptr;    /** Geometry track SOA */
  Basket *fBvector = nullptr;              /** Buffer basket used for vector API */
  Basket *fShuttleBasket = nullptr;        /** Shuttle basket from selectors to follow-up simulation stage */
  vector_t<Basket *> fStageBuffers;        /** Buffers for tracks at input of simulation stages */
  StackLikeBuffer *fStackBuffer = nullptr; /** Stack buffer tor this thread */
  TrackStat *fStat = nullptr;              /** Track statictics */
  NumaTrackBlock_t *fBlock = nullptr;      /** Current track block */
  BasketCounters *fCounters[kNstages];     /** Counters for stage handlers */

#ifdef VECCORE_CUDA
  char fBPool[sizeof(std::deque<Basket *>)]; /** Pool of empty baskets */
#else
  std::deque<Basket *> fBPool; /** Pool of empty baskets */
#endif
  vector_t<Track *> fTransported1;     // Transported tracks in current step
  int fNkeepvol = 0;     /** Number of tracks keeping the same volume */
  int fNsteps = 0;       /** Total number of steps per thread */
  int fNsnext = 0;       /** Total number of calls to getting distance to next boundary */
  int fNphys = 0;        /** Total number of steps to physics processes */
  int fNmag = 0;         /** Total number of partial steps in magnetic field */
  int fNpart = 0;        /** Total number of particles transported by the thread */
  int fNsmall = 0;       /** Total number of small steps taken */
  int fNcross = 0;       /** Total number of boundary crossings */
  int fNpushed = 0;      /** Total number of pushes with 1.E-3 */
  int fNkilled = 0;      /** Total number of tracks killed */

  geantphysics::PhysicsData  *fPhysicsData = nullptr; /** Physics data per thread */
  WorkspaceForFieldPropagation* fSpace4FieldProp = nullptr; /** Thread scratch for Field Propagation Stage */
  GUFieldPropagator       *fFieldPropagator; // For RK integration of charged particle propagation

private:
  UserDataVect_t fUserData;                /** User-defined data pointers */

private:
   // a helper function checking internal arrays and allocating more space if necessary
  template <typename T> static
  VECCORE_ATT_HOST_DEVICE
  void CheckSizeAndAlloc(T *&array, int &currentsize, size_t wantedsize) {
     if (wantedsize <= (size_t) currentsize)
      return;
    T *newarray = new T[wantedsize];
    memcpy(newarray,array,currentsize*sizeof(T));
    delete[] array;
    array = newarray;
    currentsize = wantedsize;
  }

  /**
   * @brief TaskData constructor based on a provided single buffer.
   */
  VECCORE_ATT_DEVICE
  TaskData(void *addr, size_t nTracks, int maxPerBasket, Propagator *prop = nullptr);

public:
  /** @brief TaskData constructor */
  TaskData(size_t nthreads, int maxPerBasket);

  /** @brief TaskData destructor */
  ~TaskData();

  /** @brief Attach a propagator on a numa node. */
  VECCORE_ATT_HOST_DEVICE
  void AttachPropagator(Propagator *prop, int node);
  
  /**
   * @brief Track MakeInstance based on a provided single buffer.
   */
  VECCORE_ATT_DEVICE
  static TaskData *MakeInstanceAt(void *addr, size_t nTracks, int maxPerBasket, Propagator *prop);

  /** @brief return the contiguous memory size needed to hold a Track_v */
  VECCORE_ATT_DEVICE
  static size_t SizeOfInstance(size_t nthreads, int maxPerBasket);

  /**
   * @brief Function that return double array
   *
   * @param size Size of double array
   */
  VECCORE_ATT_HOST_DEVICE
  double *GetDblArray(int size) {
    CheckSizeAndAlloc<double>(fDblArray, fSizeDbl, size);
    return fDblArray;
  }

  /**
   * @brief Function that return boolean array
   *
   * @param size Size of boolean array
   */
  VECCORE_ATT_HOST_DEVICE
  bool *GetBoolArray(int size) {
    CheckSizeAndAlloc<bool>(fBoolArray, fSizeBool, size);
    return fBoolArray;
  }

  /**
   * @brief Function that returns int array
   *
   * @param size Size of int array
   */
  VECCORE_ATT_HOST_DEVICE
  int *GetIntArray(int size) {
    CheckSizeAndAlloc<int>(fIntArray, fSizeInt, size);
    return fIntArray;
  }

  /**
   * @brief Function that returns a (per thread/task) preallocated NavigationState object
   *
   */
  VECCORE_ATT_HOST_DEVICE
  VolumePath_t *GetPath() {
    return fPath;
  }

  /** @brief Get new track from track manager */
  VECCORE_ATT_HOST_DEVICE
  Track &GetNewTrack();

  /** @brief Release the temporary track */
  VECCORE_ATT_HOST_DEVICE
  void ReleaseTrack(Track &track);

#ifndef VECCORE_CUDA
  /**
   * @brief Recycles a given basket
   *
   * @param b Pointer to current basket for recycling
   */
  void RecycleBasket(Basket *b);

#endif

  /** @brief Setter for the toclean flag */
  void SetToClean(bool flag) { fToClean = flag; }

  /** @brief Getter for the toclean flag */
  bool NeedsToClean() const { return fToClean; }

  /** @brief  Inspect simulation stages */
  void InspectStages(int istage);

  /** @brief  Set user data */
  bool SetUserData(void *data, size_t index) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    if (index >= fUserData.size())
      fUserData.resize(index + 1);
    else if (fUserData[index]) return false;
    fUserData[index] = (char*)data;
#endif
    return true;
  }

  /** @brief  Get user data */
  void *GetUserData(size_t index) { return fUserData[index]; }

private:
  /**
   * @brief Constructor TaskData
   * @todo Still not implemented
   */
  TaskData(const TaskData &);

  /**
   * @brief Operator &operator=
   * @todo Still not implemented
   */
  TaskData &operator=(const TaskData &);

  // ClassDef(TaskData, 1) // Stateful data organized per thread
};


/** &brief Class representing a user data handle to be user with user task data */
template <typename T>
class TaskDataHandle {
 friend class TDManager;
private:
  size_t fIndex;   // Data index
  char fName[50];  // Name of the handle (thanks CUDA#@!)
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  std::atomic_flag fMergeLock;        // Lock for merging user data
#endif

  TaskDataHandle(const char *name, size_t index) : fIndex(index)
  { strcpy(fName, name); }
public:
  TaskDataHandle(const TaskDataHandle &other) : fIndex(other.fIndex) {
    memcpy(fName, other.fName, 50);
  }

  TaskDataHandle &operator=(const TaskDataHandle &other) {
    fIndex = other.fIndex;
    memcpy(fName, other.fName, 50);
  }

  GEANT_FORCE_INLINE
  T *GetUserData(TaskData *td) {
    return (T*)td->GetUserData(fIndex);
  }
  
  GEANT_FORCE_INLINE
  T &operator()(TaskData *td) {
    return *(T*)td->GetUserData(fIndex);
  }

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  GEANT_FORCE_INLINE
  bool TryLock() { return !fMergeLock.test_and_set(std::memory_order_acquire); }

  GEANT_FORCE_INLINE
  void ClearLock() { fMergeLock.clear(std::memory_order_release); }
#endif

  GEANT_FORCE_INLINE
  const char *GetName() { return fName; }


/** @brief User callable, allowing to attach per-thread data of the handle type
  * @details User data corresponding to all pre-defined tokens can be allocated
  *          in MyApplication::AttachUserData, to avoid run-time checks */
  bool AttachUserData(T *data, TaskData *td) {
    return ( td->SetUserData(data, fIndex) );
  }
};

/** &brief The task data manager, distributing task data objects concurrently */
class TDManager {

using queue_t = mpmc_bounded_queue<TaskData*>;

private:
  int fMaxThreads = 0;   // Maximum number of threads
  int fMaxPerBasket = 0; // Maximum number of tracks per basket
  queue_t fQueue; // Task data queue
  vector_t<TaskData*> fTaskData;
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  std::atomic<size_t> fUserDataIndex; // Single registry point for user data
#else
  int fUserDataIndex; // CUDA does not need this
#endif
public:
  TDManager(int maxthreads, int maxperbasket) : fMaxThreads(maxthreads), fMaxPerBasket(maxperbasket), fQueue(1024), fUserDataIndex(0) {
    // Create initial task data
    for (int tid=0; tid < maxthreads; ++tid) {
      TaskData *td = new TaskData(fMaxThreads, fMaxPerBasket);
      td->fTid = fTaskData.size();
      fTaskData.push_back(td);
      fQueue.enqueue(td);
    }
  }

  ~TDManager() {
    // User data deleted by user in MyApplication::DeleteTaskData()
    for (auto td : fTaskData) delete td;
  }

  GEANT_FORCE_INLINE
  size_t GetNtaskData() const { return fTaskData.size(); }

  TaskData *GetTaskData() {
    TaskData *td;
    if (fQueue.dequeue(td)) return td;
    td =  new TaskData(fMaxThreads, fMaxPerBasket);
    td->fTid = fTaskData.size();
    fTaskData.push_back(td);
    return td;
  }

  GEANT_FORCE_INLINE
  TaskData *GetTaskData(int index) { return fTaskData[index]; }

  void ReleaseTaskData(TaskData *td) {
    while (!fQueue.enqueue(td)) {}
  }

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  template <typename T>
  TaskDataHandle<T> *RegisterUserData(const char *name) {
    return ( new TaskDataHandle<T>(name, fUserDataIndex.fetch_add(1)) );
  }

  template <typename T>
  void DeleteUserData(TaskDataHandle<T> &handle) {
    for (auto td : fTaskData) delete handle(td);
  }
#endif

  template <typename T>
  T* MergeUserData(int evslot, TaskDataHandle<T> &handle) {
    // if (handle.TryLock()) return nullptr;
    TaskData *base = fTaskData.front();
    for (auto td : fTaskData) {
      if (td == base) continue;
      // This will require the method T::Merge(int evslot, T const &other) to exist
      // The user data should book as many event slots as GeantConfig::fNbuff
      if (!handle(base).Merge(evslot, handle(td))) {
        Error("MergeUserData", "Cannot recursively merge %s data", handle.GetName());
        return nullptr;
      }
      handle(td).Clear(evslot);
    }
    handle.ClearLock();
    return &handle(base);
  }

};

} // GEANT_IMPL_NAMESPACE
} // Geant
#endif
