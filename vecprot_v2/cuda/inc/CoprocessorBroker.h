#ifndef GEANT_COPROCESSORBROKER
#define GEANT_COPROCESSORBROKER

// Implementation of TaskBroker for CUDA devices.

#ifndef GEANT_CONFIG_H
#include "Geant/Config.h"
#endif

#ifndef GEANT_TASKBROKER
#include "TaskBroker.h"
#endif

// This should be part of a global (configure time generated) header.
#ifndef VECGEOM_CUDA
#define VECGEOM_CUDA
#endif
#include "backend/cuda/Interface.h"

#ifndef __CINT__
#include <cuda.h>
#include <curand.h>
#include "random_kernel.h"
#else
class cudaStream_t;
class curandState;
#endif

class DevicePtrBase
{
   void *fPtr;
   
   DevicePtrBase(const DevicePtrBase&); // not implemented
   DevicePtrBase &operator=(const DevicePtrBase&); // not implemented

protected:
   void MemcpyToDevice(const void* what, unsigned long nbytes);
   void MemcpyToHostAsync(void* where, unsigned long nbytes, cudaStream_t stream);
   void *GetPtr() const { return fPtr; }
public:
   DevicePtrBase() : fPtr(0) {}

   virtual ~DevicePtrBase();

   void Malloc(unsigned long size);
};

template <typename T>
class DevicePtr : public DevicePtrBase
{
public:
   void Alloc(unsigned long nelems = 1) {
      Malloc(nelems*sizeof(T));
   }
   void ToDevice(const T* what, unsigned long nelems = 1) {
      MemcpyToDevice(what,nelems*sizeof(T));
   }
   void FromDevice(T* where,cudaStream_t stream) {
      // Async since we past a stream.
      MemcpyToHostAsync(where,sizeof(T),stream);
   }
   void FromDevice(T* where, unsigned long nelems , cudaStream_t stream) {
      // Async since we past a stream.
      MemcpyToHostAsync(where,nelems*sizeof(T),stream);
   }
#ifndef __CINT__
   operator T*() const { return reinterpret_cast<T*>(GetPtr()); }
#endif
};

class GeantBasketMgr;

class GPGeomManager;
class GXFieldMap;
class GPPhysicsTable;
class GPVGeometry;
class GXFieldMap;
class GPPhysicsTable;
class GPPhysics2DVector;
struct GXTrack;
class GXTrackLiason;

class TaskWorkspace;
class GeantTrack_v;

const unsigned int kMaxNumStep = 1;
const unsigned int kMaxSecondaryPerStep = 2; // maxSecondaryPerStep;

class SecondariesTable
{
public:
   DevicePtr<int> fDevStackSize;
   //DevicePtr<int> fDevOffset;
   DevicePtr<GeantTrack_v> fDevTracks;

   void Alloc(size_t maxTracks);
   void ToDevice() {
      cudaMemset(fDevStackSize,0,sizeof(int));
      //cudaMemset(fDevOffset,0,sizeof(int));
   }
};

#ifdef GXTRACKING_KERNELS
typedef
int (*kernelFunc_t)(curandState* devStates,
                    size_t nSteps,
                    size_t nTrackSize,
                    GXTrack *track, GXTrack *altTrack,
                    int *logVolumeIndices,
                    int *physVolumeIndices,
                    GXTrack *secondaries, int *secStackSize,

                    int *scratch, // array of 10
                    GXTrackLiason *scratchTrack,

                    GPGeomManager* geomManager,
                    GXFieldMap* magMap,
                    GPPhysicsTable* PhysicsTable,
                    GPPhysics2DVector* SeltzerBergerTable,

                    int nBlocks, int nThreads,
                    cudaStream_t stream);
#else
typedef
int (*kernelFunc_t)(DevicePtr<TaskWorkspace> &workSpace,
                    size_t ntracks,
                    DevicePtr<GeantTrack_v> &input,
                    DevicePtr<GeantTrack_v> &output,

                    int nBlocks, int nThreads,
                    cudaStream_t stream);
#endif

class GeantTrack;
class GeantTrack_v;
class GeantBasket;

#include "GeantTrack.h"
#include "TObject.h"
#include "sync_objects.h"

class CoprocessorBroker : public TaskBroker
{
public:
   struct Task;

   struct TaskData : public TaskBroker::TaskData, TObject {
   private:
      TaskData(const TaskData&); // not implemented
      TaskData& operator=(const TaskData&); // not implemented
   public:
      
      TaskData();
      ~TaskData();

      bool CudaSetup(unsigned int streamid, int nblocks, int nthreads, int maxTrackPerThread);

      unsigned int AddTrack(Task *task,
                            GeantBasket &basket,
                            unsigned int hostIdx);
      unsigned int TrackToDevice(Task *task, int tid,
                                 GeantBasket &basket,
                                 unsigned int startIdx);

      unsigned int TrackToHost();

      GeantBasketMgr        *fBasketMgr;
      GeantBasket           *fInputBasket;
      GeantBasket           *fOutputBasket; // Work manager track (pre)-queue

      unsigned int           fChunkSize; // Max number of track for this stream
      unsigned int           fNStaged;   // How many track have been copied to the scratch area so far.
      GeantBasketMgr        *fPrioritizer;

      int           fThreadId;
      unsigned int  fStreamId;
      cudaStream_t  fStream;

      DevicePtr<TaskWorkspace> fDevTaskWorkspace;
      DevicePtr<GeantTrack_v>  fDevTrackInput;
      vecgeom::cxx::DevicePtr<char> GetDevTrackInputBuf() {
         char *basket = (char*)&(*fDevTrackInput);
         return vecgeom::DevicePtr<char>( basket+sizeof(GeantTrack_v) );
      }
      DevicePtr<GeantTrack_v>  fDevTrackOutput;
      vecgeom::cxx::DevicePtr<char> GetDevTrackOutputtBuf() {
         char *basket = (char*)&(*fDevTrackOutput);
         return vecgeom::DevicePtr<char>( basket+sizeof(GeantTrack_v) );
      }
      SecondariesTable         fDevSecondaries;

      concurrent_queue      *fQueue;  // Queue recording whether this helper is available or not.

      void ResetNStaged() { fNStaged = 0; }
      operator cudaStream_t() { return fStream; }

      void Reset();
      void Push(concurrent_queue *q = 0);

      ClassDef(TaskData,0);
   };

public:
   typedef TaskData *Stream;
   Stream GetNextStream();

   CoprocessorBroker(const CoprocessorBroker&); // not implemented
   CoprocessorBroker& operator=(const CoprocessorBroker&); // not implemented
public:
   CoprocessorBroker();
   ~CoprocessorBroker();

   bool UploadGeometry(vecgeom::VPlacedVolume const *const volume = nullptr);
   bool UploadMagneticField(GXFieldMap** fieldMap);

   bool UploadPhysicsTable(const GPPhysicsTable *table, unsigned int nTables, GPPhysics2DVector* sbData, size_t maxZ);

   bool CudaSetup(int nblocks, int nthreads, int maxTrackPerThread);

   bool IsValid() { return fNthreads > 0; }
   bool AddTrack(GeantBasket &input, unsigned int trkid);
   void runTask(int threadid, GeantBasket &basket);
   Stream launchTask(bool wait = false);
   Stream launchTask(Task *task, bool wait = false);
   void waitForTasks();

   long GetTotalWork() { return fTotalWork; }

   int GetNstream() { return fTaskData.size(); }

private:
   char                  *fdGeometry; // Point to a GPGeomManager in GPU land.
   DevicePtr<GXFieldMap>  fdFieldMap;

   DevicePtr<GPPhysicsTable>    fd_PhysicsTable;
   DevicePtr<GPPhysics2DVector> fd_SeltzerBergerData;

public:
   struct Task {
   private:
      Task(const Task&); // not implemented
      Task& operator=(const Task&); // not implemented
   public:
      Task(kernelFunc_t kernel) : fCurrent(0), fPrevNStaged(0), fKernel(kernel), fCycles(0), fIdles(0) {}
      virtual ~Task() {}

      bool IsReadyForLaunch(unsigned int ntasks);

      TaskData     *fCurrent;  // Holder of the data to be sent to the GPU, not owned.
      unsigned int  fPrevNStaged; // Value of fCurrent->fNStaged at the time of the last call to IsReadyForLaunch
      kernelFunc_t  fKernel;   // wrapper around the cuda call to the kernel.
      unsigned int  fCycles;   // Number of times we put track in the taskData without launching (age of the data).
      unsigned int  fIdles;   // Number of times we processed basket without putting track in the taskData.

      virtual const char *Name() = 0; // Id of the task.
      virtual bool Select(GeantTrack_v &host_track, int track) = 0; // Selection routine.
   };

private:

   friend void StreamReset(cudaStream_t /* stream */, cudaError_t status, void *userData);
   friend void TrackToHost(cudaStream_t /* stream */, cudaError_t status, void *userData);

   typedef std::vector<Task*> TaskColl_t;
   typedef std::vector<TaskData*> TaskDataColl_t;
   TaskColl_t     fTasks;
   TaskDataColl_t fTaskData;

   TaskData *fNextTaskData;
   concurrent_queue fHelpers;

   int fNblocks;  // Number of cuda blocks
   int fNthreads; // Number of cuda threads
   int fMaxTrackPerThread; // Max number of tracks per cuda threads
   int fTotalWork;
};

#endif

