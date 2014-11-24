#ifndef GEANT_COPROCESSORBROKER
#define GEANT_COPROCESSORBROKER

// Implementation of TaskBroker for CUDA devices.


#ifndef GEANT_TASKBROKER
#include "TaskBroker.h"
#endif

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

class GPGeomManager;
class GXFieldMap;
class GPPhysicsTable;
class GPVGeometry;
class GXFieldMap;
class GPPhysicsTable;
class GPPhysics2DVector;
struct GXTrack;
class GXTrackLiason;

const unsigned int kMaxNumStep = 1;
const unsigned int kMaxSecondaryPerStep = 2; // maxSecondaryPerStep;

class SecondariesTable
{
public:
   DevicePtr<int> fDevStackSize;
   //DevicePtr<int> fDevOffset;
   DevicePtr<GXTrack> fDevTracks;

   void Alloc(size_t maxTracks);
   void ToDevice() {
      cudaMemset(fDevStackSize,0,sizeof(int));
      //cudaMemset(fDevOffset,0,sizeof(int));
   }
};

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

class GeantTrack;
class GeantTrack_v;
class GeantBasket;

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

      unsigned int TrackToDevice(Task *task, int tid,
                                 GeantBasket &basket,
                                 unsigned int startIdx);

      unsigned int TrackToHost();

      GXTrack               *fTrack;
      int                   *fTrackId;   // unique indentifier of the track on the CPU
      int                   *fPhysIndex;
      int                   *fLogIndex;
      unsigned int           fChunkSize; // Max number of track for this stream
      unsigned int           fNStaged;   // How many track have been copied to the scratch area so far.

      unsigned int  fStreamId;
      cudaStream_t  fStream;
      DevicePtr<curandState> fdRandStates;
      DevicePtr<GXTrack>     fDevTrack;
      DevicePtr<GXTrack>     fDevAltTrack;
      DevicePtr<int>         fDevTrackPhysIndex;
      DevicePtr<int>         fDevTrackLogIndex;
      DevicePtr<int[10]>     fDevScratch;       // For communication internal to the kernel
      DevicePtr<GXTrackLiason> fDevTrackTempSpace;// For communication internal to the kernel
      SecondariesTable       fDevSecondaries;

      int          fThreadId;

      GeantBasket           *fBasket; // Work manager track (pre)-queue
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

   bool UploadGeometry(GPVGeometry *geom);
   bool UploadMagneticField(GXFieldMap** fieldMap);

   bool UploadPhysicsTable(const GPPhysicsTable *table, unsigned int nTables, GPPhysics2DVector* sbData, size_t maxZ);

   bool CudaSetup(int nblocks, int nthreads, int maxTrackPerThread);

   bool IsValid() { return fNthreads > 0; }
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
      Task(kernelFunc_t kernel) : fCurrent(0), fKernel(kernel), fCycles(0), fIdles(0) {}
      virtual ~Task() {}

      TaskData     *fCurrent;  // Holder of the data to be sent to the GPU, not owned.
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

