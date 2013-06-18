#ifndef GEANT_COPROCESSORBROKER
#define GEANT_COPROCESSORBROKER

// Implementation of TaskBroker for CUDA devices.


#ifndef GEANT_TASKBROKER
#include "TaskBroker.h"
#endif

class DevicePtrBase
{
   void *fPtr;
protected:
   void MemcpyToDevice(const void* what, unsigned long nbytes);
   void *GetPtr() const { return fPtr; }
public:
   DevicePtrBase() : fPtr(0) {}
   
   ~DevicePtrBase();

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
   operator T*() const { return reinterpret_cast<T*>(GetPtr()); }
};

class GPFieldMap;
class GPPhysicsTable;
class GPVGeometry;
class GPFieldMap;
class GPPhysicsTable;
struct GXTrack;

#ifndef __CINT__
#include <cuda.h>
#include <curand.h>
#include "random_kernel.h"
#else
class cudaStream_t;
class curandState;
#endif

class GeantTrack;
class GeantTrackCollection;

#include "TObject.h"
#include "sync_objects.h"

class CoprocessorBroker : public TaskBroker
{
private:
   struct StreamHelper;
public:
   typedef StreamHelper *Stream;
   Stream GetNextStream();
   
public:
   CoprocessorBroker();
   ~CoprocessorBroker();
   
   bool UploadGeometry(GPVGeometry *geom);
   bool UploadMagneticField(GPFieldMap** fieldMap);
   
   bool Upload_eBremTable(const GPPhysicsTable &eBrem_table);
   bool Upload_eIoniTable(const GPPhysicsTable &eIoni_table);
   bool UploadMscTable(const GPPhysicsTable &msc_table);
   
   bool CudaSetup(int nblocks, int nthreads, int maxTrackPerThread);

   void runTask(int threadid, int nTracks, int volumeIndex, GeantTrack **tracks, int *trackin);
   Stream launchTask(bool wait = false);
   void waitForTasks();
   
   long GetTotalWork() { return fTotalWork; }

private:
   char       *fdGeometry; // Point to a GPGeomManager in GPU land.
   DevicePtr<GPFieldMap> fdFieldMap;
   
   DevicePtr<GPPhysicsTable> fd_eBremTable;
   DevicePtr<GPPhysicsTable> fd_eIoniTable;
   DevicePtr<GPPhysicsTable> fd_mscTable;
   
   struct StreamHelper : public TaskBroker::StreamHelper, TObject {
      StreamHelper();
      ~StreamHelper();

      bool CudaSetup(unsigned int streamid, int nblocks, int nthreads, int maxTrackPerThread);

      unsigned int TrackToDevice(int tid,
                                 GeantTrack **host_track, int *trackin,
                                 unsigned int startIdx, unsigned int basketSize);

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
      DevicePtr<int>         fDevTrackPhysIndex;
      DevicePtr<int>         fDevTrackLogIndex;

      int          fThreadId;
      GeantTrack **fHostTracks;

      GeantTrackCollection  *fTrackCollection; // Work manager track (pre)-queue
      concurrent_queue      *fQueue;    // Queue recording whether this helper is available or not.

      GeantTrackCollection *GetNextTrackCollection();
      void ResetNStaged() { fNStaged = 0; }
      operator cudaStream_t() { return fStream; }
   
      void Reset();
      void Push(concurrent_queue *q = 0);
   
      
      ClassDef(StreamHelper,0);
   };

   friend void StreamReset(cudaStream_t /* stream */, cudaError_t status, void *userData);
   friend void TrackToHost(cudaStream_t /* stream */, cudaError_t status, void *userData);

   StreamHelper fStream[3];
   StreamHelper *fCurrentHelper;
   concurrent_queue fHelpers;
   
   int fNblocks;
   int fNthreads;
   int fMaxTrackPerThread;
   int fKernelType;
   int fTotalWork;
   
  
};

#endif

