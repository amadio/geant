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

class GeantBasketMgr;
#include "GeantFwd.h"

namespace Geant {
#ifdef GEANT__NVCC
inline
#endif
namespace cuda {
   class GeantTrack_v;
   class GeantTaskData;
}
}

const unsigned int kMaxNumStep = 1;
const unsigned int kMaxSecondaryPerStep = 2; // maxSecondaryPerStep;


typedef
int (*kernelFunc_t)(vecgeom::cxx::DevicePtr<Geant::cuda::GeantTaskData> &workSpace, size_t workspaceSizeOf,
                    size_t ntracks,
                    vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrack_v> &input,
                    vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrack_v> &output,

                    int nBlocks, int nThreads,
                    cudaStream_t stream);
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

      Geant::GeantTaskData  *fGeantTaskData;
      GeantBasket           *fInputBasket;
      GeantBasket           *fOutputBasket; // Work manager track (pre)-queue

      unsigned int           fChunkSize; // Max number of track for this stream
      unsigned int           fNStaged;   // How many track have been copied to the scratch area so far.
      GeantBasketMgr        *fPrioritizer;

      int           fThreadId;
      unsigned int  fStreamId;
      cudaStream_t  fStream;

      size_t fDevTaskWorkspaceSizeOf;
      vecgeom::cxx::DevicePtr<Geant::cuda::GeantTaskData> fDevTaskWorkspace;
      vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrack_v>  fDevTrackInput;

      concurrent_queue      *fQueue;  // Queue recording whether this helper is available or not.

      vecgeom::cxx::DevicePtr<char> GetDevTrackInputBuf() {
         char *basket = (char*)&(*fDevTrackInput);
         return vecgeom::DevicePtr<char>( basket+ vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrack_v>::SizeOf() );
      }
      vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrack_v>  fDevTrackOutput;
      vecgeom::cxx::DevicePtr<char> GetDevTrackOutputtBuf() {
         char *basket = (char*)&(*fDevTrackOutput);
         return vecgeom::DevicePtr<char>( basket+  vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrack_v>::SizeOf() );
      }

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

   bool addTrack(int /* itr */, GeantBasket & /* basket */) { return false; }
   int SetPrioritizer() { return 0; }

   /** @brief Create the baskets for each stream */
   void CreateBaskets();

   bool UploadGeometry(vecgeom::VPlacedVolume const *const volume = nullptr);
   // bool UploadMagneticField(FieldMap** fieldMap);

   //bool UploadPhysicsTable(const PhysicsTable *table, unsigned int nTables, Physics2DVector* sbData, size_t maxZ);

   bool CudaSetup(int nblocks, int nthreads, int maxTrackPerThread);

   bool IsValid() { return fNthreads > 0; }
   bool AddTrack(GeantBasket &input, unsigned int trkid);
   void runTask(Geant::GeantTaskData &td, GeantBasket &basket);
   Stream launchTask(bool wait = false);
   Stream launchTask(Task *task, bool wait = false);
   void waitForTasks();

   long GetTotalWork() { return fTotalWork; }

   int GetNstream() { return fTaskData.size(); }

private:
   char                  *fdGeometry; // Point to a GPGeomManager in GPU land.

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
      virtual bool Select(Geant::cxx::GeantTrack_v &host_track, int track) = 0; // Selection routine.
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

