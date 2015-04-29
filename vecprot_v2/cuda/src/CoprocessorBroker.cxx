#include "CoprocessorBroker.h"


// Needed by GXEMPhysicsUtils.h for now
#include <iostream>
#include <fstream>
typedef double G4double;

#include "GPConstants.h"
#include "GXFieldMap.h"
#include "GXFieldMapData.h"

#include "GXTrack.h"
#include "GPThreeVector.h"

//Geometry
#include "GPVGeometry.h"
#include "GPUserGeometry.h"
#include "GPSimpleEcal.h"
#include "GPVPhysicalVolume.h"

//EMPhysics
#include "GPConstants.h"
#include "GPPhysics2DVector.h"

#include "GPPhysicsTable.h"
//#include "GPPhysicsTableType.h"

#include "GeantTrack.h"

#include "GXTrack.h"
#include "GXTrackLiason.h"

#include "GPEMPhysicsUtils.h"
//#include "GXRunManager.h"
//#include "GPEMPhysicsUtils.h"

#include <vector>
#include <math.h>

// #include <cuda_runtime.h>

#include "CoprocessorBrokerKernel.h"

#include "TGeoBranchArray.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

// To access the list of baskets.
#include "WorkloadManager.h"
#include "GeantScheduler.h"
#include "GeantBasket.h"
#include "GeantThreadData.h"
#include "globals.h"

#include "management/CudaManager.h"

static void HandleCudaError( cudaError_t err,
                             const char *file,
                             int line ) {
    if (err != cudaSuccess) {
       Fatal("Cuda","%s (%d) in %s at line %d\n", cudaGetErrorString( err ), err,
             file, line );
       exit( EXIT_FAILURE );
    }
}

#define HANDLE_CUDA_ERROR( err ) (HandleCudaError( err, __FILE__, __LINE__ ))

struct GeneralTask : public CoprocessorBroker::Task {
   GeneralTask() : Task( PropagateGeantTrack_gpu ) {}

   const char *Name() { return "GeneralTask"; }
   bool Select(GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle electron, which we pretend are the only
      // particle to have charge -1.

      return -1 == host_track.fChargeV[track];
   }
};

#ifdef GXTRACKING_KERNELS
struct GeneralTask : public CoprocessorBroker::Task {
   GeneralTask() : Task( tracking_gpu ) {}

   const char *Name() { return "GeneralTask"; }
   bool Select(GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle electron, which we pretend are the only
      // particle to have charge -1.

      return -1 == host_track.fChargeV[track];
   }
};

struct EnergyTask : public CoprocessorBroker::Task {
   double  fThresHold;
   TString fName;

   EnergyTask(double threshold) : Task( tracking_gpu ), fThresHold(threshold)
   {
      fName.Form("EnergyTask - %3g",threshold);
   }

   const char *Name() { return fName; }
   bool Select(GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle electron, which we pretend are the only
      // particle to have charge -1.

      return -1 == host_track.fChargeV[track] && host_track.fEV[track] > fThresHold;
   }
};

struct EnergyElectronTask : public CoprocessorBroker::Task {
   double  fThresHold;
   TString fName;

   EnergyElectronTask(double threshold) : Task( electron_multistage_gpu ), fThresHold(threshold)
   {
      fName.Form("EnergyElectronTask - %3g",threshold);
   }

   const char *Name() { return fName; }
   bool Select(GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle electron, which we pretend are the only
      // particle to have charge -1.

      return -1 == host_track.fChargeV[track] && host_track.fEV[track] > fThresHold;
   }
};

struct EnergyElectronSingleTask : public CoprocessorBroker::Task {
   double  fThresHold;
   TString fName;

   EnergyElectronSingleTask(double threshold) : Task( electron_gpu ), fThresHold(threshold)
   {
      fName.Form("EnergyElectronSingleTask - %3g",threshold);
   }

   const char *Name() { return fName; }
   bool Select(GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle electron, which we pretend are the only
      // particle to have charge -1.

      return -1 == host_track.fChargeV[track] && host_track.fEV[track] > fThresHold;
   }
};
#endif


#if 0
struct ElectronTask : public CoprocessorBroker::Task {
   ElectronTask() : Task( tracking_electron_gpu ) {}

   const char *Name() { return "ElectronTask"; }
   bool Select(GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle electron, which we pretend are the only
      // particle to have charge -1.

      return -1 == host_track.fChargeV[track];
   }
};

struct PhotonTask : public CoprocessorBroker::Task {
   PhotonTask() : Task( tracking_photon_gpu ) {}

   const char *Name() { return "PhotonTask"; }
   bool Select(GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle photon, which we pretend are the only
      // particle to have charge 0.

      return 0 == host_track.fChargeV[track];
   }
};
#endif

void DevicePtrBase::Malloc(unsigned long size) {
   HANDLE_CUDA_ERROR( cudaMalloc((void**)&fPtr,size) );
}

DevicePtrBase::~DevicePtrBase() {
   if (fPtr) HANDLE_CUDA_ERROR( cudaFree(fPtr) );
}

void DevicePtrBase::MemcpyToDevice(const void* what, unsigned long nbytes)
{
   HANDLE_CUDA_ERROR(cudaMemcpy(fPtr,what,nbytes,
                                cudaMemcpyHostToDevice) );
}

void DevicePtrBase::MemcpyToHostAsync(void* where, unsigned long nbytes, cudaStream_t stream)
{
   HANDLE_CUDA_ERROR(cudaMemcpyAsync(where, fPtr, nbytes, cudaMemcpyDeviceToHost, stream));
}

void SecondariesTable::Alloc(size_t maxTracks)
{
   fDevStackSize.Alloc();
   //fDevOffset.Alloc();
   fDevTracks.Alloc(maxTracks * kMaxNumStep * kMaxSecondaryPerStep);
}

GXFieldMap **ReadFieldMap(const char *fieldMapFile)
{
   GXFieldMap** fieldMap;
   fieldMap = (GXFieldMap **) malloc (nbinZ*sizeof(GXFieldMap *));
   for (int j = 0 ; j < nbinZ ; j++) {
      fieldMap[j] = (GXFieldMap *) malloc (nbinR*sizeof(GXFieldMap));
   }
   std::ifstream ifile(fieldMapFile, std::ios::in | std::ios::binary | std::ios::ate);
   if (ifile.is_open()) {

      //field map structure
      GXFieldMapData fd;

      std::ifstream::pos_type fsize = ifile.tellg();
      size_t dsize = sizeof(GXFieldMapData);

      long int ngrid = fsize/dsize;
      ifile.seekg (0,  std::ios::beg);

      std::cout << "... transportation ... Loading magnetic field map: "
      << fieldMapFile << std::endl;

      for(int i = 0 ; i < ngrid ; i++) {
         ifile.read((char *)&fd, sizeof(GXFieldMapData));

         //check validity of input data
         if(abs(fd.iz) > noffZ || fd.ir > nbinR) {
            std::cout << " Field Map Array Out of Range" << std::endl;
         }
         else {
            fieldMap[noffZ+fd.iz][fd.ir].Bz = fd.Bz;
            fieldMap[noffZ+fd.iz][fd.ir].Br = fd.Br;
         }
      }
      ifile.close();
   } else {
      Error("ReadFieldMap","Could not open the file %s\n",fieldMapFile);
   }
   return fieldMap;
}

bool ReadPhysicsTable(GPPhysicsTable &table, const char *filename, bool useSpline)
{

   readTable(&table,filename);
   unsigned int nv = table.nPhysicsVector;
   for(unsigned int j=0; j < nv; j++){
      table.physicsVectors[j].SetSpline(useSpline);
   }
   //  GPPhysicsTable_Print(&table);
   return true;
}


CoprocessorBroker::TaskData::TaskData() : fBasketMgr(new GeantBasketMgr(WorkloadManager::Instance()->GetScheduler(), 0, 0, true)),
                                          fInputBasket(0),fOutputBasket(0),
                                          fChunkSize(0),fNStaged(0),
                                          fPrioritizer(0),
                                          fThreadId(-1),
                                          fStreamId(0),
                                          fQueue(0)
{
   // Default constructor.
}

CoprocessorBroker::TaskData::~TaskData() {
   // Destructor.

   // We do not own fQueue.
   // We do not own the fBasket(s) (or do we?)
   // We do not own fPrioritizer

   delete fBasketMgr;
   HANDLE_CUDA_ERROR( cudaStreamDestroy(fStream) );
}

bool CoprocessorBroker::TaskData::CudaSetup(unsigned int streamid, int nblocks, int nthreads, int maxTrackPerThread)
{

   fStreamId = streamid;
   HANDLE_CUDA_ERROR( cudaStreamCreate(&fStream) );

   //prepare random engines on the device
   // fdRandStates.Alloc( nblocks*nthreads );
   // curand_setup_gpu(fdRandStates, time(NULL), nblocks, nthreads);

   unsigned int maxTrackPerKernel = nblocks*nthreads*maxTrackPerThread;
   fChunkSize = maxTrackPerKernel;

   unsigned int maxThreads = nblocks*nthreads;

   fDevTaskWorkspace.Alloc(maxThreads);
   fDevTrackInput.Alloc(maxTrackPerKernel);
   fDevTrackOutput.Alloc(maxTrackPerKernel);
   // fDevSecondaries.Alloc(maxTrackPerKernel);
   fBasketMgr->SetBcap(maxTrackPerKernel);

   fInputBasket = fBasketMgr->GetNextBasket();
   fOutputBasket = fBasketMgr->GetNextBasket();
   return true;
}

void CoprocessorBroker::TaskData::Push(concurrent_queue *q)
{
   // Add this helper to the queue and record it as the
   // default queue.

   if (q) {
      fQueue = q;
      fQueue->push(this);
   } else if (fQueue) {
      fQueue->push(this);
   }
}

void CoprocessorBroker::TaskData::Reset()
{
   // Get the stream helper ready for re-use.
   // by setting fNStaged to zero and add the helper
   // to the list of available helpers.

   ResetNStaged();
   Push();
}

CoprocessorBroker::CoprocessorBroker() : fdGeometry(0)
   ,fNextTaskData(0)
   ,fNblocks(0),fNthreads(0)
   ,fMaxTrackPerThread(1)
   ,fTotalWork(0)
//,fdFieldMap(0)
//,fdStates0(0),fdStates1(0),fdRandStates0(0),fdRandStates1(0)
{
   // Default constructor.

}

CoprocessorBroker::~CoprocessorBroker()
{
   HANDLE_CUDA_ERROR( cudaDeviceSynchronize() );
   for(unsigned int i = 0; i < fTasks.size(); ++i) {
      delete fTasks[i];
   }
   fTasks.clear();
   for(unsigned int i = 0; i < fTaskData.size(); ++i) {
      delete fTaskData[i];
   }
   fTaskData.clear();

   HANDLE_CUDA_ERROR( cudaFree(fdGeometry) );
//   cudaFree(fdFieldMap);
//   cudaFree(fd_eBremTable);
//   cudaFree(fd_eIoniTable);
//   cudaFree(fd_mscTable);
//   cudaFree(fdRandStates0);
//   cudaFree(fdRandStates1);

}

bool CoprocessorBroker::UploadGeometry(vecgeom::VPlacedVolume const *const volume)
{
   // Prepare the geometry for the device and upload it to the device's memory.

   //CudaManager::Instance().set_verbose(3);

   if (volume) vecgeom::CudaManager::Instance().LoadGeometry(volume);
   else        vecgeom::CudaManager::Instance().LoadGeometry();
   vecgeom::CudaManager::Instance().Synchronize();

   //CudaManager::Instance().PrintGeometry();
   //std::cout << std::flush;

   return true;
}

bool CoprocessorBroker::UploadMagneticField(GXFieldMap** fieldMap)
{

   fprintf(stderr, "Magnetic field uploading not implemented yet.\n");

//    printf("Creating magnetic field map on the GPU device\n");
//    //prepare fiesldMap array: fieldMap[nbinZ][nbinR];
//    GXFieldMap *bmap_h = (GXFieldMap *) malloc (nbinZ*nbinR*sizeof(GXFieldMap));
//    for (int i = 0 ; i < nbinZ ; i++) {
//       for (int j = 0 ; j < nbinR ; j++) {
//          bmap_h[i+j*nbinZ].Bz = fieldMap[i][j].Bz;
//          bmap_h[i+j*nbinZ].Br = fieldMap[i][j].Br;
//       }
//    }

//    printf("Copying data from host to device\n");

// //   cudaMalloc((void**)&fdFieldMap, nbinZ*nbinR*sizeof(GXFieldMap)) ;
// //   cudaMemcpy(fdFieldMap,bmap_h,nbinZ*nbinR*sizeof(GXFieldMap),
// //              cudaMemcpyHostToDevice);
//    fdFieldMap.Alloc(nbinZ*nbinR);
//    fdFieldMap.ToDevice(bmap_h,nbinZ*nbinR);
//    free(bmap_h);

   return true;
}

bool CoprocessorBroker::UploadPhysicsTable(const GPPhysicsTable *table, unsigned int nTables, GPPhysics2DVector* sbData, size_t maxZ)
{
   fprintf(stderr, "Physics table uploading not implemented yet.\n");

   // fd_PhysicsTable.Alloc(nTables);
   // fd_PhysicsTable.ToDevice(table,nTables);

   // fd_SeltzerBergerData.Alloc(maxZ);
   // fd_SeltzerBergerData.ToDevice(sbData,maxZ);

   return true;
}

void setup(CoprocessorBroker *broker,
           int nphi = 4,
           int nz   = 3,
           double density = 8.28)
{

   // //2. Read magnetic field map
   // const char* fieldMapFile = getenv("GP_BFIELD_MAP");
   // fieldMapFile = (fieldMapFile) ? fieldMapFile : "data/cmsExp.mag.3_8T";
   // GXFieldMap** fieldMap = ReadFieldMap(fieldMapFile);

   // //3. Create magnetic field on the device
   // broker->UploadMagneticField(fieldMap);

   // // 4. Prepare EM physics tables
   // GPPhysicsTable physicsTable[kNumberPhysicsTable];

   // TString filename;
   // for(int it = 0 ; it < kNumberPhysicsTable ; ++it) {
   //    filename.Form("data/%s",GPPhysicsTableName[it]);
   //    readTableAndSetSpline(&physicsTable[it],filename);
   // }

   // //G4SeltzerBergerModel data
   // unsigned int maxZ = 92;
   // GPPhysics2DVector* sbData = new GPPhysics2DVector[maxZ];

   // char sbDataFile[256];
   // for(unsigned int iZ = 0 ; iZ < maxZ ; ++iZ) {
   //    sprintf(sbDataFile,"data/brem_SB/br%d",iZ+1);
   //    std::ifstream fin(sbDataFile);
   //    G4bool check = RetrieveSeltzerBergerData(fin, &sbData[iZ]);
   //    if(!check) {
   //       printf("Failed To open SeltzerBerger Data file for Z= %d\n",iZ+1);
   //    }
   // }
   // printf("Copying GPPhysicsTable data from host to device\n");
   // broker->UploadPhysicsTable(physicsTable,kNumberPhysicsTable,sbData,maxZ);

   // delete [] sbData;

}

bool CoprocessorBroker::CudaSetup(int nblocks, int nthreads, int maxTrackPerThread)
{
   int deviceCount = 0;
   cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

   if (error_id != cudaSuccess) {
     Printf("Cuda CoprocessorBroker disabled because of issue #%d:%s\n", (int)error_id, cudaGetErrorString(error_id));
     return false;
   }

   setup(this);

   fNblocks = nblocks;
   fNthreads = nthreads;
   fMaxTrackPerThread = maxTrackPerThread;

   fTasks.push_back(new GeneralTask());
   /*
   fTasks.push_back(new EnergyElectronTask(6));
   fTasks.push_back(new EnergyElectronTask(4));
   fTasks.push_back(new EnergyElectronTask(2));
   fTasks.push_back(new EnergyElectronTask(0));
   */

   // Initialize global constants.
   HANDLE_CUDA_ERROR( cudaMemcpyToSymbol("gPropagator_fBmag", &(gPropagator->fBmag), sizeof(double), size_t(0), cudaMemcpyHostToDevice) );

   double tolerance = TGeoShape::Tolerance();
   HANDLE_CUDA_ERROR( cudaMemcpyToSymbol("gTolerance", &(tolerance), sizeof(double), size_t(0), cudaMemcpyHostToDevice) );

   //initialize the stream
   for(unsigned int i=0; i < 2+fTasks.size(); ++i) {
      TaskData *data = new TaskData();
      data->CudaSetup(i,nblocks,nthreads,maxTrackPerThread);
      data->Push(&fHelpers);
      fTaskData.push_back(data);
   }
   return true;
}

// bool CanHandleTrack(GeantTrack &track)
// {
//    // Currently we can only handle electron, which we pretend are the only
//    // particle to have charge -1.

//    return -1 == track.charge;
// }


unsigned int CoprocessorBroker::TaskData::AddTrack(CoprocessorBroker::Task *task,
                                                   GeantBasket &basket,
                                                   unsigned int hostIdx)
{

   if (!fInputBasket) {
      GeantBasketMgr *bmgr = basket.GetBasketMgr();
      fInputBasket = bmgr->GetNextBasket();
   }

   GeantTrack_v &input = basket.GetInputTracks();

   //  if (input.fHoles->TestBitNumber(hostIdx)) return 0;
   if (task->Select(input,hostIdx)) {
      fInputBasket->AddTrack(input,hostIdx);
      ++fNStaged;
      return 1;
   }
   return 0;
}

unsigned int CoprocessorBroker::TaskData::TrackToDevice(CoprocessorBroker::Task *task,
                                                        int tid,
                                                        GeantBasket &basket,
                                                        unsigned int startIdx)
{
   if (fThreadId == -1) {
      fThreadId = tid;
   } else if (fThreadId != tid) {
      Fatal("TrackToDevice","The stream %p is already assigned to thread %d when requested to do work for thread %d",
            this, fThreadId, tid);
      return 0;
   }

   if (!fInputBasket) {
      GeantBasketMgr *bmgr = basket.GetBasketMgr();
      fInputBasket = bmgr->GetNextBasket();
   }

   //unsigned int start = fNStaged;
   unsigned int count = 0;
   GeantTrack_v &input = basket.GetInputTracks();
   unsigned int basketSize = input.GetNtracks();
   for(unsigned int hostIdx = startIdx;
       fNStaged < fChunkSize && hostIdx < basketSize;
       ++hostIdx ) {

      ++count;

      if (input.fHoles->TestBitNumber(hostIdx))
         continue;

      if (task->Select(input,hostIdx)) {

         fInputBasket->AddTrack(input,hostIdx);

         ++fNStaged;
      }
   }

   // Track sub-ranges are not consecutive in memory.
   // // count = min(fChunkSize,basketSize-startIdx);
   // int ntrack = fNStaged - start;
   // HANDLE_CUDA_ERROR(cudaMemcpyAsync(fDevInputTrack+start, input+start, ntrack*sizeof(GXTrack),
   //                                   cudaMemcpyHostToDevice, fStream));
   // HANDLE_CUDA_ERROR(cudaMemcpyAsync(fDevTrackLogIndex+start, fLogIndex+start, ntrack*sizeof(int),
   //                                   cudaMemcpyHostToDevice, fStream));
   // HANDLE_CUDA_ERROR(cudaMemcpyAsync(fDevTrackPhysIndex+start, fPhysIndex+start, ntrack*sizeof(int),
   //                                   cudaMemcpyHostToDevice, fStream));

   return count;
}

unsigned int CoprocessorBroker::TaskData::TrackToHost()
{
   WorkloadManager *mgr = WorkloadManager::Instance();
   GeantScheduler *sch = mgr->GetScheduler();
   condition_locker &sched_locker = mgr->GetSchLocker();


   Int_t ntot = 0;
   Int_t nnew = 0;
   Int_t nkilled = 0;
   /* Int_t ninjected = */ sch->AddTracks(fOutputBasket, ntot, nnew, nkilled, nullptr /* fPrioritizer */);
   (void)ntot;
   (void)nnew;
   (void)nkilled;
   //mgr->TransportedQueue()->push(fBasket);
   sched_locker.StartOne();
   // fOutputBasket->Recycle();
   // fOutputBasket = 0;
   fOutputBasket->Clear();
   fThreadId = -1;
   return fNStaged;
}

bool CoprocessorBroker::Task::IsReadyForLaunch(unsigned int ntasks)
{
   // Return true if the stream ought to be launch because it is either full or
   // filling slowly, etc.

   if (fCurrent->fNStaged == fCurrent->fChunkSize)
      return true;

   // We do not have enough tracks

   if ( fPrevNStaged == fCurrent->fNStaged ) {
      ++( fIdles );
   } else {
      fIdles = 0;
   }
   ++( fCycles );
   fPrevNStaged = fCurrent->fNStaged;

   unsigned int idle = fIdles;
   unsigned int cycle = fCycles;
   if (fCurrent->fNStaged                // There is something
       && idle > (ntasks-1)     // We are beyond the expected/normal number of idles cycles
       && cycle > (fCurrent->fChunkSize / fCurrent->fNStaged)  // To fill we need at least that many cycles
       && 2*idle > cycle               // Our input rate has drop in half
       )
   {
      // Printf("(%d - GPU) ================= Launching idle Task %s Stream %d Idle=%d cycle=%d accumulated=%d", threadid, (*task)->Name(), stream->fStreamId, idle, cycle, stream->fNStaged);
      // if we did not make any progress in a while, assume there is no 'interesting' track left and schedule the kernel.
      return true;
   }

   // Continue to wait for more data ...
   return false;
}


// typedef void(CUDART_CB * cudaStreamCallback_t)(cudaStream_t stream, cudaError_t status, void
//                                                *userData)

void StreamReset(cudaStream_t /* stream */, cudaError_t status, void *userData)
{
   CoprocessorBroker::TaskData *helper = (CoprocessorBroker::TaskData*)userData;
   helper->Reset();
}

void TrackToHost(cudaStream_t /* stream */, cudaError_t status, void *userData)
{
   CoprocessorBroker::TaskData *helper = (CoprocessorBroker::TaskData*)userData;
   helper->TrackToHost();
}


CoprocessorBroker::Stream CoprocessorBroker::GetNextStream()
{
   // Return the current stream (one that we can still add input to)
   // or return new one.

   if (!fNextTaskData) {
      fNextTaskData = dynamic_cast<TaskData*>(fHelpers.wait_and_pop());
      if (!fNextTaskData) {
         // nothing we can do at the moment
         return 0;
      }
   }
   return fNextTaskData;
}

CoprocessorBroker::Stream CoprocessorBroker::launchTask(Task *task, bool wait /* = false */)
{

   Stream stream = task->fCurrent;
   task->fIdles = 0;
   task->fCycles = 0;
   task->fCurrent = 0;

   //Printf("(%d - GPU) == Starting kernel for task %s using stream %d with %d tracks\n",
   //       stream->fThreadId, task->Name(), stream->fStreamId, stream->fNStaged );

   //HANDLE_CUDA_ERROR(cudaMemcpyAsync(stream->fDevTrackInput->fBuf, stream->fInputBasket->fBuf, stream->fInputBasket->fBufSize,
   //                                  cudaMemcpyHostToDevice, *stream));

   fTotalWork += stream->fNStaged;
   DevicePtr<TaskWorkspace> nullDevicePtr;
   int result = task->fKernel(nullDevicePtr,
                              stream->fNStaged,
                              stream->fDevTrackInput,
                              stream->fDevTrackOutput,

                              fNblocks,fNthreads,*stream);

   if (!result) return stream;

   // Bring back the number of secondaries created.
   // int stackSize;
   // stream->fDevSecondaries.fDevStackSize.FromDevice(&stackSize, *stream);
   // stream->fDevSecondaries.fTrack.FromDevice( stream->fSecondaries, stackSize, *stream);

   // Bring back the modified tracks.
   // HANDLE_CUDA_ERROR(cudaMemcpyAsync(stream->fDevTrackOutput->fBuf, stream->fOutputBasket->fBuf,
   //                                   stream->fOutputBasket->fBufSize,
   //                                   cudaMemcpyHostToDevice, *stream));

   HANDLE_CUDA_ERROR(cudaStreamAddCallback(stream->fStream, TrackToHost, stream, 0 ));
   HANDLE_CUDA_ERROR(cudaStreamAddCallback(stream->fStream, StreamReset, stream, 0 ));

   if (wait) {
      // Use this when you need to insure the printf are actually printed.
      HANDLE_CUDA_ERROR(cudaStreamSynchronize(*stream));
   }
   // cudaDeviceSynchronize

   return stream;
}

CoprocessorBroker::Stream CoprocessorBroker::launchTask(bool wait /* = false */)
{
   Stream stream = 0;
   TaskColl_t::iterator task = fTasks.begin();
   TaskColl_t::iterator end = fTasks.end();
   while( task != end ) {

      if ((*task)->fCurrent != 0 && (*task)->fCurrent->fNStaged != 0) {

         stream = launchTask((*task),wait);
      }
      ++task;
   }
   // Return the last stream used ....
   return stream;
}


bool CoprocessorBroker::AddTrack(GeantBasket &input, unsigned int trkid)
{
   // return true if the track has been add to one of the coprocessor's task.

   bool force = false;
   //unsigned int nTracks = 1;
   unsigned int trackUsed = 0;
   //unsigned int trackStart = 0;

   TaskColl_t::iterator task = fTasks.begin();
   TaskColl_t::iterator end = fTasks.end();
   while( task != end ) {

      if ((*task)->fCurrent == 0) {
         if (fNextTaskData) {
            (*task)->fCurrent = fNextTaskData;
         } else {
            // If we do not yet have a TaskData, wait for one.
            // Consequence: if we are very busy we hang on to this
            // track until one of the task finishes.
            // It is likely better to returned that we could not
            // handle the track.
            (*task)->fCurrent = GetNextStream();
         }
         fNextTaskData = 0;
         if (!(*task)->fCurrent) break;
      }

      TaskData *stream = (*task)->fCurrent;

      unsigned int before = stream->fNStaged;
      /* unsigned int count = */ stream->AddTrack(*task,input,trkid);
      trackUsed += stream->fNStaged-before;
      // unsigned int rejected = nTracks-trackStart - (stream->fNStaged-before);
      if ( !force && !(*task)->IsReadyForLaunch(fTasks.size()) ) {
         // Continue to wait for more data ...
         break;
      }

      launchTask(*task);

      ++task;
   }

   // Missing compared to runTask si the scheduling of the most loaded task
   // if nothing new at all in a while.

   if (trackUsed) return true;
   else return false;
}


void CoprocessorBroker::runTask(int threadid, GeantBasket &basket)
                                // unsigned int nTracks, int volumeIndex, GeantTrack **tracks, int *trackin)
{
   bool force = false;

   unsigned int nTracks = basket.GetNinput();
   //<<<---------------------------------------------------------->>>
//         fprintf(stderr,"round = %d count1 = %d fNchunk=%d nTracks=%d\n",
//                 round,count1,fNchunk,nTracks);
//         cudaMemcpyAsync(track_d0, track_h+i, count1*sizeof(GXTrack),
//                         cudaMemcpyHostToDevice, stream0);
//

   unsigned int trackUsed = 0;
   TaskColl_t::iterator task = fTasks.begin();
   TaskColl_t::iterator end = fTasks.end();
   while( task != end ) {
      unsigned int trackLeft  = nTracks;
      unsigned int trackStart = 0;
      while (trackLeft) {
         if ((*task)->fCurrent == 0) {
            if (fNextTaskData) {
               (*task)->fCurrent = fNextTaskData;
            } else {
               // If we do not yet have a TaskData, wait for one.
               // Consequence: if we are very busy we hang on to this
               // basket until one of the task finishes.
               (*task)->fCurrent = GetNextStream();
            }
            fNextTaskData = 0;
            if (!(*task)->fCurrent) break;
         }
         TaskData *stream = (*task)->fCurrent;

         unsigned int before = stream->fNStaged;
         unsigned int count = stream->TrackToDevice(*task,
                                                    threadid,
                                                    basket,
                                                    trackStart);
         trackUsed += stream->fNStaged-before;
         // unsigned int rejected = nTracks-trackStart - (stream->fNStaged-before);

         //if (((*task)->fCycles % 10000) == 1)
         //   Printf("(%d - GPU) ================= Task %s Stream %d Tracks: %d seen %d skipped %d accumulated %d idles %d cycles", threadid, (*task)->Name(), stream->fStreamId, count, rejected, stream->fNStaged, (*task)->fIdles, (*task)->fCycles);

         if ( !force && !(*task)->IsReadyForLaunch(fTasks.size()) ) {
            // Continue to wait for more data ...
            break;
         }

         launchTask(*task);
         trackLeft  -= count;
         trackStart += count;
      }
      ++task;
   }

   if (trackUsed == 0) {
      // if we did not make any progress in a while, assume there is no 'interesting' track left
      // and schedule the most loaded task.

      Task *heavy = 0;
      TaskColl_t::iterator task = fTasks.begin();
      TaskColl_t::iterator end = fTasks.end();
      while( task != end ) {
         if ((*task)->fCurrent && (*task)->fCurrent->fNStaged) {
            if (heavy == 0 || (*task)->fCurrent->fNStaged > heavy->fCurrent->fNStaged) {
               heavy = *task;
            }
         }
         ++task;
      }
      if (heavy) {
         // Printf("(%d - GPU) ================= Launching heavy Task %s Stream %d Idle=%d cycle=%d accumulated=%d", threadid, heavy->Name(), heavy->fCurrent->fStreamId, heavy->fIdles, heavy->fCycles, heavy->fCurrent->fNStaged);
         launchTask(heavy);
      }
   }

   if (trackUsed != nTracks) {
      GeantTrack_v &input  = basket.GetInputTracks();
      GeantTrack_v &output = basket.GetOutputTracks();
      for(unsigned int hostIdx = 0 ;hostIdx < nTracks;
          ++hostIdx ) {

         if (input.fHoles->TestBitNumber(hostIdx))
            continue;

         input.PostponeTrack(hostIdx,output);
      }
   }
   if (!basket.GetInputTracks().IsCompact()) basket.GetInputTracks().Compact();
}

void CoprocessorBroker::waitForTasks()
{
   // Make sure all the tasks are finished.

   for(unsigned int i=0; i < fTaskData.size(); ++i) {
      if (fTaskData[i]->fNStaged) {
         HANDLE_CUDA_ERROR(cudaStreamSynchronize(*fTaskData[i]));
      }
   }
}
