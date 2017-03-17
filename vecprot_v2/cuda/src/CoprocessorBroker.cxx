#include "CoprocessorBroker.h"

#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>

#include "CoprocessorBrokerKernel.h"

// To access the list of baskets.
#include "GeantBasket.h"
#include "GeantScheduler.h"
#include "GeantTaskData.h"
#include "WorkloadManager.h"
#include "globals.h"

#include "GeantCudaUtils.h"
#include "management/CudaManager.h"

#include "GeantBasket.h"
#include "GeantEvent.h"
#include "GeantOutput.h"
#include "GeantScheduler.h"
#include "GeantTaskData.h"
#include "GeantTrack.h"
#include "GeantVApplication.h"
#include "PhysicsProcessOld.h"

using namespace Geant;

struct GeneralTask : public CoprocessorBroker::Task {
  GeneralTask() : Task(PropagateGeantTrack_gpu) {}

  const char *Name() { return "GeneralTask"; }
  bool Select(Geant::GeantTrack_v & /* host_track */, int /* track */)
  {
    // Currently we can only handle electron, which we pretend are the only
    // particle to have charge -1.

    return true;
  }
};

struct GeneralChargedTask : public CoprocessorBroker::Task {
  GeneralChargedTask(short charge) : Task(PropagateGeantTrack_gpu), fCharge(charge) {}
  short fCharge;

  const char *Name() { return "GeneralChargedTask"; }
  bool Select(Geant::GeantTrack_v &host_track, int track)
  {
    // Currently we can only handle electron, which we pretend are the only
    // particle to have charge -1.

    return fCharge == host_track.fChargeV[track];
  }
};

#if 0
struct ElectronTask : public CoprocessorBroker::Task {
   ElectronTask() : Task( tracking_electron_gpu ) {}

   const char *Name() { return "ElectronTask"; }
   bool Select(Geant::GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle electron, which we pretend are the only
      // particle to have charge -1.

      return -1 == host_track.fChargeV[track];
   }
};

struct PhotonTask : public CoprocessorBroker::Task {
   PhotonTask() : Task( tracking_photon_gpu ) {}

   const char *Name() { return "PhotonTask"; }
   bool Select(Geant::GeantTrack_v &host_track, int track)
   {
      // Currently we can only handle photon, which we pretend are the only
      // particle to have charge 0.

      return 0 == host_track.fChargeV[track];
   }
};
#endif

CoprocessorBroker::TaskData::TaskData()
    : fGeantTaskData(0), fInputBasket(0), fOutputBasket(0), fDevMaxTracks(0), fChunkSize(0), fNStaged(0),
      fPrioritizer(0), fThreadId(-1), fStreamId(0), fStream(nullptr), fDevTaskWorkspaceSizeOf(0), fQueue(0)
{
  // Default constructor.
}

CoprocessorBroker::TaskData::~TaskData()
{
  // Destructor.

  // We do not own fQueue.
  // We do not own the fBasket(s) (or do we?)
  // We do not own fPrioritizer

  delete fGeantTaskData;
  GEANT_CUDA_ERROR(cudaStreamDestroy(fStream));
}

bool CoprocessorBroker::TaskData::CudaSetup(unsigned int streamid, int nblocks, int nthreads, int maxTrackPerThread, GeantPropagator *propagator, const vecgeom::DevicePtr<Geant::cuda::GeantPropagator > &devPropagator)
{
  int maxdepth = propagator->fConfig->fMaxDepth;
  unsigned int maxTrackPerKernel = nblocks * nthreads * maxTrackPerThread;
  fChunkSize = maxTrackPerKernel;
  fDevMaxTracks = 2 * fChunkSize;

  fGeantTaskData = new Geant::GeantTaskData(nthreads, maxdepth, fDevMaxTracks);
  fGeantTaskData->fPropagator = new GeantPropagator(*propagator);
  fGeantTaskData->fTid = streamid; // NOTE: not quite the same ...
  fGeantTaskData->fBmgr = new GeantBasketMgr(propagator, 0, 0, true);
  fPrioritizer = fGeantTaskData->fBmgr;
  fPrioritizer->SetFeederQueue(propagator->fWMgr->FeederQueue());

  fStreamId = streamid;
  GEANT_CUDA_ERROR(cudaStreamCreate(&fStream));

  // prepare random engines on the device
  // fdRandStates.Alloc( nblocks*nthreads );
  // curand_setup_gpu(fdRandStates, time(NULL), nblocks, nthreads);

  unsigned int maxThreads = nblocks * nthreads;

  unsigned long size_needed = GeantTaskData::SizeOfInstance(
      maxThreads, maxdepth,
      (unsigned int)5 // maxTrackPerKernel is to much, maxTrackPerKernel/maxThreads might make more sense (i.e.
                      // maxTrackPerThread) unless we need space for extras/new tracks ...
      );
  fDevTaskWorkspaceSizeOf = size_needed;

  fDevTaskWorkspace.Malloc(maxThreads * size_needed); // Allocate(maxThreads);
  Geant::cuda::MakeInstanceArrayAt(
      fDevTaskWorkspace.GetPtr(), 4096 /* maxThreads */, size_needed, (size_t)maxThreads, maxdepth,
      (int)5 // maxTrackPerKernel is to much, maxTrackPerKernel/maxThreads might make more sense (i.e.
             // maxTrackPerThread) unless we need space for extras/new tracks ...
      , devPropagator.GetPtr()
      );

  // need to allocate enough for one object containing many tracks ...
  fDevTrackInput.Malloc(Geant::GeantTrack_v::SizeOfInstance(fDevMaxTracks, maxdepth));
  Geant::cuda::MakeInstanceAt(fDevTrackInput.GetPtr(), fDevMaxTracks, maxdepth);

  fDevTrackOutput.Malloc(Geant::GeantTrack_v::SizeOfInstance(fDevMaxTracks, maxdepth));
  Geant::cuda::MakeInstanceAt(fDevTrackOutput.GetPtr(), fDevMaxTracks, maxdepth);

  fInputBasket = new GeantBasket(propagator, fDevMaxTracks, fGeantTaskData->fBmgr);
  fInputBasket->SetMixed(true);
  // fInputBasket->SetThreshold(fThreshold.load());
  fOutputBasket = new GeantBasket(propagator, fDevMaxTracks, fGeantTaskData->fBmgr);
  fOutputBasket->SetMixed(true);
  // fOutputBasket->SetThreshold(fThreshold.load());

  return true;
}

void CoprocessorBroker::TaskData::Push(dcqueue<CoprocessorBroker::TaskData*> *q)
{
  // Add this helper to the queue and record it as the
  // default queue.

  if (q) {
    fQueue = q;
    fQueue->push(this);
  }
  else if (fQueue) {
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

CoprocessorBroker::CoprocessorBroker()
    : fdGeometry(0), fNextTaskData(0), fNblocks(0), fNthreads(0), fMaxTrackPerThread(1), fTotalWork(0),
      fNConcurrentStream(1), fIsSelective(false)
//,fdFieldMap(0)
//,fdStates0(0),fdStates1(0),fdRandStates0(0),fdRandStates1(0)
{
  // Default constructor.

  // fIsSelective = false;
  // fTasks.push_back(new GeneralTask());

  fIsSelective = true;
  fTasks.push_back(new GeneralChargedTask(-1));

  /*
  fTasks.push_back(new EnergyElectronTask(6));
  fTasks.push_back(new EnergyElectronTask(4));
  fTasks.push_back(new EnergyElectronTask(2));
  fTasks.push_back(new EnergyElectronTask(0));
  */
}

CoprocessorBroker::~CoprocessorBroker()
{
  GEANT_CUDA_ERROR(cudaDeviceSynchronize());
  for (unsigned int i = 0; i < fTasks.size(); ++i) {
    delete fTasks[i];
  }
  fTasks.clear();
  for (unsigned int i = 0; i < fTaskData.size(); ++i) {
    delete fTaskData[i];
  }
  fTaskData.clear();

  GEANT_CUDA_ERROR(cudaFree(fdGeometry));
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

  // CudaManager::Instance().set_verbose(3);

  if (volume)
    vecgeom::CudaManager::Instance().LoadGeometry(volume);
  else
    vecgeom::CudaManager::Instance().LoadGeometry();
  vecgeom::CudaManager::Instance().Synchronize();

  // CudaManager::Instance().PrintGeometry();
  // std::cout << std::flush;

  return true;
}

void setup(CoprocessorBroker * /* broker */, int /* nphi */ = 4, int /* nz */ = 3, double /* density */ = 8.28)
{

  // 2. Read magnetic field map

  // 3. Create magnetic field on the device

  // 4. Prepare EM physics tables
}

namespace Geant {
namespace cuda {
void CoprocessorBrokerInitConstant();
}
}

/** @brief Create the baskets for each stream */
void CoprocessorBroker::CreateBaskets(GeantPropagator* propagator)
{
  // We must be called after the geometry has been loaded
  // so we can use the proper maximum depth
  // and over-ride the cudaLimitStackSize.

  fDevConfig.Allocate();
  if (cudaGetLastError() != cudaSuccess) {
     printf(" ERROR ALLOC MAP\n");
     return;
  }
  fDevConfig.Construct();
  // It is save to cast here as the cuda version is 'just' missing
  // data member at the end but is otherwise the same.
  cuda::GeantConfig *conf = (cuda::GeantConfig *)propagator->fConfig;
  fDevConfig.ToDevice(conf);

  fDevPropagator.Allocate();
  if (cudaGetLastError() != cudaSuccess) {
    printf(" ERROR ALLOC MAP\n");
    return;
  }
  fDevPropagator.Construct(fNblocks * fNthreads);

   // initialize the stream
  for (unsigned int i = 0; i < GetNstream(); ++i) {
    TaskData *data = new TaskData();
    data->CudaSetup(i, fNblocks, fNthreads, fMaxTrackPerThread, propagator, fDevPropagator);
    data->Push(&fHelpers);
    fTaskData.push_back(data);
  }

  cudaDeviceSetLimit(cudaLimitStackSize, 3 * 4096);

  // Initialize global constants.
  Geant::cuda::CoprocessorBrokerInitConstant();
}

bool CoprocessorBroker::CudaSetup(int nblocks, int nthreads, int maxTrackPerThread)
{
  int deviceCount = 0;
  cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

  if (error_id != cudaSuccess) {
    Printf("Cuda CoprocessorBroker disabled because of issue #%d:%s\n", (int)error_id, cudaGetErrorString(error_id));
    return false;
  }

  fNConcurrentStream = 2; // NOTE: Is there a way to know the number really supported?

  setup(this);

  fNblocks = nblocks;
  fNthreads = nthreads;
  fMaxTrackPerThread = maxTrackPerThread;

  return true;
}

// bool CanHandleTrack(GeantTrack &track)
// {
//    // Currently we can only handle electron, which we pretend are the only
//    // particle to have charge -1.

//    return -1 == track.charge;
// }

unsigned int CoprocessorBroker::TaskData::AddTrack(CoprocessorBroker::Task *task, GeantBasket &basket,
                                                   unsigned int hostIdx)
{

  if (!fInputBasket) {
    GeantBasketMgr *bmgr = basket.GetBasketMgr();
    fInputBasket = bmgr->GetNextBasket(fGeantTaskData);
  }

  GeantTrack_v &input = basket.GetInputTracks();

  //  if (input.fHoles->TestBitNumber(hostIdx)) return 0;
  if (task->Select(input, hostIdx)) {
    fInputBasket->AddTrack(input, hostIdx);
    ++fNStaged;
    return 1;
  }
  return 0;
}

unsigned int CoprocessorBroker::TaskData::TrackToDevice(CoprocessorBroker::Task *task, int tid, GeantBasket &basket,
                                                        unsigned int startIdx)
{
  if (fThreadId == -1) {
    fThreadId = tid;
  }
  else if (fThreadId != tid) {
    Fatal("TrackToDevice", "The stream %p is already assigned to thread %d when requested to do work for thread %d",
          this, fThreadId, tid);
    return 0;
  }

  if (!fInputBasket) {
    GeantBasketMgr *bmgr = basket.GetBasketMgr();
    fInputBasket = bmgr->GetNextBasket(fGeantTaskData);
  }

  // unsigned int start = fNStaged;
  unsigned int count = 0;
  GeantTrack_v &input = basket.GetInputTracks();
  // input.PrintTracks("TrackToDevice");
  // GeantTrack_v &gputracks(fInputBasket->GetInputTracks());
  unsigned int basketSize = input.GetNtracks();

  // Well ... do physics .. just because.
  {
    GeantPropagator *propagator = fGeantTaskData->fPropagator;
    propagator->ProposeStep(basketSize, input, fGeantTaskData);
    // Apply msc for charged tracks
    propagator->ApplyMsc(basketSize, input, fGeantTaskData);
  }

  for (unsigned int hostIdx = startIdx; fNStaged < fChunkSize && hostIdx < basketSize; ++hostIdx) {

    ++count;

    if (input.fHoles->TestBitNumber(hostIdx))
      continue;

    if (task->Select(input, hostIdx)) {

      /* int t = */ fInputBasket->GetInputTracks().AddTrack(input, hostIdx);

      // Prepare the navigation state pointers in the basket
      // gputracks.fPathV[t]->ConvertToGPUPointers();
      // gputracks.fNextpathV[t]->ConvertToGPUPointers();

      // fTrackId[fNStaged] = input.PostponeTrack(hostIdx,fBasket->GetOutputTracks());
      input.MarkRemoved(hostIdx);

      ++fNStaged;
    }
  }

  // Track sub-ranges are not consecutive in memory.
  // // count = min(fChunkSize,basketSize-startIdx);
  // int ntrack = fNStaged - start;
  // GEANT_CUDA_ERROR(cudaMemcpyAsync(fDevInputTrack+start, input+start, ntrack*sizeof(GXTrack),
  //                                   cudaMemcpyHostToDevice, fStream));
  // GEANT_CUDA_ERROR(cudaMemcpyAsync(fDevTrackLogIndex+start, fLogIndex+start, ntrack*sizeof(int),
  //                                   cudaMemcpyHostToDevice, fStream));
  // GEANT_CUDA_ERROR(cudaMemcpyAsync(fDevTrackPhysIndex+start, fPhysIndex+start, ntrack*sizeof(int),
  //                                   cudaMemcpyHostToDevice, fStream));

  return count;
}

unsigned int CoprocessorBroker::TaskData::TrackToHost()
{
  WorkloadManager *mgr = fGeantTaskData->fPropagator->fWMgr;
  GeantScheduler *sch = mgr->GetScheduler();
  condition_locker &sched_locker = mgr->GetSchLocker();
  GeantTrack_v &output = *fGeantTaskData->fTransported;
  // GeantTrack_v &output = *fOutputBasket->GetInputTracks();
  GeantPropagator *propagator = fGeantTaskData->fPropagator;
  auto td = fGeantTaskData;

  if (output.GetNtracks() > output.Capacity())
    Fatal("CoprocessorBroker::TaskData::TrackToHost", "Request size in output track buffer is too large ( %d > %d )",
          output.GetNtracks(), output.Capacity());

  GeantTrack_v &transferTo(output);
  FromDeviceConversion(&transferTo, fDevTrackOutput);

  // Fix the navigation state pointers in the output basket
  for (int t = 0; t < transferTo.fMaxtracks; ++t) {
    // if (output.fHoles->TestBitNumber(t) continue;
    // if (transferTo.fPathV[t]) transferTo.fPathV[t]->ConvertToCPUPointers();
    // if (transferTo.fNextpathV[t]) transferTo.fNextpathV[t]->ConvertToCPUPointers();
  }

  // output.PrintTracks("TrackToHost");

  // Waste a lot of time ... by doing physics ...

  // Post-step actions by continuous processes for all particles. There are no
  // new generated particles at this point.
  {
    // GeantTrack_v &output = *td->fTransported;
    // GeantTrack_v &output = *fOutputBasket->GetInputTracks();
    auto ntotnext = 0;
    Material_t *mat = nullptr;

    if (1) {

      int nphys = 0;
      int nextra_at_rest = 0;
      // count phyics steps here
      auto nout = output.GetNtracks();
      for (auto itr = 0; itr < nout; ++itr) {
        if (output.fStatusV[itr] == kPhysics)
          nphys++;
        // Update track time for all output tracks (this should vectorize)
        output.fTimeV[itr] += output.TimeStep(itr, output.fStepV[itr]);
      }

      propagator->Process()->Eloss(mat, output.GetNtracks(), output, nextra_at_rest, td);
      //         if (nextra_at_rest) Printf("Extra particles: %d", nextra_at_rest);
      // Now we may also have particles killed by energy threshold
      // Do post-step actions on remaining particles
      // to do: group particles per process

      if (propagator->fConfig->fUsePhysics) {
        // Discrete processes only
        nphys = output.SortByLimitingDiscreteProcess(); // only those that kPhysics and not continous limit
        if (nphys) {
          // propagator->fNphysSteps += nphys;  dont do it here because dont count those killed in eloss
          // Do post step actions for particles suffering a given process.
          // Surviving particles are added to the output array

          // first: sample target and type of interaction for each primary tracks
          propagator->Process()->PostStepTypeOfIntrActSampling(mat, nphys, output, td);

//
// TODO: vectorized final state sampling can be inserted here through
//       a proper interface code that will:
//         1. select the necessary primary tracks based on the alrady
//            sampled interaction type
//         2. take all the member of the selected primary tracks that
//            necessary for sampling the final states
//         3. call the appropriate vector physics code that will
//            perform the physics interaction itself
//         4. update properties of the selected primary tracks,
//            insert secondary tracks into the track vector and set
//            'ntotnext' to the value of the number of secondary tracks
//            inserted to the track vector
//
#if USE_VECPHYS == 1
          propagator->fVectorPhysicsProcess->PostStepFinalStateSampling(mat, nphys, output, ntotnext, td);
#endif
          // second: sample final states (based on the inf. regarding sampled
          //         target and type of interaction above), insert them into
          //         the track vector, update primary tracks;
          propagator->Process()->PostStepFinalStateSampling(mat, nphys, output, ntotnext, td);

          if (0 /*ntotnext*/) {
            Geant::Printf("============= Basket.\n");
            output.PrintTracks();
          }
        }
      }
    }
    if (propagator->fStdApplication)
      propagator->fStdApplication->StepManager(output.GetNtracks(), output, td);
    propagator->fApplication->StepManager(output.GetNtracks(), output, td);

    // Update geometry path for crossing tracks
    ntotnext = output.GetNtracks();

    for (auto itr = 0; itr < ntotnext; ++itr) {
      output.fNstepsV[itr]++;
      if (output.fStatusV[itr] == kBoundary)
        *output.fPathV[itr] = *output.fNextpathV[itr];
    }
  }

  int ntot = 0;
  int nnew = 0;
  int nkilled = 0;
  // Printf("(%d - GPU) ================= Returning from Stream %d accumulated=%d outputNtracks=%d holes#=%lu "
  //        "basketHoles#=%lu ",
  //        fThreadId, fStreamId, fNStaged, output.GetNtracks(), transferTo.fHoles->CountBits(),
  //        output.fHoles->CountBits());
  /* int ninjected = */ sch->AddTracks(output, ntot, nnew, nkilled, fGeantTaskData);
  (void)ntot;
  (void)nnew;
  (void)nkilled;
  // mgr->TransportedQueue()->push(fBasket);
  sched_locker.StartOne();
  // fOutputBasket->Recycle(td);
  // fOutputBasket = 0;
  // fOutputBasket->Clear();
  mgr->CheckFeederAndExit();
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

  if (fPrevNStaged == fCurrent->fNStaged) {
    ++(fIdles);
  }
  else {
    fIdles = 0;
  }
  ++(fCycles);
  fPrevNStaged = fCurrent->fNStaged;

  unsigned int idle = fIdles;
  unsigned int cycle = fCycles;
  if (fCurrent->fNStaged                                     // There is something
      && idle > (ntasks - 1)                                 // We are beyond the expected/normal number of idles cycles
      && cycle > (fCurrent->fChunkSize / fCurrent->fNStaged) // To fill we need at least that many cycles
      && 2 * idle > cycle                                    // Our input rate has drop in half
      ) {
    // Printf("(%d - GPU) ================= Launching idle Task %s Stream %d Idle=%d cycle=%d accumulated=%d", threadid,
    // (*task)->Name(), stream->fStreamId, idle, cycle, stream->fNStaged);
    // if we did not make any progress in a while, assume there is no 'interesting' track left and schedule the kernel.
    return true;
  }

  // Continue to wait for more data ...
  return false;
}

// typedef void(CUDART_CB * cudaStreamCallback_t)(cudaStream_t stream, cudaError_t status, void
//                                                *userData)

void StreamReset(cudaStream_t /* stream */, cudaError_t /* status */, void *userData)
{
  CoprocessorBroker::TaskData *helper = (CoprocessorBroker::TaskData *)userData;
  helper->Reset();
}

void TrackToHost(cudaStream_t /* stream */, cudaError_t /* status */, void *userData)
{
  CoprocessorBroker::TaskData *helper = (CoprocessorBroker::TaskData *)userData;
  helper->TrackToHost();
}

void ClearTrack_v(cudaStream_t /* stream */, cudaError_t /* status */, void *userData)
{
  Geant::GeantTrack_v *tracks = (Geant::GeantTrack_v *)userData;
  tracks->Clear();
}

void FromDeviceTrackConversion(cudaStream_t /* stream */, cudaError_t /* status */, void *userData)
{
  CoprocessorBroker::TaskData *helper = (CoprocessorBroker::TaskData *)userData;
  FromDeviceConversion(helper->fGeantTaskData->fTransported, helper->fDevTrackOutput);
}

CoprocessorBroker::Stream CoprocessorBroker::GetNextStream()
{
  // Return the current stream (one that we can still add input to)
  // or return new one.

  if (!fNextTaskData) {
    fHelpers.wait_and_pop(fNextTaskData);
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

  int outstanding = 0;
  for (auto td : fTaskData)
    outstanding += td->fNStaged;
  // Printf("(%d - GPU) == Starting kernel for task %s using stream %d with %d tracks (%d total oustanding tracks)\n",
  //        stream->fThreadId, task->Name(), stream->fStreamId, stream->fNStaged, outstanding);

  if ((unsigned int)stream->fInputBasket->GetInputTracks().GetNtracks() == stream->fDevMaxTracks) {
    Fatal("CoprocessorBroker::launchTask",
          "Number of tracks allocated in input basket (%d) must be the same as on device (%d)\n",
          stream->fInputBasket->GetInputTracks().GetNtracks(), stream->fDevMaxTracks);
    return stream;
  }
  ToDevice(stream->fDevTrackInput, &(stream->fInputBasket->GetInputTracks()), *stream);
  GEANT_CUDA_ERROR(cudaStreamAddCallback(stream->fStream, ClearTrack_v, &(stream->fInputBasket->GetInputTracks()), 0));

  fTotalWork += stream->fNStaged;
  int result = task->fKernel(stream->fDevTaskWorkspace, stream->fDevTaskWorkspaceSizeOf, stream->fNStaged,
                             stream->fDevTrackInput, stream->fDevTrackOutput,

                             fNblocks, fNthreads, *stream);

  if (!result)
    return stream;

  // Bring back the number of secondaries created.
  // int stackSize;
  // stream->fDevSecondaries.fDevStackSize.FromDevice(&stackSize, *stream);
  // stream->fDevSecondaries.fTrack.FromDevice( stream->fSecondaries, stackSize, *stream);

  // Bring back the modified tracks.
  if ((unsigned int)stream->fGeantTaskData->fTransported->GetNtracks() == stream->fDevMaxTracks) {
    Fatal("CoprocessorBroker::launchTask",
          "Number of track allocated in output basket (%d) must be the same as on device (%d)\n",
          stream->fGeantTaskData->fTransported->GetNtracks(), stream->fDevMaxTracks);
    return stream;
  }
  FromDevice(stream->fGeantTaskData->fTransported, stream->fDevTrackOutput, *stream);

  Clear_gpu(stream->fDevTrackOutput, 1, 1, *stream);
  GEANT_CUDA_ERROR(cudaStreamAddCallback(stream->fStream, TrackToHost, stream, 0));
  GEANT_CUDA_ERROR(cudaStreamAddCallback(stream->fStream, StreamReset, stream, 0));

  if (wait) {
    // Use this when you need to insure the printf are actually printed.
    GEANT_CUDA_ERROR(cudaStreamSynchronize(*stream));
  }
  // cudaDeviceSynchronize

  return stream;
}

/** @brief If the coprocessor has outstanding work, return it */
GeantBasket *CoprocessorBroker::GetBasketForTransport(Geant::GeantTaskData &maintd)
{
  GeantBasketMgr *prioritizer = maintd.fBmgr;
  if (prioritizer && prioritizer->HasTracks()) {
    return prioritizer->GetBasketForTransport(&maintd);
  }
  for (auto task : fTasks) {
    if (task->fCurrent) {
      prioritizer = task->fCurrent->fPrioritizer;
      if (prioritizer && prioritizer->HasTracks()) {
        return prioritizer->GetBasketForTransport(task->fCurrent->fGeantTaskData);
      }
    }
  }
  return nullptr;
}

CoprocessorBroker::Stream CoprocessorBroker::launchTask(bool wait /* = false */)
{
  Stream stream = 0;
  TaskColl_t::iterator task = fTasks.begin();
  TaskColl_t::iterator end = fTasks.end();
  while (task != end) {

    if ((*task)->fCurrent != 0 && (*task)->fCurrent->fNStaged != 0) {

      stream = launchTask((*task), wait);
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
  // unsigned int nTracks = 1;
  unsigned int trackUsed = 0;
  // unsigned int trackStart = 0;

  TaskColl_t::iterator task = fTasks.begin();
  TaskColl_t::iterator end = fTasks.end();
  while (task != end) {

    if ((*task)->fCurrent == 0) {
      if (fNextTaskData) {
        (*task)->fCurrent = fNextTaskData;
      }
      else {
        // If we do not yet have a TaskData, wait for one.
        // Consequence: if we are very busy we hang on to this
        // track until one of the task finishes.
        // It is likely better to returned that we could not
        // handle the track.
        (*task)->fCurrent = GetNextStream();
      }
      fNextTaskData = 0;
      if (!(*task)->fCurrent)
        break;
    }

    TaskData *stream = (*task)->fCurrent;

    unsigned int before = stream->fNStaged;
    /* unsigned int count = */ stream->AddTrack(*task, input, trkid);
    trackUsed += stream->fNStaged - before;
    // unsigned int rejected = nTracks-trackStart - (stream->fNStaged-before);
    if (!force && !(*task)->IsReadyForLaunch(fTasks.size())) {
      // Continue to wait for more data ...
      break;
    }

    launchTask(*task);

    ++task;
  }

  // Missing compared to runTask si the scheduling of the most loaded task
  // if nothing new at all in a while.

  if (trackUsed)
    return true;
  else
    return false;
}

void CoprocessorBroker::runTask(GeantTaskData &td, GeantBasket &basket)
// unsigned int nTracks, int volumeIndex, GeantTrack **tracks, int *trackin)
{
  bool force = false;

  unsigned int nTracks = basket.GetNinput();

  unsigned int trackUsed = 0;
  TaskColl_t::iterator task = fTasks.begin();
  TaskColl_t::iterator end = fTasks.end();
  while (task != end) {
    unsigned int trackLeft = nTracks;
    unsigned int trackStart = 0;
    while (trackLeft) {
      if ((*task)->fCurrent == 0) {
        if (fNextTaskData) {
          (*task)->fCurrent = fNextTaskData;
        }
        else {
          // If we do not yet have a TaskData, wait for one.
          // Consequence: if we are very busy we hang on to this
          // basket until one of the task finishes.
          (*task)->fCurrent = GetNextStream();
        }
        fNextTaskData = 0;
        if (!(*task)->fCurrent)
          break;
      }
      TaskData *stream = (*task)->fCurrent;

      unsigned int before = stream->fNStaged;
      unsigned int count = stream->TrackToDevice(*task, td.fTid, basket, trackStart);
      trackUsed += stream->fNStaged - before;
      // unsigned int rejected = nTracks-trackStart - (stream->fNStaged-before);

      // if (((*task)->fCycles % 10000) == 1)
      // Printf("(%d - GPU) ================= Task %s Stream %d Tracks: %d seen, %d skipped, %d accumulated, %d idles,
      // %d cycles, %d used this iteration", threadid, (*task)->Name(), stream->fStreamId, count, rejected,
      // stream->fNStaged, (*task)->fIdles, (*task)->fCycles, stream->fNStaged-before);

      if (!force && !(*task)->IsReadyForLaunch(fTasks.size())) {
        // Continue to wait for more data ...
        break;
      }

      launchTask(*task);
      trackLeft -= count;
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
    while (task != end) {
      if ((*task)->fCurrent && (*task)->fCurrent->fNStaged) {
        if (heavy == 0 || (*task)->fCurrent->fNStaged > heavy->fCurrent->fNStaged) {
          heavy = *task;
        }
      }
      ++task;
    }
    if (heavy) {
      // Printf("(%d - GPU) ================= Launching heavy Task %s Stream %d Idle=%d cycle=%d accumulated=%d",
      // threadid, heavy->Name(), heavy->fCurrent->fStreamId, heavy->fIdles, heavy->fCycles, heavy->fCurrent->fNStaged);
      launchTask(heavy);
    }
  }

  if (trackUsed != nTracks) {
    GeantTrack_v &input = basket.GetInputTracks();
    GeantTrack_v &output = *td.fTransported;
    for (unsigned int hostIdx = 0; hostIdx < nTracks; ++hostIdx) {

      if (input.fHoles->TestBitNumber(hostIdx))
        continue;

      input.PostponeTrack(hostIdx, output);
    }
  }
  if (!basket.GetInputTracks().IsCompact())
    basket.GetInputTracks().Compact();

  size_t poolSize = td.GetBasketPoolSize();
  if (poolSize >= fTasks.size() + 20) { // MaybeCleanupBaskets requires at least 10 per thread.
    // Redistribute the baskets.
    // This is needed because the task's TaskData only uses (new) Baskets
    // while the TaskData td is the only always recycling them.
    size_t share = (poolSize - 20) / fTasks.size();
    TaskColl_t::iterator task = fTasks.begin();
    TaskColl_t::iterator end = fTasks.end();
    while (task != end) {
      if ((*task)->fCurrent) {
        for (size_t i = share; i != 0; --i) {
          auto b = td.GetNextBasket();
          if (b)
            (*task)->fCurrent->fGeantTaskData->RecycleBasket(b);
        }
      }
      ++task;
    }
  }
}

void CoprocessorBroker::waitForTasks()
{
  // Make sure all the tasks are finished.

  for (unsigned int i = 0; i < fTaskData.size(); ++i) {
    if (fTaskData[i]->fNStaged) {
      GEANT_CUDA_ERROR(cudaStreamSynchronize(*fTaskData[i]));
    }
  }
}
