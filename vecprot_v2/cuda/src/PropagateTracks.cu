#include "GeantTaskData.h"
#include "GeantTrack.h"
#include "backend/cuda/Interface.h"

#ifdef GEANT_HAVE_PHYSICS_IN_KERNEL
#include "TEFstate.h"
#include "TEXsec.h"
#include "TFinState.h"
#include "TMXsec.h"
#include "TPDecay.h"
#include "TPFstate.h"
#include "TPXsec.h"
#include "TPartIndex.h"
#include "TPrimaryGenerator.h"
#include "TTabPhysMgr.h"
#include "TTabPhysProcess.h"
#endif

/*
  C++ function gets in addition

  int nBlocks, int nThreads,
  cudaStream_t stream
 */

__global__ void PropagateGeantTrack(Geant::GeantTaskData *workSpace, size_t workspaceSizeOf, size_t ntracks,
                                    Geant::GeantTrack_v *input, Geant::GeantTrack_v *output)
{

  /* All at once would be:
  input->ComputeTransportLength(ntracks);
  input->PropagateTracks(*output,tid);
  */
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  Geant::GeantTaskData *td = (Geant::GeantTaskData *)(((char *)workSpace) + workspaceSizeOf * tid);
  td->fTransported->Clear();

#if 0
   // Test whether we use up too much memory already
   char *ptr1 = (char*)vecgeom::cuda::NavigationState::MakeInstance(3);
   char *ptr2 = new char[48];
   if (ptr1==0 || ptr2==0) {
         printf("DEBUG-GPU-4: tid=%d ptr1=%p ptr2=%p\n",tid,ptr1,ptr2);
         return;
   } else { delete [] ptr1; delete [] ptr2; }
#endif

  unsigned int itr = tid;
  while (itr < ntracks) {
    // input->ComputeTransportLengthSingle(itr,td);
    input->PropagateSingleTrack(itr, td, 0);

    // PropagateSingleTrack marks the track slot as a 'hole' and
    // Compact (called byPropagateTracks) would have moved the
    // track to td->fTransported, so we would have to move it now to output ...
    // output->AddTrackSyncAt(itr,*td->fTransported,0);

    // Move hole into output.
    output->AddTrackSyncAt(itr, *input, itr);

    itr += blockDim.x * gridDim.x;
  }
  if (output->fSelected->GetNbits() != 2 * 4096)
    printf("output bitset %ld\n", output->fSelected->GetNbits());
}

int PropagateGeantTrack_gpu(vecgeom::cxx::DevicePtr<Geant::cuda::GeantTaskData> &workSpace, size_t workspaceSizeOf,
                            size_t ntracks, vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrack_v> &input,
                            vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrack_v> &output,

                            int nBlocks, int nThreads, cudaStream_t stream)
{
  int threadsPerBlock = nThreads;
  int blocksPerGrid = nBlocks;

  // fprintf(stderr,"DEBUG-GPU-0: About to schedule the PropagateGeantTrack kernel\n");
  PropagateGeantTrack<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(workSpace, workspaceSizeOf, ntracks, input,
                                                                     output);

  return 1;
}

#if 0
// A more complete version (copied from an older revision of WorkloadManager::TransportTracks
__global__
void transport_kernel()
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  Geant::GeantTaskData *td = (Geant::GeantTaskData * )(((char *) workSpace) + workspaceSizeOf * tid);
  td->fTransported->Clear();

  GeantPropagator *propagator = td->fPropagator;

  unsigned int itr = tid;
  while (itr < ntracks) {
    // Select the discrete physics process for all particles in the basket
    if (propagator->fUsePhysics) propagator->ProposeStep(ntotransport, input, tid);
    // Apply msc for charged tracks
    propagator->ApplyMsc(ntotransport, input, tid);

    if (basket->IsMixed())
       ncross += input.PropagateTracksScalar(output,tid);
    else
       ncross += input.PropagateTracks(output,tid);

    // All tracks are now in the output track vector. Possible statuses:
    // kCrossing - particles crossing boundaries
    // kPhysics - particles reaching the point where the discrete physics process
    //            will happen.
    // kExitingSetup - particles exiting the geometry
    // kKilled - particles that could not advance in geometry after several tries

    // Post-step actions by continuous processes for all particles. There are no
    // new generated particles at this point.
    if (propagator->fUsePhysics) {
         nextra_at_rest = 0;
         gPropagator->Process()->Eloss(td->fVolume->GetMaterial(), output.GetNtracks(), output, nextra_at_rest, tid);
     }

     // Now we may also have particles killed by energy threshold
     // Do post-step actions on remaining particles
     // Loop all processes to group particles per process

     if (propagator->fUsePhysics) {
       // Discrete processes only
        int nphys = output.SortByStatus(kPhysics);
        if (nphys) {
           // Do post step actions for particles suffering a given process.
           // Surviving particles are added to the output array
           propagator->Process()->PostStep(td->fVolume->GetMaterial(), nphys, output, ntotnext, tid);
       }
    }
}
#endif
