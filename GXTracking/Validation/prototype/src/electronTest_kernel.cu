#include "electronTest_kernel.h"

#include <cuda.h>
#include <curand_kernel.h>

//Common
#include "GXTrack.h"
#include "GPPhysicsTable.h"
#include "GPPhysics2DVector.h"

//Material
#include "GPElement.h"
#include "GPMaterial.h"
#include "GPVParticleChange.h"

//EMPhysics
#include "GPeBremsstrahlung.h"
#include "GPSeltzerBergerRelModel.h"
#include "GPeIonisation.h"
#include "GPMollerBhabhaModel.h"
#include "GPIonisParamMat.h"
#include "GPUniversalFluctuation.h"
#include "GPeMultipleScattering.h"
#include "GPUrbanMscModel95.h"

#include "stdio.h"

#ifndef __CUDA_ARCH__
//#warning electronTest_kernel.cu: CUDA_ARCH undefined: Loading Root stuff
#include "GPHistoManager.hh"
#include "TH1F.h"
#else
//#warning electronTest_kernel.cu: CUDA_ARCH defined: Skipping Root stuff
#endif // __CUDA_ARCH__

FQUALIFIER
G4double GetKineticEnergy(GXTrack *track) 
{
  G4double p = sqrt(track->px*track->px + track->py*track->py 
		  + track->pz*track->pz);
  G4double mass = electron_mass_c2*track->q*track->q;
  //G4double ekin = p*p/(sqrt(p*p + mass*mass) + mass);
  G4double ekin = sqrt(p*p + mass*mass) - mass;  // simpler
  return ekin;
}

FQUALIFIER
void DefineMaterial(GPMaterial *aMat) 
{
  //PbWO4 aka CMS crystal
  GPMaterial_Constructor_ByElement(aMat,8.28*g/cm3);

  GPElement ele_Pb;
  GPElement ele_W;
  GPElement ele_O;

  GPElement_Constructor(&ele_Pb,82,207.2*g/mole);
  GPElement_Constructor(&ele_W,74,183.84*g/mole);
  GPElement_Constructor(&ele_O,8,15.9994*g/mole);

  // GPMaterial_AddElement(aMat,ele_Pb,0.45532661);
  // GPMaterial_AddElement(aMat,ele_W, 0.40403397);
  // GPMaterial_AddElement(aMat,ele_O4,0.14063942);

  GPMaterial_AddElement(aMat,ele_Pb,0.4553445);
  GPMaterial_AddElement(aMat,ele_W, 0.4040084);
  GPMaterial_AddElement(aMat,ele_O, 0.1406470);

  //update radiation length
  GPMaterial_ComputeDerivedQuantities(aMat);
}

FQUALIFIER
void StackSecondaries(GPVParticleChange *particleChange,
                      GXTrack *secTracks,
                      G4int *stackSize,
                      G4int *offset) 
{
  //store secondaries in a global stack (fixed size memory)
  G4int nsec = particleChange->GetNumberOfSecondaries();

  for(int isec = 0 ; isec < nsec ; ++isec) {
    //offset is a global counter for the last array position of secondaries 
    GXTrack* aSecTrack = particleChange->GetSecondary(isec);
#ifndef __CUDA_ARCH__
    secTracks[*stackSize] = *aSecTrack;
    ++(*stackSize);
#else
    *offset = atomicAdd(stackSize,1);
    secTracks[*offset] = *aSecTrack;
#endif
  }
}

FQUALIFIER
void Print(GXTrack *track, GPVParticleChange *particleChange)
{
  G4int nsec = particleChange->GetNumberOfSecondaries();
#ifndef __CUDA_ARCH__
  printf("CPU (step,ekin,nsec) = (%f,%f,%d)\n",track->s,track->E,nsec);
#else
  printf("GPU (step,ekin,nsec,partChange) = (%f, %f, %d, %e)\n",track->s,track->E,nsec,particleChange->GetSecondary(0)->E);
#endif
}

//-----------------------------------------------------------------------------
// Standard Electron Processes
//-----------------------------------------------------------------------------

GLOBALFUNC
void electron_kernel(curandState* devStates,
		     GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* eBrem_table, 
		     GPPhysicsTable* eIoni_table, 
		     GPPhysicsTable* eIoni_range, 
		     GPPhysicsTable* eIoni_dedx, 
		     GPPhysicsTable* eIoni_invr, 
		     GPPhysicsTable* msc_table,
		     GPPhysics2DVector* sbData,
		     GXTrack *secTracks, G4int *stackSize, G4int *offset,
		     G4int isStack, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GPMaterial aMat;
  DefineMaterial(&aMat);

  GPForceCondition condition;
  GPGPILSelection  selection;

  G4double step  = 0;
  G4double step1  = 0;
  G4double step_brem = 0;
  G4double step_ioni = 0;
  G4double step_msc  = 0;
  G4double ekin  = 0;
  G4double proposedStep = 0;

  //EM processes and models
  GPSeltzerBergerRelModel aModel(devStates,tid,sbData);
  GPeBremsstrahlung brem(devStates,tid, eBrem_table);
  brem.InitialiseProcess(&aModel);

  GPMollerBhabhaModel bModel(devStates,tid);
  GPeIonisation ioni(devStates,tid, eIoni_table, eIoni_range,
                             eIoni_dedx, eIoni_invr);
  ioni.InitialiseProcess(&bModel);

  GPUrbanMscModel95 cModel(devStates,tid);
  GPeMultipleScattering msc(devStates,tid,msc_table);
  msc.InitialiseProcess(&cModel);

  GPVParticleChange aPartChange;

  while (tid < nTrackSize) {

    ekin = GetKineticEnergy(&track[tid]);
    step = track[tid].s;

    brem.StartTracking();
    step_brem = brem.PostStepGetPhysicalInteractionLength(&track[tid],step,&condition);

    ioni.StartTracking();
    step_ioni = ioni.PostStepGetPhysicalInteractionLength(&track[tid],step,&condition);

    msc.StartTracking();
    step_msc = msc.PostStepGetPhysicalInteractionLength(&track[tid],&condition);

    //physics model defining the current step
    unsigned int model = 0;

    if(step_brem > step_ioni) model = 1;   
    if(step_ioni > step_msc) model = 2;  

    switch (model) {
    case 0 : 
      proposedStep = step_brem;
      step1 = brem.AlongStepGetPhysicalInteractionLength(&selection);
      if(step1 < proposedStep) proposedStep = step1;
      aPartChange = brem.AlongStepDoIt(&track[tid],&aMat,proposedStep);
      aPartChange = brem.PostStepDoIt(&track[tid],&aMat);
      break;
    case 1 :
      proposedStep = step_ioni;
      step1 = ioni.AlongStepGetPhysicalInteractionLength(&selection);
      if(step1 < proposedStep) proposedStep = step1;
      aPartChange = ioni.AlongStepDoIt(&track[tid],&aMat,proposedStep);
      aPartChange = ioni.PostStepDoIt(&track[tid],&aMat);
      break;
    case 2 :
      proposedStep = step_msc;                         
      if(step1 < proposedStep) proposedStep = step1;
      step1 = msc.AlongStepGetPhysicalInteractionLength(&aMat,ekin,step,
                                                        &selection);
      aPartChange = msc.AlongStepDoIt(&aMat,&track[tid]);
      aPartChange = msc.PostStepDoIt(&track[tid]);

      break;
    default :
      break;
    }

    track[tid].s = proposedStep;
    track[tid].E = aPartChange.GetProposedKineticEnergy();

    if(isStack==0) StackSecondaries(&aPartChange,secTracks,stackSize,offset);
    //    if(runType==1) Print(&track[tid], &aPartChange);

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void electron_gpu(curandState* devStates,
		  GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* eBrem_table, 
		  GPPhysicsTable* eIoni_table, 
		  GPPhysicsTable* eIoni_range, 
		  GPPhysicsTable* eIoni_dedx, 
		  GPPhysicsTable* eIoni_invr, 
		  GPPhysicsTable* msc_table,
		  GPPhysics2DVector* sbData,
		  GXTrack *secTracks, G4int *stackSize, G4int *offset,
		  G4int isStack, G4int runType, 
		  int blocksPerGrid, int threadsPerBlock) 
{
  electron_kernel<<< blocksPerGrid, threadsPerBlock >>>(devStates,track,
                  nTrackSize,eBrem_table,eIoni_table,
                  eIoni_range,eIoni_dedx,eIoni_invr, msc_table,
		  sbData,secTracks,stackSize,offset,isStack,runType);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void electron_cpu(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* eBrem_table, 
		  GPPhysicsTable* eIoni_table, 
		  GPPhysicsTable* eIoni_range, 
		  GPPhysicsTable* eIoni_dedx, 
		  GPPhysicsTable* eIoni_invr, 
		  GPPhysicsTable* msc_table,
		  GPPhysics2DVector* sbData,
		  GXTrack *secTracks, G4int *stackSize,
		  G4int isStack, G4int runType)
{
  GPMaterial aMat;
  DefineMaterial(&aMat);

  GPForceCondition condition;
  GPGPILSelection  selection;

  G4double step  = 0;
  G4double step1  = 0;
  G4double step_brem = 0;
  G4double step_ioni = 0;
  G4double step_msc  = 0;
  G4double ekin  = 0;
  G4double proposedStep = 0;

  //EM processes and models
  GPSeltzerBergerRelModel aModel(0,-1,sbData);
  GPeBremsstrahlung brem(0,-1, eBrem_table);
  brem.InitialiseProcess(&aModel);

  GPMollerBhabhaModel bModel(0,-1);
  GPeIonisation ioni(0,-1, eIoni_table, eIoni_range,
                             eIoni_dedx, eIoni_invr);
  ioni.InitialiseProcess(&bModel);

  GPUrbanMscModel95 cModel(0,-1);
  GPeMultipleScattering msc(0,-1,msc_table);
  msc.InitialiseProcess(&cModel);

  GPVParticleChange aPartChange;

  for (int tid = 0; tid < nTrackSize; tid++) {

    ekin = GetKineticEnergy(&track[tid]);
    step = track[tid].s;

    brem.StartTracking();
    step_brem = brem.PostStepGetPhysicalInteractionLength(&track[tid],step,
                                                          &condition);
    ioni.StartTracking();
    step_ioni = ioni.PostStepGetPhysicalInteractionLength(&track[tid],step,
                                                          &condition);
    msc.StartTracking();
    step_msc = msc.PostStepGetPhysicalInteractionLength(&track[tid],&condition);

    //physics model defining the current step
    unsigned int model = 0;

    if(step_brem > step_ioni) model = 1;   
    if(step_ioni > step_msc) model = 2;  

    switch (model) {
    case 0 : 
      proposedStep = step_brem;
      step1 = brem.AlongStepGetPhysicalInteractionLength(&selection);
      if(step1 < proposedStep) proposedStep = step1;
      aPartChange = brem.AlongStepDoIt(&track[tid],&aMat,proposedStep);
      aPartChange = brem.PostStepDoIt(&track[tid],&aMat);
      break;
    case 1 :
      proposedStep = step_ioni;
      step1 = ioni.AlongStepGetPhysicalInteractionLength(&selection);
      if(step1 < proposedStep) proposedStep = step1;
      aPartChange = ioni.AlongStepDoIt(&track[tid],&aMat,proposedStep);
      aPartChange = ioni.PostStepDoIt(&track[tid],&aMat);
      break;
    case 2 :
      proposedStep = step_msc;                         
      step1 = msc.AlongStepGetPhysicalInteractionLength(&aMat,ekin,step,
                                                        &selection);
      if(step1 < proposedStep) proposedStep = step1;
      aPartChange = msc.AlongStepDoIt(&aMat,&track[tid]);
      aPartChange = msc.PostStepDoIt(&track[tid]);
      break;
    default :
      break;
    }

    track[tid].s = proposedStep;
    track[tid].E = aPartChange.GetProposedKineticEnergy();

    if(isStack==0) StackSecondaries(&aPartChange,secTracks,stackSize,0);
    //    if(runType==1) Print(&track[tid], &aPartChange);
  }
}


//-----------------------------------------------------------------------------
// GPeBremsstrahlungProcess + Fixed size memory for secondaries
//-----------------------------------------------------------------------------

GLOBALFUNC
void brem_kernel(curandState* devStates, 
		 GXTrack *track, size_t nTracks,
		 GPPhysicsTable* eBrem_table, 
		 GPPhysicsTable* eIoni_table, 
		 GPPhysicsTable* msc_table,
		 GPPhysics2DVector* sbData,
		 GXTrack *secTracks, G4int *stackSize, G4int *offset, 
		 G4int isStack, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GPMaterial aMat;
  DefineMaterial(&aMat);

  //EM models and processes
  GPSeltzerBergerRelModel aModel(devStates,tid,sbData);
  GPeBremsstrahlung brem(devStates,tid,eBrem_table);
  brem.InitialiseProcess(&aModel);

  //local varialbles
  G4double step  = 0;
  G4double step1 = 0;
  G4double step2 = 0;
  GPForceCondition condition;
  GPGPILSelection  selection;
  GPVParticleChange aPartChange;

  while (tid < nTracks) {

    step = track[tid].s;
    brem.StartTracking();

    step1 = brem.PostStepGetPhysicalInteractionLength(&track[tid],step,&condition);
    step2 = brem.AlongStepGetPhysicalInteractionLength(&selection);
    track[tid].s = (step1 < step2) ? step1 : step2;

    aPartChange = brem.AlongStepDoIt(&track[tid],&aMat,track[tid].s);
    aPartChange = brem.PostStepDoIt(&track[tid],&aMat);
    track[tid].E = aPartChange.GetProposedKineticEnergy();

    if(isStack==0) StackSecondaries(&aPartChange,secTracks,stackSize,offset);
    //if(runType==1) Print(&track[tid], &aPartChange);

    tid += blockDim.x * gridDim.x;
  }

  //ensure that all threads complete before copying
  __syncthreads();

}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void brem_gpu(curandState* devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* eBrem_table, 
	      GPPhysicsTable* eIoni_table, 
	      GPPhysicsTable* msc_table,
	      GPPhysics2DVector* sbData,
	      GXTrack *secTrack, G4int *stackSize, G4int *offset, 
	      G4int isStack, G4int runType,
	      int blocksPerGrid, int threadsPerBlock) 
{
  brem_kernel<<< blocksPerGrid, threadsPerBlock >>>(devStates,track,
		       nTrackSize,eBrem_table,eIoni_table,msc_table,
		       sbData,secTrack,stackSize,offset,isStack,runType);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void brem_cpu(GXTrack *track, size_t nTracks,
	      GPPhysicsTable* eBrem_table, 
	      GPPhysicsTable* eIoni_table, 
	      GPPhysicsTable* msc_table,
	      GPPhysics2DVector* sbData,
	      GXTrack *secTracks, G4int *stackSize,
	      G4int isStack, G4int runType)
{
  GPMaterial aMat;
  DefineMaterial(&aMat);


  //EM models and processes
  GPSeltzerBergerRelModel aModel(0,-1,sbData);
  GPeBremsstrahlung brem(0,-1, eBrem_table);
  brem.InitialiseProcess(&aModel);

  //local variables
  G4double step  = 0;
  G4double step1 = 0;
  G4double step2 = 0;
  GPForceCondition condition;
  GPGPILSelection  selection;
  GPVParticleChange aPartChange;

#ifndef __CUDA_ARCH__
  G4double inppx,inppy,inppz,outpx,outpy,outpz,secpx,secpy,secpz;
  GPHistoManager& hmgr = GPHistoManager::getInstance();
#endif

  for (int tid = 0; tid < nTracks; ++tid) {

    step = track[tid].s;
    brem.StartTracking();
    G4double mass = electron_mass_c2;

#ifndef __CUDA_ARCH__
    inppx = track[tid].px;
    inppy = track[tid].py;
    inppz = track[tid].pz;
    G4double inpMomentum = sqrt(inppx*inppx+inppy*inppy+inppz*inppz);
    G4double inpTotEnergy = sqrt(inpMomentum*inpMomentum+mass*mass);
#ifdef GPUDEBUG
    printf("\nkernel - Before brem: tid=%i Etot=%f Ekin=%f  P=(%f;%f;%f)\n",tid,
	   (track[tid].E+mass)/MeV, track[tid].E/MeV, inppx/MeV, inppy/MeV, inppz/MeV);
#endif // GPUDEBUG
#endif // __CUDA_ARCH__

    step1 = brem.PostStepGetPhysicalInteractionLength(&track[tid],step,&condition);
    step2 = brem.AlongStepGetPhysicalInteractionLength(&selection);
    track[tid].s = (step1 < step2) ? step1 : step2;

    aPartChange = brem.AlongStepDoIt(&track[tid],&aMat,track[tid].s); // fParticleChange.ParticleChangeForLoss_InitializeForAlongStep(track);
    aPartChange = brem.PostStepDoIt(&track[tid],&aMat);       // fParticleChange.ParticleChangeForLoss_InitializeForPostStep(track);
    // then currentModel->SampleSecondaries(), cM->GetSec(),
    // followed by fPC->AddSec(), fPC->GetProposedKE() (and fPC->ProposeTrackStatus(stopKill) if needed)

    G4double outKinEnergy = aPartChange.GetProposedKineticEnergy();
    G4double outTotEnergy = outKinEnergy + mass;
    G4double outMomentum = sqrt(outTotEnergy*outTotEnergy-mass*mass);
    GPThreeVector outMomVec = GPThreeVector_mult( aPartChange.GetProposedMomentumDirection(), outMomentum );
    track[tid].E = outKinEnergy;
    track[tid].px = GPThreeVector_x(outMomVec);
    track[tid].py = GPThreeVector_y(outMomVec);
    track[tid].pz = GPThreeVector_z(outMomVec);

#ifndef __CUDA_ARCH__
    //== make some plots: electron
    outpx = track[tid].px;
    outpy = track[tid].py;
    outpz = track[tid].pz;
#ifdef GPUDEBUG
    printf(" outElectron: tid=%i Etot=%f, P=(%e;%e;%e)\n",tid,outpx,outpy,outpz);
    printf(" kernel - After brem: tid=%i E=%f Ekin=%f  P=(%f;%f;%f)\n",tid,
	   (track[tid].E+mass)/MeV,
	   track[tid].E/MeV,
	   track[tid].px/MeV,
	   track[tid].py/MeV,
	   track[tid].pz/MeV);
#endif // GPUDEBUG

    hmgr.getHisto("hbremStepLengthPost").Fill(step1/mm);
    hmgr.getHisto("hbremStepLengthAlong").Fill(step2/mm);
    hmgr.getHisto("hbremStepLengthDiff").Fill((step2-step1)/mm);
    hmgr.getHisto("hbremEnergyLoss").Fill((inpTotEnergy-outTotEnergy)/MeV);

    //== make some plots: electron change due to physics process
    G4double outEleAngle = acos((inppx*outpx+inppy*outpy+inppz*outpz)/(inpMomentum*outMomentum));
    hmgr.getHisto("hbremEnergy").Fill(outTotEnergy-mass);
    hmgr.getHisto("hbremAngle").Fill(outEleAngle);

    G4int nsec = aPartChange.GetNumberOfSecondaries();
    G4double secmass, secmom, secKinE, secTotE;
    for(int isec = 0 ; isec < nsec ; ++isec) {
      //offset is a global counter for the last array position of secondaries 
      GXTrack* aSecTrack = aPartChange.GetSecondary(isec);
      secpx = aSecTrack->px;
      secpy = aSecTrack->py;
      secpz = aSecTrack->pz;
      secmom = sqrt(secpx*secpx+secpy*secpy+secpz*secpz);
      secmass = electron_mass_c2 * aSecTrack->q*aSecTrack->q;
      secTotE = sqrt(secmom*secmom+secmass*secmass);
      secKinE = secTotE - secmass;
      hmgr.getHisto("hbremSecEnergy").Fill(secKinE/MeV);

      G4double secAngle = acos((inppx*secpx+inppy*secpy+inppz*secpz)/(inpMomentum*secmom));
      hmgr.getHisto("hbremSecAngle").Fill(secAngle);

#ifdef GPUDEBUG
      printf("Gamma: tid=%i, E=%e, kinE=%e,  P=(%e;%e;%e)\n",tid,secTotE/MeV,secKinE/MeV,secpx/MeV,secpy/MeV,secpz/MeV);
      printf("Rel.difference: tid=%i, P=(%e;%e;%e)\n",tid,
	     (inppx-outpx-secpx)/inppx*100,
       	     (inppy-outpy-secpy)/inppy*100,
	     (inppz-outpz-secpz)/inppz*100);
      printf("Angles: pc=%e, gamma=%e\n\n",outEleAngle,secAngle);
#endif
    }
#endif // __CUDA_ARCH__

    // store secondaries in a global stack (fixed size memory)
    if(isStack==0) StackSecondaries(&aPartChange,secTracks,stackSize,0);
    //if(runType==1) Print(&track[tid], &aPartChange);
  }
}

//-----------------------------------------------------------------------------
// GPeBremsstrahlungProcess + Dynamic memory for secondaries
//-----------------------------------------------------------------------------

GLOBALFUNC
void brem_kernel_dma(curandState* devStates, GXTrack *track, 
		     size_t nTracks,
		     GPPhysicsTable* eBrem_table, 
		     GPPhysicsTable* eIoni_table, 
		     GPPhysicsTable* msc_table,
		     GPPhysics2DVector* sbData,
		     GXSecContainer *container,
		     G4int *stackSize,
		     G4int *offset) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GPMaterial aMat;
  DefineMaterial(&aMat);

  //EM models and processes
  GPSeltzerBergerRelModel aModel(devStates,tid,sbData);
  GPeBremsstrahlung brem(devStates,tid,eBrem_table);
  brem.InitialiseProcess(&aModel);

  //local varialbles
  G4double step  = 0;
  G4double step1 = 0;
  G4double step2 = 0;
  GPForceCondition condition;
  GPGPILSelection  selection;
  GPVParticleChange aPartChange;

  //secondary
  __shared__ GXTrack* secBuffer;
  __shared__ int* tsec;

  G4int maxSec = (nTracks%gridDim.x !=0 ) ?
    (nTracks/gridDim.x + 1)*maxSecondaryPerStep : 
    (nTracks/gridDim.x)*maxSecondaryPerStep ;

  if(threadIdx.x ==0 ) {
    secBuffer = (GXTrack*)malloc(maxSec*sizeof(GXTrack));
    tsec = (int*)malloc(sizeof(int));
    *tsec = 0;
  }

  __syncthreads();

  while (tid < nTracks) {

    step = track[tid].s;
    brem.StartTracking();

    step1 = brem.PostStepGetPhysicalInteractionLength(&track[tid],step,
						      &condition);
    step2 = brem.AlongStepGetPhysicalInteractionLength(&selection);
    track[tid].s = (step1 < step2) ? step1 : step2;

    aPartChange = brem.AlongStepDoIt(&track[tid],&aMat,track[tid].s);
    aPartChange = brem.PostStepDoIt(&track[tid],&aMat);
    track[tid].E = aPartChange.GetProposedKineticEnergy();

    //secondary collection
    GXTrack aGamma = aModel.GetSecondary();
    secBuffer[atomicAdd(tsec,1)] = aGamma;

    tid += blockDim.x * gridDim.x;
  }

  //ensure that all threads complete before copying
  __syncthreads();

  //store the size/address of secondaries and copy them to the global memory
  if(threadIdx.x ==0 ) {

    G4int blockIndex = blockIdx.x;
    atomicAdd(stackSize,*tsec);

    container[blockIndex].size = *tsec;
    container[blockIndex].addr = (GXTrack*)malloc((*tsec)*sizeof(GXTrack));
    memcpy(container[blockIndex].addr,secBuffer,(*tsec)*sizeof(GXTrack));
 
    //free the shared memory
    free(secBuffer);
    free(tsec);
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void brem_gpu_dma(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* eBrem_table, 
		  GPPhysicsTable* eIoni_table, 
		  GPPhysicsTable* msc_table,
		  GPPhysics2DVector* sbData,
		  GXSecContainer *secContainer_d,
		  G4int *stackSize,
		  G4int *offset,
		  int blocksPerGrid, int threadsPerBlock) 
{
  brem_kernel_dma<<< blocksPerGrid, threadsPerBlock >>>(devStates,track,
		     nTrackSize,eBrem_table,eIoni_table,msc_table,
		     sbData,secContainer_d,stackSize,offset);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void brem_cpu_dma(GXTrack *track, size_t nTracks,
		  GPPhysicsTable* eBrem_table, 
		  GPPhysicsTable* eIoni_table, 
		  GPPhysicsTable* msc_table,
		  GPPhysics2DVector* sbData)
{

  GPMaterial aMat;
  DefineMaterial(&aMat);

  GPForceCondition condition;
  GPGPILSelection  selection;

  G4double step  = 0;
  G4double step1 = 0;
  G4double step2 = 0;

  //EM models and processes
  GPSeltzerBergerRelModel aModel(0,-1,sbData);
  GPeBremsstrahlung brem(0,-1, eBrem_table);
  brem.InitialiseProcess(&aModel);

  GPVParticleChange aPartChange;

  //secondary
  G4int stackSize = 0;

  GXSecContainer *container 
    = (GXSecContainer*) malloc(nTracks*sizeof(GXSecContainer));

  for (int tid = 0; tid < nTracks; ++tid) {

    step = track[tid].s;
    brem.StartTracking();

    step1 = brem.PostStepGetPhysicalInteractionLength(&track[tid],step,
						      &condition);
    step2 = brem.AlongStepGetPhysicalInteractionLength(&selection);
    track[tid].s = (step1 < step2) ? step1 : step2;

    aPartChange = brem.AlongStepDoIt(&track[tid],&aMat,track[tid].s);
    aPartChange = brem.PostStepDoIt(&track[tid],&aMat);
    track[tid].E = aPartChange.GetProposedKineticEnergy();

    //secondary collection
    G4int nsec = aPartChange.GetNumberOfSecondaries();
    GXTrack* secondaries = (GXTrack*)malloc(nsec*sizeof(GXTrack));

    for(int isec = 0 ; isec < nsec ; ++isec) {
      GXTrack* aSecTrack = aPartChange.GetSecondary(isec);
      secondaries[isec] = *aSecTrack;
    }

    container[tid].size = nsec;
    container[tid].addr = secondaries;
    stackSize += nsec;

  }

  GXTrack *secTracks = (GXTrack*)malloc(stackSize*sizeof(GXTrack));
  G4int offset = 0;

  for(int it = 0 ; it < nTracks ; ++it) {
    memcpy(secTracks+offset,container[it].addr,
           container[it].size*sizeof(GXTrack));
    free(container[it].addr);
    offset += container[it].size; 
  }

  free(secTracks);
}

//-----------------------------------------------------------------------------
// GPeIonisationProcess + Preallocated fixed memory for secondaries
//-----------------------------------------------------------------------------

GLOBALFUNC
void ioni_kernel(curandState* devStates, GXTrack *track, 
		 size_t nTracks,
		 GPPhysicsTable* eBrem_table, 
		 GPPhysicsTable* eIoni_table, 
		 GPPhysicsTable* eIoni_range, 
		 GPPhysicsTable* eIoni_dedx, 
		 GPPhysicsTable* eIoni_invr, 
		 GPPhysicsTable* msc_table,
		 GPPhysics2DVector* sbData,
		 GXTrack *secTracks, G4int *stackSize, G4int *offset, 
		 G4int isStack, G4int runType)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GPMaterial aMat;
  DefineMaterial(&aMat);

  //EM models and processes
  GPMollerBhabhaModel aModel(devStates,tid);
  GPeIonisation ioni(devStates,tid, eIoni_table, eIoni_range,
		     eIoni_dedx, eIoni_invr);
  ioni.InitialiseProcess(&aModel);

  //local variables
  G4double step  = 0;
  G4double step1 = 0;
  G4double step2 = 0;
  GPForceCondition condition;
  GPGPILSelection  selection;
  GPVParticleChange aPartChange;

  while (tid < nTracks) {

    step = track[tid].s;
    ioni.StartTracking();

    step1 = ioni.PostStepGetPhysicalInteractionLength(&track[tid],step,
						      &condition);
    step2 = ioni.AlongStepGetPhysicalInteractionLength(&selection); 
    track[tid].s = (step1 < step2) ? step1 : step2;

    aPartChange = ioni.AlongStepDoIt(&track[tid],&aMat,track[tid].s);
    aPartChange = ioni.PostStepDoIt(&track[tid],&aMat);
    track[tid].E = aPartChange.GetProposedKineticEnergy();

    if(isStack==0) StackSecondaries(&aPartChange,secTracks,stackSize,offset);
    //    if(runType==1) Print(&track[tid], &aPartChange);

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void ioni_gpu(curandState* devStates, GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* eBrem_table, 
	      GPPhysicsTable* eIoni_table, 
	      GPPhysicsTable* eIoni_range, 
	      GPPhysicsTable* eIoni_dedx, 
	      GPPhysicsTable* eIoni_invr, 
	      GPPhysicsTable* msc_table,
	      GPPhysics2DVector* sbData,
	      GXTrack *secTracks, G4int *stackSize, G4int *offset, 
	      G4int isStack, G4int runType, 
	      int blocksPerGrid, int threadsPerBlock) 
{
  ioni_kernel<<< blocksPerGrid, threadsPerBlock >>>(devStates,track,
		       nTrackSize,eBrem_table,eIoni_table,
	               eIoni_range,eIoni_dedx,eIoni_invr,msc_table,
		       sbData,secTracks,stackSize,offset,isStack,runType);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void ioni_cpu(GXTrack *track, size_t nTracks,
	      GPPhysicsTable* eBrem_table, 
	      GPPhysicsTable* eIoni_table, 
	      GPPhysicsTable* eIoni_range, 
	      GPPhysicsTable* eIoni_dedx, 
	      GPPhysicsTable* eIoni_invr, 
	      GPPhysicsTable* msc_table,
	      GPPhysics2DVector* sbData,
	      GXTrack *secTracks, G4int *stackSize,
	      G4int isStack, G4int runType)
{
  GPMaterial aMat;
  DefineMaterial(&aMat);
  G4double mass = electron_mass_c2;

  //EM models and processes
  GPMollerBhabhaModel aModel(0,-1);
  GPeIonisation ioni(0,-1, eIoni_table, eIoni_range,
		     eIoni_dedx, eIoni_invr);
  ioni.InitialiseProcess(&aModel);

  //local variables
  G4double step  = 0;
  G4double step1 = 0;
  G4double step2 = 0;
  GPForceCondition condition;
  GPGPILSelection  selection;
  GPVParticleChange aPartChange;

#ifndef __CUDA_ARCH__
  GPHistoManager& hmgr = GPHistoManager::getInstance();
#endif

  for (int tid = 0; tid < nTracks; ++tid) {

    step = track[tid].s;
    ioni.StartTracking();

#ifndef __CUDA_ARCH__
    G4double inppx = track[tid].px;
    G4double inppy = track[tid].py;
    G4double inppz = track[tid].pz;
    G4double inpMomentum = sqrt(inppx*inppx+inppy*inppy+inppz*inppz);
    G4double inpTotEnergy = sqrt(inpMomentum*inpMomentum+mass*mass);
#ifdef GPUDEBUG
    printf("\nkernel - Before ioni: tid=%i Etot=%f Ekin=%f  P=(%f;%f;%f)\n",tid,
           (track[tid].E+mass)/MeV, track[tid].E/MeV, inppx/MeV, inppy/MeV, inppz/MeV);
#endif // GPUDEBUG
#endif // __CUDA_ARCH__

    step1 = ioni.PostStepGetPhysicalInteractionLength(&track[tid],step,&condition);
    step2 = ioni.AlongStepGetPhysicalInteractionLength(&selection);
    track[tid].s = (step1 < step2) ? step1 : step2;

    aPartChange = ioni.AlongStepDoIt(&track[tid],&aMat,track[tid].s);
    aPartChange = ioni.PostStepDoIt(&track[tid],&aMat);

    G4double outKinEnergy = aPartChange.GetProposedKineticEnergy();
    G4double outTotEnergy = outKinEnergy + mass;
    G4double outMomentum = sqrt(outTotEnergy*outTotEnergy - mass*mass);
    GPThreeVector outMomVector = GPThreeVector_mult(aPartChange.GetProposedMomentumDirection(), outMomentum );
    track[tid].s = (step1 < step2) ? step1 : step2;
    track[tid].E = outKinEnergy;
    track[tid].px = GPThreeVector_x(outMomVector);
    track[tid].py = GPThreeVector_y(outMomVector);
    track[tid].pz = GPThreeVector_z(outMomVector);

#ifndef __CUDA_ARCH__
    //== make some plots: electron
    G4double outpx = track[tid].px;
    G4double outpy = track[tid].py;
    G4double outpz = track[tid].pz;
#ifdef GPUDEBUG
    printf(" Ioni: outElectron: tid=%i Etot=%f MeV, P[MeV]=(%e;%e;%e)\n",tid,outpx,outpy,outpz);
    printf(" kernel - After ioni: tid=%i E=%f MeV, Ekin=%f MeV, P[MeV]=(%f;%f;%f)\n",tid,
           (track[tid].E+mass)/MeV,
           track[tid].E/MeV,
           track[tid].px/MeV,
           track[tid].py/MeV,
           track[tid].pz/MeV);
#endif // GPUDEBUG

    hmgr.getHisto("hioniStepLengthPost").Fill(step1/mm);
    hmgr.getHisto("hioniStepLengthAlong").Fill(step2/mm);
    hmgr.getHisto("hioniStepLengthDiff").Fill((step2-step1)/mm);
    hmgr.getHisto("hioniEnergyLoss").Fill((inpTotEnergy-outTotEnergy)/MeV);

    //== make some plots: electron change due to physics process
    G4double outEleAngle = acos((inppx*outpx+inppy*outpy+inppz*outpz)/(inpMomentum*outMomentum));
    hmgr.getHisto("hioniEnergy").Fill(outTotEnergy);
    hmgr.getHisto("hioniAngle").Fill(outEleAngle);

    G4int nsec = aPartChange.GetNumberOfSecondaries();
    for(int isec = 0 ; isec < nsec ; ++isec) {
      //offset is a global counter for the last array position of secondaries 
      GXTrack* aSecTrack = aPartChange.GetSecondary(isec);
      G4double secpx = aSecTrack->px;
      G4double secpy = aSecTrack->py;
      G4double secpz = aSecTrack->pz;
      G4double secMomentum = sqrt(secpx*secpx+secpy*secpy+secpz*secpz);
      hmgr.getHisto("hioniSecEnergy").Fill(secMomentum/MeV);

      G4double secAngle = acos((inppx*secpx+inppy*secpy+inppz*secpz)/(inpMomentum*secMomentum));
      hmgr.getHisto("hioniSecAngle").Fill(secAngle);

#ifdef GPUDEBUG
      printf("Gamma: tid=%i, E=%e, P=(%e;%e;%e)\n",tid,secMomentum/MeV,secpx/MeV,secpy/MeV,secpz/MeV);
      printf("Rel.difference: tid=%i, P=(%e;%e;%e)\n",tid,
             (inppx-outpx-secpx)/inppx*100,
             (inppy-outpy-secpy)/inppy*100,
             (inppz-outpz-secpz)/inppz*100);
      printf("Angles: pc=%e, gamma=%e\n\n",outEleAngle/MeV,secAngle);
#endif
    }
#endif // __CUDA_ARCH__

    if(isStack==0) StackSecondaries(&aPartChange,secTracks,stackSize,0);
    //    if(runType==1) Print(&track[tid], &aPartChange);
  }
}

//-----------------------------------------------------------------------------
// GPeIonisationProcess + Dynamic memory for secondaries
//-----------------------------------------------------------------------------

GLOBALFUNC
void ioni_kernel_dma(curandState* devStates, GXTrack *track, 
		     size_t nTracks,
		     GPPhysicsTable* eBrem_table, 
		     GPPhysicsTable* eIoni_table, 
		     GPPhysicsTable* eIoni_range, 
		     GPPhysicsTable* eIoni_dedx, 
		     GPPhysicsTable* eIoni_invr, 
		     GPPhysicsTable* msc_table,
		     GPPhysics2DVector* sbData,
		     GXSecContainer *container,
		     G4int *stackSize,
		     G4int *offset) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GPMaterial aMat;
  DefineMaterial(&aMat);

  //EM models and processes
  GPMollerBhabhaModel aModel(devStates,tid);
  GPeIonisation ioni(devStates,tid, eIoni_table, eIoni_range,
		     eIoni_dedx, eIoni_invr);
  ioni.InitialiseProcess(&aModel);

  //local variables
  G4double step  = 0;
  G4double step1 = 0;
  G4double step2 = 0;
  GPForceCondition condition;
  GPGPILSelection  selection;
  GPVParticleChange aPartChange;

  //secondary
  __shared__ GXTrack* secBuffer;
  __shared__ int* tsec;

  G4int maxSec = (nTracks%gridDim.x !=0 ) ?
    (nTracks/gridDim.x + 1)*maxSecondaryPerStep : 
    (nTracks/gridDim.x)*maxSecondaryPerStep ;

  if(threadIdx.x ==0 ) {
    secBuffer = (GXTrack*)malloc(maxSec*sizeof(GXTrack));
    tsec = (int*)malloc(sizeof(int));
    *tsec = 0;
  }

  __syncthreads();

  //  G4int nsec = 0;
  //  GXTrack aTrack;

  while (tid < nTracks) {

    step = track[tid].s;
    ioni.StartTracking();

    step1 = ioni.PostStepGetPhysicalInteractionLength(&track[tid],step,
						      &condition);
    step2 = ioni.AlongStepGetPhysicalInteractionLength(&selection);
    track[tid].s = (step1 < step2) ? step1 : step2;

    aPartChange = ioni.AlongStepDoIt(&track[tid],&aMat,track[tid].s);
    aPartChange = ioni.PostStepDoIt(&track[tid],&aMat);

    track[tid].E = aPartChange.GetProposedKineticEnergy();

    //secondary
    GXTrack aSecTrack = aModel.GetSecondary();
    secBuffer[atomicAdd(tsec,1)] = aSecTrack;

    tid += blockDim.x * gridDim.x;
  }

  //ensure that all threads complete before copying
  __syncthreads();

  //store the size/address of secondaries and copy them to the global memory
  if(threadIdx.x ==0 ) {

    G4int blockIndex = blockIdx.x;
    atomicAdd(stackSize,*tsec);

    container[blockIndex].size = *tsec;
    container[blockIndex].addr = (GXTrack*)malloc((*tsec)*sizeof(GXTrack));
    memcpy(container[blockIndex].addr,secBuffer,(*tsec)*sizeof(GXTrack));
 
    //free the shared memory
    free(secBuffer);
    free(tsec);
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void ioni_gpu_dma(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* eBrem_table, 
		  GPPhysicsTable* eIoni_table, 
		  GPPhysicsTable* eIoni_range, 
		  GPPhysicsTable* eIoni_dedx, 
		  GPPhysicsTable* eIoni_invr, 
		  GPPhysicsTable* msc_table,
		  GPPhysics2DVector* sbData,
		  GXSecContainer *secContainer_d,
		  G4int *stackSize,
		  G4int *offset,
		  int blocksPerGrid, int threadsPerBlock) 
{
  ioni_kernel_dma<<< blocksPerGrid, threadsPerBlock >>>(devStates,track,
		       nTrackSize,eBrem_table,eIoni_table,
	               eIoni_range,eIoni_dedx,eIoni_invr,
                       msc_table,
		       sbData,secContainer_d,stackSize,offset);
}

void ioni_cpu_dma(GXTrack *track, size_t nTracks,
		  GPPhysicsTable* eBrem_table, 
		  GPPhysicsTable* eIoni_table, 
		  GPPhysicsTable* eIoni_range, 
		  GPPhysicsTable* eIoni_dedx, 
		  GPPhysicsTable* eIoni_invr, 
		  GPPhysicsTable* msc_table,
		  GPPhysics2DVector* sbData)
{

  GPMaterial aMat;
  DefineMaterial(&aMat);

  //EM models and processes
  GPMollerBhabhaModel aModel(0,-1);
  GPeIonisation ioni(0,-1, eIoni_table, eIoni_range,
		     eIoni_dedx, eIoni_invr);
  ioni.InitialiseProcess(&aModel);

  //local variables
  G4double step  = 0;
  G4double step1 = 0;
  G4double step2 = 0;
  GPForceCondition condition;
  GPGPILSelection  selection;
  GPVParticleChange aPartChange;

  //secondary
  G4int stackSize = 0;

  GXSecContainer *container 
    = (GXSecContainer*) malloc(nTracks*sizeof(GXSecContainer));

  for (int tid = 0; tid < nTracks; ++tid) {

    step = track[tid].s;
    ioni.StartTracking();

    step1 = ioni.PostStepGetPhysicalInteractionLength(&track[tid],step,
						      &condition);
    step2 = ioni.AlongStepGetPhysicalInteractionLength(&selection);
    track[tid].s = (step1 < step2) ? step1 : step2;

    aPartChange = ioni.AlongStepDoIt(&track[tid],&aMat,track[tid].s);
    aPartChange = ioni.PostStepDoIt(&track[tid],&aMat);
    track[tid].E = aPartChange.GetProposedKineticEnergy();

    //secondary collection
    G4int nsec = aPartChange.GetNumberOfSecondaries();;
    GXTrack* secondaries = (GXTrack*)malloc(nsec*sizeof(GXTrack));

    for(int isec = 0 ; isec < nsec ; ++isec) {
      GXTrack* aSecTrack = aPartChange.GetSecondary(isec);
      secondaries[isec] = *aSecTrack;
    }

    container[tid].size = nsec;
    container[tid].addr = secondaries;
    stackSize += nsec;

  }

  GXTrack *secTracks = (GXTrack*)malloc(stackSize*sizeof(GXTrack));
  G4int offset = 0;

  for(int it = 0 ; it < nTracks ; ++it) {
    memcpy(secTracks+offset,container[it].addr,
           container[it].size*sizeof(GXTrack));
    free(container[it].addr);
    offset += container[it].size; 
  }

  free(secTracks);
}

//-----------------------------------------------------------------------------
// GPeMultipleScatteringProcess
//-----------------------------------------------------------------------------

GLOBALFUNC
void msc_kernel(curandState* devStates, GXTrack *track, 
		size_t nTrackSize,
		GPPhysicsTable* eBrem_table, 
		GPPhysicsTable* eIoni_table, 
		GPPhysicsTable* msc_table,
		GPPhysics2DVector* sbData,
		GXTrack *secTracks, G4int *stackSize, G4int *offset, 
		G4int isStack, G4int runType)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GPMaterial aMat;
  DefineMaterial(&aMat);

  //EM models and processes
  GPUrbanMscModel95 cModel(devStates,tid);
  GPeMultipleScattering msc(devStates,tid,msc_table);
  msc.InitialiseProcess(&cModel);

  //local vairables
  G4double step  = 0;
  G4double ekin  = 0;
  G4double step1 = 0;
  G4double step2 = 0;
  GPForceCondition condition;
  GPGPILSelection  selection;
  GPVParticleChange aPartChange;

  while (tid < nTrackSize) {

    ekin = GetKineticEnergy(&track[tid]);
    step = track[tid].s;
    msc.StartTracking();

    step1 = msc.PostStepGetPhysicalInteractionLength(&track[tid],&condition);
    step2 = msc.AlongStepGetPhysicalInteractionLength(&aMat,ekin,step,
    						      &selection);
    aPartChange = msc.AlongStepDoIt(&aMat,&track[tid]);
    aPartChange = msc.PostStepDoIt(&track[tid]);

    track[tid].s = (step1 < step2) ? step1 : step2;
    track[tid].E = aPartChange.GetProposedKineticEnergy();

    if(isStack==0) StackSecondaries(&aPartChange,secTracks,stackSize,offset);
    //    if(runType==1) Print(&track[tid], &aPartChange);

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void msc_gpu(curandState* devStates,
	     GXTrack *track, size_t nTrackSize,
	     GPPhysicsTable* eBrem_table, 
	     GPPhysicsTable* eIoni_table, 
	     GPPhysicsTable* msc_table,
	     GPPhysics2DVector* sbData,
	     GXTrack *secTracks, G4int *stackSize, G4int *offset, 
	     G4int isStack, G4int runType,
	     int blocksPerGrid, int threadsPerBlock) 
{
  msc_kernel<<< blocksPerGrid, threadsPerBlock >>>(devStates,track,
		       nTrackSize,eBrem_table,eIoni_table,msc_table,
		       sbData,secTracks,stackSize,offset,isStack,runType);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void msc_cpu(GXTrack *track, size_t nTrackSize,
	     GPPhysicsTable* eBrem_table, 
	     GPPhysicsTable* eIoni_table, 
	     GPPhysicsTable* msc_table,
	     GPPhysics2DVector* sbData,
	     GXTrack *secTracks, G4int *stackSize,
	     G4int isStack, G4int runType)
{
  GPMaterial aMat;
  DefineMaterial(&aMat);

  //EM models and processes
  GPUrbanMscModel95 cModel(0,-1);
  GPeMultipleScattering msc(0,-1,msc_table);
  msc.InitialiseProcess(&cModel);

  //local vairables
  G4double step  = 0;
  G4double ekin  = 0;
  G4double step1 = 0;
  G4double step2 = 0;
  GPForceCondition condition;
  GPGPILSelection  selection;
  GPVParticleChange aPartChange;

  for (int tid = 0; tid < nTrackSize; ++tid) {

    ekin = GetKineticEnergy(&track[tid]);
    step = track[tid].s;
    msc.StartTracking();

    step1 = msc.PostStepGetPhysicalInteractionLength(&track[tid],&condition);
    step2 = msc.AlongStepGetPhysicalInteractionLength(&aMat,ekin,step,
    						      &selection);
    aPartChange = msc.AlongStepDoIt(&aMat,&track[tid]);
    aPartChange = msc.PostStepDoIt(&track[tid]);

    track[tid].s = (step1 < step2) ? step1 : step2;
    track[tid].E = aPartChange.GetProposedKineticEnergy();

    if(isStack==0) StackSecondaries(&aPartChange,secTracks,stackSize,0);
    //    if(runType==1) Print(&track[tid], &aPartChange);

  }
}

//-----------------------------------------------------------------------------
// Other kernels
//-----------------------------------------------------------------------------
#include "dma_kernel.cu"

//-----------------------------------------------------------------------------
// Common
//-----------------------------------------------------------------------------
#include "GPStep.cu"
#include "GPStepPoint.cu"
#include "GPThreeVector.cu"
#include "GPPhysicsTable.cu"
#include "GPVParticleChange.cu"

//-----------------------------------------------------------------------------
// GPMaterial
//-----------------------------------------------------------------------------
#include "GPElement.cu"
#include "GPMaterial.cu"
#include "GPIonisParamMat.cu"

//EMPhysics
#include "GPeBremsstrahlung.cu"
#include "GPSeltzerBergerRelModel.cu"
#include "GPPhysics2DVector.cu"
#include "GPeIonisation.cu"
#include "GPMollerBhabhaModel.cu"
#include "GPUniversalFluctuation.cu"
#include "GPeMultipleScattering.cu"
#include "GPUrbanMscModel95.cu"
