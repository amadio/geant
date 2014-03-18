#include <iostream>

#include <cuda.h>
#include <curand_kernel.h>

//Common
#include "GXTrack.h"
#include "GPPhysicsTable.h"
#include "GPForceCondition.h"
#include "GPVParticleChange.h"

//Material
#include "GPElement.h"
#include "GPMaterial.h"
#include "GPSandiaTable.h"

//Photon Processes/Models
#include "GPPhotonModel.h"
#include "GPPhotonProcess.h"
#include "GPEmProcessType.h"

#include "stdio.h"

FQUALIFIER
G4double GetKineticEnergy(GXTrack *track) {
  G4double p = sqrt(track->px*track->px + track->py*track->py 
		  + track->pz*track->pz);
  G4double mass = electron_mass_c2*track->q*track->q;
  G4double ekin = p*p/(sqrt(p*p + mass*mass) + mass);
  return ekin;
}

FQUALIFIER
void DefineMaterial(GPMaterial *aMat) 
{
  //PbWO4 aka CMS crystal
  GPMaterial_Constructor_ByElement(aMat,8.28*g/cm3);

  GPElement ele_Pb;
  GPElement ele_W;
  GPElement ele_O4;
  GPElement_Constructor(&ele_Pb,82,207.2*g/mole);
  GPElement_Constructor(&ele_W,74,183.84*g/mole);
  GPElement_Constructor(&ele_O4,8*4,15.9994*4*g/mole);

  GPMaterial_AddElement(aMat,ele_Pb,0.45532661);
  GPMaterial_AddElement(aMat,ele_W,0.40403397);
  GPMaterial_AddElement(aMat,ele_O4,0.14063942);
}

FQUALIFIER
G4double GetGPIL(GPPhotonProcess* process, GPMaterial *aMat, 
	     GXTrack *track, GPForceCondition* condition) 
{
  G4double step = track->s;
  process->StartStepping();
  G4double proposedStep = 
    process->PostStepGetPhysicalInteractionLength(track,aMat,step,condition);
  return proposedStep;
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
  for(int isec = 0 ; isec < nsec ; ++isec) {
    GXTrack* strack = particleChange->GetSecondary(isec);
    printf("    (isec,ekin,px) = (%d,%f,%f)\n",isec,strack->E,strack->px);
  }
#else
  printf("GPU (step,ekin,nsec) = (%f,%f,%d)\n",track->s,track->E,nsec);
  for(int isec = 0 ; isec < nsec ; ++isec) {
    GXTrack* strack = particleChange->GetSecondary(isec);
    printf("    (isec,ekin,px) = (%d,%f,%f)\n",isec,strack->E,strack->px);
  }
#endif
}

//-----------------------------------------------------------------------------
// GPComptonScattering + G4KleinNishinaCompton
//-----------------------------------------------------------------------------

GLOBALFUNC
void compt_kernel(curandState* devStates,
		  GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* gCompt_table, 
		  GPPhysicsTable* gConv_table, 
		  GPPhysicsTable* gPhot_table,
		  GXTrack *secTracks, G4int *stackSize, G4int *offset,
		  G4int isStack, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GPMaterial aMat;
  DefineMaterial(&aMat);

  //  G4double step  = 0;
  G4double proposedStep = 0;
  GPForceCondition condition;
  GPVParticleChange particleChange;

  //EM models
  GPPhotonModel comptModel(devStates,tid);
  comptModel.SetProcessType(kCompton);

  //EM processes
  GPPhotonProcess compt(devStates,tid, gCompt_table);
  compt.InitialiseProcess(kCompton,&comptModel);

  while (tid < nTrackSize) {

    compt.StartTracking(&track[tid]);

    proposedStep = GetGPIL(&compt,&aMat,&track[tid],&condition);
    particleChange = compt.PostStepDoIt(&track[tid],&aMat);

    track[tid].s = proposedStep;
    track[tid].E = particleChange.GetProposedKineticEnergy();

    //store secondary particles
    if(isStack==0) 
      StackSecondaries(&particleChange,secTracks,stackSize,offset);    
    if(runType==1) Print(&track[tid],&particleChange);

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void compt_gpu(curandState* devStates,
	       GXTrack *track, size_t nTrackSize,
	       GPPhysicsTable* gCompt_table, 
	       GPPhysicsTable* gConv_table, 
	       GPPhysicsTable* gPhot_table,
	       GXTrack *secTracks, G4int *stackSize, G4int *offset,
	       G4int isStack, G4int runType,
	       int blocksPerGrid, int threadsPerBlock) 
{
  compt_kernel<<< blocksPerGrid, threadsPerBlock >>>(devStates,track,
		       nTrackSize,gCompt_table,gConv_table,gPhot_table,
		       secTracks,stackSize,offset,isStack, runType);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void compt_cpu(GXTrack *track, size_t nTrackSize,
	       GPPhysicsTable* gCompt_table, 
	       GPPhysicsTable* gConv_table, 
	       GPPhysicsTable* gPhot_table,
	       GXTrack *secTracks, G4int *stackSize,
	       G4int isStack, G4int runType)
{
  GPMaterial aMat;
  DefineMaterial(&aMat);

  G4double proposedStep = 0;
  GPForceCondition condition;
  GPVParticleChange particleChange;

  //EM models
  GPPhotonModel comptModel(0,-1);
  comptModel.SetProcessType(kCompton);

  //EM processes
  GPPhotonProcess compt(0,-1,gCompt_table);
  compt.InitialiseProcess(kCompton,&comptModel);

  for (int tid = 0; tid < nTrackSize; tid++) {

    compt.StartTracking(&track[tid]);

    proposedStep = GetGPIL(&compt,&aMat,&track[tid],&condition);
    particleChange = compt.PostStepDoIt(&track[tid],&aMat);

    track[tid].s = proposedStep;
    track[tid].E = particleChange.GetProposedKineticEnergy();
  
    //store secondary particles
    if(isStack==0) 
      StackSecondaries(&particleChange,secTracks,stackSize,0);    
    if(runType==1) Print(&track[tid],&particleChange);

  }

}

//-----------------------------------------------------------------------------
// GPGammaConversion + G4BetheHeitlerModel
//-----------------------------------------------------------------------------

GLOBALFUNC
void conv_kernel(curandState* devStates,
		 GXTrack *track, size_t nTrackSize,
		 GPPhysicsTable* gCompt_table, 
		 GPPhysicsTable* gConv_table, 
		 GPPhysicsTable* gPhot_table,
		 GXTrack *secTracks, G4int *stackSize, G4int *offset,
		 G4int isStack, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GPMaterial aMat;
  DefineMaterial(&aMat);

  G4double proposedStep = 0;
  GPForceCondition condition;
  GPVParticleChange particleChange;

  //EM models

  GPPhotonModel convModel(devStates,tid);
  convModel.SetProcessType(kConversion);

  //EM processes
  GPPhotonProcess conv(devStates,tid, gConv_table);
  conv.InitialiseProcess(kConversion,&convModel);

  while (tid < nTrackSize) {

    conv.StartTracking(&track[tid]);

    proposedStep = GetGPIL(&conv,&aMat,&track[tid],&condition);
    particleChange = conv.PostStepDoIt(&track[tid],&aMat);

    track[tid].s = proposedStep;
    track[tid].E = particleChange.GetProposedKineticEnergy();

    //store secondary particles
    if(isStack==0) 
      StackSecondaries(&particleChange,secTracks,stackSize,offset);    
    if(runType==1) Print(&track[tid],&particleChange);

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void conv_gpu(curandState* devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* gCompt_table, 
	      GPPhysicsTable* gConv_table, 
	      GPPhysicsTable* gPhot_table,
	      GXTrack *secTracks, G4int *stackSize, G4int *offset,
	      G4int isStack, G4int runType,
	      int blocksPerGrid, int threadsPerBlock) 
{
  conv_kernel<<< blocksPerGrid, threadsPerBlock >>>(devStates,track,
		       nTrackSize,gCompt_table,gConv_table,gPhot_table,
		       secTracks,stackSize,offset,isStack,runType);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void conv_cpu(GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* gCompt_table, 
	      GPPhysicsTable* gConv_table, 
	      GPPhysicsTable* gPhot_table,
	      GXTrack *secTracks, G4int *stackSize,
	      G4int isStack, G4int runType)
{
  GPMaterial aMat;
  DefineMaterial(&aMat);

  G4double proposedStep = 0;
  GPForceCondition condition;
  GPVParticleChange particleChange;

  //EM models
  GPPhotonModel convModel(0,-1);
  convModel.SetProcessType(kConversion);

  //EM processes
  GPPhotonProcess conv(0,-1, gConv_table);
  conv.InitialiseProcess(kConversion,&convModel);

  for (int tid = 0; tid < nTrackSize; tid++) {

    conv.StartTracking(&track[tid]);

    proposedStep = GetGPIL(&conv,&aMat,&track[tid],&condition);
    particleChange = conv.PostStepDoIt(&track[tid],&aMat);

    track[tid].s = proposedStep;
    track[tid].E = particleChange.GetProposedKineticEnergy();

    GXTrack* eminus = particleChange.GetSecondary(0);
    GXTrack* eplus  = particleChange.GetSecondary(1);

    //store secondary particles
    if(isStack==0) 
      StackSecondaries(&particleChange,secTracks,stackSize,0);    
    if(runType==1) Print(&track[tid],&particleChange);
  }
}

//-----------------------------------------------------------------------------
// G4PhotoElectricEffect + G4PEEffectFluoModel
//-----------------------------------------------------------------------------

GLOBALFUNC
void phot_kernel(curandState* devStates,
		 GXTrack *track, size_t nTrackSize,
		 GPPhysicsTable* gCompt_table, 
		 GPPhysicsTable* gConv_table, 
		 GPPhysicsTable* gPhot_table,
		 GXTrack *secTracks, G4int *stackSize, G4int *offset,
		 G4int isStack, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GPMaterial aMat;
  DefineMaterial(&aMat);

  G4double proposedStep = 0;
  GPForceCondition condition;
  GPVParticleChange particleChange;

  //EM models
  GPPhotonModel photModel(devStates,tid);
  photModel.SetProcessType(kPhotoElectric);

  GPSandiaTable aSandiaTable;
  GPSandiaTable_Constructor(&aSandiaTable);
  photModel.SetaSandiaTable(&aSandiaTable);

  //EM processes
  GPPhotonProcess phot(devStates,tid, gPhot_table);
  phot.InitialiseProcess(kPhotoElectric,&photModel);

  while (tid < nTrackSize) {

    phot.StartTracking(&track[tid]);

    proposedStep = GetGPIL(&phot,&aMat,&track[tid],&condition);
    particleChange = phot.PostStepDoIt(&track[tid],&aMat);

    track[tid].s = proposedStep;
    track[tid].E = particleChange.GetProposedKineticEnergy();

    //store secondary particles
    if(isStack==0) 
      StackSecondaries(&particleChange,secTracks,stackSize,offset);    
    if(runType==1) Print(&track[tid],&particleChange);

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void phot_gpu(curandState* devStates, 
	      GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* gCompt_table, 
	      GPPhysicsTable* gConv_table, 
	      GPPhysicsTable* gPhot_table,
	      GXTrack *secTracks, G4int *stackSize, G4int *offset,
	      G4int isStack, G4int runType,
	      int blocksPerGrid, int threadsPerBlock) 
{
  phot_kernel<<< blocksPerGrid, threadsPerBlock >>>(devStates,track,
		       nTrackSize,gCompt_table,gConv_table,gPhot_table,
		       secTracks,stackSize,offset,isStack,runType);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void phot_cpu(GXTrack *track, size_t nTrackSize,
	      GPPhysicsTable* gCompt_table, 
	      GPPhysicsTable* gConv_table, 
	      GPPhysicsTable* gPhot_table,
	      GXTrack *secTracks,
	      G4int *stackSize,
	      G4int isStack, G4int runType)
{
  GPMaterial aMat;
  DefineMaterial(&aMat);

  G4double proposedStep = 0;
  GPForceCondition condition;
  GPVParticleChange particleChange;

  //EM models
  GPPhotonModel photModel(0,-1);
  photModel.SetProcessType(kPhotoElectric);

  GPSandiaTable aSandiaTable;
  GPSandiaTable_Constructor(&aSandiaTable);
  photModel.SetaSandiaTable(&aSandiaTable);

  //EM processes
  GPPhotonProcess phot(0,-1, gPhot_table);
  phot.InitialiseProcess(kPhotoElectric,&photModel);

  for (int tid = 0; tid < nTrackSize; tid++) {

    phot.StartTracking(&track[tid]);

    proposedStep = GetGPIL(&phot,&aMat,&track[tid],&condition);
    particleChange = phot.PostStepDoIt(&track[tid],&aMat);

    track[tid].s = proposedStep;
    track[tid].E = particleChange.GetProposedKineticEnergy();

    //store secondary particles
    if(isStack==0) 
      StackSecondaries(&particleChange,secTracks,stackSize,0);    
    if(runType==1) Print(&track[tid],&particleChange);
  }
}
//-----------------------------------------------------------------------------
// GPPhotonProcess + GPPhotonModel
//-----------------------------------------------------------------------------

GLOBALFUNC
void gamma_kernel(curandState* devStates, 
		  GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* gCompt_table, 
		  GPPhysicsTable* gConv_table, 
		  GPPhysicsTable* gPhot_table,
		  GXTrack *secTracks, G4int *stackSize, G4int *offset,
		  G4int isStack, G4int runType) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GPMaterial aMat;
  DefineMaterial(&aMat);

  G4double step_compt = 0;
  G4double step_conv = 0;
  G4double step_phot = 0;
  G4double proposedStep = 0;
  GPForceCondition condition;
  GPVParticleChange particleChange;

  //Photon models
  GPPhotonModel comptModel(devStates,tid);
  comptModel.SetProcessType(kCompton);

  GPPhotonModel convModel(devStates,tid);
  convModel.SetProcessType(kConversion);

  GPPhotonModel photModel(devStates,tid);
  photModel.SetProcessType(kPhotoElectric);

  GPSandiaTable aSandiaTable;
  GPSandiaTable_Constructor(&aSandiaTable);
  photModel.SetaSandiaTable(&aSandiaTable);

  //Photon processes
  GPPhotonProcess compt(devStates,tid, gCompt_table);
  compt.InitialiseProcess(kCompton,&comptModel);

  GPPhotonProcess conv(devStates,tid, gConv_table);
  conv.InitialiseProcess(kConversion,&convModel);

  GPPhotonProcess phot(devStates,tid, gPhot_table);
  phot.InitialiseProcess(kPhotoElectric,&photModel);

  while (tid < nTrackSize) {

    compt.StartTracking(&track[tid]);
    conv.StartTracking(&track[tid]);
    phot.StartTracking(&track[tid]);

    step_compt = GetGPIL(&compt,&aMat,&track[tid],&condition);
    step_conv  = GetGPIL(&conv,&aMat,&track[tid],&condition);
    step_phot  = GetGPIL(&phot,&aMat,&track[tid],&condition);

    //physics process defining the current step
    GPEmProcessType process = kCompton;

    if(step_compt > step_conv) process = kConversion;   
    if(step_conv > step_phot) process = kPhotoElectric;      

    switch (process) {
    case kCompton : 
      proposedStep = step_compt;
      particleChange = compt.PostStepDoIt(&track[tid],&aMat);
      break;
    case kConversion :
      proposedStep = step_conv;
      particleChange = conv.PostStepDoIt(&track[tid],&aMat);
      break;
    case kPhotoElectric :
      proposedStep = step_phot;
      particleChange = phot.PostStepDoIt(&track[tid],&aMat);
      break;
    default :
      break;
    }

    track[tid].s = proposedStep;
    track[tid].E = particleChange.GetProposedKineticEnergy();

    //store secondary particles
    if(isStack==0) 
      StackSecondaries(&particleChange,secTracks,stackSize,offset);    
    if(runType==1) Print(&track[tid],&particleChange);

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void gamma_gpu(curandState* devStates, 
	       GXTrack *track, size_t nTrackSize,
	       GPPhysicsTable* gCompt_table, 
	       GPPhysicsTable* gConv_table, 
	       GPPhysicsTable* gPhot_table,
	       GXTrack *secTracks, G4int *stackSize, G4int *offset,
	       G4int isStack, G4int runType,
	       int blocksPerGrid, int threadsPerBlock) 
{
  gamma_kernel<<< blocksPerGrid, threadsPerBlock >>>(devStates,track,
		       nTrackSize,gCompt_table,gConv_table,gPhot_table,
		       secTracks,stackSize,offset,isStack, runType);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void gamma_cpu(GXTrack *track, size_t nTrackSize,
	       GPPhysicsTable* gCompt_table, 
	       GPPhysicsTable* gConv_table, 
	       GPPhysicsTable* gPhot_table,
	       GXTrack *secTracks, G4int *stackSize,
	       G4int isStack, G4int runType)
{
  GPMaterial aMat;
  DefineMaterial(&aMat);

  G4double step_compt = 0;
  G4double step_conv = 0;
  G4double step_phot = 0;
  G4double proposedStep = 0;
  GPForceCondition condition;
  GPVParticleChange particleChange;

  //Photon models
  GPPhotonModel comptModel(0,-1);
  comptModel.SetProcessType(kCompton);

  GPPhotonModel convModel(0,-1);
  convModel.SetProcessType(kConversion);

  GPPhotonModel photModel(0,-1);
  photModel.SetProcessType(kPhotoElectric);

  GPSandiaTable aSandiaTable;
  GPSandiaTable_Constructor(&aSandiaTable);
  photModel.SetaSandiaTable(&aSandiaTable);

  //Photon processes
  GPPhotonProcess compt(0,-1, gCompt_table);
  compt.InitialiseProcess(kCompton,&comptModel);

  GPPhotonProcess conv(0,-1, gConv_table);
  conv.InitialiseProcess(kConversion,&convModel);

  GPPhotonProcess phot(0,-1, gPhot_table);
  phot.InitialiseProcess(kPhotoElectric,&photModel);

  for (int tid = 0; tid < nTrackSize; tid++) {

    compt.StartTracking(&track[tid]);
    conv.StartTracking(&track[tid]);
    phot.StartTracking(&track[tid]);

    step_compt = GetGPIL(&compt,&aMat,&track[tid],&condition);
    step_conv  = GetGPIL(&conv,&aMat,&track[tid],&condition);
    step_phot  = GetGPIL(&phot,&aMat,&track[tid],&condition);

    //physics process defining the current step
    GPEmProcessType process = kCompton;

    if(step_compt > step_conv) process = kConversion;   
    if(step_conv > step_phot) process = kPhotoElectric;      

    switch (process) {
    case kCompton : 
      proposedStep = step_compt;
      particleChange = compt.PostStepDoIt(&track[tid],&aMat);
      break;
    case kConversion :
      proposedStep = step_conv;
      particleChange = conv.PostStepDoIt(&track[tid],&aMat);
      break;
    case kPhotoElectric :
      proposedStep = step_phot;
      particleChange = phot.PostStepDoIt(&track[tid],&aMat);
      break;
    default :
      break;
    }

    track[tid].s = proposedStep;
    track[tid].E = particleChange.GetProposedKineticEnergy();

    //store secondary particles
    if(isStack==0) 
      StackSecondaries(&particleChange,secTracks,stackSize,0);    
    if(runType==1) Print(&track[tid],&particleChange);
  }
}

//-----------------------------------------------------------------------------
// Other kernels
//-----------------------------------------------------------------------------
#include "random_kernel.cu"

//-----------------------------------------------------------------------------
// Common
//-----------------------------------------------------------------------------
#include "GPThreeVector.cu"
#include "GPPhysicsTable.cu"
#include "GPPhysics2DVector.cu"
#include "GPVParticleChange.cu"

//-----------------------------------------------------------------------------
// Material
//-----------------------------------------------------------------------------
#include "GPElement.cu"
#include "GPMaterial.cu"
#include "GPStep.cu"
#include "GPStepPoint.cu"
#include "GPAtomicShells.cu"
#include "GPSandiaTable.cu"

//-----------------------------------------------------------------------------
// Photon processes/models
//-----------------------------------------------------------------------------
#include "GPPhotonModel.cu"
#include "GPPhotonProcess.cu"
