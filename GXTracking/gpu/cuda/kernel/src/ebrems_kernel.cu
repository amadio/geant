#include <iostream>

#include <cuda.h>
#include <curand_kernel.h>

//Common
#include "GXTrack.h"
#include "GPPhysicsTable.h"
#include "GPPhysicsTableType.h"
#include "GPPhysics2DVector.h"

#include "GXVDiscreteSampling.h"
#include "GXSeltzerBerger.h"
#include "GXSamplingTexture.h"
#include "GXHistogramManager.h"

void SetupSamplingTexture(int size,
			  G4double* PDFX_d, G4double* PDFY_d)
{
  //texture
  GXSamplingTexture samplingTexture;
  samplingTexture.Load(size, PDFX_d, PDFY_d);
}

FQUALIFIER
void StoreSecondary(GXTrack* aTrack, GXTrack *aSecondaryStack, 
		    G4int *aStackSize) 
{
#ifndef __CUDA_ARCH__
  aSecondaryStack[*aStackSize] = *aTrack;
  ++(*aStackSize);
#else
  int theOffset = atomicAdd(aStackSize,1);
  aSecondaryStack[theOffset] = *aTrack;
#endif
}

//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_k0(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  
  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;

  while (tid < nTrackSize) {
    energy = track[tid].E;
    track[tid].px = energy + 1.0*tid;
    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_k0(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_k0<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ);
}

//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_k1(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  
  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;
  
  while (tid < nTrackSize) {
    energy =  track[tid].E;
    model.SampleByRandom(energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();
    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_k1(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_k1<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ);
}

//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_k2(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  
  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;
  G4double Z = 82;

  while (tid < nTrackSize) {

    energy = track[tid].E;
    model.SampleByInversePDF(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_k2(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_k2<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ);
}

//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_k3(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  
  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;
  G4double Z = 82;

  while (tid < nTrackSize) {
    energy = track[tid].E;
    model.SampleByInversePDFTexture(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_k3(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_k3<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ);
}

//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_k4(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  
  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;
  G4double Z = 82;

  while (tid < nTrackSize) {
    energy = track[tid].E;
    model.SampleByInversePDFLinearInterpolation(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_k4(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_k4<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ);
}

//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_k5(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;
  G4double Z = 82;

  while (tid < nTrackSize) {
    energy = track[tid].E;
    model.SampleByInversePDFTextureLinearInterpolation(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_k5(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_k5<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ);
}

//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_k6(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFA,PDFX,PDFQ);
  
  G4double energy = -999.;
  G4double Z = 82;

  while (tid < nTrackSize) {
    energy = track[tid].E;
    model.SampleByAlias(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_k6(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_k6<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ);
}


//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_k7(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFA,PDFX,PDFQ);
  
  G4double energy = -999.;
  G4double Z = 82;

  while (tid < nTrackSize) {
    energy = track[tid].E;
    model.SampleByAliasLinearInterpolation(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_k7(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_k7<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ);
}

//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_k8(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  
  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4int ntrial;
  G4double energy = -999.;
  G4double Z = 82;
  
  while (tid < nTrackSize) {
    ntrial = 0;
    energy = track[tid].E;
    model.SampleByCompositionRejection(Z,energy,ntrial);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_k8(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_k8<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ);
}

//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_k9(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  
  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4int ntrial;
  G4double energy = -999.;
  G4double Z = 82;
  
  while (tid < nTrackSize) {
    ntrial = 0;
    energy = track[tid].E;
    model.SampleByAverageTrials(Z,2.0,energy,ntrial);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_k9(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_k9<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ);
}

//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_v0(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ,
		     G4int *stackSize, GXTrack *secTrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  
  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;
  G4double Z = 82;

  while (tid < nTrackSize) {

    energy = track[tid].E;
    model.SampleByInversePDF(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    GXTrack aSec;
    model.CreateSecondary(&aSec,&track[tid],-1);
    StoreSecondary(&aSec,secTrack,stackSize);

    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_v0(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_v0<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ,
						      stackSize, secTrack);
}

//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_v1(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ,
		     G4int *stackSize, GXTrack *secTrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;
  G4double Z = 82;

  while (tid < nTrackSize) {
    energy = track[tid].E;
    model.SampleByInversePDFLinearInterpolation(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    GXTrack aSec;
    model.CreateSecondary(&aSec,&track[tid],-1);
    StoreSecondary(&aSec,secTrack,stackSize);

    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_v1(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_v1<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ,
						      stackSize, secTrack);
}

//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_v2(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ,
		     G4int *stackSize, GXTrack *secTrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFA,PDFX,PDFQ);
  
  G4double energy = -999.;
  G4double Z = 82;

  while (tid < nTrackSize) {
    energy = track[tid].E;
    model.SampleByAlias(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    GXTrack aSec;
    model.CreateSecondary(&aSec,&track[tid],-1);
    StoreSecondary(&aSec,secTrack,stackSize);

    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_v2(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_v2<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ,
						      stackSize, secTrack);
}


//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_v3(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ,
		     G4int *stackSize, GXTrack *secTrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFA,PDFX,PDFQ);
  
  G4double energy = -999.;
  G4double Z = 82;

  while (tid < nTrackSize) {
    energy = track[tid].E;
    model.SampleByAliasLinearInterpolation(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    GXTrack aSec;
    model.CreateSecondary(&aSec,&track[tid],-1);
    StoreSecondary(&aSec,secTrack,stackSize);

    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_v3(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_v3<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						      track, nTrackSize, 
						      physicsTable, sbData,
						      PDFX, PDFY, PDFA, PDFQ,
						      stackSize, secTrack);
}

//------------------------------------------------------------------------------
GLOBALFUNC
void model_kernel_v4(curandState* devStates, GXTrack *track, size_t nTrackSize,
		     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		     G4double *PDFX,G4double *PDFY, G4int *PDFA,G4double *PDFQ,
		     G4int *stackSize, GXTrack *secTrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4int ntrial;
  G4double energy = -999.;
  G4double Z = 82;

  while (tid < nTrackSize) {
    ntrial = 0;
    energy = track[tid].E;
    model.SampleByCompositionRejection(Z,energy,ntrial);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    GXTrack aSec;
    model.CreateSecondary(&aSec,&track[tid],-1);
    StoreSecondary(&aSec,secTrack,stackSize);

    tid += blockDim.x * gridDim.x;
  }
}

void model_gpu_v4(curandState* devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack,
		  int blocksPerGrid, int threadsPerBlock) 
{
  model_kernel_v4<<<blocksPerGrid, threadsPerBlock>>>(devStates,
						       track, nTrackSize, 
						       physicsTable, sbData,
						       PDFX, PDFY, PDFA, PDFQ,
						       stackSize, secTrack);
}

//CPU
//------------------------------------------------------------------------------
void model_cpu_k0(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;

  for (int tid = 0; tid < nTrackSize; tid++) {
    energy = track[tid].E;
    track[tid].px = energy + 1.0*tid;
  }
}

//------------------------------------------------------------------------------
void model_cpu_k1(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;

  for (int tid = 0; tid < nTrackSize; tid++) {
    energy =  track[tid].E;
    model.SampleByRandom(energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();
  }
}

//------------------------------------------------------------------------------
void model_cpu_k2(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    energy = track[tid].E;
    model.SampleByInversePDF(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();
  }
}

//------------------------------------------------------------------------------
void model_cpu_k3(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    energy = track[tid].E;
    model.SampleByInversePDFTexture(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();
  }
}

//------------------------------------------------------------------------------
void model_cpu_k4(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    energy = track[tid].E;
    model.SampleByInversePDFLinearInterpolation(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();
  }
}

//------------------------------------------------------------------------------
void model_cpu_k5(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    energy = track[tid].E;
    model.SampleByInversePDFTextureLinearInterpolation(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();
  }
}

//------------------------------------------------------------------------------
void model_cpu_k6(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFA,PDFX,PDFQ);
  
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    energy = track[tid].E;
    model.SampleByAlias(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();
  }
}

//------------------------------------------------------------------------------
void model_cpu_k7(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFA,PDFX,PDFQ);
  
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    energy = track[tid].E;
    model.SampleByAliasLinearInterpolation(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();
  }
}

//------------------------------------------------------------------------------
void model_cpu_k8(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4int ntrial;
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    ntrial = 0;
    energy = track[tid].E;
    model.SampleByCompositionRejection(Z,energy,ntrial);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();
  }
}

//------------------------------------------------------------------------------
void model_cpu_k9(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4int ntrial;
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    ntrial = 0;
    energy = track[tid].E;
    model.SampleByAverageTrials(Z,2.0,energy,ntrial);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();
  }
}

//------------------------------------------------------------------------------
void model_cpu_v0(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    energy = track[tid].E;
    model.SampleByInversePDF(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    GXTrack aSec;
    model.CreateSecondary(&aSec,&track[tid],-1);
    StoreSecondary(&aSec,secTrack,stackSize);
  }
}

//------------------------------------------------------------------------------
void model_cpu_v1(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    energy = track[tid].E;
    model.SampleByInversePDFLinearInterpolation(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    GXTrack aSec;
    model.CreateSecondary(&aSec,&track[tid],-1);
    StoreSecondary(&aSec,secTrack,stackSize);
  }
}

//------------------------------------------------------------------------------
void model_cpu_v2(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFA,PDFX,PDFQ);
  
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    energy = track[tid].E;
    model.SampleByAlias(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    GXTrack aSec;
    model.CreateSecondary(&aSec,&track[tid],-1);
    StoreSecondary(&aSec,secTrack,stackSize);
  }
}

//------------------------------------------------------------------------------
void model_cpu_v3(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFA,PDFX,PDFQ);
  
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    energy = track[tid].E;
    model.SampleByAliasLinearInterpolation(Z,energy);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    GXTrack aSec;
    model.CreateSecondary(&aSec,&track[tid],-1);
    StoreSecondary(&aSec,secTrack,stackSize);
  }
}

//------------------------------------------------------------------------------
void model_cpu_v4(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  
  G4int ntrial;
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    ntrial = 0;
    energy = track[tid].E;
    model.SampleByCompositionRejection(Z,energy,ntrial);
    track[tid].px = model.GetSecondaryEnergy()*model.GetSecondarySinTheta();

    GXTrack aSec;
    model.CreateSecondary(&aSec,&track[tid],-1);
    StoreSecondary(&aSec,secTrack,stackSize);
  }
}

//------------------------------------------------------------------------------
void cpu_val(GXTrack *track, size_t nTrackSize,
	     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
	     G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
	     GXHistogramManager* theHisto)
{
  int tid = -1;
  curandState* devStates = 0;

  GXSeltzerBerger model(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFX,PDFY);
  GXSeltzerBerger model2(devStates,tid,sbData,1.0,1001.,NROW,NCOL,PDFA,PDFX,PDFQ);
  
  G4int ntrial;
  G4double energy = -999.;
  G4double Z = 82;

  for (int tid = 0; tid < nTrackSize; tid++) {
    ntrial = 0;
    energy = track[tid].E;
    model.SampleByCompositionRejection(Z,energy,ntrial);

    theHisto->h_einc->Fill(energy);
    theHisto->h_ntrial->Fill(ntrial);

    theHisto->h_energy[0]->Fill(model.GetSecondaryEnergy());
    theHisto->h_angle[0]->Fill(model.GetSecondarySinTheta());

    model.SampleByInversePDF(Z,energy);
    theHisto->h_energy[1]->Fill(model.GetSecondaryEnergy());
    theHisto->h_angle[1]->Fill(model.GetSecondarySinTheta());

    model.SampleByInversePDFLinearInterpolation(Z,energy);
    theHisto->h_energy[2]->Fill(model.GetSecondaryEnergy());
    theHisto->h_angle[2]->Fill(model.GetSecondarySinTheta());

    model2.SampleByAlias(Z,energy);
    theHisto->h_energy[3]->Fill(model2.GetSecondaryEnergy());
    theHisto->h_angle[3]->Fill(model2.GetSecondarySinTheta());

    model2.SampleByAliasLinearInterpolation(Z,energy);
    theHisto->h_energy[4]->Fill(model2.GetSecondaryEnergy());
    theHisto->h_angle[4]->Fill(model2.GetSecondarySinTheta());
  }
}

#include "random_kernel.cu"

//-----------------------------------------------------------------------------
// Common
//-----------------------------------------------------------------------------
#include "GPThreeVector.cu"
#include "GPPhysicsTable.cu"
#include "GPPhysics2DVector.cu"
#include "GXVDiscreteSampling.cu"
#include "GXSeltzerBerger.cu"
#include "GXHistogramManager.cu"

//-----------------------------------------------------------------------------
// FrameWork and Application
//-----------------------------------------------------------------------------
#include "GXTrackHandler.cc"
