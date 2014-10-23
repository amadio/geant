#ifndef MODEL_KERNEL_H
#define MODEL_KERNEL_H 1

#include "GXTrack.h"
#include "GXFieldMap.h"
#include "GPPhysicsTable.h"
#include "GPPhysics2DVector.h"
#include "GXHistogramManager.h"

#include <curand.h>
#include <curand_kernel.h>

void SetupSamplingTexture(int size, G4double *X, G4double *Y);

const G4int nKernels = 10;
const G4int nModels  =  5;

// GPU

void model_gpu_k0(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int nBlocks, int nThreads); 

void model_gpu_k1(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int nBlocks, int nThreads); 

void model_gpu_k2(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int nBlocks, int nThreads); 

void model_gpu_k3(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int nBlocks, int nThreads); 

void model_gpu_k4(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int nBlocks, int nThreads); 

void model_gpu_k5(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int nBlocks, int nThreads); 

void model_gpu_k6(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int nBlocks, int nThreads); 

void model_gpu_k7(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int nBlocks, int nThreads); 

void model_gpu_k8(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int nBlocks, int nThreads); 

void model_gpu_k9(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  int nBlocks, int nThreads); 

typedef void (*KernelFunc_t)(curandState *devStates, 
			     GXTrack *track, size_t nTrackSize,
			     GPPhysicsTable* physicsTable, 
			     GPPhysics2DVector* sbData,
			     G4double *PDFX, G4double *PDFY, 
			     G4int *PDFA, G4double *PDFQ,
                             int nBlocks, int nThreads);

KernelFunc_t KernelFunc[] = {model_gpu_k0, 
                             model_gpu_k1, 
                             model_gpu_k2,
                             model_gpu_k3,
                             model_gpu_k4,
                             model_gpu_k5,
                             model_gpu_k6, 
                             model_gpu_k7, 
                             model_gpu_k8, 
                             model_gpu_k9}; 

//GPU output validation

void model_gpu_v0(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack,
		  int nBlocks, int nThreads); 

void model_gpu_v1(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack,
		  int nBlocks, int nThreads); 

void model_gpu_v2(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack,
		  int nBlocks, int nThreads); 

void model_gpu_v3(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack,
		  int nBlocks, int nThreads); 

void model_gpu_v4(curandState *devStates, GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack,
		  int nBlocks, int nThreads); 

typedef void (*GpuVal_t)(curandState *devStates, 
			 GXTrack *track, size_t nTrackSize,
			 GPPhysicsTable* physicsTable, 
			 GPPhysics2DVector* sbData,
			 G4double *PDFX, G4double *PDFY, 
			 G4int *PDFA, G4double *PDFQ,
			 G4int *stackSize, GXTrack *secTrack,
			 int nBlocks, int nThreads);

GpuVal_t GpuVal[] = {model_gpu_v0, 
		     model_gpu_v1,
		     model_gpu_v2,
		     model_gpu_v3,
		     model_gpu_v4};

//CPU

void model_cpu_k0(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ);

void model_cpu_k1(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ);

void model_cpu_k2(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ);

void model_cpu_k3(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ);

void model_cpu_k4(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ);

void model_cpu_k5(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ);

void model_cpu_k6(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ);

void model_cpu_k7(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ);

void model_cpu_k8(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ);

void model_cpu_k9(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ);

typedef void (*CpuFunc_t)(GXTrack *track, size_t nTrackSize,
                          GPPhysicsTable* physicsTable,
			  GPPhysics2DVector* sbData,
                          G4double *PDFX, G4double *PDFY,
			  G4int *PDFA, G4double *PDFQ);

CpuFunc_t CpuFunc[] = {model_cpu_k0,
		       model_cpu_k1,
		       model_cpu_k2,
                       model_cpu_k3,
                       model_cpu_k4,
                       model_cpu_k5,
                       model_cpu_k6,
                       model_cpu_k7,
                       model_cpu_k8,
		       model_cpu_k9}; 

//CPU output validation

void model_cpu_v0(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack);

void model_cpu_v1(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack);

void model_cpu_v2(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack);

void model_cpu_v3(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack);

void model_cpu_v4(GXTrack *track, size_t nTrackSize,
		  GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
		  G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
		  G4int *stackSize, GXTrack *secTrack);

typedef void (*CpuVal_t)(GXTrack *track, size_t nTrackSize,
                          GPPhysicsTable* physicsTable,
			  GPPhysics2DVector* sbData,
                          G4double *PDFX, G4double *PDFY,
			  G4int *PDFA, G4double *PDFQ,
                          G4int *stackSize, GXTrack *secTrack);

CpuVal_t CpuVal[] = {model_cpu_v0,
		     model_cpu_v1,
		     model_cpu_v2,
		     model_cpu_v3,
		     model_cpu_v4};

//validation of MC methods on CPU

void cpu_val(GXTrack *track, size_t nTrackSize,
	     GPPhysicsTable* physicsTable, GPPhysics2DVector* sbData,
	     G4double *PDFX, G4double *PDFY, G4int *PDFA, G4double *PDFQ,
	     GXHistogramManager* theHisto);

#endif
