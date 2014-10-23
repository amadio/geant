#include <iostream>
#include <fstream>
#include "stdio.h"

#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>

#include "GXTrack.h"
#include "GPEMPhysicsUtils.h"
#include "GXDiscreteSampling.h"
#include "GXSamplingTexture.h"
#include "GPPhysicsTableType.h"
#include "moller_kernel.h"
#include "random_kernel.h"

//Magnetic Field
#include "GPConstants.h"
#include "GXFieldMap.h"
#include "GXFieldMapData.h"

#include <fstream>
#include <iomanip>

#include "GXHistogramManager.h"
#include "GXTrackHandler.h"

//CUDA
#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>

using namespace std;

int main (int argc, char* argv[]) 
{
  // arguments: [runType]
  int runType     = 1; // 0 test 1 normal 2 validation

  if(argc >= 2) runType = atoi(argv[1]);
  if( runType <0 || runType > 2) {
    std::cout << "Usage: runType=[0-2]... " << std::endl;
    return 0;
  }

  int nDevice;
  bool cudaEnabled = false;

  cudaGetDeviceCount(&nDevice);

  if(nDevice > 0) {
    cudaDeviceReset();
    cudaEnabled = true;
    std::cout << "CUDA Enabled with " << nDevice << " Device(s)" << std::endl;
  }
  else {
    std::cout << "Waning: No Cuda Capable Device ... " << std::endl;
  }

  //set the default number of threads and thread blocks
  int theNBlocks  =  26;
  int theNThreads = 192;

  if(runType == 0 ) { 
    theNBlocks  =  2;
    theNThreads =  4;
  }
  char* cudaNBlocks;
  char* cudaNThreads;

  cudaNBlocks = getenv("GP_CUDA_NBLOCKS");
  if(cudaNBlocks) theNBlocks = atoi(cudaNBlocks);
  cudaNThreads = getenv("GP_CUDA_NTHREADS");
  if(cudaNThreads) theNThreads = atoi(cudaNThreads);
  
  std::cout << "...moller kernels <<<" << theNBlocks << "," 
            << theNThreads <<">>> (...) ..." << std::endl;
  
  // prepare physics table
  GPPhysicsTable physicsTable[kNumberPhysicsTable];

  char filename[256];
  for(int it = 0 ; it < kNumberPhysicsTable ; ++it) {
    sprintf(filename,"data/%s",GPPhysicsTableName[it]);
    readTableAndSetSpline(&physicsTable[it],filename);
  }

  //G4SeltzerBergerModel data
  G4int maxZ = 92;
  GPPhysics2DVector* sbData =
    (GPPhysics2DVector*) malloc(maxZ*sizeof(GPPhysics2DVector));

  char sbDataFile[256];
  for(G4int iZ = 0 ; iZ < maxZ ; iZ++) {  
    sprintf(sbDataFile,"data/brem_SB/br%d",iZ+1);
    std::ifstream fin(sbDataFile);
    G4bool check = RetrieveSeltzerBergerData(fin, &sbData[iZ]);
    if(!check) {
      printf("Failed To open SeltzerBerger Data file for Z= %d\n",iZ+1);
    }
  }

  printf("Copying Physics tables from host to device\n");
  
  GPPhysicsTable* physicsTable_d;
  GPPhysics2DVector* sbData_d;
  
  if(cudaEnabled) {
    cudaMalloc((void**)&physicsTable_d, 
               kNumberPhysicsTable*sizeof(GPPhysicsTable));
    cudaMemcpy(physicsTable_d, &physicsTable, 
               kNumberPhysicsTable*sizeof(GPPhysicsTable),
               cudaMemcpyHostToDevice);

    cudaMalloc((void**)&sbData_d,maxZ*sizeof(GPPhysics2DVector));
    cudaMemcpy(sbData_d, sbData, maxZ*sizeof(GPPhysics2DVector),
               cudaMemcpyHostToDevice);
  }

  //start time memory allocation on GPU
  cudaEvent_t start, stop;
  cudaEventCreate (&start);
  cudaEventCreate (&stop);

  float elapsedTimeInit_PDF = 0.0;
  float elapsedTimeInit_Alias = 0.0;
  float elapsedTimeAlloc_Tables = 0.0;

  //prepare PDF
  const G4int Z = 82;
  const G4double minP = 1.0;
  const G4double maxP = 1001.;

  G4double *PDFX;
  G4double *PDFY;
  G4int    *PDFA;
  G4double *PDFQ;

  G4double **pdfx = new G4double *[NROW];
  G4double **pdfy = new G4double *[NROW];
  G4int    **pdfa = new G4int    *[NROW];
  G4double **pdfq = new G4double *[NROW];
  
  for(int j=0;j<NROW;j++) {
    pdfx[j] = new G4double [NCOL];
    pdfy[j] = new G4double [NCOL];
    pdfa[j] = new G4int    [NCOL];
    pdfq[j] = new G4double [NCOL];
  }
  
  cudaEventRecord (start,0);
  //---------------------------------------------------------------------------
  BuildInversePDF_MollerBhabha(Z, true, minP,maxP,
			       NROW,NCOL, pdfx,pdfy); 
  //---------------------------------------------------------------------------
  cudaEventRecord (stop,0);
  cudaEventSynchronize (stop);
  cudaEventElapsedTime (&elapsedTimeInit_PDF,start,stop);
  printf("Time to build Moller PDF Tables = %f\n",elapsedTimeInit_PDF);

  cudaEventRecord (start,0);
  //---------------------------------------------------------------------------
  BuildAliasSamplingTable_MollerBhabha(Z, true, minP,maxP,
				       NROW, NCOL, pdfx, pdfa, pdfq); 
  //---------------------------------------------------------------------------
  cudaEventRecord (stop,0);
  cudaEventSynchronize (stop);
  cudaEventElapsedTime (&elapsedTimeInit_Alias,start,stop);
  printf("Time to build Moller Alias Tables = %f\n",elapsedTimeInit_Alias);
  
  PDFX = (G4double*) malloc(NROW*NCOL*sizeof(G4double));
  PDFY = (G4double*) malloc(NROW*NCOL*sizeof(G4double));
  PDFA = (G4int*)    malloc(NROW*NCOL*sizeof(G4int));
  PDFQ = (G4double*) malloc(NROW*NCOL*sizeof(G4double));
  
  for(int i = 0; i < NROW; ++i) {
    for(int j = 0; j < NCOL ; ++j) {
      PDFX[i*NCOL+j] = pdfx[i][j]; 
      PDFY[i*NCOL+j] = pdfy[i][j]; 
      PDFA[i*NCOL+j] = pdfa[i][j]; 
      PDFQ[i*NCOL+j] = pdfq[i][j]; 
    }
  }
  for(int j=0;j<NROW;j++) {
    delete [] pdfx[j];
    delete [] pdfy[j];
    delete [] pdfa[j];
    delete [] pdfq[j];
  }
  delete [] pdfx;
  delete [] pdfy;
  delete [] pdfa;
  delete [] pdfq;
  
  cudaEventRecord (start,0);

  G4double *PDFX_d;
  G4double *PDFY_d;
  G4int    *PDFA_d;
  G4double *PDFQ_d;
    
  cudaMalloc((void**)&PDFX_d, NROW*NCOL*sizeof(G4double));
  cudaMalloc((void**)&PDFY_d, NROW*NCOL*sizeof(G4double));
  cudaMalloc((void**)&PDFA_d, NROW*NCOL*sizeof(G4int));
  cudaMalloc((void**)&PDFQ_d, NROW*NCOL*sizeof(G4double));
  
  cudaMemcpy(PDFX_d, PDFX, NROW*(NCOL)*sizeof(G4double), 
               cudaMemcpyHostToDevice);
  cudaMemcpy(PDFY_d, PDFY, NROW*(NCOL)*sizeof(G4double), 
	     cudaMemcpyHostToDevice);
  cudaMemcpy(PDFA_d, PDFA, NROW*NCOL*sizeof(G4int), 
	     cudaMemcpyHostToDevice);
  cudaMemcpy(PDFQ_d, PDFQ, NROW*(NCOL)*sizeof(G4double), 
	     cudaMemcpyHostToDevice);
  
  cudaEventRecord (stop,0);
  cudaEventSynchronize (stop);
  cudaEventElapsedTime (&elapsedTimeAlloc_Tables,start,stop);
  
  printf("Time to allocate Sampling Tables = %f\n",elapsedTimeAlloc_Tables);

  //Binding texture memory for sampling (PDF) tables 
  SetupSamplingTexture(NROW*NCOL*sizeof(G4double), PDFX_d, PDFY_d);
  
  // get the list of tracks from the stack
  printf("Preparing input tracks at host\n");

  GXHistogramManager* theHisto = GXHistogramManager::Instance();
  theHisto->BookHistograms("moller.root");

  GXTrackHandler *theTrackHandler = new GXTrackHandler();

  //event-loop
  int nevent = 10;

  //allocate memory for host tracks
  int nTracks = 26*192*16;
  if(runType == 0 ) nTracks = 16;
  GXTrack *track_h;
  cudaHostAlloc((void**) &track_h, nTracks*sizeof(GXTrack),
		cudaHostAllocDefault);
  GXTrack *track_c = (GXTrack *) malloc (nTracks*sizeof(GXTrack));
  
  for (int i = 0 ; i < nevent ; i++) {

    printf("Event %d> Transporation in kernel for %d tracks\n",i,nTracks);

    //generate tracks with random position and momentum
    theTrackHandler->GenerateRandomTracks(nTracks,1.0,minP,maxP);

    track_c = theTrackHandler->GetTracks();
    for(size_t i = 0 ; i < nTracks ; i++){
       theTrackHandler->CopyTrack(&track_c[i],&track_h[i]);
    }

    float elapsedTimeGPU_Up = 0.0;
    float elapsedTimeGPU_Down = 0.0;

    //GPU
    //start time for gpu g4 kernel
    cudaEventRecord (start,0);

    GXTrack *track_d;
    cudaMalloc((void**)&track_d, nTracks*sizeof(GXTrack));
    cudaMemcpy(track_d, track_h, nTracks*sizeof(GXTrack), 
               cudaMemcpyHostToDevice);

    //stop time for kernel
    cudaEventRecord (stop,0);
    cudaEventSynchronize (stop);
    cudaEventElapsedTime (&elapsedTimeGPU_Up,start,stop);

    //prepare random engines on the device
    curandState* devStates = 0;
    cudaMalloc(&devStates, theNBlocks*theNThreads* sizeof(curandState));
    curand_setup_gpu(devStates, time(NULL), theNBlocks, theNThreads);

    //GPU kernel
    float elapsedTimeGPU[nKernels];
    printf("MC Model:");
    printf(" Dummy  Random PDF    IPDFT  IPDF2  IPDF2T Alias  Alias2 Geant4 G4Avg\n");
    printf("GPU Time:");
    for(unsigned int k = 0 ; k < nKernels ; ++k) {
      cudaEventRecord (start,0);
      if(cudaEnabled) {
        //-----------------------------------------------------------------
        KernelFunc[k](devStates, track_d, nTracks, physicsTable_d, sbData_d, 
                      PDFX_d, PDFY_d, PDFA_d, PDFQ_d, 
		      theNBlocks, theNThreads);
        //-----------------------------------------------------------------
      }
      cudaEventRecord (stop,0);
      cudaEventSynchronize (stop);
      cudaEventElapsedTime (&elapsedTimeGPU[k],start,stop);
      printf(" %6.3f",elapsedTimeGPU[k]);
    }     
    printf("\n");

    //CPU 
    float elapsedTimeCPU[nKernels];

    printf("CPU Time:");
    for(unsigned int k = 0 ; k < nKernels ; ++k) {
       cudaEventRecord (start,0);
      //--------------------------------------------------------------
      CpuFunc[k](track_c, nTracks, physicsTable, sbData, 
                 PDFX, PDFY, PDFA, PDFQ);
      //--------------------------------------------------------------
      // This does the work on the CPU

      cudaEventRecord (stop,0);
      cudaEventSynchronize (stop);
      cudaEventElapsedTime (&elapsedTimeCPU[k],start,stop);

      printf(" %6.3f",elapsedTimeCPU[k]);
    }   
    printf("\n");

    //Ratio CPU/GPU
    printf("CPU/GPU :");
    for(unsigned int k = 0 ; k < nKernels ; ++k) {
      printf(" %6.3f",elapsedTimeCPU[k]/elapsedTimeGPU[k]);
    }
    printf("\n");

    //Ratio G4/Model GPU
    printf("GPU G4/M:");
    for(unsigned int k = 0 ; k < nKernels ; ++k) {
      printf(" %6.3f",elapsedTimeGPU[8]/elapsedTimeGPU[k]);
    }
    printf("\n");

    //Ratio G4/Model CPU
    printf("CPU G4/M:");
    for(unsigned int k = 0 ; k < nKernels ; ++k) {
      printf(" %6.3f",elapsedTimeCPU[8]/elapsedTimeCPU[k]);
    }
    printf("\n");

    //Validation of MC methods
    cpu_val(track_c, nTracks, physicsTable, sbData, 
	    PDFX, PDFY, PDFA, PDFQ, theHisto);

    //Output (secondary particle) validation on GPU 
    GXTrack *secTracks_d;
    cudaMalloc((void**)&secTracks_d, nTracks*sizeof(GXTrack));
    GXTrack *secTracks_h = (GXTrack*) malloc(nTracks*sizeof(GXTrack));
    GXTrack *secTracks_c = (GXTrack*) malloc(nTracks*sizeof(GXTrack));

    //@@@counters of secondaries
    G4int stackSize = 0;
    //atomic counter for the total number of secondaries
    G4int *stackSize_d;
    cudaMalloc((void**)&stackSize_d, sizeof(G4int));
    cudaMemset(stackSize_d,0,sizeof(G4int));

    if(runType == 2) {
      float elapsedTimeGPUVal[nModels];
      printf("MC Model:");
      printf(" PDF    IPDF2  Alias  Alias2 Geant4\n");
      printf("GPU Val :");
      for(unsigned int k = 0 ; k < nModels ; ++k) {
	cudaMemset(stackSize_d,0,sizeof(G4int));
	cudaEventRecord (start,0);
	if(cudaEnabled) {
	  //-----------------------------------------------------------------
	  GpuVal[k](devStates, track_d, nTracks, physicsTable_d, sbData_d, 
		    PDFX_d, PDFY_d, PDFA_d, PDFQ_d, 
		    stackSize_d, secTracks_d, theNBlocks, theNThreads);
	  //-----------------------------------------------------------------
	}

	cudaEventRecord (stop,0);
	cudaEventSynchronize (stop);
	cudaEventElapsedTime (&elapsedTimeGPUVal[k],start,stop);
	printf(" %6.3f",elapsedTimeGPUVal[k]);
	
	cudaMemcpyAsync(secTracks_h, secTracks_d, nTracks*sizeof(GXTrack),
			cudaMemcpyDeviceToHost);

	for(int i = 0 ; i < nTracks ; ++i) {
	  theHisto->g_energy[k]->Fill(secTracks_h[i].E);
	  G4double px = secTracks_h[i].px;
	  G4double py = secTracks_h[i].py;
	  G4double pz = secTracks_h[i].pz;
	  G4double sint = sqrt(px*px+py*py)/sqrt(px*px+py*py+pz*pz);
	  theHisto->g_angle[k]->Fill(sint);
	}
      }
      printf("\n");

      printf("CPU Val :");
      float elapsedTimeCPUVal[nModels];
      for(unsigned int k = 0 ; k < nModels ; ++k) {
	stackSize = 0;
	cudaEventRecord (start,0);
	//--------------------------------------------------------------
	CpuVal[k](track_c, nTracks, physicsTable, sbData, 
		  PDFX, PDFY, PDFA, PDFQ, &stackSize,secTracks_c);
	//--------------------------------------------------------------
	cudaEventRecord (stop,0);
	cudaEventSynchronize (stop);
	cudaEventElapsedTime (&elapsedTimeCPUVal[k],start,stop);
	
	for(int i = 0 ; i < nTracks ; ++i) {
	  theHisto->c_energy[k]->Fill(secTracks_c[i].E);
	  G4double px = secTracks_c[i].px;
	  G4double py = secTracks_c[i].py;
	  G4double pz = secTracks_c[i].pz;
	  G4double sint = sqrt(px*px+py*py)/sqrt(px*px+py*py+pz*pz);
	  theHisto->c_angle[k]->Fill(sint);
	}
	printf(" %6.3f",elapsedTimeCPUVal[k]);
      }   
      printf("\n");

     //Ratio CPU/GPU
      printf("Val C/G :");
      for(unsigned int k = 0 ; k < nModels ; ++k) {
	printf(" %6.3f",elapsedTimeCPUVal[k]/elapsedTimeGPUVal[k]);
      }
      printf("\n");
 
      //Ratio Model/G4 GPU
      printf("Val G/G4:");
      for(unsigned int k = 0 ; k < nModels ; ++k) {
	printf(" %6.3f",elapsedTimeGPUVal[4]/elapsedTimeGPUVal[k]);
      }
      printf("\n");

      //Ratio Model/G4 CPU
      printf("Val C/G4:");
      for(unsigned int k = 0 ; k < nModels ; ++k) {
	printf(" %6.3f",elapsedTimeCPUVal[4]/elapsedTimeCPUVal[k]);
      }
      printf("\n");
    }

    if(cudaEnabled) {
      cudaFree(track_d);
      cudaFree(devStates);
    }

    cudaFree(secTracks_d);
    cudaFree(stackSize_d);
    free(secTracks_h);
    free(secTracks_c);
  }

  cudaFreeHost(track_h);

  //clean up: free memory on device and host
  free(PDFX);
  free(PDFY);
  free(PDFA);
  free(PDFQ);

  //end of event-loop
  if(cudaEnabled) {
    cudaFree(physicsTable_d);
    cudaFree(sbData_d);
    cudaFree(PDFX_d);
    cudaFree(PDFY_d);
    cudaFree(PDFA_d);
    cudaFree(PDFQ_d);
  }

  free(sbData);
  delete theHisto;
  delete theTrackHandler;
}
