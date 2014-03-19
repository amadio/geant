#include "stdio.h"

#include "random_kernel.h"

#include "GXRunManager.h"
#include "GXTrackHandler.h"

//Application
#include "GXSimpleEcal.h"
#include "GXGPUManager.h"
#include "GXMICManager.h"

using namespace std;

int main (int argc, char* argv[]) 
{
  int NBlocks   =  32;
  int NThreads  = 128;
  if(argc >= 2) NBlocks  = atoi(argv[1]);
  if(argc >= 3) NThreads = atoi(argv[2]);

  GXRunManager* runManager = GXRunManager::Instance();

  //mandatory - Geometry
  runManager->ConstructGeometry(new GXSimpleEcal());

  //mandatory - CoprocessorManager
  runManager->RegisterCoprocessorManager(new GXGPUManager(NBlocks,NThreads));

  //mandatory - Initialization
  runManager->Initialization();

  //optinal - provides any external scheduler
  GXGPUManager* gpuManager = 
    dynamic_cast<GXGPUManager*> (runManager->GetCoprocessorManager());
  GXTrackHandler *theHandler = gpuManager->GetTrackHandler();

  //artificial task-loop: should be replaced by a set of concurrent tasks 
  //                      processing a basket of tracks

  int ntasks = 10;

  for (int i = 0 ; i < ntasks ; ++i) {

    size_t nTracks = 4096*16;

    // populate tracks with track position, momentum, energy, steplength
    theHandler->GenerateRandomTracks(nTracks);
    printf("Task %d> N(e-,gamma)=(%d, %d)\n",i,
	  theHandler->GetNumberOfElectrons(),theHandler->GetNumberOfPhotons());

    //H2D
    gpuManager->UploadTaskData();

    //GPU Task
    gpuManager->LaunchTask();

    //D2H
    gpuManager->DownloadTaskData();

    //CPU Task
    gpuManager->LaunchCPUTask();

    //deallocate memory of task data
    gpuManager->DeallocateTaskData();

    //Performance summary
    gpuManager->PrintPerformance(i);

  }
  //end of event-loop

  gpuManager->DeallocateDeviceMemory();
  delete theHandler;
}
