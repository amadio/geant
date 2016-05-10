#include <iostream>

#include "ProcessBenchmarker.h"
#include "base/SystemOfUnits.h"
#include "SamplingMethod.h"
#include "GUPhysicsProcessIndex.h"
//#include "GUPhysicsProcessName.h"

using namespace vecphys;

int main(int argc, char* argv[])
{
  //default run
  int ntracks = 4992;
  int nrepetitions = 100;
  double minEnergy =  500.*MeV;
  double maxEnergy =  minEnergy; 
  int emProcess = GUPhysicsProcessIndex::kNullProcess ; //all models
  int materialMode = 0; //selected materials
  int runMode = -1; 
  
  if(argc >= 2) ntracks =      atoi(argv[1]);
  if(argc >= 3) nrepetitions = atoi(argv[2]);
  if(argc >= 4) {
    minEnergy  =  atof(argv[3]) * MeV;
    maxEnergy = minEnergy;  // Assume mono-energetic if no max is defined
    if(argc >= 5) {
       maxEnergy  =  atof(argv[4]) * MeV;
    }
  }

  if(argc >= 6)  {
    emProcess = atoi(argv[5]);   
    if(emProcess < -1 || emProcess >= kNumberPhysicsProcess ) {
      std::cout << "  Illegal vector physics process " << emProcess 
		<< "! Should be [-1:" << kNumberPhysicsProcess-1 << "]" << std::endl;
      exit(0);
    }
  }

  if(argc >= 7)  {
    runMode = atoi(argv[6]);   
    if(runMode < -1 || runMode >= 4) {
      std::cout << "  Illegal run mode " << runMode 
		<< "! Should be [-1:" << 3 << "] [-1,0,1,2,3]=[all,Scalar,Vector,Geant3,GeantV]" << std::endl;
      exit(0);
    }
  }

  if(argc >= 8)  {
    materialMode = atoi(argv[7]);   
    if(materialMode < 0 || materialMode > 1) {
      std::cout << "  Illegal material mode " << materialMode 
		<< "! Should be [0:1] (see " << std::endl;
      exit(0);
    }
  }

  ProcessBenchmarker tester;
  tester.SetNTracks(ntracks);
  tester.SetRepetitions(nrepetitions);

  tester.SetMinP( minEnergy );
  tester.SetMaxP( maxEnergy );
  tester.SetEmProcess( emProcess );
  tester.SetMaterialMode( materialMode );
  tester.SetRunMode( runMode );
  
  int status = tester.RunBenchmark();

  if(status==1) std::cout << "RunProcessBenchmark Failed" << std::endl;
  return status;
}
