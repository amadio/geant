#include <iostream>

#include "GUBenchmarker.h"
#include "base/SystemOfUnits.h"
#include "SamplingMethod.h"
#include "GUPhysicsModelName.h"

using namespace vecphys;

int main(int argc, char* argv[])
{
  //default run
  int ntracks = 4992;
  int nrepetitions = 100;
  double minEnergy =  500.*MeV;
  double maxEnergy =  minEnergy; 
  SamplingMethod sampleType = SamplingMethod::kAlias;
  int emModel = GUPhysicsModelIndex::kNullModel ; //all models
  int materialMode = 0; //32 selected elements
  
  if(argc >= 2) ntracks =      atoi(argv[1]);
  if(argc >= 3) nrepetitions = atoi(argv[2]);
  if(argc >= 4) {
     minEnergy  =  atof(argv[3]);
     std::cout << "  Min energy (MeV) = " << minEnergy << std::endl;

     maxEnergy = minEnergy;  // Assume mono-energetic if no max is defined
     if(argc >= 5) {
        maxEnergy  =  atof(argv[4]);
        std::cout << "  Max energy (MeV) = " << maxEnergy << std::endl;
     } else {
        std::cout << "  Mono-energetic> max energy (MeV) = " << maxEnergy << std::endl;
     }
  } else {
     std::cout << " Using defaults:  min energy (MeV) = " << minEnergy << std::endl;
     std::cout << "                  max energy (MeV) = " << maxEnergy << std::endl;
  }

  if(argc >= 6)  {
    int type = atoi(argv[5]);   
    if(type==0) sampleType = SamplingMethod::kAlias;
    if(type==1) sampleType = SamplingMethod::kRejection;
    if(type==2) sampleType = SamplingMethod::kUnpack;
    if(sampleType>=0 && sampleType < SamplingMethod::kNumberSamplingMethod ) {
      std::cout << "  Using sampling method = " << sampleType << std::endl;
    }
    else {
      std::cout << "  Illegal sampling method " << sampleType << std::endl;
      exit(0);
    }
  } else {
      std::cout << "  Using the default sampling method = " << sampleType << std::endl;
  }

  if(argc >= 7)  {
    emModel = atoi(argv[6]);   
    if(emModel == -1 ) {
      std::cout << "  Validation for all available vector EM physics models" << std::endl;
    }
    else if(emModel >= 0 && emModel < kNumberPhysicsModel ) {
      std::cout << "  Validation for vector EM physics model  = " 
		<< GUPhysicsModelName[emModel] << std::endl;
    }
    else {
      std::cout << "  Illegal vector physics model " << emModel 
		<< "! Should be [-1:" << kNumberPhysicsModel-1 << "]" << std::endl;
      exit(0);
    }
  } else {
      std::cout << "  Validation for all available vector EM physics models" << std::endl;
  }

  if(argc >= 8)  {
    materialMode = atoi(argv[7]);   
    if(materialMode < 0 || materialMode > 1) {
      std::cout << "  Illegal material mode " << materialMode 
		<< "! Should be [0:1] (see " << std::endl;
      exit(0);
    }
  }
  else {
     std::cout << "  Validation with 16 selected elements" << std::endl;
  }
  //todo: refine setting arguments with a bash style or with a file

  GUBenchmarker tester;
  tester.SetNTracks(ntracks);
  tester.SetRepetitions(nrepetitions);

  // P = E for gammas (only)
  tester.SetMinP( minEnergy );
  tester.SetMaxP( maxEnergy );
  tester.SetSampleType( sampleType );
  tester.SetEmModel( emModel );
  tester.SetMaterialMode( materialMode );
  
  int status = tester.RunBenchmark();

  if(status==1) std::cout << "RunBenchmark Failed" << std::endl;
  return status;
}
