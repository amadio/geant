#include <iostream>

#include "GUBenchmarker.h"
#include "base/SystemOfUnits.h"


using namespace vecphys;

int main(int argc, char* argv[])
{

  //default run
  int ntracks = 4992;
  int nrepetitions = 100;

  // double minEnergy =  10.0 * KeV;
  // double maxEnergy = 500.0 * MeV;

  double minEnergy =  500.*MeV;
  double maxEnergy =  minEnergy; 
  
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

  GUBenchmarker tester;
  tester.SetNTracks(ntracks);
  tester.SetRepetitions(nrepetitions);

  // P = E for gammas (only)
  tester.SetMinP( minEnergy );
  tester.SetMaxP( maxEnergy );
  
  int status = tester.RunBenchmark();

  if(status==1) std::cout << "RunBenchmark Failed" << std::endl;
  return status;
}
