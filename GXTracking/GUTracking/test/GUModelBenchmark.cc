#include <iostream>

#include "GUBenchmarker.h"

using namespace vecphys;

int main(int argc, char* argv[]) {

  //default run
  int ntracks = 4992;
  int nrepetitions = 100;
  constexpr double MeV = 0.001;
  double energy = 0.5 * MeV;
  
  if(argc >= 2) ntracks =      atoi(argv[1]);
  if(argc >= 3) nrepetitions = atoi(argv[2]);
  if(argc >= 4) {
     energy  =  atof(argv[2]);  
     std::cout << " Main: argument energy = " << energy << std::endl;
  } else {
     std::cout << " Main: default  energy = " << energy << std::endl;
  }

  GUBenchmarker tester;
  tester.SetNTracks(ntracks);
  tester.SetRepetitions(nrepetitions);
  int status = tester.RunBenchmark();

  if(status==1) std::cout << "RunBenchmark Failed" << std::endl;
  return status;
}
