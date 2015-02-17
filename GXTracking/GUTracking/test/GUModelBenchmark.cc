#include <iostream>

#include "GUBenchmarker.h"

using namespace vecphys;

int main(int argc, char* argv[]) {

  //default run
  int ntracks = 4992;
  int nrepetitions = 100;

  if(argc >= 2) ntracks = atoi(argv[1]);
  if(argc >= 3) nrepetitions = atoi(argv[2]);

  GUBenchmarker tester;
  tester.SetNTracks(ntracks);
  tester.SetRepetitions(nrepetitions);
  int status = tester.RunBenchmark();

  if(status==1) std::cout << "RunBenchmark Failed" << std::endl;
  return status;
}
