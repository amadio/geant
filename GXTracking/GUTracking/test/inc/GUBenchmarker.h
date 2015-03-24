#ifndef GUBenchmarker_H
#define GUBenchmarker_H 1

#include "base/Global.h"

namespace vecphys {

class GUTrackHandler;

class GUBenchmarker {

public:

  GUBenchmarker();
  ~GUBenchmarker();

  int RunBenchmark();

  void SetNTracks(const int ntracks) { fNtracks = ntracks; }
  void SetRepetitions(const unsigned repetitions) { 
    fRepetitions = repetitions; 
  }

private:
    
  int  RunBenchmarkInteract();

  void PrepareTargetElements(int *targetElements, int ntracks);
  
  void RunGeant4();
  void RunScalar();
  void RunVector();
  void CheckTimer();
  void CheckRandom();

#ifdef VECPHYS_CUDA
  void RunCuda();
#endif

protected:
  void
  ReportInteraction( double incomingEnergy, int    targetElement,
                     double GammaOut_E,     double GammaOut_Pz,
                     double electron_En,    double electron_Pz,
                     bool   print_Uz= false
     );
  // Print results - one per line, to enable debugging

private:
  GUTrackHandler *fTrackHandler;

  int fNtracks;
  unsigned fRepetitions;
  int fVerbosity;
};

} // end namespace vecphys

#endif
