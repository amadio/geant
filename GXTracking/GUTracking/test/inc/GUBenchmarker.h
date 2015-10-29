#ifndef GUBenchmarker_H
#define GUBenchmarker_H 1

#include "base/Global.h"
#include "SamplingMethod.h"

namespace vecphys {

  //class GUTrack;
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

  void SetMinP(double pMin) { fMinP= pMin; }
  void SetMaxP(double pMax) { fMaxP= pMax; }
  void SetSampleType(SamplingMethod type) { fSampleType = type ; }

  void SetMonoEnergeticBeam(double E){ SetMinP(E); SetMaxP(E); } // For gamma only now!

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
  ///phihome/syjun/devel/sbTest/test/src/GUBenchmarker.cc:1:  GUComptonKleinNishina *fmodel;

  GUTrackHandler *fTrackHandler;

  int fNtracks;
  unsigned fRepetitions;
  int fVerbosity;

  double fMinP, fMaxP;  // Minimum and Maximum momentum of primaries
  SamplingMethod fSampleType;
};

} // end namespace vecphys

#endif
