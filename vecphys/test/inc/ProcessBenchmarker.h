#ifndef ProcessBenchmarker_H
#define ProcessBenchmarker_H 1

#include "MaterialHandler.h"
#include "SamplingMethod.h"
#include "base/VecPhys.h"

namespace vecphys {

class GUTrackHandler;

class ProcessBenchmarker {

public:
  ProcessBenchmarker();
  ~ProcessBenchmarker();

  int RunBenchmark();

  void SetNTracks(const int ntracks) { fNtracks = ntracks; }
  void SetRepetitions(const unsigned repetitions) { fRepetitions = repetitions; }

  void SetMinP(double pMin) { fMinP = pMin; }
  void SetMaxP(double pMax) { fMaxP = pMax; }
  void SetEmProcess(int process) { fEmProcess = process; }
  void SetMaterialMode(int materialMode) { fMaterialMode = materialMode; }
  void SetRunMode(int runMode) { fRunMode = runMode; }

private:
  int RunBenchmarkProcess();
  void PrepareTargetElements(int *targetElements, int ntracks);

  void RunScalar();
  void RunVector();
  void RunGeant3();
  void RunGeantV();

private:
  GUTrackHandler *fTrackHandler;
  MaterialHandler *fMaterialHandler;

  int fNtracks;
  unsigned fRepetitions;
  int fVerbosity;

  double fMinP;
  double fMaxP;
  int fEmProcess;
  int fMaterialMode;
  int fRunMode;
};

} // end namespace vecphys

#endif
