#include "ProcessBenchmarker.h"
#include "ProcessBenchmarker_cpu.h"

#include "GUPhysicsProcessIndex.h"
#include "GUPhysicsProcessName.h"

#include "base/Stopwatch.h"
using vecgeom::Stopwatch;
#include "base/SystemOfUnits.h"

#include "GUTrackHandler.h"

#ifdef VECPHYS_ROOT
#include "GUHistogram.h"
#endif

namespace vecphys {

// Mono-energetic test
constexpr double defaultMinP = 500. * MeV;
constexpr double defaultMaxP = defaultMinP;

ProcessBenchmarker::ProcessBenchmarker()
    : fNtracks(4992), fRepetitions(1), fVerbosity(1), fMinP(defaultMinP), fMaxP(defaultMaxP),
      fEmProcess(GUPhysicsProcessIndex::kNullProcess), fMaterialMode(0), fRunMode(-1)
{
  fTrackHandler    = new GUTrackHandler();
  fMaterialHandler = vecphys::MaterialHandler::Instance();
  fMaterialHandler->BuildMaterialTable();
}

ProcessBenchmarker::~ProcessBenchmarker()
{
  delete fTrackHandler;
}

int ProcessBenchmarker::RunBenchmark()
{
  printf(" Run ProcessBenchmarker with arguments: %d %d %f %f %d %d %d\n", fNtracks, fRepetitions, fMinP, fMaxP,
         fEmProcess, fRunMode, fMaterialMode);
  printf(" Ntracks       = %d\n", fNtracks);
  printf(" NRepetitions  = %d\n", fRepetitions);
  printf(" MinP (MeV)    = %f\n", fMinP);
  printf(" MaxP (MeV)    = %f\n", fMaxP);
  if (fEmProcess == -1) {
    printf(" EM Process    = -1 ( Using all available vector EM physics processes)\n");
  } else {
    printf(" EM Process    =  %d ( Process = %s )\n", fEmProcess, GUPhysicsProcessName[fEmProcess]);
  }

  if (fRunMode == -1) {
    printf(" Run Mode      = -1 ( Using all available test modes)\n");
  } else {
    printf(" Run Mode      =  %d ( [-1,0,1,2,3,4]=[all,Scalar,Vector,Geant3,GeantV,Cuda] )\n", fRunMode);
  }

  int errorcode = 0;
  errorcode += RunBenchmarkProcess();
  return (errorcode) ? 1 : 0;
}

int ProcessBenchmarker::RunBenchmarkProcess()
{
  if (fVerbosity > 0) {
  }

  int mismatches = 0;

  if (fRunMode == 0 || fRunMode == -1) RunScalar();
  if (fRunMode == 1 || fRunMode == -1) RunVector();
  if (fRunMode == 2 || fRunMode == -1) RunGeant3();
  if (fRunMode == 3 || fRunMode == -1) RunGeantV();
  if (fRunMode == 4 || fRunMode == -1) {
#ifdef VECPHYS_CUDA
    RunCuda();
#endif
  }

  return mismatches;
}

void ProcessBenchmarker::RunScalar()
{
#ifdef VECPHYS_ROOT
  GUHistogram *histogram = new GUHistogram("process_scalar.root", fMaxP);
#endif

  int *targetElements = new int[fNtracks];
  GUTrack *itrack_aos = (GUTrack *)malloc(fNtracks * sizeof(GUTrack));

  Real_t elapsedTotal[kNumberPhysicsProcess];
  Real_t elapsedT[kNumberPhysicsProcess];
  for (int k        = 0; k < kNumberPhysicsProcess; ++k)
    elapsedTotal[k] = 0.;

  for (unsigned r = 0; r < fRepetitions; ++r) {
    fMaterialHandler->PrepareMaterialIndex(targetElements, fNtracks, fMaterialMode);
    fTrackHandler->GenerateRandomTracks(fNtracks, fMinP, fMaxP);

    GUTrack *track_aos = fTrackHandler->GetAoSTracks();
    fTrackHandler->SortAoSTracksByEnergy(track_aos, fNtracks);

    for (int k = 0; k < kNumberPhysicsProcess; ++k) {
      if (fEmProcess == GUPhysicsProcessIndex::kNullProcess || fEmProcess == k) {
        fTrackHandler->CopyAoSTracks(track_aos, itrack_aos, fNtracks);
        elapsedT[k] = 0.0;
        elapsedT[k] = ScalarKernelFunc[k](fNtracks, itrack_aos, targetElements);
        elapsedTotal[k] += elapsedT[k];

#ifdef VECPHYS_ROOT
        histogram->RecordTime(k, elapsedT[k]);
        for (int i = 0; i < fNtracks; ++i) {
          histogram->RecordHistosProc(k, itrack_aos[i].E, itrack_aos[i].nint, itrack_aos[i].s, itrack_aos[i].lambda);
        }
#endif
      }
    }
  }

  for (int k = 0; k < kNumberPhysicsProcess; ++k) {
    if (fEmProcess == GUPhysicsProcessIndex::kNullProcess || fEmProcess == k) {
      printf("%s  Scalar Total time of %3d reps = %6.3f sec\n", GUPhysicsProcessName[k], fRepetitions, elapsedTotal[k]);
    }
  }

  free(itrack_aos);
  delete[] targetElements;

#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}

void ProcessBenchmarker::RunVector()
{
#ifdef VECPHYS_ROOT
  GUHistogram *histogram = new GUHistogram("process_vector.root", fMaxP);
#endif
  // SOA tracks
  GUTrackHandler *handler_in = new GUTrackHandler(fNtracks);
  GUTrack_v itrack_soa       = handler_in->GetSoATracks();

  int *targetElements = new int[fNtracks];

  Real_t elapsedTotal[kNumberPhysicsProcess];
  Real_t elapsedT[kNumberPhysicsProcess];

  for (int k        = 0; k < kNumberPhysicsProcess; ++k)
    elapsedTotal[k] = 0.;

  for (unsigned r = 0; r < fRepetitions; ++r) {

    fMaterialHandler->PrepareMaterialIndex(targetElements, fNtracks, fMaterialMode);
    fTrackHandler->GenerateRandomTracks(fNtracks, fMinP, fMaxP);

    GUTrack_v &track_soa = fTrackHandler->GetSoATracks();
    fTrackHandler->SortSoATracksByEnergy(track_soa, fNtracks);

    for (int k = 0; k < kNumberPhysicsProcess; ++k) {
      if (fEmProcess == GUPhysicsProcessIndex::kNullProcess || fEmProcess == k) {
        fTrackHandler->CopySoATracks(track_soa, itrack_soa, fNtracks);

        elapsedT[k] = 0.0;
        elapsedT[k] = VectorKernelFunc[k](itrack_soa, targetElements);
        elapsedTotal[k] += elapsedT[k];

#ifdef VECPHYS_ROOT
        histogram->RecordTime(k, elapsedT[k]);
        for (int i = 0; i < fNtracks; ++i) {
          histogram->RecordHistosProc(k, itrack_soa.E[i], itrack_soa.nint[i], itrack_soa.s[i], itrack_soa.lambda[i]);
        }
#endif
      }
    }
  }

  for (int k = 0; k < kNumberPhysicsProcess; ++k) {
    if (fEmProcess == GUPhysicsProcessIndex::kNullProcess || fEmProcess == k) {
      printf("%s  Vector Total time of %3d reps = %6.3f sec\n", GUPhysicsProcessName[k], fRepetitions, elapsedTotal[k]);
    }
  }
  delete handler_in;
  delete[] targetElements;

#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}

void ProcessBenchmarker::RunGeantV()
{
#ifdef VECPHYS_ROOT
  GUHistogram *histogram = new GUHistogram("process_geantv.root", fMaxP);
#endif
  // SOA tracks
  GUTrackHandler *handler_in = new GUTrackHandler(fNtracks);
  GUTrack_v itrack_soa       = handler_in->GetSoATracks();

  int *targetElements = new int[fNtracks];

  Real_t elapsedTotal[kNumberPhysicsProcess];
  Real_t elapsedT[kNumberPhysicsProcess];

  for (int k        = 0; k < kNumberPhysicsProcess; ++k)
    elapsedTotal[k] = 0.;

  for (unsigned r = 0; r < fRepetitions; ++r) {

    fMaterialHandler->PrepareMaterialIndex(targetElements, fNtracks, fMaterialMode);
    fTrackHandler->GenerateRandomTracks(fNtracks, fMinP, fMaxP);

    GUTrack_v &track_soa = fTrackHandler->GetSoATracks();
    fTrackHandler->SortSoATracksByEnergy(track_soa, fNtracks);

    for (int k = 0; k < kNumberPhysicsProcess; ++k) {
      if (fEmProcess == GUPhysicsProcessIndex::kNullProcess || fEmProcess == k) {
        fTrackHandler->CopySoATracks(track_soa, itrack_soa, fNtracks);

        elapsedT[k] = 0.0;
        elapsedT[k] = GeantVKernelFunc[k](itrack_soa, targetElements);
        elapsedTotal[k] += elapsedT[k];

#ifdef VECPHYS_ROOT
        histogram->RecordTime(k, elapsedT[k]);
        for (int i = 0; i < fNtracks; ++i) {
          histogram->RecordHistosProc(k, itrack_soa.E[i], itrack_soa.nint[i], itrack_soa.s[i], itrack_soa.lambda[i]);
        }
#endif
      }
    }
  }

  for (int k = 0; k < kNumberPhysicsProcess; ++k) {
    if (fEmProcess == GUPhysicsProcessIndex::kNullProcess || fEmProcess == k) {
      printf("%s  GeantV Total time of %3d reps = %6.3f sec\n", GUPhysicsProcessName[k], fRepetitions, elapsedTotal[k]);
    }
  }
  delete handler_in;
  delete[] targetElements;

#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}

void ProcessBenchmarker::RunGeant3()
{
#ifdef VECPHYS_ROOT
  GUHistogram *histogram = new GUHistogram("process_geant3.root", fMaxP);
#endif
  int *targetElements = new int[fNtracks];
  GUTrack *itrack_aos = (GUTrack *)malloc(fNtracks * sizeof(GUTrack));

  Real_t elapsedTotal[kNumberPhysicsProcess];
  Real_t elapsedT[kNumberPhysicsProcess];
  for (int k        = 0; k < kNumberPhysicsProcess; ++k)
    elapsedTotal[k] = 0.;

  for (unsigned r = 0; r < fRepetitions; ++r) {
    fMaterialHandler->PrepareMaterialIndex(targetElements, fNtracks, fMaterialMode);
    fTrackHandler->GenerateRandomTracks(fNtracks, fMinP, fMaxP);

    GUTrack *track_aos = fTrackHandler->GetAoSTracks();
    fTrackHandler->SortAoSTracksByEnergy(track_aos, fNtracks);

    for (int k = 0; k < kNumberPhysicsProcess; ++k) {
      if (fEmProcess == GUPhysicsProcessIndex::kNullProcess || fEmProcess == k) {
        fTrackHandler->CopyAoSTracks(track_aos, itrack_aos, fNtracks);
        elapsedT[k] = 0.0;
        elapsedT[k] = Geant3KernelFunc[k](fNtracks, itrack_aos, targetElements);
        elapsedTotal[k] += elapsedT[k];

#ifdef VECPHYS_ROOT
        histogram->RecordTime(k, elapsedT[k]);
        for (int i = 0; i < fNtracks; ++i) {
          histogram->RecordHistosProc(k, itrack_aos[i].E, itrack_aos[i].nint, itrack_aos[i].s, itrack_aos[i].lambda);
        }
#endif
      }
    }
  }

  for (int k = 0; k < kNumberPhysicsProcess; ++k) {
    if (fEmProcess == GUPhysicsProcessIndex::kNullProcess || fEmProcess == k) {
      printf("%s  a.k.a-G3 Total time of %3d reps = %6.3f sec\n", GUPhysicsProcessName[k], fRepetitions,
             elapsedTotal[k]);
    }
  }

  free(itrack_aos);
  delete[] targetElements;

#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}

} // end namespace vecphys
