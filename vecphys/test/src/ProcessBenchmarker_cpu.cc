#include "GUHistogram.h"
#include "base/Stopwatch.h"
using vecgeom::Stopwatch;

//#include "GUAuxFunctions.h"

#include "ElectronProcess.h"
#include "PhotonProcess.h"

namespace vecphys {

// Scalar

Real_t ScalarPhotonProcess(int ntracks, GUTrack *itrack_aos, int *materialIndex)
{
  static vecphys::cxx::PhotonProcess process(0, -1);

  static Stopwatch timer;
  Real_t elapsedTime = 0.;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    process.GetStepLengthAndProcess<ScalarBackend>(itrack_aos[i], materialIndex[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t ScalarElectronProcess(int ntracks, GUTrack *itrack_aos, int *materialIndex)
{
  //  dummy for now
  static vecphys::cxx::ElectronProcess process(0, -1);

  static Stopwatch timer;
  Real_t elapsedTime = 0.;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    process.GetStepLengthAndProcess<ScalarBackend>(itrack_aos[i], materialIndex[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

// Vector

Real_t VectorPhotonProcess(GUTrack_v &itrack_soa, int *materialIndex)
{
  static vecphys::cxx::PhotonProcess process(0, -1);

  static Stopwatch timer;
  Real_t elapsedTime = 0.;

  timer.Start();

  process.GetStepLengthAndProcess<VectorBackend>(itrack_soa, materialIndex);

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t VectorElectronProcess(GUTrack_v &itrack_soa, int *materialIndex)
{
  static vecphys::cxx::ElectronProcess process(0, -1);

  static Stopwatch timer;
  Real_t elapsedTime = 0.;

  timer.Start();

  process.GetStepLengthAndProcess<VectorBackend>(itrack_soa, materialIndex);

  elapsedTime = timer.Stop();

  return elapsedTime;
}

// a.k.a Geant3

Real_t Geant3PhotonProcess(int ntracks, GUTrack *itrack_aos, int *materialIndex)
{
  static vecphys::cxx::PhotonProcess process(0, -1);

  static Stopwatch timer;
  Real_t elapsedTime = 0.;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    process.G3StepLengthAndProcess<ScalarBackend>(itrack_aos[i], materialIndex[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t Geant3ElectronProcess(int ntracks, GUTrack *itrack_aos, int *materialIndex)
{
  static vecphys::cxx::ElectronProcess process(0, -1);

  static Stopwatch timer;
  Real_t elapsedTime = 0.;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    process.G3StepLengthAndProcess<ScalarBackend>(itrack_aos[i], materialIndex[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

// GeantV - Hybrid

Real_t GeantVPhotonProcess(GUTrack_v &itrack_soa, int *materialIndex)
{
  static vecphys::cxx::PhotonProcess process(0, -1);

  static Stopwatch timer;
  Real_t elapsedTime = 0.;

  timer.Start();

  process.GVStepLengthAndProcess<VectorBackend>(itrack_soa, materialIndex);

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t GeantVElectronProcess(GUTrack_v &itrack_soa, int *materialIndex)
{
  static vecphys::cxx::ElectronProcess process(0, -1);

  static Stopwatch timer;
  Real_t elapsedTime = 0.;

  timer.Start();

  process.GVStepLengthAndProcess<VectorBackend>(itrack_soa, materialIndex);

  elapsedTime = timer.Stop();

  return elapsedTime;
}

} // end namespace vecphys
