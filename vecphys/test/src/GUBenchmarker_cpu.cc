#include "base/Stopwatch.h"
using vecgeom::Stopwatch;
#include "GUHistogram.h"

#include "ComptonKleinNishina.h"
#include "ConversionBetheHeitler.h"
#include "PhotoElectronSauterGavrila.h"
#include "IonisationMoller.h"
#include "BremSeltzerBerger.h"
#include "SamplingMethod.h"

namespace vecphys {

// Scalar

Real_t ScalarKleinNishina(int ntracks, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos,
                          SamplingMethod sampleType)
{
  static vecphys::cxx::ComptonKleinNishina model(0, -1);
  model.SetSamplingMethod(sampleType);
  model.Initialization();

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    //    printf("this elec[%d] = %d\n",i,targetElements[i]);
    model.Interact<ScalarBackend>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  // validation for the total cross section
  /*
  double sigma = 0;
  for(int i = 0 ; i < ntracks ; ++i) {
    model.AtomicCrossSection<ScalarBackend>(itrack_aos[i], targetElements[i],sigma);
  }
  */
  return elapsedTime;
}

Real_t ScalarHybridCompton(int ntracks, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos,
                           SamplingMethod sampleType)
{
  static vecphys::cxx::ComptonKleinNishina model(0, -1);
  model.SetSamplingMethod(sampleType);
  model.Initialization();

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    model.ModelInteract<ScalarBackend>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t ScalarBetheHeitler(int ntracks, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos,
                          SamplingMethod sampleType)
{
  static vecphys::cxx::ConversionBetheHeitler model(0, -1);
  model.SetSamplingMethod(sampleType);
  model.Initialization();

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    model.Interact<ScalarBackend>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t ScalarSauterGavrila(int ntracks, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos,
                           SamplingMethod sampleType)
{
  static vecphys::cxx::PhotoElectronSauterGavrila model(0, -1);
  model.SetSamplingMethod(sampleType);
  model.Initialization();

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    model.ModelInteract<ScalarBackend>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t ScalarMollerBhabha(int ntracks, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos,
                          SamplingMethod sampleType)
{
  static vecphys::cxx::IonisationMoller model(0, -1);
  model.SetSamplingMethod(sampleType);
  model.Initialization();

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    model.Interact<ScalarBackend>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t ScalarSeltzerBerger(int ntracks, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos,
                           SamplingMethod sampleType)
{
  static vecphys::cxx::BremSeltzerBerger model(0, -1);
  model.SetSamplingMethod(sampleType);
  model.Initialization();

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    model.Interact<ScalarBackend>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

// Vector

Real_t VectorKleinNishina(GUTrack_v &itrack_soa, int *targetElements, GUTrack_v &otrack_soa, SamplingMethod sampleType)
{
  static vecphys::cxx::ComptonKleinNishina model(0, -1);
  model.SetSamplingMethod(sampleType);
  model.Initialization();

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  //  int ntracks = itrack_soa.numTracks;

  if (sampleType == SamplingMethod::kUnpack) {
    timer.Start();

    model.InteractUnpack<VectorBackend>(itrack_soa, targetElements, otrack_soa);

    elapsedTime = timer.Stop();
  } else {
    timer.Start();

    model.Interact<VectorBackend>(itrack_soa, targetElements, otrack_soa);

    elapsedTime = timer.Stop();
  }

  // validation for the total cross section
  /*
  double* sigma  = new double [ntracks];
  for(int i = 0 ; i < ntracks ; ++i) {
    sigma[i] = 0;
  }

  model.AtomicCrossSection<VectorBackend>(itrack_soa, targetElements,sigma);

  delete [] sigma;
  */
  return elapsedTime;
}

Real_t VectorHybridCompton(GUTrack_v &itrack_soa, int *targetElements, GUTrack_v &otrack_soa, SamplingMethod sampleType)
{
  static vecphys::cxx::ComptonKleinNishina model(0, -1);
  model.SetSamplingMethod(sampleType);
  model.Initialization();

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  if (sampleType == SamplingMethod::kUnpack) {
    timer.Start();

    model.InteractUnpack<VectorBackend>(itrack_soa, targetElements, otrack_soa);

    elapsedTime = timer.Stop();
  } else {
    timer.Start();

    model.ModelInteract<VectorBackend>(itrack_soa, targetElements, otrack_soa);

    elapsedTime = timer.Stop();
  }

  return elapsedTime;
}

Real_t VectorBetheHeitler(GUTrack_v &itrack_soa, int *targetElements, GUTrack_v &otrack_soa, SamplingMethod sampleType)
{
  static vecphys::cxx::ConversionBetheHeitler model(0, -1);
  model.SetSamplingMethod(sampleType);
  model.Initialization();

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  model.Interact<VectorBackend>(itrack_soa, targetElements, otrack_soa);

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t VectorSauterGavrila(GUTrack_v &itrack_soa, int *targetElements, GUTrack_v &otrack_soa, SamplingMethod sampleType)

{
  static vecphys::cxx::PhotoElectronSauterGavrila model(0, -1);
  model.SetSamplingMethod(sampleType);
  model.Initialization();

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  model.ModelInteract<VectorBackend>(itrack_soa, targetElements, otrack_soa);

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t VectorMollerBhabha(GUTrack_v &itrack_soa, int *targetElements, GUTrack_v &otrack_soa, SamplingMethod sampleType)
{
  static vecphys::cxx::IonisationMoller model(0, -1);
  model.SetSamplingMethod(sampleType);
  model.Initialization();

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  model.Interact<VectorBackend>(itrack_soa, targetElements, otrack_soa);

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t VectorSeltzerBerger(GUTrack_v &itrack_soa, int *targetElements, GUTrack_v &otrack_soa, SamplingMethod sampleType)
{
  static vecphys::cxx::BremSeltzerBerger model(0, -1);
  model.SetSamplingMethod(sampleType);
  model.Initialization();

  Real_t elapsedTime = 0.0;
  static Stopwatch timer;

  timer.Start();

  model.Interact<VectorBackend>(itrack_soa, targetElements, otrack_soa);

  elapsedTime = timer.Stop();

  return elapsedTime;
}

// Geant4 composition and rejection

Real_t G4KleinNishina(int ntracks, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos)
{
  static vecphys::cxx::ComptonKleinNishina model(0, -1);

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    model.InteractG4<ScalarBackend>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  // validation for the total cross section
  /*
  double sigma = 0.0;
  for(int i = 0 ; i < ntracks ; ++i) {
    model.AtomicCrossSectionG4<ScalarBackend>(itrack_aos[i], targetElements[i], sigma);
  }
  */
  return elapsedTime;
}

Real_t G4HybridCompton(int ntracks, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos)
{
  static vecphys::cxx::ComptonKleinNishina model(0, -1);

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    model.InteractG4<ScalarBackend>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t G4BetheHeitler(int ntracks, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos)
{
  static vecphys::cxx::ConversionBetheHeitler model(0, -1);

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    model.InteractG4<ScalarBackend>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t G4SauterGavrila(int ntracks, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos)
{
  static vecphys::cxx::PhotoElectronSauterGavrila model(0, -1);

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    model.InteractG4<ScalarBackend>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t G4MollerBhabha(int ntracks, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos)
{
  static vecphys::cxx::IonisationMoller model(0, -1);

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    model.InteractG4<ScalarBackend>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Real_t G4SeltzerBerger(int ntracks, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos)
{
  static vecphys::cxx::BremSeltzerBerger model(0, -2);

  // fThreadId = -2 is used not to build the alias table for testing the Geant4
  // composition and rejection method. This convention should be temporary
  // and is only used for the purpose of validation or other testings.
  // see BremSeltzerBerger::Initialization().

  static Stopwatch timer;
  Real_t elapsedTime = 0.0;

  timer.Start();

  for (int i = 0; i < ntracks; ++i) {
    model.InteractG4<ScalarBackend>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

} // end namespace vecphys
