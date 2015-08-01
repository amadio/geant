#include "base/Stopwatch.h"
#include "GUHistogram.h"

#include "ComptonKleinNishina.h"
#include "ConversionBetheHeitler.h"
#include "PhotoElectronSauterGavrila.h"
#include "IonisationMoller.h"
#include "BremSeltzerBerger.h"

namespace vecphys {

// Scalar

Precision ScalarKleinNishina(int ntracks, 
	                     GUTrack* itrack_aos,
			     int *targetElements,
			     GUTrack* otrack_aos)
{
  static vecphys::cxx::ComptonKleinNishina model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {
    model.Interact<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  //validation for the total cross section
  double sigma = 0;
  for(int i = 0 ; i < ntracks ; ++i) {
    model.AtomicCrossSection<kScalar>(itrack_aos[i], targetElements[i],sigma);
  }

  return elapsedTime;
}

Precision ScalarBetheHeitler(int ntracks, 
			     GUTrack* itrack_aos,
			     int *targetElements,
			     GUTrack* otrack_aos)
{
  static vecphys::cxx::ConversionBetheHeitler model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {
    model.Interact<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision ScalarSauterGavrila(int ntracks, 
			      GUTrack* itrack_aos,
			      int *targetElements,
			      GUTrack* otrack_aos)
{
  static vecphys::cxx::PhotoElectronSauterGavrila model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {
    model.Interact<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision ScalarMollerBhabha(int ntracks, 
			     GUTrack* itrack_aos,
			     int *targetElements,
			     GUTrack* otrack_aos)
{
  static vecphys::cxx::IonisationMoller model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {
    model.Interact<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision ScalarSeltzerBerger(int ntracks, 
			      GUTrack* itrack_aos,
			      int *targetElements,
			      GUTrack* otrack_aos)
{
  static vecphys::cxx::BremSeltzerBerger model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {
    model.Interact<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

// Vector

Precision VectorKleinNishina(GUTrack_v& itrack_soa,
			     int *targetElements,
			     GUTrack_v& otrack_soa)
{
  static vecphys::cxx::ComptonKleinNishina model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  int ntracks = itrack_soa.numTracks;

  timer.Start();

  model.Interact<kVc>(itrack_soa, targetElements, otrack_soa);

  elapsedTime = timer.Stop();

  //validation for the total cross section
  double* sigma  = new double [ntracks];
  for(int i = 0 ; i < ntracks ; ++i) {
    sigma[i] = 0;  
  }

  model.AtomicCrossSection<kVc>(itrack_soa, targetElements,sigma);

  delete [] sigma;

  return elapsedTime;

}

Precision VectorBetheHeitler(GUTrack_v& itrack_soa,
			     int *targetElements,
			     GUTrack_v& otrack_soa)
{
  static vecphys::cxx::ConversionBetheHeitler model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  model.Interact<kVc>(itrack_soa, targetElements, otrack_soa);

  elapsedTime = timer.Stop();

  return elapsedTime;

}

Precision VectorSauterGavrila(GUTrack_v& itrack_soa,
			      int *targetElements,
			      GUTrack_v& otrack_soa)
{
  static vecphys::cxx::PhotoElectronSauterGavrila model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  model.Interact<kVc>(itrack_soa, targetElements, otrack_soa);

  elapsedTime = timer.Stop();

  return elapsedTime;

}

Precision VectorMollerBhabha(GUTrack_v& itrack_soa,
			     int *targetElements,
			     GUTrack_v& otrack_soa)
{
  static vecphys::cxx::IonisationMoller model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  model.Interact<kVc>(itrack_soa, targetElements, otrack_soa);

  elapsedTime = timer.Stop();

  return elapsedTime;

}

Precision VectorSeltzerBerger(GUTrack_v& itrack_soa,
			      int *targetElements,
			      GUTrack_v& otrack_soa)
{
  static vecphys::cxx::BremSeltzerBerger model(0,-1);

  Precision elapsedTime = 0.0;
  static Stopwatch timer;

  timer.Start();

  model.Interact<kVc>(itrack_soa, targetElements, otrack_soa);

  elapsedTime = timer.Stop();

  return elapsedTime;
}

// Geant4 composition and rejection

Precision G4KleinNishina(int ntracks, 
	                 GUTrack* itrack_aos,
			 int *targetElements,
			 GUTrack* otrack_aos)
{
  static vecphys::cxx::ComptonKleinNishina model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {
    model.InteractG4<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  //validation for the total cross section
  double sigma = 0.0;
  for(int i = 0 ; i < ntracks ; ++i) {
    model.AtomicCrossSectionG4<kScalar>(itrack_aos[i], targetElements[i], sigma);
  }
  return elapsedTime;
}

Precision G4BetheHeitler(int ntracks, 
			 GUTrack* itrack_aos,
			 int *targetElements,
			 GUTrack* otrack_aos)
{
  static vecphys::cxx::ConversionBetheHeitler model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {
    model.InteractG4<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision G4SauterGavrila(int ntracks, 
			  GUTrack* itrack_aos,
			  int *targetElements,
			  GUTrack* otrack_aos)
{
  static vecphys::cxx::PhotoElectronSauterGavrila model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {
    model.InteractG4<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision G4MollerBhabha(int ntracks, 
	                 GUTrack* itrack_aos,
			 int *targetElements,
			 GUTrack* otrack_aos)
{
  static vecphys::cxx::IonisationMoller model(0,-1);
  //  static vecphys::cxx::GUMollerBhabha model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {
    model.InteractG4<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision G4SeltzerBerger(int ntracks, 
			  GUTrack* itrack_aos,
			  int *targetElements,
			  GUTrack* otrack_aos)
{
  static vecphys::cxx::BremSeltzerBerger model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {
    model.InteractG4<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

} // end namespace vecphys
