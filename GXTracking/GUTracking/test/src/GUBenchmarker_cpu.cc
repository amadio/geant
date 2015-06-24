#include "base/Stopwatch.h"
#include "GUHistogram.h"

#include "GUComptonKleinNishina.h"
#include "GVComptonKleinNishina.h"
#include "GUConversionBetheHeitler.h"
#include "GUPhotoElectronSauterGavrila.h"
#include "GUMollerBhabha.h"
#include "GUSeltzerBerger.h"

namespace vecphys {

// Scalar

Precision ScalarKleinNishina(int ntracks, 
	                     GUTrack* itrack_aos,
			     int *targetElements,
			     GUTrack* otrack_aos)
{
  static vecphys::cxx::GUComptonKleinNishina model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {

    model.Interact<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision ScalarVKleinNishina(int ntracks, 
	                      GUTrack* itrack_aos,
			      int *targetElements,
			      GUTrack* otrack_aos)
{
  static vecphys::cxx::GVComptonKleinNishina model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {

    model.Interact<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision ScalarBetheHeitler(int ntracks, 
			     GUTrack* itrack_aos,
			     int *targetElements,
			     GUTrack* otrack_aos)
{
  static vecphys::cxx::GUConversionBetheHeitler model(0,-1);

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
  static vecphys::cxx::GUPhotoElectronSauterGavrila model(0,-1);

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
  static vecphys::cxx::GUMollerBhabha model(0,-1);

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
  static vecphys::cxx::GUSeltzerBerger model(0,-1);

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
  static vecphys::cxx::GUComptonKleinNishina model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  model.Interact<kVc>(itrack_soa, targetElements, otrack_soa);

  elapsedTime = timer.Stop();

  return elapsedTime;

}

Precision VectorVKleinNishina(GUTrack_v& itrack_soa,
			      int *targetElements,
			      GUTrack_v& otrack_soa)
{
  static vecphys::cxx::GVComptonKleinNishina model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  model.Interact<kVc>(itrack_soa, targetElements, otrack_soa);

  elapsedTime = timer.Stop();

  return elapsedTime;

}

Precision VectorBetheHeitler(GUTrack_v& itrack_soa,
			     int *targetElements,
			     GUTrack_v& otrack_soa)
{
  static vecphys::cxx::GUConversionBetheHeitler model(0,-1);

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
  static vecphys::cxx::GUPhotoElectronSauterGavrila model(0,-1);

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
  static vecphys::cxx::GUMollerBhabha model(0,-1);

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
  static vecphys::cxx::GUSeltzerBerger model(0,-1);

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
  static vecphys::cxx::GUComptonKleinNishina model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {
    model.InteractG4<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision G4VKleinNishina(int ntracks, 
	                  GUTrack* itrack_aos,
                          int *targetElements,
			  GUTrack* otrack_aos)
{
  static vecphys::cxx::GVComptonKleinNishina model(0,-1);

  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  for(int i = 0 ; i < ntracks ; ++i) {
    model.InteractG4<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
  }

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision G4BetheHeitler(int ntracks, 
			 GUTrack* itrack_aos,
			 int *targetElements,
			 GUTrack* otrack_aos)
{
  static vecphys::cxx::GUConversionBetheHeitler model(0,-1);

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
  static vecphys::cxx::GUPhotoElectronSauterGavrila model(0,-1);

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
  static vecphys::cxx::GUMollerBhabha model(0,-1);

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
  static vecphys::cxx::GUSeltzerBerger model(0,-1);

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
