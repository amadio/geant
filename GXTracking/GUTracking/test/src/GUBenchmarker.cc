#include "GUBenchmarker.h"

#include "base/Stopwatch.h"

#include "GUAliasSampler.h"
#include "GUComptonKleinNishina.h"
#include "GUTrackHandler.h"

namespace vecphys {

GUBenchmarker::GUBenchmarker()
    : fNtracks(4992),
      fRepetitions(1),
      fVerbosity(1)
{
  fTrackHandler = new GUTrackHandler();
}

GUBenchmarker::~GUBenchmarker() 
{
  delete fTrackHandler;
}

int GUBenchmarker::RunBenchmark() {
  int errorcode=0;
  errorcode+=RunBenchmarkInteract();
  return (errorcode)? 1 : 0;
}

int GUBenchmarker::RunBenchmarkInteract() {
  if (fVerbosity > 0) {
    printf("RunBenchmarkInteract: Ntracks = %d, fRepetitions = %d\n",
	   fNtracks,fRepetitions);
  }

  int mismatches = 0;

  fTrackHandler->GenerateRandomTracks(fNtracks);
  // Run all benchmarks
  RunGeant4();
  RunScalar();
  RunVector();

#ifdef VECPHYS_CUDA
  RunCuda();
#endif
  return mismatches;
}

void GUBenchmarker::RunScalar() 
{
  //prepare input tracks
  GUTrack* itrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));
  GUTrack* otrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));
  itrack_aos = fTrackHandler->GetAoSTracks();

  int *targetElements = new int [fNtracks];
  for(int i = 0 ; i < fNtracks ; ++i) {
    targetElements[i] = i ;
  }

  GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);

  Stopwatch timer;

  for (unsigned r = 0; r < fRepetitions; ++r) {

    timer.Start();

    //AOS
    for(int i = 0 ; i < fNtracks ; ++i) {
      model->Interact<kScalar>(itrack_aos[i],targetElements[i],otrack_aos[i]);
    }

    Precision elapsedScalar = timer.Stop();

    if (fVerbosity > 0) {
      printf("Scalar Task %d >: %6.3fs\n",r,elapsedScalar);
    }

    if (fVerbosity > 1) {
      for(unsigned i = 0; i < 8 ; ++i) {
        printf(" E[%d]= %f\n",i,otrack_aos[i].E);
      }
    }
  }

  free(itrack_aos);
  free(otrack_aos);
}

void GUBenchmarker::RunGeant4() 
{
  GUTrack* itrack_aos = fTrackHandler->GetAoSTracks();
  GUTrack* otrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));
  itrack_aos = fTrackHandler->GetAoSTracks();

  int *targetElements = new int [fNtracks];
  for(int i = 0 ; i < fNtracks ; ++i) {
    targetElements[i] = i ;
  }

  vecphys::cxx::GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);

  Stopwatch timer;
  for (unsigned r = 0; r < fRepetitions; ++r) {
    timer.Start();
    //AOS
    for(int i = 0 ; i < fNtracks ; ++i) {
      model->InteractG4<kScalar>(itrack_aos[i],targetElements[i],otrack_aos[i]);
    }

    Precision elapsedScalar = timer.Stop();
    if (fVerbosity > 0) {
      printf("Geant4 Task %d >: %6.3fs\n",r,elapsedScalar);
    }
  }

  if (fVerbosity > 1) {
    for(unsigned i = 0; i < 4 ; ++i) printf(" E[%d]= %f\n",i,otrack_aos[i].E);
  }
}

void GUBenchmarker::RunVector()
{
  //input SOA tracks
  GUTrack_v track_in = fTrackHandler->GetSoATracks();
  GUTrackHandler *handler_out = new GUTrackHandler(fNtracks);
  GUTrack_v track_out = handler_out->GetSoATracks();

  int *targetElements = new int[fNtracks];
  for(int i = 0 ; i < fNtracks ; ++i) {
    targetElements[i] = i ;
  }

  vecphys::cxx::GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);

  Stopwatch timer;
  for (unsigned r = 0; r < fRepetitions; ++r) {
    timer.Start();

    model->Interact(track_in,targetElements,track_out);

    Precision elapsedVector = timer.Stop();
    if (fVerbosity > 0) {
      printf("Vector Task %d >: %6.3fs\n",r,elapsedVector);
    }
  }

  if (fVerbosity > 1) {
    for(unsigned i = 0; i < 4 ; ++i) printf(" E[%d]= %f\n",i,(track_out.E)[i]);
  }
}

} // end namespace vecphys
