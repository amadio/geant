#include "GUBenchmarker.h"

#include "base/Stopwatch.h"

#include "GUAliasSampler.h"
#include "GUComptonKleinNishina.h"
#include "GUTrackHandler.h"

namespace vecphys {

GUBenchmarker::GUBenchmarker()
    : fNtracks(4992),
      fRepetitions(10),
      fVerbosity(1)
{
  ;
}

GUBenchmarker::~GUBenchmarker() 
{
  ;
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

  // Run all benchmarks
  RunScalar();
  RunVector();

#ifdef VECPHYS_CUDA
  RunCuda();
#endif
  return mismatches;
}

void GUBenchmarker::RunScalar() 
{
  GUTrackHandler *handler_in = new GUTrackHandler();
  handler_in->GenerateRandomTracks(fNtracks);

  GUTrack* track_aos = handler_in->GetAoSTracks();
  GUTrack* track_aos_out = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));

  int *targetElements = new int [fNtracks];
  for(int i = 0 ; i < fNtracks ; ++i) {
    targetElements[i] = i ;
  }

  Stopwatch timer;
  timer.Start();

  vecphys::cxx::GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);

  for (unsigned r = 0; r < fRepetitions; ++r) {
    //AOS
    for(int i = 0 ; i < fNtracks ; ++i) {
      model->Interact(track_aos[i],targetElements[i],&track_aos_out[i]);
    }
  }
  Precision elapsedScalar = timer.Stop();
  if (fVerbosity > 0) {
    printf("Scalar: %6.3fs\n",elapsedScalar);
  }
}

void GUBenchmarker::RunVector()
{
  GUTrackHandler *handler_in = new GUTrackHandler();
  handler_in->GenerateRandomTracks(fNtracks);

  GUTrack_v track_in = handler_in->GetSoATracks();

  GUTrackHandler *handler_out = new GUTrackHandler(fNtracks);
  GUTrack_v track_out = handler_out->GetSoATracks();

  //put GUComptonProcess here which will call
  vecphys::cxx::GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);

  int *targetElements = new int[fNtracks];
  for(int i = 0 ; i < fNtracks ; ++i) {
    targetElements[i] = i ;
  }

  Stopwatch timer;
  timer.Start();
  for (unsigned r = 0; r < fRepetitions; ++r) {
    //SOA
    model->Interact(track_in,targetElements,&track_out);
  }
  Precision elapsedVector = timer.Stop();
  timer.Start();
  if (fVerbosity > 0) {
    printf("Vector: %6.3fs\n",elapsedVector);
  }
}

} // end namespace vecphys
