#include "GUBenchmarker.h"

#include "base/Stopwatch.h"

#include "GUAliasSampler.h"
#include "GUComptonKleinNishina.h"
#include "GUTrackHandler.h"

#ifdef VECPHYS_ROOT
#include "GUHistogram.h"
#endif

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
#ifdef VECPHYS_ROOT
  GUHistogram *histogram = new GUHistogram("scalar.root");
#endif

  int *targetElements = new int [fNtracks];
  for(int i = 0 ; i < fNtracks ; ++i) {
    targetElements[i] = i ;
  }

  GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);
  GUTrack* otrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));

  Stopwatch timer;

  for (unsigned r = 0; r < fRepetitions; ++r) {
    //prepare input tracks
    fTrackHandler->GenerateRandomTracks(fNtracks);
    GUTrack* itrack_aos = fTrackHandler->GetAoSTracks();

    timer.Start();

    //AOS
    for(int i = 0 ; i < fNtracks ; ++i) {
      model->Interact<kScalar>(itrack_aos[i],targetElements[i],otrack_aos[i]);
    }

    Precision elapsedScalar = timer.Stop();
    if (fVerbosity > 0) {
      printf("Scalar Task %d > %6.3f sec\n",r,elapsedScalar);
    }

#ifdef VECPHYS_ROOT
    histogram->ftime->Fill(elapsedScalar);
    for(int i = 0 ; i < fNtracks ; ++i) {
      histogram->fenergy->Fill(otrack_aos[i].E);
      histogram->fangle->Fill(otrack_aos[i].pz/otrack_aos[i].E);
    }
#endif    

  }
  free(otrack_aos);

  delete model;

#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}

void GUBenchmarker::RunGeant4() 
{
#ifdef VECPHYS_ROOT
  GUHistogram *histogram = new GUHistogram("geant4.root");
#endif

  int *targetElements = new int [fNtracks];
  for(int i = 0 ; i < fNtracks ; ++i) {
    targetElements[i] = i ;
  }

  vecphys::cxx::GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);
  GUTrack* otrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));

  Stopwatch timer;
  for (unsigned r = 0; r < fRepetitions; ++r) {

    fTrackHandler->GenerateRandomTracks(fNtracks);
    GUTrack* itrack_aos = fTrackHandler->GetAoSTracks();

    timer.Start();
    //AOS
    for(int i = 0 ; i < fNtracks ; ++i) {
      model->InteractG4<kScalar>(itrack_aos[i],targetElements[i],otrack_aos[i]);
    }

    Precision elapsedG4 = timer.Stop();
    if (fVerbosity > 0) {
      printf("Geant4 Task %d > %6.3f sec\n",r,elapsedG4);
    }

#ifdef VECPHYS_ROOT
    histogram->ftime->Fill(elapsedG4);
    for(int i = 0 ; i < fNtracks ; ++i) {
      histogram->fenergy->Fill(otrack_aos[i].E);
      histogram->fangle->Fill(otrack_aos[i].pz/otrack_aos[i].E);
    }
#endif    

  }

  free(otrack_aos);
  delete model;
#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}

void GUBenchmarker::RunVector()
{
#ifdef VECPHYS_ROOT
  GUHistogram *histogram = new GUHistogram("vector.root");
#endif
  //output SOA tracks
  GUTrackHandler *handler_out = new GUTrackHandler(fNtracks);
  GUTrack_v otrack_soa = handler_out->GetSoATracks();

  int *targetElements = new int[fNtracks];
  for(int i = 0 ; i < fNtracks ; ++i) {
    targetElements[i] = i ;
  }

  vecphys::cxx::GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);

  Stopwatch timer;
  for (unsigned r = 0; r < fRepetitions; ++r) {

    fTrackHandler->GenerateRandomTracks(fNtracks);
    GUTrack_v itrack_soa = fTrackHandler->GetSoATracks();

    timer.Start();

    model->Interact<kVc>(itrack_soa,targetElements,otrack_soa);

    Precision elapsedVector = timer.Stop();
    if (fVerbosity > 0) {
      printf("Vector Task %d > %6.3f sec\n",r,elapsedVector);
    }

#ifdef VECPHYS_ROOT
    histogram->ftime->Fill(elapsedVector);
    for(int i = 0 ; i < fNtracks ; ++i) {
      histogram->fenergy->Fill(otrack_soa.E[i]);
      histogram->fangle->Fill(otrack_soa.pz[i]/otrack_soa.E[i]);
    }
#endif

  }

  delete model;
#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}

} // end namespace vecphys
