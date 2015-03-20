#include "GUBenchmarker.h"

#include "base/Stopwatch.h"

#include "GUAliasSampler.h"
#include "GUComptonKleinNishina.h"
#include "GUTrackHandler.h"

#define  VECPHYS_ROOT 1

#ifdef VECPHYS_ROOT
#include "GUHistogram.h"
#endif

static int fElementMode = 1;
//   Values:
//   1  - Single material
//   

static double minP= 1.0, maxP=1.0;
static bool   verbose = false;

namespace vecphys {

GUBenchmarker::GUBenchmarker()
    : fNtracks(4992),
      fRepetitions(1),
      fVerbosity(0)
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
  // if (fVerbosity > 0) {
  printf("RunBenchmarkInteract: Ntracks = %d, fRepetitions = %d\n",
	   fNtracks,fRepetitions);
  // }

  int mismatches = 0;

  // Run all benchmarks
  CheckTimer();
  CheckRandom();

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
  // PrepareTargetElements( targetElements, fNtracks);

  GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);
  GUTrack* otrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));

  if( verbose ) std::cout << "Scalar Run" << std::endl;
  
  Stopwatch timer;
  Precision elapsedScalarTotal= 0.0;
  Precision elapsedScalarTotal2= 0.0;

  for (unsigned r = 0; r < fRepetitions; ++r)
  {
    PrepareTargetElements(targetElements, fNtracks);
    // In 'random' mode, it should change for every iteration
     
    //prepare input tracks
    fTrackHandler->GenerateRandomTracks(fNtracks, minP, maxP);    
    GUTrack* itrack_aos = fTrackHandler->GetAoSTracks();

    timer.Start();

    //AOS
    for(int i = 0 ; i < fNtracks ; ++i) {
       if( verbose ) std::cout << " E_in = "  << itrack_aos[i].E
                               << " elem = "  << targetElements[i];
       model->Interact<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
       if( verbose )  std::cout << " E_out = "  << itrack_aos[i].E
                                << " E_elec = " << otrack_aos[i].E  << std::endl;
    }

    Precision elapsedScalar = timer.Stop();
    elapsedScalarTotal  += elapsedScalar;
    elapsedScalarTotal2 += elapsedScalar*elapsedScalar;
    
#ifdef VERBOSE
    if (fVerbosity > 0) {
      printf("Scalar Task %d > %6.3f sec\n",r,elapsedScalar);
    }
#endif
    
#ifdef VECPHYS_ROOT
    histogram->ftime->Fill(elapsedScalar);
    for(int i = 0 ; i < fNtracks ; ++i) {
      histogram->fenergy->Fill(otrack_aos[i].E);
      histogram->fangle->Fill(otrack_aos[i].pz/otrack_aos[i].E);
    }
#endif    

  }
  free(otrack_aos);

  printf("Scalar Task Total time of %3d reps = %6.3f sec\n", fRepetitions, elapsedScalarTotal);
  delete model;

#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}



void GUBenchmarker::PrepareTargetElements(int *targetElements, int ntracks)
{
   int mainElement = 8;   // Oxygen
   constexpr int NumFx = 16;
   constexpr int maxElement = 100;  // Should obtain it from elsewhere ...
   
   int Element[NumFx] = { 82, 74, 8, 7, 6, 13, 18, 22, 26, 27, 30, 48, 54, 64, 79, 92 };  
                      //  Pb   W  O  N  C  Al  Ar, Ti  Fe  Cu  Zn  Cd  Xe  Gd  Au   U

   static int noCalls=0 ;
   // static int lastIndex= 1;

   noCalls++;
   
   bool report = (noCalls == 1 );
   
   if( (fElementMode == 0) || (fElementMode > maxElement ) )     //  All Elements
   {
      if( report ) 
         std::cout << " Generating Target Elements with Random mode - mode # = " << fElementMode
                   << " numFx= " << NumFx << std::endl;
      
      for(int i = 0 ; i < ntracks ; ++i) {
         targetElements[i] = ( i % maxElement) + 1 ;
      }
   }
   else if( fElementMode == 1 )
   {
      if( report ) 
         std::cout << " Using *Constant* Target Element - mode # = " << fElementMode
                   << " Element used is Z= " << mainElement << std::endl;
      
      for(int i = 0 ; i < ntracks ; ++i) {
         targetElements[i] = mainElement;
      }
   }
   else if( fElementMode <= NumFx )
   {
      int numElements = fElementMode;
      if( report )
         std::cout << " Generating Target Elements from table of elements "
                   << " - mode # = " << fElementMode << std::endl
                   << " Using " << numElements << " number of elements. " << std::endl;
      int indEl;
      for(int i = 0 ; i < ntracks ; ++i) {
         indEl = ( i % numElements ) ;
         targetElements[i] = Element[ indEl ]; 
      }
      // lastIndex= indEl; 
      
   }else{
      // Cycle through different numbers of elements
      int numElements = fElementMode;
      if( report )      
         std::cout << " Generating Target Elements cycling through "
                   << " - mode # = " << fElementMode << std::endl
                   << " Using " << numElements << " number of elements. " << std::endl;
      for(int i = 0 ; i < ntracks ; ++i) {
         targetElements[i] = 
            ( (i * 101 + (noCalls%numElements) * 59 ) % numElements ) + 1;

         assert ( (targetElements[i] > 0)  &&  (targetElements[i] <= numElements) );
      }
      // lastIndex= indEl; 
   }
}
   
void GUBenchmarker::RunGeant4() 
{
#ifdef VECPHYS_ROOT
  GUHistogram *histogram = new GUHistogram("geant4.root");
#endif

  int *targetElements = new int [fNtracks];
  // PrepareTargetElements(targetElements, fNtracks);

  vecphys::cxx::GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);
  GUTrack* otrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));

  Stopwatch timer;
  Precision  elapsedG4Total = 0.0;

  if( verbose ) 
     std::cout << " Geant4 Run Output " << std::endl;
  
  for (unsigned r = 0; r < fRepetitions; ++r) {
    PrepareTargetElements(targetElements, fNtracks);
    // In 'random' mode, it should change for every iteration

    fTrackHandler->GenerateRandomTracks(fNtracks, minP, maxP);
    GUTrack* itrack_aos = fTrackHandler->GetAoSTracks();

    timer.Start();
    //AOS
    for(int i = 0 ; i < fNtracks ; ++i) {
       if( verbose ) std::cout << " E_in = "  << itrack_aos[i].E
                               << " elem = "  << targetElements[i];
      model->InteractG4<kScalar>(itrack_aos[i],targetElements[i],otrack_aos[i]);
      if( verbose )  std::cout << " E_out = "  << itrack_aos[i].E
                               <<  " E_elec = " << otrack_aos[i].E  << std::endl;
    }

    Precision elapsedG4 = timer.Stop();
    elapsedG4Total += elapsedG4;
#ifdef VERBOSE
    if (fVerbosity > 0) {
      printf("Geant4 Task %d > %6.3f sec\n",r,elapsedG4);
    }
#endif
    
#ifdef VECPHYS_ROOT
    histogram->ftime->Fill(elapsedG4);
    for(int i = 0 ; i < fNtracks ; ++i) {
      histogram->fenergy->Fill(otrack_aos[i].E);
      histogram->fangle->Fill(otrack_aos[i].pz/otrack_aos[i].E);
    }
#endif    

  }
  printf("Geant4 Task: Total time of %3d reps = %6.3f sec , per iter/per track: %8.4g sec\n",
    fRepetitions ,elapsedG4Total, elapsedG4Total/fRepetitions );

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
  // PrepareTargetElements( targetElements, fNtracks);
  
  vecphys::cxx::GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);

  Stopwatch timer;
  Precision elapsedVectorTotal= 0.0;

  for (unsigned r = 0; r < fRepetitions; ++r) {
    PrepareTargetElements(targetElements, fNtracks);
    // In 'random' mode, it should change for every iteration

    fTrackHandler->GenerateRandomTracks(fNtracks, minP, maxP);
    
    GUTrack_v itrack_soa = fTrackHandler->GetSoATracks();

    timer.Start();

    model->Interact<kVc>(itrack_soa,targetElements,otrack_soa);

    Precision elapsedVector = timer.Stop();
    elapsedVectorTotal += elapsedVector;
#ifdef VERBOSE
    if (fVerbosity > 0) {
      printf("Vector Task %d > %6.3f sec\n",r,elapsedVector);
    }
#endif
#ifdef VECPHYS_ROOT
    histogram->ftime->Fill(elapsedVector);
    for(int i = 0 ; i < fNtracks ; ++i) {
      histogram->fenergy->Fill(otrack_soa.E[i]);
      histogram->fangle->Fill(otrack_soa.pz[i]/otrack_soa.E[i]);
    }
#endif

  }

  printf("Vector Task Total time of %3d reps = %6.3f sec\n",fRepetitions ,elapsedVectorTotal);

  delete model;
#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}

void GUBenchmarker::CheckTimer()
{
  Stopwatch timer;
  Precision elapsedNullTotal= 0.0;

  for (unsigned r = 0; r < fRepetitions; ++r) {
    timer.Start();
    Precision elapsedNull = timer.Stop();
    elapsedNullTotal += elapsedNull;
  }
  printf("Null Task Timer: Total time of %3d reps = %6.3f sec\n",fRepetitions ,elapsedNullTotal);
} 

#include <Vc/Vc>
  using Vc::double_v;

void GUBenchmarker::CheckRandom() 
{
  size_t nsample = 1000000;

  //test Vc random
  std::cout << "Vc::Size  = " <<  Vc::Vector<Precision>::Size << std::endl;

  Stopwatch timer; 
  timer.Start();

  double_v a, b, c;

  for (size_t i = 0; i < nsample ; ++i) {
    a = double_v::Random();
    b = double_v::Random();
    c += a + b;
  }

  timer.Stop();

  double vctime =  timer.Elapsed();

  std::cout << "vector sum  = " <<  c[0] + c[1] << std::endl;

  timer.Start();

  //test scalar random
  double sa = 0;
  double sb = 0;
  double sc = 0;

  for (size_t i = 0; i < nsample*2; ++i) {
    sa = (double)rand()/RAND_MAX;
    sb = (double)rand()/RAND_MAX;
    sc += sa + sb;
  }

  timer.Stop();

  double srtime =  timer.Elapsed();

  std::cout << "scalar sum  = " <<  sc << std::endl;
  std::cout << "time for sampling " << nsample << " random numbers: "
    << "(Vc,rand) = (" << vctime << ":" << srtime << ")" <<std::endl;
}

} // end namespace vecphys
