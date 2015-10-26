#include "GUBenchmarker.h"
#include "GUBenchmarker_cpu.h"
#include "GUPhysicsModelName.h"

#include "base/Stopwatch.h"
#include "base/SystemOfUnits.h"

#include "GUTrackHandler.h"

#ifdef VECPHYS_ROOT
#include "GUHistogram.h"
#endif

static int fElementMode = 1;
//   Values:
//   1  - Single material
//   

static bool   verbose = false; // true;

// #define UNIT_CORRECTION 1000

#ifdef UNIT_CORRECTION
//  Use this option to check whether the 'G4' option has a unit error
const double unitFactor = UNIT_CORRECTION;
//  Trial: Can use value of 1000 case conversion from GeV to MeV is needed.
#endif

namespace vecphys {

// Mono-energetic test
constexpr double defaultMinP = 500.*MeV;
constexpr double defaultMaxP = defaultMinP;
//  Equal maxp & minP ==> monoenergetic beam

// 'Broad spectrum' test 
// static double defaultMinP=   5.0 * KeV;
// static double defaultMaxP= 500.0 * MeV;

GUBenchmarker::GUBenchmarker()
    : fNtracks(4992),
      fRepetitions(1),
      fVerbosity(1),
      fMinP(defaultMinP),  
      fMaxP(defaultMaxP)   
{
  fTrackHandler = new GUTrackHandler();
}

GUBenchmarker::~GUBenchmarker() 
{
  delete fTrackHandler;
}

int GUBenchmarker::RunBenchmark()
{
  printf("fMinP = %f, fMaxP = %f \n",fMinP, fMaxP);
  int errorcode=0;
  errorcode+=RunBenchmarkInteract();
  return (errorcode)? 1 : 0;
}

int GUBenchmarker::RunBenchmarkInteract()
{
  if (fVerbosity > 0) {
    printf("RunBenchmarkInteract: Ntracks = %d, fRepetitions = %d\n",
  	   fNtracks,fRepetitions);
  }

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
  GUHistogram *histogram = new GUHistogram("scalar.root", fMaxP); //maxE= maxP
#endif

  if( verbose ) std::cout << "Scalar Run" << std::endl;
  
  int *targetElements = new int [fNtracks];
  GUTrack* itrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));
  GUTrack* otrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));

  Precision elapsedTotal[kNumberPhysicsModel];
  Precision elapsedT[kNumberPhysicsModel];
  for(int k = 0 ; k < kNumberPhysicsModel ; ++k) elapsedTotal[k] = 0.; 

  for (unsigned r = 0; r < fRepetitions; ++r)
  {
    PrepareTargetElements(targetElements, fNtracks);
    // In 'random' mode, it should change for every iteration
     
    //prepare input tracks
    fTrackHandler->GenerateRandomTracks(fNtracks, fMinP, fMaxP);
    GUTrack* track_aos = fTrackHandler->GetAoSTracks();

    for(unsigned int k = 0 ; k < kNumberPhysicsModel ; ++k) {
      fTrackHandler->CopyAoSTracks(track_aos,itrack_aos);
      elapsedT[k] = 0.0;
      elapsedT[k] = ScalarKernelFunc[k](fNtracks,itrack_aos,targetElements,
                                        otrack_aos);
      elapsedTotal[k] += elapsedT[k];

#ifdef VECPHYS_ROOT
      histogram->RecordTime(k,elapsedT[k]);
      for(int i = 0 ; i < fNtracks ; ++i) {
        histogram->RecordHistos(k,track_aos[i].E,
	  		        itrack_aos[i].E,
			        itrack_aos[i].pz/itrack_aos[i].E,
			        otrack_aos[i].E,
			        otrack_aos[i].pz/otrack_aos[i].E);
      } 
#endif    
    }
  }

  for(int k = 0 ; k < kNumberPhysicsModel ; ++k) {
    printf("%s  Scalar Total time of %3d reps = %6.3f sec\n",
           GUPhysicsModelName[k], fRepetitions, elapsedTotal[k]);
  }

  free(itrack_aos);
  free(otrack_aos);
  delete [] targetElements;
#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}

void GUBenchmarker::RunGeant4() 
{
#ifdef VECPHYS_ROOT
  GUHistogram *histogram = new GUHistogram("geant4.root", fMaxP); //maxE= maxP
#endif

  if( verbose ) std::cout << " Geant4 Run Output " << std::endl;

  int *targetElements = new int [fNtracks];
  GUTrack* itrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));
  GUTrack* otrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));

  Precision elapsedTotal[kNumberPhysicsModel];
  Precision elapsedT[kNumberPhysicsModel];
  for(int k = 0 ; k < kNumberPhysicsModel ; ++k) elapsedTotal[k] = 0.; 

  for (unsigned r = 0; r < fRepetitions; ++r) {

    PrepareTargetElements(targetElements, fNtracks);
    // In 'random' mode, it should change for every iteration

    fTrackHandler->GenerateRandomTracks(fNtracks, fMinP, fMaxP);
    GUTrack* track_aos = fTrackHandler->GetAoSTracks();

    for(unsigned int k = 0 ; k < kNumberPhysicsModel ; ++k) {
      fTrackHandler->CopyAoSTracks(track_aos,itrack_aos);

      elapsedT[k] = 0.0;
      elapsedT[k] = G4KernelFunc[k](fNtracks,itrack_aos,targetElements,
                                    otrack_aos);
      elapsedTotal[k] += elapsedT[k];

#ifdef VECPHYS_ROOT
      histogram->RecordTime(k,elapsedT[k]);
      for(int i = 0 ; i < fNtracks ; ++i) {
        histogram->RecordHistos(k,track_aos[i].E,
	  		        itrack_aos[i].E,
			        itrack_aos[i].pz/itrack_aos[i].E,
			        otrack_aos[i].E,
			        otrack_aos[i].pz/otrack_aos[i].E);
      } 
#endif    
    }
  }    

  for(int k = 0 ; k < kNumberPhysicsModel ; ++k) {
    printf("%s  Geant4 Total time of %3d reps = %6.3f sec\n",
           GUPhysicsModelName[k], fRepetitions, elapsedTotal[k]);
  }

  free(itrack_aos);
  free(otrack_aos);
  delete [] targetElements;
#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}

void GUBenchmarker::RunVector()
{
#ifdef VECPHYS_ROOT
  GUHistogram *histogram = new GUHistogram("vector.root", fMaxP);//maxE = fMaxP
#endif

  //output SOA tracks
  GUTrackHandler *handler_out = new GUTrackHandler(fNtracks);
  GUTrack_v otrack_soa = handler_out->GetSoATracks();

  GUTrackHandler *handler_in = new GUTrackHandler(fNtracks);
  GUTrack_v itrack_soa = handler_in->GetSoATracks();

  int *targetElements = new int[fNtracks];
  
  Precision elapsedTotal[kNumberPhysicsModel];
  Precision elapsedT[kNumberPhysicsModel];
  for(int k = 0 ; k < kNumberPhysicsModel ; ++k) elapsedTotal[k] = 0.; 

  for (unsigned r = 0; r < fRepetitions; ++r) {

    PrepareTargetElements(targetElements, fNtracks);
    // In 'random' mode, it should change for every iteration

    fTrackHandler->GenerateRandomTracks(fNtracks, fMinP, fMaxP);
    GUTrack_v& track_soa = fTrackHandler->GetSoATracks();
    
    for(unsigned int k = 0 ; k < kNumberPhysicsModel ; ++k) {
      fTrackHandler->CopySoATracks(track_soa,itrack_soa);

      elapsedT[k] = 0.0;
      elapsedT[k] = VectorKernelFunc[k](itrack_soa,targetElements,otrack_soa);
      elapsedTotal[k] += elapsedT[k];

#ifdef VECPHYS_ROOT
      histogram->RecordTime(k,elapsedT[k]);
      for(int i = 0 ; i < fNtracks ; ++i) {
        histogram->RecordHistos(k,track_soa.E[i],
	  		        itrack_soa.E[i],
			        itrack_soa.pz[i]/itrack_soa.E[i],
			        otrack_soa.E[i],
			        otrack_soa.pz[i]/otrack_soa.E[i]);
      } 
#endif    
    }
  }    

  for(int k = 0 ; k < kNumberPhysicsModel ; ++k) {
    printf("%s  Vector Total time of %3d reps = %6.3f sec\n",
           GUPhysicsModelName[k], fRepetitions, elapsedTotal[k]);
  }
  delete handler_in;
  delete handler_out;
  delete [] targetElements;
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
  printf("Null Task Timer: Total time of %3d reps = %6.3f sec\n", fRepetitions,
         elapsedNullTotal);
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

  double sumC= 0.0;
  for (size_t i = 0; i <  Vc::Vector<Precision>::Size ; ++i) {
     sumC += c[i]; 
  }
  
  std::cout << "vector sum  = " <<  sumC << std::endl;

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
   
void
GUBenchmarker::ReportInteraction( double incomingEnergy,
                                  int    targetElement,
                                  double GammaOut_E,
                                  double GammaOut_Pz,
                                  double electron_En,
                                  double electron_Pz,
                                  bool   print_Uz
   )
{
   double GammaOut_Uz = GammaOut_E  > 0.0 ? GammaOut_Pz / GammaOut_E : -2.0  ;
   double electron_Uz = electron_En > 0.0 ? electron_Pz / electron_En : -2.0  ;
   double SumEout=      GammaOut_E + electron_En;
   
   printf( " E_in = %8.5f  elem = %4d E_out = %8.4g  E_elec = %8.4g  Sum(E_out)= %8.4g  Diff(In-Out)= %8.3g   ",
           incomingEnergy,
           targetElement,
           GammaOut_E,
           electron_En,
           SumEout,
           incomingEnergy - SumEout
      );
   if( print_Uz )
      printf(" |  g.uz= %8.4g e-.uz = %8.4g",  GammaOut_Uz, electron_Uz );
   printf("\n"); 
}

} // end namespace vecphys
