#include "GUBenchmarker.h"

#include "base/Stopwatch.h"

#include "GUAliasSampler.h"
#include "GUComptonKleinNishina.h"
#include "GUTrackHandler.h"
#include "SystemOfUnits.h"

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
      fVerbosity(0),
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
  GUHistogram *histogram = new GUHistogram("scalar.root", fMaxP);  // maxE= maxP (photon)    
#endif

  int *targetElements = new int [fNtracks];

  GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);
  GUTrack* otrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));

  if( verbose ) std::cout << "Scalar Run" << std::endl;
  
  Stopwatch timer;
  Precision elapsedScalarTotal= 0.0;
  Precision elapsedScalarTotal2= 0.0;

  Precision* incomingEn = new Precision[fNtracks];

  for (unsigned r = 0; r < fRepetitions; ++r)
  {
    PrepareTargetElements(targetElements, fNtracks);
    // In 'random' mode, it should change for every iteration
     
    //prepare input tracks
    fTrackHandler->GenerateRandomTracks(fNtracks, fMinP, fMaxP);
    GUTrack* itrack_aos = fTrackHandler->GetAoSTracks();

    for(int i = 0 ; i < fNtracks ; ++i) {
       incomingEn[i] = itrack_aos[i].E;
    }

    timer.Start();

    //AOS
    for(int i = 0 ; i < fNtracks ; ++i) {
       model->Interact<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
    }

    Precision elapsedScalar = timer.Stop();
    elapsedScalarTotal  += elapsedScalar;
    elapsedScalarTotal2 += elapsedScalar*elapsedScalar;
    
    if (fVerbosity > 0) {
      printf("Scalar Task %d > %6.3f sec\n",r,elapsedScalar);
    }
    
#ifdef VECPHYS_ROOT
    histogram->RecordTime(elapsedScalar);
    for(int i = 0 ; i < fNtracks ; ++i) {
       // histogram->fenergy->Fill(otrack_aos[i].E);
       // histogram->fangle->Fill(otrack_aos[i].pz/otrack_aos[i].E);
       histogram->RecordHistos( incomingEn[i],
                                itrack_aos[i].E,                  // Outgoing gamma  energy
                                itrack_aos[i].pz/itrack_aos[i].E, // Outgoing gamma  u_z 
                                otrack_aos[i].E,                  // Electron  energy
                                otrack_aos[i].pz/otrack_aos[i].E); // Electron  u_z
    }    
#endif    

  }
  free(otrack_aos);

  // double averageScalar = elapsedScalarTotal  / fRepetitions;
  // double averSqrScalar = elapsedScalarTotal2 / fRepetitions;  
  printf("Scalar Task Total time of %3d reps = %6.3f sec\n", fRepetitions, elapsedScalarTotal);
  delete model;

  delete[] incomingEn; 
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
  GUHistogram *histogram = new GUHistogram("geant4.root", fMaxP);  // maxE= maxP (photon)
#endif

  int *targetElements = new int [fNtracks];
  // PrepareTargetElements(targetElements, fNtracks);

  vecphys::cxx::GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);
  GUTrack* otrack_aos = (GUTrack*) malloc(fNtracks*sizeof(GUTrack));

  Stopwatch timer;
  Precision  elapsedG4Total = 0.0;

  if( verbose ) 
     std::cout << " Geant4 Run Output " << std::endl;

  Precision* incomingEn = new Precision[fNtracks];


#ifdef UNIT_CORRECTION
  const double  invFactor = 1.0 / unitFactor;

  if( unitFactor != 1.0 ) 
     std::cout << " Using TRIAL Unit factor =  " << unitFactor << std::endl;
#endif
  
  for (unsigned r = 0; r < fRepetitions; ++r) {
    PrepareTargetElements(targetElements, fNtracks);
    // In 'random' mode, it should change for every iteration

    fTrackHandler->GenerateRandomTracks(fNtracks, fMinP, fMaxP);
    GUTrack* itrack_aos = fTrackHandler->GetAoSTracks();

    for(int i = 0 ; i < fNtracks ; ++i) {
       incomingEn[i] = itrack_aos[i].E;   // Record the unaltered value
#ifdef UNIT_CORRECTION
       itrack_aos[i].E  *= unitFactor;
       itrack_aos[i].px *= unitFactor;
       itrack_aos[i].py *= unitFactor;
       itrack_aos[i].pz *= unitFactor;
#endif
       // incomingEn[i] = itrack_aos[i].E;  // Record the altered value
    }
    
    timer.Start();
    //AOS
    for(int i = 0 ; i < fNtracks ; ++i) {

      model->InteractG4<kScalar>(itrack_aos[i], targetElements[i], otrack_aos[i]);
      //     **********
    }

    Precision elapsedG4 = timer.Stop();
    elapsedG4Total += elapsedG4;
    if (fVerbosity > 0) {
      printf("Geant4 Task %d > %6.3f sec\n",r,elapsedG4);
    }
    
#ifdef VECPHYS_ROOT
    histogram->RecordTime(elapsedG4);
    for(int i = 0 ; i < fNtracks ; ++i)
    {
#ifdef UNIT_CORRECTION
       // incomingEn[i]    *= invFactor; // (if) Recorded the altered value
       itrack_aos[i].E  *= invFactor;
       itrack_aos[i].px *= invFactor;
       itrack_aos[i].py *= invFactor;
       itrack_aos[i].pz *= invFactor;

       otrack_aos[i].E  *= invFactor;
       otrack_aos[i].px *= invFactor;
       otrack_aos[i].py *= invFactor;
       otrack_aos[i].pz *= invFactor;
#endif
       histogram->RecordHistos( incomingEn[i],
                                itrack_aos[i].E,                  // Outgoing gamma  energy
                                itrack_aos[i].pz/itrack_aos[i].E, // Outgoing gamma  u_z    
                                otrack_aos[i].E,                  // Electron  energy
                                otrack_aos[i].pz/otrack_aos[i].E); // Electron  u_z

       if( verbose )
       {
          ReportInteraction( incomingEn[i],
                             targetElements[i],                             
                             itrack_aos[i].E,                  // Outgoing gamma  energy
                             itrack_aos[i].pz/itrack_aos[i].E, // Outgoing gamma  u_z    
                             otrack_aos[i].E,                  // Electron  energy
                             otrack_aos[i].pz/otrack_aos[i].E); // Electron  u_z
          
          printf( " E_in = %8.5f  elem = %4d E_out = %8.4g  E_elec = %8.4g  Sum(E_out)= %8.4g  Diff(In-Out)= %8.3g \n",
                  incomingEn[i], 
                  targetElements[i], 
                  itrack_aos[i].E,
                  otrack_aos[i].E,
                  itrack_aos[i].E + otrack_aos[i].E, 
                  incomingEn[i] - ( itrack_aos[i].E + otrack_aos[i].E )
                ); 
       }
    }    
#endif    

  }
  printf("Geant4 Task: Total time of %3d reps = %6.3f sec , per iter/per track: %8.4g sec\n",
    fRepetitions ,elapsedG4Total, elapsedG4Total/fRepetitions );

  free(otrack_aos);

  delete[] incomingEn;
  delete model;
#ifdef VECPHYS_ROOT
  delete histogram;
#endif
}

void GUBenchmarker::RunVector()
{
#ifdef VECPHYS_ROOT
  GUHistogram *histogram = new GUHistogram("vector.root", fMaxP);  // maxE= fMaxP (photon)  
#endif
  //output SOA tracks
  GUTrackHandler *handler_out = new GUTrackHandler(fNtracks);
  GUTrack_v otrack_soa = handler_out->GetSoATracks();

  int *targetElements = new int[fNtracks];
  // PrepareTargetElements( targetElements, fNtracks);
  Precision* incomingEn = new Precision[fNtracks];
  
  vecphys::cxx::GUComptonKleinNishina *model = new GUComptonKleinNishina(0,-1);

   Stopwatch timer;
  Precision elapsedVectorTotal= 0.0;

  for (unsigned r = 0; r < fRepetitions; ++r) {
    PrepareTargetElements(targetElements, fNtracks);
    // In 'random' mode, it should change for every iteration

    fTrackHandler->GenerateRandomTracks(fNtracks, fMinP, fMaxP);
    
    GUTrack_v itrack_soa = fTrackHandler->GetSoATracks();

    for(int i = 0 ; i < fNtracks ; ++i) {
       incomingEn[i] = itrack_soa.E[i];
    }
    
    timer.Start();

    model->Interact<kVc>(itrack_soa,targetElements,otrack_soa);

    Precision elapsedVector = timer.Stop();
    elapsedVectorTotal += elapsedVector;
    if (fVerbosity > 0) {
      printf("Vector Task %d > %6.3f sec\n",r,elapsedVector);
    }
#ifdef VECPHYS_ROOT
    // histogram->ftime->Fill(elapsedVector);
    for(int i = 0 ; i < fNtracks ; ++i) {
       // histogram->fenergy->Fill(otrack_soa.E[i]);
       // histogram->fangle->Fill(otrack_soa.pz[i]/otrack_soa.E[i]);
    }

    histogram->RecordTime(elapsedVector);
    for(int i = 0 ; i < fNtracks ; ++i) {
      // histogram->fenergy->Fill(otrack_aos[i].E);
      // histogram->fangle->Fill(otrack_aos[i].pz/otrack_aos[i].E);

       histogram->RecordHistos( incomingEn[i],
                                itrack_soa.E[i],                  // Outgoing gamma  energy
                                itrack_soa.pz[i]/itrack_soa.E[i], // Outgoing gamma  u_z    
                                otrack_soa.E[i],                  // Electron  energy
                                otrack_soa.pz[i]/otrack_soa.E[i]); // Electron  u_z

       if( verbose )
          GUBenchmarker::ReportInteraction( incomingEn[i],
                                            targetElements[i], 
                                itrack_soa.E[i],                  // Outgoing gamma  energy
                                itrack_soa.pz[i]/itrack_soa.E[i], // Outgoing gamma  u_z    
                                otrack_soa.E[i],                  // Electron  energy
                                otrack_soa.pz[i]/otrack_soa.E[i]); // Electron  u_z          
    } 
#endif

  }

  printf("Vector Task Total time of %3d reps = %6.3f sec\n",fRepetitions ,elapsedVectorTotal);

  delete[] incomingEn;
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

  double sumC= 0.0;
  for (size_t i = 0; i <  Vc::Vector<Precision>::Size ; ++i) {
     sumC += c[i]; 
  }
  
  // std::cout << "vector sum  = " <<  c[0] + c[1] << std::endl;
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
