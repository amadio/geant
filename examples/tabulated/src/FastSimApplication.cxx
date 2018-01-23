#ifdef USE_ROOT
#include "FastSimApplication.h"

#ifdef USE_VECGEOM_NAVIGATOR
  #include "management/GeoManager.h"
  using vecgeom::GeoManager;
#else 
#include "TGeoNode.h"
#endif

#include "GeantRunManager.h"
#include "GeantFactoryStore.h"
#include "GeantTrackVec.h"
#include "GeantPropagator.h"
#include "GeantTaskData.h"
#include "globals.h"

#ifdef USE_ROOT
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#endif

#include <cassert>
#include <string>

using std::min;
using std::max;
using namespace Geant;

//______________________________________________________________________________
FastSimApplication::FastSimApplication(GeantRunManager *runmgr)
  : GeantVApplication(runmgr), fInitialized(false), fMHist(), 
    fRatioMomentumInTracker( nullptr ), fRatioEnergyInEcal( nullptr ),
    fRatioEnergyInHcal( nullptr )
{
  //std::cout << "APPLICATION : FastSimApplication::FastSimApplication" << std::endl;  // Debug
  #ifdef USE_ROOT
  fRatioMomentumInTracker = new TH1F( "TrackerRatioP", "Momentum smeared in tracker", 100, 0.8, 1.2 );
  fRatioMomentumInTracker->GetXaxis()->SetTitle( "p_smeared/p_true" );
  fRatioMomentumInTracker->GetYaxis()->SetTitle( "Entries" );
  fRatioEnergyInEcal = new TH1F( "EcalRatioE", "Energy smeared in ECAL", 100, 0.8, 1.2 );
  fRatioEnergyInEcal->GetXaxis()->SetTitle( "E_smeared/E_true" );
  fRatioEnergyInEcal->GetYaxis()->SetTitle( "Entries" );
  fRatioEnergyInHcal = new TH1F( "HcalRatioE", "Energy smeared in HCAL", 100, 0.0, 2.0 );
  fRatioEnergyInHcal->GetXaxis()->SetTitle( "E_smeared/E_true" );
  fRatioEnergyInHcal->GetYaxis()->SetTitle( "Entries" );
  #endif 
}

//______________________________________________________________________________
bool FastSimApplication::Initialize() {
  if ( fInitialized ) return true;
  // Fill the vectors isTrackerVolume, isEcalVolume, isHcalVolume
  std::cout << "APPLICATION : FastSimApplication::Initialize" << std::endl;  // Debug
  int numberOfVolumes = fRunMgr->GetNvolumes();
  vector_t<Volume_t const *> &lvolumes = fRunMgr->GetVolumes();
  printf("Found %d logical volumes", numberOfVolumes);
  const Volume_t *aVolume;
  isTrackerVolume.reserve( numberOfVolumes+1 );
  isEcalVolume.reserve( numberOfVolumes+1 );
  isHcalVolume.reserve( numberOfVolumes+1 );
  for ( auto ivol = 0; ivol < numberOfVolumes; ivol++ ) {
    aVolume = lvolumes[ivol];
    if (!aVolume) break;
#ifdef USE_VECGEOM_NAVIGATOR
    int idVol = aVolume->id();
#else
    int idVol = aVolume->GetNumber();
#endif
    std::string nameVolume;
    if ( aVolume ) nameVolume = aVolume->GetName();
    std::transform(nameVolume.begin(), nameVolume.end(), nameVolume.begin(), ::tolower);
    int flagTracker = 0;
    int flagEcal = 0;
    int flagHcal = 0;
    if ( nameVolume.find(std::string("tracker")) != std::string::npos ) {
      flagTracker = 1;
    } else if ( nameVolume.find(std::string("ecal")) != std::string::npos ) {
      flagEcal = 1;
    } else if ( nameVolume.find(std::string("hcal")) != std::string::npos ) {
      flagHcal = 1;
    }
    isTrackerVolume[ idVol ] = flagTracker;
    isEcalVolume[ idVol ] = flagEcal;
    isHcalVolume[ idVol ] = flagHcal;
    std::cout << "\t name=" << nameVolume << " ; id=" << idVol
              << " ; (my-counter=" << ivol << ")"
              << " ; flagTracker=" << flagTracker << " ; flagEcal=" << flagEcal 
              << " ; flagHcal=" << flagHcal << std::endl;
  }
  fInitialized = true;  
  return true;
}

//______________________________________________________________________________
void FastSimApplication::StepManager( int npart, const GeantTrack_v &tracks, GeantTaskData *td ) {
  // Application stepping manager. 
  //std::cout << "APPLICATION : FastSimApplication::StepManager" << std::endl;  // Debug
  GeantPropagator *propagator = td->fPropagator;
  // Loop all tracks, and for those with fast-simulation process, fill the 
  // proper histograms (which depend on the track position in the detector).
  // Note: there is no simple and clean way to get the initial (i.e. unsmeared)
  //       momentum/energy, therefore we set it by hand.
  int ivol;
  Volume_t const *vol;
  for ( int itr = 0; itr < npart; itr++ ) {
    // Consider only tracks with fast-simulation process (identified by the number 999).
    if ( tracks.fProcessV[itr] != 999 ) continue;
    vol = tracks.GetVolume( itr );
    #ifdef USE_VECGEOM_NAVIGATOR
      ivol = vol->id();
    #else
      ivol = vol->GetNumber();
    #endif
    //std::cout << "\t Volume=" << vol->GetName() << " ; ivol=" << ivol;  // Debug
    if ( propagator->fNthreads > 1 ) fMHist.lock();
    // Workaround for the unknown initial energy of the projectile: fixed by hand!
    double ekin_true = 50.0;  // [GeV]
    if ( isTrackerVolume[ ivol ] ) {
      //std::cout << " --> is TRACKER volume ! ";  // Debug
      double p_true = std::sqrt( ekin_true*ekin_true + 2.0*tracks.fMassV[itr]*ekin_true );
      fRatioMomentumInTracker->Fill( tracks.fPV[itr]/p_true, 1.0 );
    } else if ( isEcalVolume[ ivol ] ) {
      //std::cout << " --> is ECAL volume ! ";  // Debug
      fRatioEnergyInEcal->Fill( tracks.fEdepV[itr]/ekin_true, 1.0 );
    } else if ( isHcalVolume[ ivol ] ) {
      //std::cout << " --> is HCAL volume ! ";  // Debug
      fRatioEnergyInHcal->Fill( tracks.fEdepV[itr]/ekin_true, 1.0 );
    }
    //std::cout << std::endl;  // Debug
    if ( propagator->fNthreads > 1 ) fMHist.unlock();      
  } 
  return;  
}

//______________________________________________________________________________
void FastSimApplication::FinishRun() {
  //std::cout << "APPLICATION : FastSimApplication::FinishRun" << std::endl;  // Debug
  #ifdef  USE_ROOT
  TFile* rootFile = TFile::Open( "output.root", "RECREATE" );
  fRatioMomentumInTracker->Write();
  fRatioEnergyInEcal->Write();
  fRatioEnergyInHcal->Write();
  rootFile->Close();
  #endif
}
#endif
