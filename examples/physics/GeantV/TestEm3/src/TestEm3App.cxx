
#include "TestEm3App.h"

// vecgeom::GeoManager
#include "management/GeoManager.h"
//
#include "GeantEvent.h"
#include "GeantTrackVec.h"
#include "GeantRunManager.h"
#include "GeantTaskData.h"
#include "globals.h"

#include "Geant/Error.h"

#include "SystemOfUnits.h"

#include "TestEm3DetectorConstruction.h"
#include "TestEm3PrimaryGenerator.h"


#include <cassert>
#include <iostream>
#include <iomanip>

/*
#ifdef USE_ROOT
  #include "TFile.h"
#endif
*/

namespace userapplication {

TestEm3App::TestEm3App(Geant::GeantRunManager *runmgr, TestEm3DetectorConstruction *det, TestEm3PrimaryGenerator *gun)
  : Geant::GeantVApplication(runmgr), fDetector(det), fPrimaryGun(gun) {
  fIsPerformance         = false;
  fInitialized           = false;
  // all these will be set properly at initialization
  fNumPrimaryPerEvent    =   -1;
  fNumBufferedEvents     =   -1;
  fPrimaryParticleCharge =   -1.;
  fDataHandlerEvents     = nullptr;
//  fDataHandlerRun        = nullptr;
  fData                  = nullptr;
}


TestEm3App::~TestEm3App() {
  if (fData) {
    delete fData;
  }
}


void TestEm3App::AttachUserData(Geant::GeantTaskData *td) {
  if (fIsPerformance)
   return;
  // Create application specific thread local data structure to collecet/handle thread local multiple per-event data
  // structure. Provide number of event-slots and number of primaries per event
  TestEm3ThreadDataEvents *eventData = new TestEm3ThreadDataEvents(fNumBufferedEvents, fNumPrimaryPerEvent, fNumAbsorbers);
  fDataHandlerEvents->AttachUserData(eventData, td);
}

bool TestEm3App::Initialize() {
  if (fIsPerformance)
    return true;
  // Initialize application. Geometry must be loaded.
  if (fInitialized)
    return true;
  // check if the detector is set and get all information that the detector should provide
  if (!fDetector) {
    Geant::Error("TestEm3App::Initialize", "Geometry not available!");
    return false;
  }

  //get detector parameters/logical volume IDs
  fNumAbsorbers = fDetector->GetNumberOfAbsorbersPerLayer();
  fAbsorberLogicalVolumeID.resize(fNumAbsorbers+3,-1);    // +3 logical volume (world,calo,layer)
  for (int k=0; k<fNumAbsorbers; k++){
    int id = fDetector->GetAbsorberLogicalVolumeID(k);
  	fAbsorberLogicalVolumeID[id] = k;  // map to absorber index
  }
  // get all information that the primary generator (simple gun) should provide
  if (!fPrimaryGun) {
    Geant::Error("TestEm3App::Initialize", "PrimaryGenerator not available!");
    return false;
  }
  fPrimaryParticleCharge = fPrimaryGun->GetPrimaryParticle()->GetPDGCharge();
  //
  // get number of primary per event and number of event-slots from Geant::GeantConfig
  fNumPrimaryPerEvent    = fRunMgr->GetConfig()->fNaverage;
  fNumBufferedEvents     = fRunMgr->GetConfig()->fNbuff;
  //
  // register thread local user data and get handler for them
  fDataHandlerEvents = fRunMgr->GetTDManager()->RegisterUserData<TestEm3ThreadDataEvents>("TestEm3ThreadDataEvents");
  //
  // create the unique, global data struture that will be used to store cumulated per-primary data during the simulation
  fData        = new TestEm3Data(fNumAbsorbers);
  //
  fInitialized = true;
  return true;
}


void TestEm3App::SteppingActions(Geant::GeantTrack &track, Geant::GeantTaskData *td) {
  if (fIsPerformance)
   return;
  // it is still a bit tricky but try to get the ID of the logical volume in which the current step was done
  Node_t const *current;
  int idvol = -1;
  int ilev = -1;
  ilev = track.Path()->GetCurrentLevel() - 1;
  if (ilev<1) {
    return;
  }
  current = track.Path()->Top();
  if (!current) {
    return;
  }
  idvol  = current->GetLogicalVolume()->id();
  // get some particle properties
  const geantphysics::Particle *part = geantphysics::Particle::GetParticleByInternalCode(track.GVcode());
  int    pdgCode = part->GetPDGCode();
  double  charge = part->GetPDGCharge();
//  double  ekin   = track.fE-track.fMass;
//  bool isPrimary = ( track.fGeneration==0 );
  //
  // get the user defined thread local data structure per-primary particle for: the event-slot index (that defines the
  // per-event data structure) and the primary index (that defines the per-primary data structure within that per-event
  // data structure). NOTE: each tracks stores the event-slot and primary partcile index that event and primary particle
  // within that event the track belongs to.
  TestEm3DataPerPrimary &dataPerPrimary =  (*fDataHandlerEvents)(td).GetDataPerEvent(track.EventSlot()).GetDataPerPrimary(track.PrimaryParticleIndex());
  // do the scoring if the current step was done in the target logical volume

  int currentAbsorber=fAbsorberLogicalVolumeID[idvol];
  if (currentAbsorber>-1) {
    // collet charged/neutral steps that were done in the target (do not count the creation step i.e. secondary tracks
    // that has just been added in this step)
    if (track.Status() != Geant::kNew) {
      if (charge==0.0) {
        dataPerPrimary.AddNeutralStep();
        dataPerPrimary.AddNeutralTrackL(track.GetStep(), currentAbsorber);
      } else {
        dataPerPrimary.AddChargedStep();
        dataPerPrimary.AddChargedTrackL(track.GetStep(), currentAbsorber);
      }
      dataPerPrimary.AddEdepInAbsorber(track.Edep(), currentAbsorber);
    }
    // collect secondary particle type statistics
    if (track.Status() == Geant::kNew) {
      switch(pdgCode) {
        // gamma
        case  22 : dataPerPrimary.AddGamma();
                   break;
        // e
        case  11 : dataPerPrimary.AddElectron();
                  break;
        // e+
        case -11 : dataPerPrimary.AddPositron();
                   break;
      }
    }
  }
}


void TestEm3App::FinishEvent(Geant::GeantEvent *event) {
  if (fIsPerformance)
   return;
  // merge the thread local data (filled in the SteppingActions() and distributed now in the different threads) that
  // belongs to the event (that occupied a given event-slot) that has been just transported
  TestEm3ThreadDataEvents *data = fRunMgr->GetTDManager()->MergeUserData(event->GetSlot(), *fDataHandlerEvents);
  // after the merge, we write the data into the user defined unique, global data structure. However, since more than
  // one thread can write into this global data structure, we need to protect the global data object by a lock:
  TestEm3DataPerEvent &dataPerEvent = data->GetDataPerEvent(event->GetSlot());
  fMutex.lock();
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fData->AddDataPerPrimary(dataPerEvent.GetDataPerPrimary(i));
  }
  fMutex.unlock();
  // clear the currently added ("master") thread local data (the event-slot where the currently finished event were)
  data->Clear(event->GetSlot());
  return;
}

void TestEm3App::FinishRun() {
  if (fIsPerformance)
   return;
  double norm = (double)fRunMgr->GetNprimaries();
  norm        = 1./norm;
  //
  // normalize mean and mean2 values with the number of primary particles transported
  double meanChSteps   = fData->GetChargedSteps()*norm;
//  double meanChSteps2  = fData->GetChargedSteps2()*norm;
  double meanNeSteps   = fData->GetNeutralSteps()*norm;
//  double meanNeSteps2  = fData->GetNeutralSteps2()*norm;

  double meanChTrackL[fNumAbsorbers];
  double meanChTrackL2[fNumAbsorbers];
//  double meanNeTrackL[fNumAbsorbers];
//  double meanNeTrackL2[fNumAbsorbers];
  double meanEdep[fNumAbsorbers];
  double meanEdep2[fNumAbsorbers];

  for (int k=0; k<fNumAbsorbers; k++){
	  meanChTrackL[k]  = fData->GetChargedTrackL(k)*norm;
	  meanChTrackL2[k] = fData->GetChargedTrackL2(k)*norm;
//	  meanNeTrackL[k]  = fData->GetNeutralTrackL(k)*norm;
//	  meanNeTrackL2[k] = fData->GetNeutralTrackL2(k)*norm;
	  meanEdep[k]      = fData->GetEdepInAbsorber(k)*norm;
	  meanEdep2[k]     = fData->GetEdepInAbsorber2(k)*norm;
  }

  double meanNGamma    = fData->GetGammas()*norm;
  double meanNElectron = fData->GetElectrons()*norm;
  double meanNPositron = fData->GetPositrons()*norm;

  // prepare for computing sigmas
//  double rmsChSteps    = std::sqrt(std::abs(meanChSteps2  - meanChSteps*meanChSteps));
//  double rmsNeSteps    = std::sqrt(std::abs(meanNeSteps2  - meanNeSteps*meanNeSteps));

  double rmsChTrackL[fNumAbsorbers];
//  double rmsNeTrackL[fNumAbsorbers];
  double rmsEdep[fNumAbsorbers];

  for (int k=0; k<fNumAbsorbers; k++){
  	rmsChTrackL[k]   = std::sqrt(std::abs(meanChTrackL2[k] - meanChTrackL[k]*meanChTrackL[k]));
//  	rmsNeTrackL[k]   = std::sqrt(std::abs(meanNeTrackL2[k] - meanNeTrackL[k]*meanNeTrackL[k]));
  	rmsEdep[k]       = std::sqrt(std::abs(meanEdep2[k]     - meanEdep[k]*meanEdep[k]));
  }
  //
  // printout
  std::cout<< " \n ==================================   Run summary   =========================================== \n" << std::endl;
  std::cout<< std::setprecision(3);
  std::cout<< "  The run was " << fRunMgr->GetNprimaries()                           << " "
                               << fPrimaryGun->GetPrimaryParticle()->GetName()       << " of with E = "
                               << fPrimaryGun->GetPrimaryParticleEnergy()            << " [GeV] "
                               << std::endl;
  std::cout<< std::endl;
  std::cout<< std::setprecision(4);
  std::cout<< " \n ---------------------------------------------------------------------------------------------- \n" << std::endl;
  std::cout<< "    material        Edep         RMS     sqrt(E0(GeV))*rmsE/Emean      track-length  \n" << std::endl;
  double eprim = fPrimaryGun->GetPrimaryParticleEnergy();
  for (int k=0; k<fNumAbsorbers; k++) {
    const std::string absMatName = fDetector->GetAbsorberMaterial(k)->GetName();
    double absMeanEdep           = meanEdep[k]/geant::MeV;
    double absEdepRMS            = rmsEdep[k]/geant::MeV;
    double absMeanChTrl          = meanChTrackL[k]/geant::cm;
    double absChTrlRMS           = rmsChTrackL[k]/geant::cm;
//    double absEvis               = meanEdep[k]/eprim;
    double absResolution         = 0.;
    if (meanEdep[k]>0.) {
      absResolution = 100.*std::sqrt(eprim/geant::GeV)*rmsEdep[k]/meanEdep[k];
    }
    double absResolutionRMS     = absResolution*std::sqrt(norm);
    // print out
    std::cout << "  " <<std::setw(20) << absMatName << ": "
      << std::setprecision(6)
      << std::setw(6)  << absMeanEdep      << " MeV " << " :  "
      << std::setprecision(6)
      << std::setw(5)  << absEdepRMS       << " MeV "
      << std::setw(5) << absResolution     << " +- "
      << std::setw(5)  << absResolutionRMS << " %  "
      << std::setprecision(6)
      << std::setw(5) << absMeanChTrl      << " cm "   << " +- "
      << std::setw(4)  << absChTrlRMS      << " cm "
      << std::endl;
  }
  std::cout<< " \n ---------------------------------------------------------------------------------------------- \n" << std::endl;

  std::cout << std::setprecision(6);
  std::cout << " Mean number of gamma          " << meanNGamma    << std::endl;
  std::cout << " Mean number of e-             " << meanNElectron << std::endl;
  std::cout << " Mean number of e+             " << meanNPositron << std::endl;
  std::cout << std::setprecision(6)
            << " Mean number of charged steps  " << meanChSteps << std::endl;
  std::cout << " Mean number of neutral steps  " << meanNeSteps << std::endl;

  std::cout<< std::endl;
  std::cout<< " \n ============================================================================================== \n" << std::endl;
/*
  //
  // print the merged histogram into file
  std::string filename(fHist1FileName);
#ifdef USE_ROOT
  //ROOT-style TH1F output histogram of energy depositions by primaries
  filename.append(".root");
  TFile *file = new TFile(filename.c_str(),"RECREATE");
  TH1F  *rootHist = runData->GetHisto1();
  rootHist->Write();
  file->Close();
#else
  //ASCII-style histogram of energy depositions by primaries
  filename.append(".dat");
  FILE  *f        = fopen(filename.c_str(),"w");
  Hist  *hist     = runData->GetHisto1();
  double dEDep   = hist->GetDelta();
  for (int i=0; i<hist->GetNumBins(); ++i) {
    double EDep  = hist->GetX()[i];
    double val    = hist->GetY()[i];
    fprintf(f,"%d\t%lg\t%lg\n",i,EDep+0.5*dEDep,val*norm); // norm = 1/nPrimaries, so the resulting plot is normalized to number of primaries
  }
  fclose(f);
#endif
*/
}

}  // namespace userapplication
