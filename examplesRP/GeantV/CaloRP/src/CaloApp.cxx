//===--- CaloApp.cxx - GeantV ---------------------------------*- C++ -*-===//
// @file CaloDetectorConstruction.h
// @brief Implementation of application initialization and scoring for GeantV calorimeter prototype
// @author Ryan Schmitz, derived from example by M Novak
// @date August 1, 2017
//
//     //===----------------------------------------------------------------------===//
//
#include "CaloApp.h"

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
#include "PhysicsData.h"
#include "MaterialCuts.h"
#include "Material.h"
#include "LightTrack.h"
#include "PhysicsProcess.h"
#include "EMPhysicsProcess.h"
#include "PhysicsManagerPerParticle.h"


#include "CaloDetectorConstruction.h"
#include "CaloPrimaryGenerator.h"


#include <cassert>
#include <iostream>
#include <iomanip>

#ifdef USE_ROOT
  #include "TFile.h"
#endif


namespace userapplication {

CaloApp::CaloApp(Geant::GeantRunManager *runmgr, CaloDetectorConstruction *det, CaloPrimaryGenerator *gun)
  : Geant::GeantVApplication(runmgr), fDetector(det), fPrimaryGun(gun) {
  fHist1FileName         = "CaloAppHist";
  fInitialized           = false;
  // all these will be set properly at initialization
  fNumPrimaryPerEvent    =   -1;
  fNumBufferedEvents     =   -1;
  fHist1NumBins          =   100;
  fHist1Min              =    0.;
  fHist1Max              =   30.;
  fPrimaryParticleCharge =   -1.;
  fDataHandlerEvents     = nullptr;
  fDataHandlerRun        = nullptr;
  fData                  = nullptr;
}


CaloApp::~CaloApp() {
  if (fData) {
    delete fData;
  }
}


void CaloApp::AttachUserData(Geant::GeantTaskData *td) {
  // Create application specific thread local data structure to collecet/handle thread local multiple per-event data
  // structure. Provide number of event-slots and number of primaries per event
  CaloAppThreadDataEvents *eventData = new CaloAppThreadDataEvents(fNumBufferedEvents, fNumPrimaryPerEvent);
  fDataHandlerEvents->AttachUserData(eventData, td);
  // Create application specific thread local data structure to collecet/handle thread local run-global data structure.
  CaloAppThreadDataRun *runData = new CaloAppThreadDataRun();
  runData->CreateHisto1(fHist1NumBins, fHist1Min, fHist1Max);
  fDataHandlerRun->AttachUserData(runData, td);
}

bool CaloApp::Initialize() {
  // Initialize application. Geometry must be loaded.
  if (fInitialized)
    return true;
  // check if the detector is set and get all information that the detector should provide
  if (!fDetector) {
    Geant::Error("CaloApp::Initialize", "Geometry not available!");
    return false;
  }

  //get detector parameters/logical volume IDs
  fNumAbsorbers = fDetector->GetNumAbsorbers();
  for (int k=0; k<fNumAbsorbers; k++){
  	fAbsorberLogicalVolumeID[k] = fDetector->GetAbsorberLogicalVolumeID(k);
  }

  // get all information that the primary generator (simple gun) should provide
  if (!fPrimaryGun) {
    Geant::Error("CaloApp::Initialize", "PrimaryGenerator not available!");
    return false;
  }
  fPrimaryParticleCharge = fPrimaryGun->GetPrimaryParticle()->GetPDGCharge();
  //
  // get number of primary per event and number of event-slots from Geant::GeantConfig
  fNumPrimaryPerEvent    = fRunMgr->GetConfig()->fNaverage;
  fNumBufferedEvents     = fRunMgr->GetConfig()->fNbuff;
  //
  // register thread local user data and get handler for them
  fDataHandlerEvents = fRunMgr->GetTDManager()->RegisterUserData<CaloAppThreadDataEvents>("CaloAppThreadDataEvents");
  fDataHandlerRun    = fRunMgr->GetTDManager()->RegisterUserData<CaloAppThreadDataRun>("CaloAppThreadDataRun");
  //
  // create the unique, global data struture that will be used to store cumulated per-primary data during the simulation
  fData        = new CaloAppData();
  //
  // CREATE PhysicsData here: should be done at the init of PhysicsProcessHandler but
  // GeantTaskData are constructed later than that call
  for (int i=0; i<fRunMgr->GetNthreadsTotal(); ++i) {
    fRunMgr->GetTDManager()->GetTaskData(i)->fPhysicsData = new geantphysics::PhysicsData();
  }
  //
  fInitialized = true;
  return true;
}


void CaloApp::SteppingActions(Geant::GeantTrack &track, Geant::GeantTaskData *td) {
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
  const geantphysics::Particle *part = geantphysics::Particle::GetParticleByInternalCode(track.fGVcode);
  int    pdgCode = part->GetPDGCode();
  double  charge = part->GetPDGCharge();
  double  ekin   = track.fE-track.fMass;
  //
  bool isTransmit  = ( (track.fXdir>0. && track.fXpos>0.0 && track.fStatus==Geant::kBoundary) && ( ekin > 0.0 ) );
  bool isReflected = ( (track.fXdir<0. && track.fXpos<0.0 && track.fStatus==Geant::kBoundary) && ( ekin > 0.0 ) );
  bool isPrimary   = ( track.fGeneration==0 );
  //
  // get the user defined thread local data structure per-primary particle for: the event-slot index (that defines the
  // per-event data structure) and the primary index (that defines the per-primary data structure within that per-event
  // data structure). NOTE: each tracks stores the event-slot and primary partcile index that event and primary particle
  // within that event the track belongs to.
  CaloAppDataPerPrimary &dataPerPrimary =  (*fDataHandlerEvents)(td).GetDataPerEvent(track.fEvslot).GetDataPerPrimary(track.PrimaryParticleIndex());
  // do the scoring if the current step was done in the target logical volume

  int currentAbsorber=0;
  bool validVolume=false;

  for (int k=0; k<fNumAbsorbers; k++){
	if (idvol==fAbsorberLogicalVolumeID[k]){
		currentAbsorber=k;
		validVolume=true;
		break;
	}
  }


  if (validVolume) {
    // collet charged/neutral steps that were done in the target (do not count the creation step i.e. secondary tracks
    // that has just been added in this step)
    if (track.fStatus!=Geant::kNew) {
      if (charge==0.0) {
        dataPerPrimary.AddNeutralStep();
        dataPerPrimary.AddNeutralTrackL(track.fStep,currentAbsorber);
      } else {
        dataPerPrimary.AddChargedStep();
        dataPerPrimary.AddChargedTrackL(track.fStep,currentAbsorber);
      }
      dataPerPrimary.AddEdepInAbsorber(track.fEdep,currentAbsorber);
    }
    // collect secondary particle type statistics
    if (track.fStatus==Geant::kNew) {
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
    // energy leakage
    if (isTransmit || isReflected) {
      double energyLeak = ekin;
      // e+ created during the simulation so add its 2 e- rest mass energy
      if (!isPrimary && pdgCode==-11) {
        energyLeak += 2.*geant::kElectronMassC2;
      }
      if (isPrimary) {
        dataPerPrimary.AddELeakPrimary(energyLeak);
      } else {
        dataPerPrimary.AddELeakSecondary(energyLeak);
      }
    }
    // add energy deposited in first absorber (per primary) to histogram
    if (isTransmit && isPrimary) {
          // get the user defined thread local data structure for the run
          CaloAppThreadDataRun  &dataRun =  (*fDataHandlerRun)(td);
          dataRun.GetHisto1()->Fill(dataPerPrimary.GetEdepInAbsorber(0));
    }
  }
}


void CaloApp::FinishEvent(Geant::GeantEvent *event) {
  // merge the thread local data (filled in the SteppingActions() and distributed now in the different threads) that
  // belongs to the event (that occupied a given event-slot) that has been just transported
  CaloAppThreadDataEvents *data = fRunMgr->GetTDManager()->MergeUserData(event->GetSlot(), *fDataHandlerEvents);
  // after the merge, we write the data into the user defined unique, global data structure. However, since more than
  // one thread can write into this global data structure, we need to protect the global data object by a lock:
  CaloAppDataPerEvent &dataPerEvent = data->GetDataPerEvent(event->GetSlot());
  fMutex.lock();
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fData->AddDataPerPrimary(dataPerEvent.GetDataPerPrimary(i));
  }
  fMutex.unlock();
  // clear the currently added ("master") thread local data (the event-slot where the currently finished event were)
  data->Clear(event->GetSlot());
  return;
}

void CaloApp::FinishRun() {
  double norm = (double)fRunMgr->GetNprimaries();
  norm        = 1./norm;
  //
  // merge the run-global thread local data from the working threads: i.e. the thread local histograms
  CaloAppThreadDataRun *runData = fRunMgr->GetTDManager()->MergeUserData(-1, *fDataHandlerRun);
  //
  // normalize mean and mean2 values with the number of primary particles transported
  double meanChSteps   = fData->GetChargedSteps()*norm;
  double meanChSteps2  = fData->GetChargedSteps2()*norm;
  double meanNeSteps   = fData->GetNeutralSteps()*norm;
  double meanNeSteps2  = fData->GetNeutralSteps2()*norm;

  double meanChTrackL[fNumAbsorbers+1];
  double meanChTrackL2[fNumAbsorbers+1];
  double meanNeTrackL[fNumAbsorbers+1];
  double meanNeTrackL2[fNumAbsorbers+1];
  double meanEdep[fNumAbsorbers+1];
  double meanEdep2[fNumAbsorbers+1];

  for (int k=1;k<=fNumAbsorbers;k++){
	  meanChTrackL[k]  = fData->GetChargedTrackL(k)*norm;
	  meanChTrackL2[k] = fData->GetChargedTrackL2(k)*norm;
	  meanNeTrackL[k]  = fData->GetNeutralTrackL(k)*norm;
	  meanNeTrackL2[k] = fData->GetNeutralTrackL2(k)*norm;
	  meanEdep[k]      = fData->GetEdepInAbsorber(k)*norm;
	  meanEdep2[k]     = fData->GetEdepInAbsorber2(k)*norm;
  }

  double meanNGamma    = fData->GetGammas()*norm;
  double meanNElectron = fData->GetElectrons()*norm;
  double meanNPositron = fData->GetPositrons()*norm;

  double meanELeakPr   = fData->GetELeakPrimary()*norm;
  double meanELeakPr2  = fData->GetELeakPrimary2()*norm;
  double meanELeakSec  = fData->GetELeakSecondary()*norm;
  double meanELeakSec2 = fData->GetELeakSecondary2()*norm;
  // prepare for computing sigmas
  double rmsChSteps    = meanChSteps2  - meanChSteps*meanChSteps;
  double rmsNeSteps    = meanNeSteps2  - meanNeSteps*meanNeSteps;

  double rmsChTrackL[fNumAbsorbers+1];
  double rmsNeTrackL[fNumAbsorbers+1];
  double rmsEdep[fNumAbsorbers+1];

  for (int k=0; k<fNumAbsorbers; k++){
	rmsChTrackL[k]   = meanChTrackL2[k] - meanChTrackL[k]*meanChTrackL[k];
	rmsNeTrackL[k]   = meanNeTrackL2[k] - meanNeTrackL[k]*meanNeTrackL[k];
	rmsEdep[k]       = meanEdep2[k]     - meanEdep[k]*meanEdep[k];
  }
  double rmsELeakPr    = meanELeakPr2  - meanELeakPr*meanELeakPr;
  double rmsELeakSec   = meanELeakSec2 - meanELeakSec*meanELeakSec;
  // compute sigmas and write it into rms..
  for (int k=0; k<fNumAbsorbers; k++){
	  if (rmsChTrackL[k]>0.) {
	    rmsChTrackL[k]  = std::sqrt(rmsChTrackL[k]*norm);
	  } else {
	    rmsChTrackL[k]  = 0.;
	  }
	  if (rmsNeTrackL[k]>0.) {
	    rmsNeTrackL[k]  = std::sqrt(rmsNeTrackL[k]*norm);
	  } else {
	    rmsNeTrackL[k]  = 0.;
	  }
	  if (rmsEdep[k]>0.) {
	    rmsEdep[k]     = std::sqrt(rmsEdep[k]*norm);
	  } else {
	    rmsEdep[k]     = 0.;
	  }
  }

  if (rmsChSteps>0.) {
    rmsChSteps  = std::sqrt(rmsChSteps*norm);
  } else {
    rmsChSteps  = 0.;
  }
  if (rmsNeSteps>0.) {
    rmsNeSteps  = std::sqrt(rmsNeSteps*norm);
  } else {
    rmsNeSteps  = 0.;
  }
  if (rmsELeakPr>0.) {
    rmsELeakPr  = std::sqrt(rmsELeakPr*norm);
  } else {
    rmsELeakPr  = 0.;
  }
  if (rmsELeakSec>0.) {
    rmsELeakSec = std::sqrt(rmsELeakSec*norm);
  } else {
    rmsELeakSec = 0.;
  }
  //
  // some additional quantities:
  // 
for(int k=0; k<fNumAbsorbers; k++){
///////////////////////for loop start, looping over absorbers///////////////////////////////////////////
  std::cout << "---------------Started data analysis for " << fDetector->GetAbsorberMaterialName(k) << " absorber---------------------" << std::endl;
  const  geantphysics::Material* absorberMaterial = fDetector->GetAbsorberMaterial(k);
  double absorberDensity     = absorberMaterial->GetDensity();
  double absorberThickness   = fDetector->GetAbsorberThickness(k)*fDetector->GetNumLayers();
  double meandEdx          = 0.;
  if (meanChSteps>0.) {
    meandEdx = meanEdep[k]/absorberThickness;
  }
  double meanStoppingPower = meandEdx/absorberDensity;
  // get the target MaterialCuts and get the physics manager for the primary particle
  const geantphysics::MaterialCuts *absorberMatCut    =
     geantphysics::MaterialCuts::GetMaterialCut(fDetector->GetDetectorRegionIndex(),absorberMaterial->GetIndex());
  geantphysics::PhysicsManagerPerParticle *pManager =
     fPrimaryGun->GetPrimaryParticle()->GetPhysicsManagerPerParticlePerRegion(absorberMatCut->GetRegionIndex());
  // get kEnergyLoss process(es) for the primary particle (active in the target region) compute restricted and full dEdX
  double dEdXRestrComputed = 0.;
  double dEdXFullComputed  = 0.;
  if (pManager->HasEnergyLossProcess()) {
    const std::vector<geantphysics::PhysicsProcess*> &procVect =  pManager->GetListProcesses();
    for (size_t ip=0; ip<procVect.size(); ++ip) {
      if (procVect[ip]->GetType()==geantphysics::ProcessType::kEnergyLoss) {
        geantphysics::EMPhysicsProcess *emProc = static_cast<geantphysics::EMPhysicsProcess*>(procVect[ip]);
        dEdXRestrComputed += emProc->ComputeDEDX(absorberMatCut, fPrimaryGun->GetPrimaryParticleEnergy(), fPrimaryGun->GetPrimaryParticle());
        dEdXFullComputed  += emProc->ComputeDEDX(absorberMatCut, fPrimaryGun->GetPrimaryParticleEnergy(), fPrimaryGun->GetPrimaryParticle(),true);
      }
    }
  }
  double stpPowerRestrComputed = dEdXRestrComputed/absorberDensity;
  double stpPowerFullComputed  = dEdXFullComputed/absorberDensity;
  //
  //
  // printout
  std::cout<< " \n ==================================   Run summary   =========================================== \n" << std::endl;
  std::cout<< std::setprecision(3);
  std::cout<< "  The run was " << fRunMgr->GetNprimaries()                           << " "
                               << fPrimaryGun->GetPrimaryParticle()->GetName()       << " of "
                               << fPrimaryGun->GetPrimaryParticleEnergy()/geant::MeV << " [MeV] through "
                               << absorberThickness/geant::um                          << " [um] of "
                               << fDetector->GetAbsorberMaterialName(k)          << " ("
                               << absorberDensity/(geant::g/geant::cm3)                << " [g/cm3])"
                               << std::endl;
  std::cout<< std::endl;
  std::cout<< std::setprecision(4);
  std::cout<< "  Total energy deposit in this absorber type per event = " << meanEdep[k]/geant::keV <<  " +- " << rmsEdep[k]/geant::keV << " [keV] " << std::endl;
  std::cout<< std::endl;
  std::cout<< "  -----> Mean dE/dx = " << meandEdx/(geant::MeV/geant::cm)                        << " [MeV/cm]     ("
                                       << meanStoppingPower/(geant::MeV*geant::cm2/geant::g)     << " [MeV*cm2/g]) "
                                       << std::endl;
  std::cout<< std::endl;
  std::cout<< "  From formulas : " << std::endl
           << "    restricted dEdx = " << dEdXRestrComputed/(geant::MeV/geant::cm)               << " [MeV/cm]     ("
                                       << stpPowerRestrComputed/(geant::MeV*geant::cm2/geant::g) << " [MeV*cm2/g]) "
                                       << std::endl;
  std::cout<< "    full dEdx       = " << dEdXFullComputed/(geant::MeV/geant::cm)                << " [MeV/cm]     ("
                                       << stpPowerFullComputed/(geant::MeV*geant::cm2/geant::g)  << " [MeV*cm2/g]) "
                                       << std::endl;
  std::cout<< std::endl;
//  std::cout<< "  Energy balance :  edep + eleak = " << (meanEdep[k] + meanELeakPr + meanELeakSec)/geant::MeV << " [MeV] " << std::endl;
//  std::cout<< std::endl;
  std::cout<< "  Total track length (charged) in absorber per event = " << meanChTrackL[k]/geant::um << " +- " << rmsChTrackL[k]/geant::um <<  " [um] "<<std::endl;
  std::cout<< "  Total track length (neutral) in absorber per event = " << meanNeTrackL[k]/geant::um << " +- " << rmsNeTrackL[k]/geant::um <<  " [um] "<< std::endl;
  std::cout<< std::endl;
  std::cout<< " \n ============================================================================================== \n" << std::endl;
  //
  
std::cout << "---------------Finished data analysis for " << fDetector->GetAbsorberMaterialName(k) << " absorber-------------------" << std::endl;
}
//////////////////////////////////////////////end for loop ///////////////////////////////////////////////
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
  std::cout<< "  Direct Energy Deposition of Primaries histogram is written into file " << filename <<  std::endl;
  std::cout << "\n==================== General Sim Results =====================================\n\n";
  std::cout<< "  Leakage :  primary = " << meanELeakPr/geant::MeV  << " +- " << rmsELeakPr/geant::MeV  << " [MeV] "
           << "  secondaries = "        << meanELeakSec/geant::MeV << " +- " << rmsELeakSec/geant::MeV << " [MeV] "
           << std::endl;
  std::cout<< "  Number of steps (charged) in absorber per event = " << meanChSteps << " +- " << rmsChSteps << std::endl;
  std::cout<< "  Number of steps (neutral) in absorber per event = " << meanNeSteps << " +- " << rmsNeSteps << std::endl;
  std::cout<< std::endl;
  std::cout<< "  Number of secondaries per event : Gammas = " << meanNGamma <<";   electrons = " << meanNElectron << ";   positrons = " << meanNPositron << std::endl;
  std::cout<< std::endl;
}

}  // namespace userapplication
