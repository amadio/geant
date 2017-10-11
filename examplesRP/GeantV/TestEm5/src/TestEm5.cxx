#include "TestEm5.h"

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


#include "UserDetectorConstruction.h"
#include "UserPrimaryGenerator.h"


#include <cassert>
#include <iostream>
#include <iomanip>


namespace userapplication {

TestEm5::TestEm5(Geant::GeantRunManager *runmgr, UserDetectorConstruction *det, UserPrimaryGenerator *gun)
  : Geant::GeantVApplication(runmgr), fDetector(det), fPrimaryGun(gun) {
  fHist1FileName         = "testEm5Hist1.dat";
  fInitialized           = false;
  // all these will be set properly at initialization
  fTargetLogicalVolumeID =   -1;
  fNumPrimaryPerEvent    =   -1;
  fNumBufferedEvents     =   -1;
  fHist1NumBins          =   61;
  fHist1Min              =    0.;
  fHist1Max              =   30.;
  fPrimaryParticleCharge =   -1.;
  fDataHandlerEvents     = nullptr;
  fDataHandlerRun        = nullptr;
  fData                  = nullptr;
}


TestEm5::~TestEm5() {
  if (fData) {
    delete fData;
  }
}


void TestEm5::AttachUserData(Geant::GeantTaskData *td) {
  // Create application specific thread local data structure to collecet/handle thread local multiple per-event data
  // structure. Provide number of event-slots and number of primaries per event
  TestEm5ThreadDataEvents *eventData = new TestEm5ThreadDataEvents(fNumBufferedEvents, fNumPrimaryPerEvent);
  fDataHandlerEvents->AttachUserData(eventData, td);
  // Create application specific thread local data structure to collecet/handle thread local run-global data structure.
  TestEm5ThreadDataRun *runData = new TestEm5ThreadDataRun();
  runData->CreateHisto1(fHist1NumBins, fHist1Min, fHist1Max);
  fDataHandlerRun->AttachUserData(runData, td);
}

bool TestEm5::Initialize() {
  // Initialize application. Geometry must be loaded.
  if (fInitialized)
    return true;
  // check if the detector is set and get all information that the detector should provide
  if (!fDetector) {
    Geant::Error("TestEm5::Initialize", "Geometry not available!");
    return false;
  }
  fTargetLogicalVolumeID = fDetector->GetTargetLogicalVolumeID();
  // get all information that the primary generator (simple gun) should provide
  if (!fPrimaryGun) {
    Geant::Error("TestEm5::Initialize", "PrimaryGenertaor not available!");
    return false;
  }
  fPrimaryParticleCharge = fPrimaryGun->GetPrimaryParticle()->GetPDGCharge();
  //
  // get number of primary per event and number of event-slots from Geant::GeantConfig
  fNumPrimaryPerEvent    = fRunMgr->GetConfig()->fNaverage;
  fNumBufferedEvents     = fRunMgr->GetConfig()->fNbuff;
  //
  // register thread local user data and get handler for them
  fDataHandlerEvents = fRunMgr->GetTDManager()->RegisterUserData<TestEm5ThreadDataEvents>("TestEm5ThreadDataEvents");
  fDataHandlerRun    = fRunMgr->GetTDManager()->RegisterUserData<TestEm5ThreadDataRun>("TestEm5ThreadDataRun");
  //
  // create the unique, global data struture that will be used to store cumulated per-primary data during the simulation
  fData        = new TestEm5Data();
  //
  fInitialized = true;
  return true;
}


void TestEm5::SteppingActions(Geant::GeantTrack &track, Geant::GeantTaskData *td) {
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
  TestEm5DataPerPrimary &dataPerPrimary =  (*fDataHandlerEvents)(td).GetDataPerEvent(track.fEvslot).GetDataPerPrimary(track.PrimaryParticleIndex());
  // do the scoring if the current step was done in the target logical volume
  if (idvol==fTargetLogicalVolumeID) {
    // collet charged/neutral steps that were done in the target (do not count the creation step i.e. secondary tracks
    // that has just been added in this step)
    if (track.fStatus!=Geant::kNew) {
      if (charge==0.0) {
        dataPerPrimary.AddNeutralStep();
        dataPerPrimary.AddNeutralTrackL(track.fStep);
      } else {
        dataPerPrimary.AddChargedStep();
        dataPerPrimary.AddChargedTrackL(track.fStep);
      }
      dataPerPrimary.AddEdepInTarget(track.fEdep);
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
    // transmitted: last step in target
    if (isTransmit) {
      if (isPrimary) {
        dataPerPrimary.AddPrimaryTransmitted();
      }
      if (charge==fPrimaryParticleCharge) {
        dataPerPrimary.SetOneTransmitted(1.);
      }
    }
    // reflected: last step in target
    if (isReflected) {
      if (isPrimary) {
        dataPerPrimary.AddPrimaryReflected();
      }
      if (charge==fPrimaryParticleCharge) {
        dataPerPrimary.SetOneReflected(1.);
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
    // angular distribution of transmitted charged particles: F(theta)[deg^-2] when written out
    if (isTransmit && charge!=0.) {
      // compute angular: angle of direction measured from x-dir
      double cost  = track.fXdir;
      if (cost>0.0) {
        double theta = std::acos(cost);
        double ww    = geant::degree*geant::degree;
        theta        = theta/geant::degree;
        if (theta>=fHist1Min && theta<fHist1Max) {
          // get the user defined thread local data structure for the run
          TestEm5ThreadDataRun  &dataRun =  (*fDataHandlerRun)(td);
          dataRun.GetHisto1()->Fill(theta,ww);
        }
      }
    }
  }
}


void TestEm5::FinishEvent(Geant::GeantEvent *event) {
  // merge the thread local data (filled in the SteppingActions() and distributed now in the different threads) that
  // belongs to the event (that occupied a given event-slot) that has been just transported
  TestEm5ThreadDataEvents *data = fRunMgr->GetTDManager()->MergeUserData(event->GetSlot(), *fDataHandlerEvents);
  // after the merge, we write the data into the user defined unique, global data structure. However, since more than
  // one thread can write into this global data structure, we need to protect the global data object by a lock:
  TestEm5DataPerEvent &dataPerEvent = data->GetDataPerEvent(event->GetSlot());
  fMutex.lock();
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fData->AddDataPerPrimary(dataPerEvent.GetDataPerPrimary(i));
  }
  fMutex.unlock();
  // clear the currently added ("master") thread local data (the event-slot where the currently finished event were)
  data->Clear(event->GetSlot());
  return;
}

void TestEm5::FinishRun() {
  double norm = (double)fRunMgr->GetNprimaries();
  norm        = 1./norm;
  //
  // merge the run-global thread local data from the working threads: i.e. the thread local histograms
  TestEm5ThreadDataRun *runData = fRunMgr->GetTDManager()->MergeUserData(-1, *fDataHandlerRun);
  //
  // normalize mean and mean2 values with the number of primary particles transported
  double meanChSteps   = fData->GetChargedSteps()*norm;
  double meanChSteps2  = fData->GetChargedSteps2()*norm;
  double meanNeSteps   = fData->GetNeutralSteps()*norm;
  double meanNeSteps2  = fData->GetNeutralSteps2()*norm;
  double meanChTrackL  = fData->GetChargedTrackL()*norm;
  double meanChTrackL2 = fData->GetChargedTrackL2()*norm;
  double meanNeTrackL  = fData->GetNeutralTrackL()*norm;
  double meanNeTrackL2 = fData->GetNeutralTrackL2()*norm;
  double meanNGamma    = fData->GetGammas()*norm;
  double meanNElectron = fData->GetElectrons()*norm;
  double meanNPositron = fData->GetPositrons()*norm;
  double meanPrimTrans = fData->GetPrimaryTransmitted()*norm;
  double meanPrimRefl  = fData->GetPrimaryReflecred()*norm;
  double meanOneTrans  = fData->GetOneTransmitted()*norm;
  double meanOneRefl   = fData->GetOneReflected()*norm;
  double meanEdep      = fData->GetEdepInTarget()*norm;
  double meanEdep2     = fData->GetEdepInTarget2()*norm;
  double meanELeakPr   = fData->GetELeakPrimary()*norm;
  double meanELeakPr2  = fData->GetELeakPrimary2()*norm;
  double meanELeakSec  = fData->GetELeakSecondary()*norm;
  double meanELeakSec2 = fData->GetELeakSecondary2()*norm;
  // prepare for computing sigmas
  double rmsChSteps    = meanChSteps2  - meanChSteps*meanChSteps;
  double rmsNeSteps    = meanNeSteps2  - meanNeSteps*meanNeSteps;
  double rmsChTrackL   = meanChTrackL2 - meanChTrackL*meanChTrackL;
  double rmsNeTrackL   = meanNeTrackL2 - meanNeTrackL*meanNeTrackL;
  double rmsEdep       = meanEdep2     - meanEdep*meanEdep;
  double rmsELeakPr    = meanELeakPr2  - meanELeakPr*meanELeakPr;
  double rmsELeakSec   = meanELeakSec2 - meanELeakSec*meanELeakSec;
  // compute sigmas and write it into rms..
  if (rmsChTrackL>0.) {
    rmsChTrackL  = std::sqrt(rmsChTrackL*norm);
  } else {
    rmsChTrackL  = 0.;
  }
  if (rmsNeTrackL>0.) {
    rmsNeTrackL  = std::sqrt(rmsNeTrackL*norm);
  } else {
    rmsNeTrackL  = 0.;
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
  if (rmsEdep>0.) {
    rmsEdep     = std::sqrt(rmsEdep*norm);
  } else {
    rmsEdep     = 0.;
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
  const  geantphysics::Material* targetMaterial = fDetector->GetTargetMaterial();
  double targetDensity     = targetMaterial->GetDensity();
  double targetThickness   = -2.*fDetector->GetTargetXStart();
  double meandEdx          = 0.;
  if (meanChSteps>0.) {
    meandEdx = meanEdep/targetThickness;
  }
  double meanStoppingPower = meandEdx/targetDensity;
  // get the target MaterialCuts and get the physics manager for the primary particle
  const geantphysics::MaterialCuts *targetMatCut    =
     geantphysics::MaterialCuts::GetMaterialCut(fDetector->GetTargetRegionIndex(),targetMaterial->GetIndex());
  geantphysics::PhysicsManagerPerParticle *pManager =
     fPrimaryGun->GetPrimaryParticle()->GetPhysicsManagerPerParticlePerRegion(targetMatCut->GetRegionIndex());
  // get kEnergyLoss process(es) for the primary particle (active in the target region) compute restricted and full dEdX
  double dEdXRestrComputed = 0.;
  double dEdXFullComputed  = 0.;
  if (pManager && pManager->HasEnergyLossProcess()) {
    const std::vector<geantphysics::PhysicsProcess*> &procVect =  pManager->GetListProcesses();
    for (size_t ip=0; ip<procVect.size(); ++ip) {
      if (procVect[ip]->GetType()==geantphysics::ProcessType::kEnergyLoss) {
        geantphysics::EMPhysicsProcess *emProc = static_cast<geantphysics::EMPhysicsProcess*>(procVect[ip]);
        dEdXRestrComputed += emProc->ComputeDEDX(targetMatCut, fPrimaryGun->GetPrimaryParticleEnergy(), fPrimaryGun->GetPrimaryParticle());
        dEdXFullComputed  += emProc->ComputeDEDX(targetMatCut, fPrimaryGun->GetPrimaryParticleEnergy(), fPrimaryGun->GetPrimaryParticle(),true);
      }
    }
  }
  double stpPowerRestrComputed = dEdXRestrComputed/targetDensity;
  double stpPowerFullComputed  = dEdXFullComputed/targetDensity;
  //
  //
  // printout
  std::cout<< " \n ==================================   Run summary   =========================================== \n" << std::endl;
  std::cout<< std::setprecision(3);
  std::cout<< "  The run was " << fRunMgr->GetNprimaries()                           << " "
                               << fPrimaryGun->GetPrimaryParticle()->GetName()       << " of "
                               << fPrimaryGun->GetPrimaryParticleEnergy()/geant::MeV << " [MeV] through "
                               << targetThickness/geant::um                          << " [um] of "
                               << fDetector->GetTargetMaterial()->GetName()          << " ("
                               << targetDensity/(geant::g/geant::cm3)                << " [g/cm3])"
                               << std::endl;
  std::cout<< std::endl;
  std::cout<< std::setprecision(4);
  std::cout<< "  Total energy deposit in target per event = " << meanEdep/geant::keV <<  " +- " << rmsEdep/geant::keV << " [keV] " << std::endl;
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
  std::cout<< "  Leakage :  primary = " << meanELeakPr/geant::MeV  << " +- " << rmsELeakPr/geant::MeV  << " [MeV] "
           << "  secondaries = "        << meanELeakSec/geant::MeV << " +- " << rmsELeakSec/geant::MeV << " [MeV] "
           << std::endl;
  std::cout<< "  Energy balance :  edep + eleak = " << (meanEdep + meanELeakPr + meanELeakSec)/geant::MeV << " [MeV] " << std::endl;
  std::cout<< std::endl;
  std::cout<< "  Total track length (charged) in absorber per event = " << meanChTrackL/geant::um << " +- " << rmsChTrackL/geant::um <<  " [um] "<<std::endl;
  std::cout<< "  Total track length (neutral) in absorber per event = " << meanNeTrackL/geant::um << " +- " << rmsNeTrackL/geant::um <<  " [um] "<< std::endl;
  std::cout<< std::endl;
  std::cout<< "  Number of steps (charged) in absorber per event = " << meanChSteps << " +- " << rmsChSteps << std::endl;
  std::cout<< "  Number of steps (neutral) in absorber per event = " << meanNeSteps << " +- " << rmsNeSteps << std::endl;
  std::cout<< std::endl;
  std::cout<< "  Number of secondaries per event : Gammas = " << meanNGamma <<";   electrons = " << meanNElectron << ";   positrons = " << meanNPositron << std::endl;
  std::cout<< std::endl;
  std::cout<< "  Number of events with the primary particle transmitted = " << meanPrimTrans/geant::perCent << " [%] " <<std::endl;
  std::cout<< "  Number of events with at least  1 particle transmitted (same charge as primary) = " << meanOneTrans/geant::perCent << " [%] " << std::endl;
  std::cout<< std::endl;
  std::cout<< "  Number of events with the primary particle reflected = " << meanPrimRefl/geant::perCent << " [%] " <<std::endl;
  std::cout<< "  Number of events with at least  1 particle reflected (same charge as primary) = " << meanOneRefl/geant::perCent << " [%] " << std::endl;
  std::cout<< std::endl;
  std::cout<< "  Normalised angular distribution histogram is written into file: " << fHist1FileName << std::endl;
  std::cout<< " \n ============================================================================================== \n" << std::endl;
  //
  // print the merged histogram into file
  FILE  *f        = fopen(fHist1FileName.c_str(),"w");
  Hist  *hist     = runData->GetHisto1();
  double dTheta   = hist->GetDelta();
  for (int i=0; i<hist->GetNumBins(); ++i) {
    double theta  = hist->GetX()[i];
    double factor = geant::kTwoPi*(std::cos(theta*geant::degree) - std::cos((theta+dTheta)*geant::degree));
    double val    = hist->GetY()[i];
    fprintf(f,"%d\t%lg\t%lg\n",i,theta+0.5*dTheta,val*norm/factor);
  }
  fclose(f);
}


}  // namespace userapplication
