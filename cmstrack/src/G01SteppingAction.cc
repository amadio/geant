//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
/// \file G01SteppingAction.cc
/// \brief Implementation of the G01SteppingAction class

#include "TGeoManager.h"

#include "time.h"
#include "VTfileio.h"
#include "G01SteppingAction.hh"
#include "G01DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4TransportationManager.hh"
#include "G4Transportation.hh"
#include "G4ProcessTable.hh"
#include "G4Geantino.hh"
#include "G4TouchableHistory.hh"
#include "G4AffineTransform.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G01SteppingAction* G01SteppingAction::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G01SteppingAction* G01SteppingAction::Instance()
{
// Static acces function via G4RunManager 

  return fgInstance;
}      

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G01SteppingAction::G01SteppingAction()
: G4UserSteppingAction(),
  fEnergy(0.)
{ 
  fgInstance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G01SteppingAction::~G01SteppingAction()
{ 
  fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G01SteppingAction::UserSteppingAction(const G4Step* step)
{
   // Get the tree to fill
  static VTfileio * io = VTfileio::I();
  static clock_t cputime;
  static clock_t cpustep=clock();
  static clock_t cputrack;
  static G4double cpsi = 1./CLOCKS_PER_SEC;
  static char volold[132]=" ";
  static int ivl;
  static int ish;

  // Get the navigators
  static G4Transportation* tr = dynamic_cast<G4Transportation*>(G4ProcessTable::GetProcessTable()->FindProcess("Transportation",G4Geantino::Geantino()));
  //  static G4EventManager *em = G4EventManager::GetEventManager();

  // We get the info at the beg of step
  G4StepPoint *startStep = step->GetPreStepPoint();
  
  //  const G4TouchableHandle touchhandle = startStep->GetTouchableHandle();
  // const G4VTouchable *touch = touchhandle();
  const G4VTouchable *touch = startStep->GetTouchableHandle()();

  G4LogicalVolume* vstart = touch->GetVolume()->GetLogicalVolume();

  G4VPhysicalVolume *v2 
     = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();

  G4bool endOfTrack = TRUE;
  G4String vendName = "OutOfWorld";
  G4LogicalVolume* vend = 0;
  if(v2) {
     vend = v2->GetLogicalVolume();
     vendName = vend->GetName();
     endOfTrack = FALSE;
  } 

  G4Track *track = step->GetTrack();
  G4ThreeVector dstart = startStep->GetMomentumDirection();
  G4ThreeVector pstart = startStep->GetMomentum();
  G4ThreeVector xstart = startStep->GetPosition();
  G4ThreeVector xend = step->GetPostStepPoint()->GetPosition();

  int trid = track->GetTrackID();
  int trpid = track->GetParentID();

  int begend = 0;
  if(track->GetCurrentStepNumber()==1) {
     cputrack = clock();
     if(io->IsNewEvent()) {
	printf("Total primaries = %d\n",io->GetPrimaries());
	cputime = clock();
	cpustep = clock();
     }
     if(trpid==0) {
	printf("#%d %-15s %g MeV ",trid,
	       track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName().c_str(),
	       track->GetTotalEnergy());
	if(trid<io->GetPrimaries()) {
	   int eta = trid*(clock()-cputime)*cpsi/(io->GetPrimaries()-trid);
	   int hours = eta/3600;
	   int mins = (eta-hours*3600)/60;
	   int secs = eta-hours*3600-mins*60;
	   printf("ETA %2dh%2.2dm%2.2ds",hours,mins,secs);
	}
	printf("\n");
     } 
#if VERBOSE
     G4cout << "Track Starting " << track->GetParticleDefinition()->GetParticleName() << G4endl;
#endif
     begend = 1;
  }

  //  if(std::abs(xend.z())>11000) track->SetTrackStatus(fStopAndKill);
  const char *gvol = vstart->GetName().c_str();
  if(strcmp(gvol,volold)) {
     ivl = io->VolumeIndex(gvol);
     ish = io->ShapeIndex(vstart->GetSolid()->GetEntityType().c_str());
     strncpy(volold,gvol,131);
  }
  //  if(ivl==2891 && track->GetKineticEnergy() < 1*CLHEP::MeV) track->SetTrackStatus(fStopAndKill);
  // if(ivl==1623 && track->GetKineticEnergy() < 1*CLHEP::MeV) track->SetTrackStatus(fStopAndKill);
  
  if(track->GetTrackStatus()==fStopAndKill || endOfTrack) {
#if VERBOSE
     G4cout << "Track Stopping" << G4endl;
#endif
     if(begend==1) begend=3;
     else begend = 2;
  }
  double snext = 0;
#if defined(SNEXTG4)
  snext = tr->GetLinearStepLength();
#endif

  int iproc = io->ProcessIndex( step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName().c_str());

  G4double safety = step->GetPreStepPoint()->GetSafety();
  io->Fill(xstart.x(), xstart.y(), xstart.z(), pstart.x(), pstart.y(), pstart.z(), 
	   track->GetParticleDefinition()->GetPDGEncoding(),
	   ivl,ish,safety, snext, step->GetStepLength(), 0, iproc, begend,
	   trid,trpid,(clock()-cputrack)*cpsi,(clock()-cpustep)*cpsi);
  cpustep=clock();

#if VERBOSE
  G4cout.setf(std::ios::scientific);
  G4cout << std::setprecision(2) << "# " << track->GetCurrentStepNumber() << " from " << vstart->GetName() 
	 << " to " << vendName << " p " << pstart << G4endl;
  G4cout << std::setprecision(2) << step->GetStepLength() << " from " << xstart << " to " << xend
	 << " status " << track->GetTrackStatus() << " " 
	 << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << " snext " << snext  << G4endl;
  G4cout.unsetf(std::ios::scientific);
#endif
#if NAVTEST
  G4double point[3];
  G4double dir[3];
  G4double norm = pstart.mag();
  if(norm>1e-10 && safety > 5e-10) {
     norm=1./norm;
     point[0]=xstart.x()*0.1;point[1]=xstart.y()*0.1;point[2]=xstart.z()*0.1;
     dir[0]=pstart.x()*norm;dir[1]=pstart.y()*norm;dir[2]=pstart.z()*norm;
     const char* rvol = gGeoManager->InitTrack(point,dir)->GetVolume()->GetName();
     //	 printf("Now in %s\n",rvol);
     if(strcmp(rvol,gvol)) {
	printf("\nROOT vol %s != Geant4 vol %s sf %9.3g sn %9.3g st %9.3g\n",rvol,gvol,safety, snext, step->GetStepLength());
	// Root part
	const char *rpath = gGeoManager->GetPath();
	printf("ROOT path %s\n",rpath);
	char rpoint[2048] = "\0";
	for ( Int_t i=0; i<gGeoManager->GetLevel(); ++i) {
	   Double_t plocal[3];
	   gGeoManager->GetMotherMatrix(gGeoManager->GetLevel()-1-i)->MasterToLocal(point,plocal);
	   strcat(rpoint,Form("/%8.03g,%8.03g,%8.03g",plocal[0],plocal[1],plocal[2]));
	}
	printf("%s\n",rpoint);
	const G4NavigationHistory *navhistory = touch->GetHistory();
	char gpath[1024] = "\0";
	char gpoint[2048] = "\0";
	for ( G4int i=0; i<=navhistory->GetDepth(); i++ ) {
	   if(navhistory->GetVolume(i)) {
	      navhistory->GetVolume(i)->CheckOverlaps(1000,1e-3,FALSE,1);
	      strcat(gpath,"/");
	      strcat(gpath,navhistory->GetVolume(i)->GetName().c_str());
	      strcat(gpath,"_");
	      switch(navhistory->GetVolumeType(i)) {
	      case kNormal:
		 strcat(gpath,"n");
		 break;
	      case kReplica:
		 strcat(gpath,"n");
		 break;
	      case kParameterised:
		 strcat(gpath,"p");
		 break;
	      }
	      strcat(gpath,Form("%d",navhistory->GetReplicaNo(i)));
	      G4AffineTransform gtr = navhistory->GetTransform(i);
	      G4ThreeVector ptr = gtr.TransformPoint(xstart);
	      strcat(gpoint,Form("/%8.03g,%8.03g,%8.03g",ptr.x()*.1,ptr.y()*.1,ptr.z()*.1));
	   } else 
	      strcat(gpath,"/<Null>");
	}
	printf("G4   path %s\n",gpath);
	printf("%s\n",gpoint);
	//	G4cout << "G4 path " << *navhistory << " Root path " << rpath << G4endl;
     } else {
	//	    printf("%f %f %f %f\n",fPx,fPy,fPz,norm);
	gGeoManager->FindNextBoundaryAndStep();
	//	hstep->Fill(gGeoManager->GetStep()*10-fStep);
     }
  }
#endif

  G4double edep = step->GetTotalEnergyDeposit();
  fEnergy += edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G01SteppingAction::Reset()
{
  fEnergy = 0.;
}

