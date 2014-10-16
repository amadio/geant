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
  static int ntrack=0;
  static int npleft=0;
  static clock_t cputime;

  // Get the navigators
  static G4Transportation* tr = dynamic_cast<G4Transportation*>(G4ProcessTable::GetProcessTable()->FindProcess("Transportation",G4Geantino::Geantino()));
  static G4EventManager *em = G4EventManager::GetEventManager();

  // We get the info at the beg of step
  G4StepPoint *startStep = step->GetPreStepPoint();

  G4LogicalVolume* vstart = startStep->GetTouchableHandle()
     ->GetVolume()->GetLogicalVolume();

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

  int begend = 0;
  if(track->GetCurrentStepNumber()==1) {
     if(io->IsNewEvent()) {
	ntrack=em->GetStackManager()->GetNTotalTrack();
	npleft=ntrack;
	cputime = clock();
     }
#if VERBOSE
     G4cout << "Track Starting " << track->GetParticleDefinition()->GetParticleName() << G4endl;
#endif
     begend = 1;
     // fill the initial step that
  }

  if(std::abs(xend.z())>11000) track->SetTrackStatus(fStopAndKill);

  if(track->GetTrackStatus()==fStopAndKill || endOfTrack) {
#if VERBOSE
     G4cout << "Track Stopping" << G4endl;
#endif
     if(begend==1) begend=3;
     else begend = 2;
     int nremain = em->GetStackManager()->GetNTotalTrack();
     if(nremain < npleft) {
	int eta = nremain*(clock()-cputime)/(ntrack-nremain)/CLOCKS_PER_SEC;
	int hours = eta/3600;
	int mins = (eta-hours*3600)/60;
	int secs = eta-hours*3600-mins*60;
	printf("ETA %2.2dh%2.2dm%2.2ds\n",hours,mins,secs);
	npleft = nremain;
     }
  }
  double snext = tr->GetLinearStepLength();

  int ivl = io->VolumeIndex(vstart->GetName().c_str());
  int iproc = io->ProcessIndex( step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName().c_str());

  io->Fill(xstart.x(), xstart.y(), xstart.z(), pstart.x(), pstart.y(), pstart.z(), 
	   track->GetParticleDefinition()->GetPDGEncoding(),
	   ivl,step->GetPreStepPoint()->GetSafety(), snext, step->GetStepLength(), 0, iproc, begend);

#if VERBOSE
  G4cout.setf(std::ios::scientific);
  G4cout << std::setprecision(2) << "# " << track->GetCurrentStepNumber() << " from " << vstart->GetName() 
	 << " to " << vendName << " p " << pstart << G4endl;
  G4cout << std::setprecision(2) << step->GetStepLength() << " from " << xstart << " to " << xend
	 << " status " << track->GetTrackStatus() << " " 
	 << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << " snext " << snext  << G4endl;
  G4cout.unsetf(std::ios::scientific);
#endif

  G4double edep = step->GetTotalEnergyDeposit();
  fEnergy += edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G01SteppingAction::Reset()
{
  fEnergy = 0.;
}

