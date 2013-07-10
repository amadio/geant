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
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"
#include "B1materials.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction* B1PrimaryGeneratorAction::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const B1PrimaryGeneratorAction* B1PrimaryGeneratorAction::Instance()
{
// Static acces function via G4RunManager 

  return fgInstance;
}      

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction() : 
   G4VUserPrimaryGeneratorAction(),
   fEnergy(1*GeV),
   fVerbose(1)
{

   //   theMessenger = new G4ParticleGunMessenger(this);

   fgInstance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  fgInstance = 0;
}

void B1PrimaryGeneratorAction::RandomDir(G4ThreeVector& direction) const {
     G4double theta = CLHEP::pi * G4UniformRand();
     G4double phi   = CLHEP::twopi * G4UniformRand();
     
     G4double sinTh=  std::sin(theta);
     G4double xdir = sinTh * std::cos(phi);
     G4double ydir = sinTh * std::sin(phi);
     G4double zdir = std::cos(theta);

     direction.set(xdir,ydir,zdir);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  static G4String nameElectron("e-");
  static G4String nameProton("proton");
  static G4ParticleDefinition *electron=particleTable->FindParticle(nameElectron);
  static G4ParticleDefinition *proton=particleTable->FindParticle(nameProton);
 
  G4PrimaryParticle *particle;
  G4PrimaryVertex *vertex;
  G4ThreeVector position;
  G4ThreeVector direction;
  for(G4int im=0; im<nmaterials; ++im) {
     position.set(MaterialPosition[im][0],MaterialPosition[im][1],MaterialPosition[im][2]);
     if(fVerbose) {
	G4cout << "Generating..." << G4endl 
	       << "vertex at " << position << G4endl;
     }
     vertex = new G4PrimaryVertex(position,0);
     
     particle = new G4PrimaryParticle(electron);
     particle->SetKineticEnergy(fEnergy);
     RandomDir(direction);
     particle->SetMomentumDirection(direction);
     if(fVerbose) {
	G4cout << particle->GetParticleDefinition()->GetParticleName() 
	       << " with momentum " << particle->GetMomentum()
	       << G4endl;
     }
     vertex->SetPrimary(particle);

     particle = new G4PrimaryParticle(proton);
     particle->SetKineticEnergy(fEnergy);
     RandomDir(direction);
     particle->SetMomentumDirection(direction);
     if(fVerbose) {
	G4cout << particle->GetParticleDefinition()->GetParticleName() 
	       << " with momentum " << particle->GetMomentum()
	       << G4endl;
     }
     vertex->SetPrimary(particle);

     anEvent->AddPrimaryVertex(vertex);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

