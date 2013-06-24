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

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(6.*MeV);
  
  fgInstance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
  fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
  
  unsigned int matNo= 1;

  // Copied from Detector Construction
  const G4double kLattice = 15*cm;
  const G4double kRadius = 5*cm;
  const G4int nballs = 10;
  G4double world_sizeXY = kLattice*(nballs-1)+2*(kLattice-kRadius);
  G4double world_sizeZ = 2*(kLattice-kRadius);
  
  G4double xplace = kLattice-kRadius+(matNo%nballs)*kLattice-0.5*world_sizeXY;
  G4double yplace = kLattice-kRadius+(matNo/nballs)*kLattice-0.5*world_sizeXY;
  G4double zplace = 0.0;
  
#if 0
  G4double envSizeR = 0;
  G4Sphere* envSphere = NULL;
  
  G4LogicalVolume* envLV
  = G4LogicalVolumeStore::GetInstance()->GetVolume(materialVec[matNo]);

  if ( envLV ) envSphere = dynamic_cast<G4Sphere*>(envLV->GetSolid());
  if ( envSphere ) {
    envSizeR = envSphere->GetRadius();  // XHalfLength()*2.;
  }
  else  {
    G4cerr << "Envelope volume of 'sphere' shape not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }

  G4double size = 0.8; 
  G4double r     = size * envSizeR * G4UniformRand();
  G4double theta = CLHEP::pi * G4UniformRand();
  G4double phi   = 2.0 * CLHEP::pi * G4UniformRand();
  
  G4double rsinTh= r * std::sin(theta);
  xplace += rsinTh * std::cos(phi);
  yplace += rsinTh * std::sin(phi);
  zplace += r*std::cos(theta);
#endif
 
  fParticleGun->SetParticlePosition(G4ThreeVector(xplace, yplace, zplace));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

