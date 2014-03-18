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
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

// #include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <iostream>
#include <string>

#define RAND G4UniformRand()


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(G4int npartsPerEvent = 1)
  : _nParticlesPerEvent(npartsPerEvent), _fixedEnergy(0)
{
  particleGun  = new G4ParticleGun();
  // Detector = (DetectorConstruction*)
             // G4RunManager::GetRunManager()->GetUserDetectorConstruction();  
  
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  particleGun->SetParticleDefinition( particleTable->FindParticle("e-") );

  rndmFlag = "off";
  _trkfile = new std::ifstream("trkfile");
  // skip one line
  std::string tmp;
  std::getline(*_trkfile,tmp);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
  _trkfile->close();
  delete _trkfile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  for(int i=0; i<_nParticlesPerEvent; ++i) {
    // GenerateASingleRandomizedPrimary(anEvent);
    GenerateASinglePrimaryFromFile(anEvent);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateASinglePrimaryFromFile(G4Event* anEvent) {
  static unsigned int icount = 0;
  G4double x,y,z,px,py,pz,Ek,q,s;
  (*_trkfile) >> q >> x >> y >> z >> s >> px >> py >> pz >> Ek;

  if(icount++ % 1000 == 0) {
    std::cout <<"PrimGenAction: track #"<< icount <<" read: "<< q <<"   "<< x <<"   "<< y <<"   "<< z <<"   "<< s
	      <<"   "<< px <<"   "<< py <<"   "<< pz <<"   "<< Ek << std::endl;
  }


  G4ThreeVector position( x, y, z );
  G4double momentum = sqrt(px*px+py*py+pz*pz);
  G4double mass = particleGun->GetParticleDefinition()->GetPDGMass();

  particleGun->SetParticleMomentumDirection(G4ThreeVector( px/momentum, py/momentum, pz/momentum ) );
  particleGun->SetParticleEnergy(sqrt(momentum*momentum+mass*mass)-mass);
  particleGun->SetParticlePosition( position );
  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateASingleRandomizedPrimary(G4Event* anEvent)
{
  //this function is called at the begining of event
  G4double kinEmin=20*MeV, kinEmax=1000*MeV;
  if(_fixedEnergy>0) {
    kinEmin = _fixedEnergy;
    kinEmax = kinEmin;
  }
  G4double ecalRmin=1290*mm, ecalRmax=1520*mm, ecalZmax=3000*mm;
  G4double mass = G4ParticleTable::GetParticleTable()->FindParticle("e-")->GetPDGMass();

  // randomize position and momentum - 
  // duplicating kinematics from GXTracking/test/electronTest.cc
  G4double rho = ecalRmin + (ecalRmax-ecalRmin)*RAND;
  G4double zpos = ecalZmax*(2*RAND-1.0);
  G4double phi = twopi*RAND;
  G4double theta = atan(rho/zpos);
  if(theta<0) theta += pi;

  G4double costheta = cos(theta);
  G4double sintheta = sin(theta);
  G4double cosphi = cos(phi);
  G4double sinphi = sin(phi);

  G4double kinE = kinEmin + (kinEmax-kinEmin)*RAND;
  G4double p = sqrt(kinE*(kinE+2*mass));
  G4ThreeVector position( rho*cosphi, rho*sinphi, zpos );
  G4ThreeVector momentum( G4ThreeVector( p*sintheta*cosphi, p*sintheta*sinphi, p*costheta ) );

  particleGun->SetParticleMomentumDirection(G4ThreeVector(sintheta*cosphi,sintheta*sinphi,costheta));
  particleGun->SetParticleEnergy(sqrt(p*p+mass*mass));
  particleGun->SetParticlePosition( position );
  particleGun->GeneratePrimaryVertex(anEvent);
}
