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
/// \file persistency/gdml/G01/src/G01PrimaryGeneratorAction.cc
/// \brief Implementation of the G01PrimaryGeneratorAction class
//
//
// $Id: G01PrimaryGeneratorAction.cc 68025 2013-03-13 13:43:46Z gcosmo $
//
//

#include "G01PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "Pythia8/Pythia.h"

#include "VTfileio.h"

using namespace Pythia8;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static Pythia pythia;

G01PrimaryGeneratorAction::G01PrimaryGeneratorAction()
   : G4VUserPrimaryGeneratorAction(), fEnergy(1*GeV), fVerbose(FALSE)
{
   pythia.readString("Beams:eCM = 7000.");
   pythia.readString("HardQCD:all = on");
   pythia.readString("PhaseSpace:pTHatMin = 20.");
   pythia.readString("Random:setSeed = on");
   pythia.readString("Random:seed = 0");
   pythia.init();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G01PrimaryGeneratorAction::~G01PrimaryGeneratorAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G01PrimaryGeneratorAction::RandomDir(G4ThreeVector& direction) const {
     G4double theta = CLHEP::pi * G4UniformRand();
     G4double phi   = CLHEP::twopi * G4UniformRand();
     
     G4double sinTh=  std::sin(theta);
     G4double xdir = sinTh * std::cos(phi);
     G4double ydir = sinTh * std::sin(phi);
     G4double zdir = std::cos(theta);

     direction.set(xdir,ydir,zdir);
}

void G01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   // generate electron and one proton
   // default particle kinematic
   G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
   G4PrimaryParticle *particle;
   G4PrimaryVertex *vertex;
   G4ThreeVector position;
   G4ThreeVector direction;
#define ONE
#if !defined(ONE)
   static G4String nameElectron("e-");
   static G4String nameProton("proton");
   static G4ParticleDefinition *electron=particleTable->FindParticle(nameElectron);
   static G4ParticleDefinition *proton=particleTable->FindParticle(nameProton);   
   position.set(0,0,0);
   
   vertex = new G4PrimaryVertex(position,0);

   for(int i=0; i<4; ++i) {
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
   }
   
   anEvent->AddPrimaryVertex(vertex); 
#endif
#if defined(ONE)
   // Generator. Process selection. LHC initialization. Histogram.

   position.set(0,0,0);
   vertex = new G4PrimaryVertex(position,0);

   int nprim = 0;
   while (!pythia.next());
   for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal()) {
	 G4ParticleDefinition *g4p = 
	    particleTable->FindParticle(pythia.event[i].id());
	 if(!g4p) {
	    cout << "EEEEK" << endl;
	    exit(1);
	 }
	 if(pythia.event[i].id()!=g4p->GetPDGEncoding()) {
	    cout << "EEEEK2!!!!!!" << endl;
	    exit(1);
	 }
	 particle = new G4PrimaryParticle(g4p);
	 particle->Set4Momentum(pythia.event[i].px()*1000,
				pythia.event[i].py()*1000,
				pythia.event[i].pz()*1000,
				pythia.event[i].e()*1000);
	 vertex->SetPrimary(particle);
#if VERBOSE
	 G4cout << g4p->GetParticleName() << " mass "
		<< pythia.event[i].m() << " energy " 
		<< pythia.event[i].e() << G4endl;
#endif
	 ++nprim;
      }
   anEvent->AddPrimaryVertex(vertex); 
   pythia.stat();
   VTfileio::I()->SetPrimaries(nprim);
   printf("Generating event with %d primaries\n",nprim);
#endif
}
