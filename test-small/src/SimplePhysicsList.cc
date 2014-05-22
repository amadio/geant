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

#include "SimplePhysicsList.hh"

#include "G4ProcessManager.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SimplePhysicsList::SimplePhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SimplePhysicsList::~SimplePhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SimplePhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle(); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SimplePhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructTotal();
  // ConstructDecay();  // Should be embedded in the 'Total' Process - tbc

  //  HadronPhysicsFTFP_BERT_WP();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4PhysicsListHelper.hh"
// #include "G4BestUnit.hh"

#include "TotalPhysicsProcess.hh"

#include "TabulatedProcess.hh"
#include "VectorizedProcess.hh"
#include "TabulatedDataManager.hh"
#include "TPartIndex.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SimplePhysicsList::ConstructTotal()
{
  // G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  
  G4cout << " SimplePhysicsList: constructing one TotalPhysicsProcess per particle " << G4endl;
  
  // Must use the old functionality for adding processes
  //  - the new one works only for recognised processes, e.g. Brem, compton ..
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    // G4String particleName = particle->GetParticleName();
    
    G4VRestContinuousDiscreteProcess* totalPhysics=
       new TotalPhysicsProcess(particle->GetParticleName());

    // ph->RegisterProcess( totalPhysics, particle);
    pmanager->AddContinuousProcess( totalPhysics );
    pmanager->AddDiscreteProcess( totalPhysics );
    pmanager->AddRestProcess( totalPhysics );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"

void SimplePhysicsList::ConstructDecay()
{
  // Add Decay Process
  G4cout << "Constructing separate Decay process(es) for each particle." << G4endl;
  
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    
    if (theDecayProcess->IsApplicable(*particle))
    {
      pmanager ->AddProcess(theDecayProcess);
      
      // set ordering for PostStepDoIt and AtRestDoIt
      
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SimplePhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "SimplePhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue // / CLHEP::Unit::mm
    << " mm" << G4endl;
    // G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  //
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
  SetCutValue(defaultCutValue, "proton");

  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SimplePhysicsList::HadronPhysicsFTFP_BERT_WP()
{
}
