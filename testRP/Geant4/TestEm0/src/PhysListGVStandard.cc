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
// $Id: PhysListGVStandard.cc
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListGVStandard.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicsListHelper.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4EmParameters.hh"

#include "G4BuilderType.hh"
#include "G4LossTableManager.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListGVStandard::PhysListGVStandard(const G4String& name) : G4VPhysicsConstructor(name)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(1);
  // inactivate energy loss fluctuations
  param->SetLossFluctuations(false);
  // set min/max energy for tables: 100 eV - 100 TeV by default
  //param->SetMinEnergy(100*eV);
  //param->SetMaxEnergy(100*TeV);
  // set lowest kinetic energy i.e. tracking cut for charged particles having energy loss process: 1 keV by default
  //param->SetLowestElectronEnergy(1*keV);
  // activate/inactivate integral approach: true by default
  // param->SetIntegral(true);
  // inactivate to use cuts as final range
  param->SetUseCutAsFinalRange(false);

  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListGVStandard::~PhysListGVStandard()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListGVStandard::ConstructProcess()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // Add standard EM Processes
  //
  auto aParticleIterator = GetParticleIterator();
  aParticleIterator->reset();
  while( (*aParticleIterator)() ){
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
//      ph->RegisterProcess(new G4PhotoElectricEffect, particle);
        
      G4double LivermoreHighEnergyLimit = 100*TeV;
      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      //LIVERMORE PE
      G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel = new G4LivermorePhotoElectricModel();
        
      theLivermorePhotoElectricModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);
      ph->RegisterProcess(thePhotoElectricEffect, particle);
      
      ph->RegisterProcess(new G4ComptonScattering, particle);
      ph->RegisterProcess(new G4GammaConversion, particle);
    } else if (particleName == "e-") {
//      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      //
      G4eIonisation* eIoni = new G4eIonisation();
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
    } else if (particleName == "e+") {
//      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      //
      G4eIonisation* eIoni = new G4eIonisation();
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      //
//      ph->RegisterProcess(new G4eplusAnnihilation(), particle);
    }
  }

  // Deexcitation
  //
//  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
//  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
