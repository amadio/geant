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
// $Id: SampleInteractions.cc  J. Apostolakis  $
//
// -------------------------------------------------------------------
//      Geant4 Sampler file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      File name:     SampleInteractions
//
//      Author:        J. Apostolakis
//      Creation date: 20 June 2013
//
//      Adapted from:  test35.cc (HARP) by V.Ivanchenko
//
//      Modifications:
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4ElementVector.hh"

#include "G4VDiscreteProcess.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4ProcessManager.hh"
#include "G4FastStep.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleChangeForDecay.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4ParticleChangeForMSC.hh"

#include "G4ParticleTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Nucleus.hh"
#include "G4IonTable.hh"
#include "G4HadronicProcess.hh"

#include "G4UnitsTable.hh"
#include "G4StateManager.hh"

#include "G4NistManager.hh"

#include "G4ForceCondition.hh"
#include "G4TouchableHistory.hh"

#include "G4VCrossSectionDataSet.hh"

#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "G4Timer.hh"

#include "SampleInteractions.hh"

#include <TFinState.h>

using namespace std;

G4double GetNuclearMass( G4int, G4int, G4int ); // G4Material* material );
const G4MaterialCutsCouple* FindMaterialCutsCouple( G4Material* mat );

int SampleInteractions(
                          G4Material* material,
			  const G4ThreeVector *pos,
                          const G4ParticleDefinition* part,
                          G4VProcess* process,
                          G4double energy,
                          G4double sigmae,
                          G4double theStep,
                          G4int    nevt,
                          G4int    verbose
   )
{
   G4int ret= -1;
  
   G4VDiscreteProcess* discProc= dynamic_cast<G4VDiscreteProcess*>(process);
   G4VContinuousDiscreteProcess* contdProc= dynamic_cast<G4VContinuousDiscreteProcess*>(process);

   if( discProc || contdProc ) {
      ret= SampleDiscreteInteractions( material, pos, part, process, energy, sigmae, theStep, nevt, verbose);
   } else {
      G4cout << " Process " << process->GetProcessName() << " for particle " << part->GetParticleName() << " has no discrete part. " << G4endl;
   }
   return ret;
}

int SampleDiscreteInteractions(
                          G4Material* material, 
			  const G4ThreeVector *pos,
                          const G4ParticleDefinition* part,
                          G4VProcess* proc,
                          G4double energy,
                          G4double sigmae,
                          G4double theStep,
                          G4int    nevt,
                          G4int    verbose
                               )
{
    std::cout << "# Start of Sample Interactions for "
     << " Material= " << material->GetName()
     << " Process = " << proc->GetProcessName()
     << "  E= " << energy
     << "  sigmaE= " << sigmae
     << "     #####" << endl;

    // -------- Projectile
    const G4ThreeVector &aPosition = *pos;
    G4ThreeVector aDirection      = G4ThreeVector(0.0,0.0,1.0);

    G4double mass = part->GetPDGMass();
    //    if(m_pmax == 0.0) { m_pmax = emax; }
    //    else              { emax   = m_pmax; }
    // energy = sqrt(m_p*m_p + mass*mass) - mass; 
    G4DynamicParticle* dParticle= new G4DynamicParticle(part,aDirection,energy);

    // ------- Define target A
    const G4Element* elm = material->GetElement(0); 
    G4int A = (G4int)(elm->GetN()+0.5);
    G4int Z = (G4int)(elm->GetZ()+0.5);

    // ------- Binning 
    cout << "The particle:  " << part->GetParticleName() << endl;
    if(verbose > 0) {
      cout << "The material:  " << material->GetName() 
	   << "  Z= " << Z << "  A= " << A << endl;
      cout << "The step:      " << theStep/mm << " mm" << endl;
      cout << "The position:  " << aPosition/mm << " mm" << endl;
      cout << "The direction: " << aDirection << endl;
      // cout << "The time:      " << aTime/ns << " ns" << endl;
    }
  
    G4double amass = GetNuclearMass( Z, A, verbose ); 
    cout << "Ekin = " << energy/GeV << " GeV" << endl;

  // G4VCrossSectionDataSet* GetCrossSectionDS(const G4ParticleDefinition*, G4Material *);

  // G4VCrossSectionDataSet* cs = GetCrossSectionDS(part, material);

  // G4double cross_sec = 0.0;

    // -------- Track
    G4double aTime(0.0);
    G4Track* gTrack;
    gTrack = new G4Track(dParticle,aTime,aPosition);
    G4TouchableHandle touchable=new G4TouchableHistory();
    G4TransportationManager::GetTransportationManager()->
       GetNavigatorForTracking()->LocateGlobalPointAndUpdateTouchableHandle(aPosition,
								       aDirection,
								       touchable,
								       false);

    gTrack->SetTouchableHandle(touchable);

    // -------- Step

    G4Step* step;
    step = new G4Step();
    step->SetTrack(gTrack);
    gTrack->SetStep(step);
    
    G4StepPoint *aPoint, *bPoint;
    aPoint = new G4StepPoint();
    aPoint->SetPosition(aPosition);
    aPoint->SetMaterial(material);
    G4double safety = 10000.*cm;
    aPoint->SetSafety(safety);
    aPoint->SetTouchableHandle(touchable);
    step->SetPreStepPoint(aPoint);

    const G4MaterialCutsCouple* mcc=
        FindMaterialCutsCouple( material );
  
    aPoint->SetMaterialCutsCouple( mcc );
    aPoint->SetMaterial( material );
  
    bPoint = new G4StepPoint(*aPoint);
    G4ThreeVector bPosition = aPosition + aDirection*theStep;
    bPoint->SetPosition(bPosition);
    step->SetPostStepPoint(bPoint);
    step->SetStepLength(theStep);

    if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
      cout << "G4StateManager PROBLEM: Not able to set it to Idle state! " << endl;

    G4Timer* timer = new G4Timer();
    timer->Start();

    G4double e;
    G4VParticleChange* aChange = 0;

    for (G4int iter=0; iter<nevt; ++iter) {

      G4int modu = 10000;
      if(verbose>1 || iter == modu*(iter/modu)) {
        G4cout << "### " << iter << "-th event start " << G4endl;
      }
      // if(saverand) { CLHEP::HepRandom::saveEngineStatus("random.txt"); }
     
      // Decide  energy  of projectile 
      //         ******
      G4double e0 = energy;
      if(sigmae > 0.0) {
        do {e0 = G4RandGauss::shoot(energy,sigmae);} while (e0 < 0.0);
      }
      dParticle->SetKineticEnergy(e0);
      gTrack->SetStep(step);
      gTrack->SetKineticEnergy(e0);
      G4cout << "Projectile " << gTrack->GetParticleDefinition()->GetParticleName() << G4endl;

      G4LorentzVector labv;
      labv = G4LorentzVector(0., 0., sqrt(e0*(e0 + 2.0*mass)), 
			     e0 + mass + amass);

      aChange = proc->PostStepDoIt(*gTrack,*step); 
      // ** Cause Discrete Interaction **

      G4int n = aChange->GetNumberOfSecondaries();

      G4int nbar = 0;
      G4int n_pr = 0;
      G4int n_nt = 0;
      G4int n_pi = 0;
      const G4DynamicParticle* sec = 0;
      G4ParticleDefinition* pd;
      G4int j;
      
      for(j=0; j<n; ++j) {
        sec = aChange->GetSecondary(j)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        G4int enc= pd->GetPDGEncoding();
        if(pd->GetPDGMass() > 100.*MeV) { ++nbar; }
        if(enc==2212) { ++n_pr; }
        if(enc==2112) { ++n_nt; }
        if(std::fabs(enc)==211 || enc==210) { ++n_pi; }
      }

      printf(" Interaction %5d:  Created %2d secondaries.  %2d hadrons (%2d protons, %2d neutrons), %2d pions\n", iter, n, nbar, n_pr, n_nt, n_pi );
      
      G4ThreeVector  mom;
      G4LorentzVector fm;
      //  Examine the secondaries
      
      for(G4int i=0; i<n; ++i) {
        G4double p, mass1, px, py, pt, theta; // , x;

        sec = aChange->GetSecondary(i)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        mom = sec->GetMomentumDirection();
        e   = sec->GetKineticEnergy();
	if (e < 0.0) { e = 0.0; }

        theta = mom.theta();
        // G4double cost  = cos(theta);

        mass1 = pd->GetPDGMass();
	p = sqrt(e*(e + 2.0*mass1));
	mom *= p;
	fm = G4LorentzVector(mom, e + mass1);
	labv -= fm;
        px = mom.x();
        py = mom.y();
	//        pz = mom.z();
        pt = sqrt(px*px +py*py);

        printf(" Sec[%2d]= %10s (%10d)", i, (const char*) pd->GetParticleName(), pd->GetPDGEncoding()); // , mom.z(), pt, theta );
        printf(" Z= %3d B= %3d ", pd->GetAtomicNumber(), pd->GetBaryonNumber() );
        printf(" pz=%12.4g, pt=%12.4g, theta=%12.6g\n", mom.z(), pt, theta );

        delete aChange->GetSecondary(i);
      }
      aChange->Clear();
      
      G4cout << " E/p balance = " << labv << G4endl;
    }
  
    timer->Stop();
    G4cout << "  "  << *timer << G4endl;

    // Clean up
    delete timer;

    // A step owns its Step points 
    //   - must ensure that they are valid or null (and not pointint to same StepPt)
    delete step; 

    // A track owns its dynamic particle - that will be deleted by it
    delete gTrack;

    return true;
}

