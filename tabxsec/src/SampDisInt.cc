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
#include "G4VEnergyLossProcess.hh"
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

#include "SampleDisInt.hh"

#include <TFinState.h>

using namespace std;

G4double GetNuclearMass( G4int, G4int, G4int ); // G4Material* material );
const G4MaterialCutsCouple* FindMaterialCutsCouple( G4Material* mat );

int SampDisInt(
               G4Material* material,
               G4ThreeVector *pos,
               G4DynamicParticle *dpart,
               G4VProcess* proc,
               G4int    nevt,
               G4int    verbose,
               TFinState& fs)
{
  
  G4VDiscreteProcess* discProc= dynamic_cast<G4VDiscreteProcess*>(proc);
  G4VContinuousDiscreteProcess* contdProc= dynamic_cast<G4VContinuousDiscreteProcess*>(proc);
  
  if( !discProc && !contdProc ) {
    G4cout << " Process " << proc->GetProcessName()
    << " for particle " << dpart->GetParticleDefinition()->GetParticleName()
    << " has no discrete part. " << G4endl;
    return -1;
  }
  
  // -------- Simple structure to hold one final states
  struct Finstat_t {
    G4int *pid;
    G4float (*mom)[3];
  };
  
  Finstat_t *fstat = new Finstat_t[nevt];
  G4int *npart = new G4int[nevt];
  G4float *kerma = new G4float[nevt];
  G4float *weight = new G4float[nevt];
  
  // -------- We chose a short step because we do not want anything
  // -------- to happen during the step
  
  const G4double theStep=1*mm;
  
  // -------- Projectile
  const G4ThreeVector &aPosition = *pos;
  G4ThreeVector aDirection      = dpart->GetMomentumDirection();
  const G4DynamicParticle *dpsave = new G4DynamicParticle(*dpart);
  
  G4double mass = dpart->GetParticleDefinition()->GetPDGMass();
  
  // ------- Define target A
  const G4Element* elm = material->GetElement(0);
  G4int A = (G4int)(elm->GetN()+0.5);
  G4int Z = (G4int)(elm->GetZ()+0.5);
  G4double amass = GetNuclearMass( Z, A, verbose );
  
  // ------- Printout
  if(verbose>1) {
    G4cout << G4endl << "Process   : " << proc->GetProcessName() << G4endl
    << "Particle  : " << dpart->GetParticleDefinition()->GetParticleName()
    << "; p " << dpart->Get4Momentum() << G4endl
    << "Material  : " << material->GetName()
    << "  Z " << Z
    << "  A " << A << " mass " << amass << G4endl
    << "Step      : " << theStep/mm << "mm Position " << aPosition/mm << "mm" << G4endl;
  }
  
  G4double e0 = dpart->GetKineticEnergy();
  
  // G4VCrossSectionDataSet* GetCrossSectionDS(const G4ParticleDefinition*, G4Material *);
  // G4VCrossSectionDataSet* cs = GetCrossSectionDS(part, material);
  
  // -------- Track
  G4Track* gTrack;
  gTrack = new G4Track(new G4DynamicParticle(*dpart),0,G4ThreeVector(aPosition));
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
  
  G4VEnergyLossProcess *elp;
  if((elp = dynamic_cast<G4VEnergyLossProcess*>(proc))) {
    G4ForceCondition fcond;
    elp->PostStepGetPhysicalInteractionLength(*gTrack,theStep,&fcond);
  }
  
  if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
    cout << "G4StateManager PROBLEM: Not able to set it to Idle state! " << G4endl;
  
  G4double e;
  G4VParticleChange* aChange = 0;
  
  
  //----------------------------- Start event loop -------------------------------
  for (G4int iter=0; iter<nevt; ++iter) {
    
    G4bool needendl = FALSE;
    G4int modu = 10000;
    if(verbose>1 || iter == modu*(iter/modu)) {
      G4cout << "### " << iter << "-th event: "
      << dpart->GetParticleDefinition()->GetParticleName() << " "
      << proc->GetProcessName()
      << " @" << dpart->Get4Momentum() << " on "
      << material->GetName() << " (" << Z <<"," << A << "," << amass << ")";
      needendl=TRUE;
    }
    
    gTrack->SetStep(step);
    gTrack->SetKineticEnergy(e0);
    
    G4LorentzVector labv;
    labv = G4LorentzVector(gTrack->GetMomentum(), gTrack->GetTotalEnergy());
    if(labv != dpsave->Get4Momentum()) {
      if(needendl) {
        G4cout << G4endl;
        needendl=FALSE;
      }
      G4cout << "The momentum of the projectile has changed from "
      << labv << " to "
      << dpart->Get4Momentum() << G4endl;
      exit(1);
    }
    G4LorentzVector pcons(labv);
    G4LorentzVector opart(0.,0.,0.,0.);
    
    if(dynamic_cast<G4HadronInelasticProcess*>(proc))
      pcons[3]+=amass;
    if(dpart->GetParticleDefinition()->GetParticleName() == G4String("e+")
       && proc->GetProcessName() == G4String("annihil"))
      pcons[3]+=mass;
    // ----------------------------------- ** Cause Discrete Interaction ** ----------------------------
    
    aChange = proc->PostStepDoIt(*gTrack,*step);
    // -- See wheter the original particle is still alive
    if((aChange->GetTrackStatus() == fAlive) || (aChange->GetTrackStatus() == fStopButAlive)) {
      // -- Original particle is alive -- let's get it
      printf("after %s on %s %s is still alive\n",
             (const char *)proc->GetProcessName(),
             (const char *)material->GetName(),
             (const char *)dpart->GetParticleDefinition()->GetParticleName());
      G4Step *tstep = aChange->UpdateStepForPostStep(step);
      opart.set(tstep->GetTrack()->GetMomentum(),tstep->GetTrack()->GetTotalEnergy());
      if((G4String("hIoni") == proc->GetProcessName()) || (G4String("ionIoni") == proc->GetProcessName())
         || (G4String("eIoni") == proc->GetProcessName())) {
        // the electron is NOT produced... ;-)
        opart[3]-=G4ParticleTable::GetParticleTable()->FindParticle("e-")->GetPDGMass();
      }
    }
    
    // ----------------------------------- Trying to find out the angle in elastic ----------------------------
    if(proc->GetProcessName() == G4String("hadElastic")) {
      G4double cost = labv.vect().cosTheta(*((G4ParticleChange*) aChange)->GetMomentumDirection());
      if(needendl) {
        G4cout << G4endl;
        needendl=FALSE;
      }
      G4cout << "Elastic scattering by cosTheta "
      << cost << " (" << 180*std::acos(cost)/std::acos(-1) << " deg)" << G4endl;
      G4double ken = ((G4ParticleChange*) aChange)->GetEnergy();
      G4double mas = ((G4ParticleChange*) aChange)->GetMass();
      G4double pmo = std::sqrt(ken*(2*ken+mas));
      opart.set(*((G4ParticleChange*) aChange)->GetMomentumDirection()*pmo,ken+mas);
      G4cout << " Output mom in elastic scatt " << opart << " mass " << mas << G4endl;
      opart[3]-=amass;
    }
    pcons -= opart;
    
    // ----------------------------------- Generated secondaries ----------------------------
    G4int n = aChange->GetNumberOfSecondaries();
    
    npart[iter] = 0;
    kerma[iter] = 0;
    weight[iter] = 1; // for the moment all events have the same prob
    if(n) {
      fstat[iter].mom = new G4float[n][3];
      fstat[iter].pid = new G4int[n];
    } else {
      fstat[iter].mom = 0;
      fstat[iter].pid = 0;
    }
    
    G4int nbar = 0;
    G4int n_pr = 0;
    G4int n_nt = 0;
    G4int n_pi = 0;
    G4int n_ga = 0;
    G4int n_el = 0;
    G4int n_po = 0;
    const G4DynamicParticle* sec = 0;
    G4ParticleDefinition* pd;
    G4int j;
    
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
      G4int enc= pd->GetPDGEncoding();
      if(enc==-2212) ++nbar;
      if(enc==2212) ++n_pr;
      if(enc==2112) ++n_nt;
      if(enc==22) ++n_ga;
      if(enc==11) ++n_el;
      if(enc==-11) ++n_po;
      if(std::fabs(enc)==211 || enc==210) { ++n_pi; }
      
      theta = mom.theta();
      
      mass1 = pd->GetPDGMass();
      p = sqrt(e*(e + 2.0*mass1));
      mom *= p;
      fm = sec->Get4Momentum();
      pcons -= fm;
      px = mom.x();
      py = mom.y();
      //        pz = mom.z();
      pt = sqrt(px*px +py*py);
      
      G4int g5pid = TPartIndex::I()->PartIndex(enc);
      if(g5pid<0) {
        kerma[iter]+=e+mass1;
      } else {
        fstat[iter].pid[npart[iter]] = g5pid;
        fstat[iter].mom[npart[iter]][0] = mom.x();
        fstat[iter].mom[npart[iter]][1] = mom.x();
        fstat[iter].mom[npart[iter]][2] = mom.x();
        ++npart[iter];
      }
      if(verbose > 1) {
        if(needendl) {
          G4cout << G4endl;
          needendl=FALSE;
        }
        G4cout << " Sec[" << i << "]="
        << pd->GetParticleName() << " (" << pd->GetPDGEncoding() << ") "
        << " Z= " << pd->GetAtomicNumber() << " B= " << pd->GetBaryonNumber()
        << " p= " << sec->Get4Momentum()
        << G4endl;
      }
    }
    if(verbose>1 || iter == modu*(iter/modu)) {
      G4cout << ": " << n << " sec (" << n_pr << " protons, "
      << n_nt << " neutrons, " << n_pi << " pi), "
      << n_el << " e-, " << n_po << " e+, " << n_ga << " g"
      << G4endl;
    }
    const G4double prec=1e-5;
    G4double err=0;
    if(n) {
      G4double ptot = labv.vect().mag();
      G4double ermax = std::abs(pcons[3]);
      if(labv[3]) ermax/=labv[3];
      for(G4int i=0; i<3; ++i) {
        if(ptot)
          err = std::abs(pcons[3])/ptot;
        if(err>ermax) ermax=err;
      }
      if(ermax>prec) {
        G4cout << setfill('-') << setw(120) << "-" << setfill(' ') << setw(0) << G4endl
        <<"Dubious E/p balance " << setiosflags(ios::scientific) << setprecision(2) << ermax
        << " p=" << setiosflags(ios::fixed) << setprecision(6) << pcons << " for "
        << proc->GetProcessName() << " of "
        << dpart->GetParticleDefinition()->GetParticleName() << " @ "
        << dpart->Get4Momentum() << " on "
        << material->GetName() << " (" << Z <<"," << A << "," << amass << ")" << G4endl
        << "Particles generated:" << G4endl;
        G4FastStep *uChange = dynamic_cast<G4FastStep*>(aChange);
        G4ParticleChange *vChange = dynamic_cast<G4ParticleChange*>(aChange);
        G4ParticleChangeForDecay *wChange = dynamic_cast<G4ParticleChangeForDecay*>(aChange);
        G4ParticleChangeForGamma *xChange = dynamic_cast<G4ParticleChangeForGamma*>(aChange);
        G4ParticleChangeForLoss *yChange = dynamic_cast<G4ParticleChangeForLoss*>(aChange);
        G4ParticleChangeForMSC *zChange = dynamic_cast<G4ParticleChangeForMSC*>(aChange);
        
        /*	    G4cout << " G4FastStep " << uChange
         << " GParticleChange " << vChange
         << " GparticleChangeForDecay " << wChange
         << " G4ParticleChangeForGamma " << xChange
         << " G4ParticleChangeForLoss " << yChange
         << " G4ParticleChangeForMSC " << zChange
         << G4endl;*/
        if(yChange) {
          /*	       const G4Track *trnew = yChange->GetCurrentTrack();
           G4cout << "Momentum " << trnew->GetMomentum()
           << " kinEnergy " << trnew->GetKineticEnergy()
           << " mass " << trnew->GetParticleDefinition()->GetPDGMass()
           << G4endl;*/
        }
        for(G4int i=0; i<n; ++i) {
          
          sec = aChange->GetSecondary(i)->GetDynamicParticle();
          pd  = sec->GetDefinition();
          G4cout << "  #" << setw(3) << i << setw(0) << ": "
		      << setiosflags(ios::left) << setw(10) << pd->GetParticleName()
		      << " (" << setiosflags(ios::right) << setw(6) << pd->GetPDGEncoding() << ") " << setw(0)
		      << " Z:" << pd->GetAtomicNumber() << " B:" << pd->GetBaryonNumber()
		      << " p:" << sec->Get4Momentum()
		      << G4endl;
        }
        G4cout << setfill('-') << setw(120) << "-" << setfill(' ') << setw(0) << G4endl;
      }
    }
    for(G4int i=0; i<n; ++i)
      delete aChange->GetSecondary(i);
    aChange->Clear();
    
  }
  //----------------------------- End event loop -------------------------------
  
  //----------------------------- Build final structure ------------------------
  G4int ntotp = 0;
  for(G4int i=0; i<nevt; ++i) ntotp+=npart[i];
  G4int *tpid = new G4int[ntotp];
  G4float (*tmom)[3] = new G4float[ntotp][3];
  
  G4int ip=0;
  for(G4int i=0; i<nevt; ++i) {
    for(G4int j=0; j<npart[i]; ++j ) {
      tmom[ip+j][0]=fstat[i].mom[j][0];
      tmom[ip+j][1]=fstat[i].mom[j][1];
      tmom[ip+j][2]=fstat[i].mom[j][2];
      tpid[ip+j]=fstat[i].pid[j];
    }
    ip+=npart[i];
  }
  
  fs.SetFinState(nevt, weight, kerma, npart, tmom, tpid);
  
  // A step owns its Step points
  //   - must ensure that they are valid or null (and not pointint to same StepPt)
  delete step;
  
  // A track owns its dynamic particle - that will be deleted by it
  delete gTrack;
  
  delete [] fstat;
  delete [] npart;
  delete [] kerma;
  delete [] weight;
  
  delete [] tmom;
  delete [] tpid;
  
  return true;
}

G4double GetNuclearMass( G4int Z, G4int N, G4int verbose )
{
  G4double mass= 0.0;
  G4Nucleus targetNucleus;
  
  if(verbose>1) G4cout << "Nucleus with N= " << N << "  Z= " << Z << G4endl;
  targetNucleus.SetParameters((G4double)N, (G4double)Z);
  mass = targetNucleus.AtomicMass((G4double)N, (G4double)Z);
  if(verbose>1) G4cout << "Mass from targetNucleus(MeV)= " << mass/MeV << G4endl;
  mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, N);
  if(verbose>1) G4cout << "Mass from IonTable(MeV)=      " << mass/MeV << G4endl;
  return mass;
}

#include "G4ProductionCutsTable.hh"
// G4MaterialCutsCouple

const G4MaterialCutsCouple* FindMaterialCutsCouple( G4Material* mat )
{
  const G4MaterialCutsCouple* couple = 0;
  
  G4ProductionCutsTable* theCoupleTable =
  G4ProductionCutsTable::GetProductionCutsTable();
  
  size_t numOfCouples = theCoupleTable->GetTableSize();
  
  for (size_t i=0; i<numOfCouples; i++) {
    couple = theCoupleTable->GetMaterialCutsCouple(i);
    if (couple->GetMaterial() == mat) break;
  }
  if(couple == 0) {
    G4cerr << "Fatal ERROR> Not found couple for material " << mat->GetName() << G4endl;
    exit(1);
  }
  return couple;
}
