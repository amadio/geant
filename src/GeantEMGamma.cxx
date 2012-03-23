// Toy physics processes for our propagator prototype. Currently including:
// - single scattering as a discrete process
// - energy loss as continuous process
// - generic interaction as discrete process, producing secondaries

#include "GeantEMGamma.h"
#include "GeantVolumeBasket.h"
#include "TMath.h"
// #include "TH1.h"
// #include "TF1.h"
#include "TRandom.h"
#include "globals.h"
#include "GeantTrack.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"             // For TGeoRotation
// #include "TGeoBranchArray.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
// #include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "ComptonCrossSection.h"    // For ComptonCrossSection - not separate
#include "TAtomicShells.h" 

#include <cstdio>

ClassImp(GammaCompton)

GammaCompton::GammaCompton(const char *name) : PhysicsProcess(name), 
   fTcut(0.0),
   fTmax(0.) 
{
   TObject::SetBit(kDiscrete);
   fComptonXS= new ComptonCrossSection(); 
}

//______________________________________________________________________________
void GammaCompton::ComputeIntLen(TGeoVolume *vol, 
                                      Int_t ntracks, 
                                      Int_t *trackin, 
                                      Double_t *lengths, 
                                      Int_t tid)
{
// Generates an interaction length for the Compton process.
//
// trackin and lengths arrays should be be correctly set by the caller
   // const Double_t kC1 = 500.;
   // const Double_t xlen = TMath::Limits<double>::Max();
   Int_t itrack;
   Double_t density = 1.0e-100;  //  finite minimum density

   TGeoMaterial *mat = vol->GetMaterial();
   if (mat) density = mat->GetDensity();

   // Write in the thread space for the current basket  ... ???
   Double_t *rndArray = &gPropagator->fDblArray[2*tid*ntracks];
   Int_t irnd = 0;
   TRandom *rngEngine= gPropagator->fRndm[tid]; 
   rngEngine->RndmArray(ntracks, rndArray);

   GeantTrack** gTracks= gPropagator->fTracks; 
   static const Int_t PdgGamma= 22; 

   //  Could use count of number of gammas which are alive ...
   for (Int_t i=0; i<ntracks; i++) {   
      Double_t sigma= 0.0, meanFreePath= DBL_MAX;
      Double_t lengthC= DBL_MAX;
      itrack = trackin[i];

      GeantTrack* ptrack= gTracks[itrack]; 
      // Check that the particle is a gamma && alive
      if( (ptrack->pdg == PdgGamma) && (ptrack->IsAlive()) ){
         Double_t     kinEnergy= ptrack->e; 
         TGeoMaterial *mat = vol->GetMaterial();

         sigma= fComptonXS->CrossSectionForMaterial( *mat, kinEnergy ); // , fTcut, fTmax); 
         if( sigma > 0.0 ) {
            // Use sigma to sample step size HERE
            meanFreePath= 1.0/sigma; 
            Double_t theNumberOfInteractionLengthsLeft= - std::log( rndArray[irnd++] ); 
            lengthC = theNumberOfInteractionLengthsLeft * meanFreePath; 
         }
      }
      if( lengthC < lengths[itrack] )
      {
         lengths[itrack] = lengthC; 
         // limited[itrack] = kComptonProcess;  // Flag that Compton limited this one
      }
   }
}

//______________________________________________________________________________
void GammaCompton::PostStep(TGeoVolume *vol,
                                 Int_t  ntracks,
                                 Int_t *trackin, 
                                 Int_t &nout, 
                                 Int_t* trackout, 
                                 Int_t  tid)
{
// Do post-step actions on particle after scattering process. Surviving tracks
// copied in trackout
   // Compute the max theta angle opening after scattering.
   // const Double_t ctmax = TMath::Cos(1.*TMath::DegToRad()); // 1 degree
   // Double_t   theta, phi, scale,thetav,phiv; 
   // Double_t   dir[3];
   // Double_t   dirnew[3];
   GeantTrack   *track = 0;
   Int_t        itrack;
   // Double_t p;
   Double_t     *rndArray = &gPropagator->fDblArray[2*tid*ntracks];
   // Int_t        irnd = 0;
   TRandom      *fRndEngine= gPropagator->fRndm[tid]; 
   TGeoMaterial *tMaterial = vol->GetMaterial();

   fRndEngine->RndmArray(2*ntracks, rndArray);
   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];
      track  = gPropagator->fTracks[itrack];
      if( 1 ) // should be ( limited[i] == kComptonProcess ) 
      {
         //  Only if Compton limited the step does it create a secondary
         TLorentzVector  gamma4m(  track->px, track->py, track->pz, track->e);
         TLorentzVector  electron4m( 0.0, 0.0, 0.0, 0.0);            //  
         // Int_t           electronOut;   -> Always creates 'out' electron
         Double_t        enDeposit; 
         SampleSecondaries( *tMaterial, 
                            gamma4m,      // electronOut, 
                            electron4m, 
                            enDeposit ); 

         if( 1 ) //  ( electron4m.E > E_threshold )
         {            
            // static TParticlePDG* electronPDG = TDatabasePDG::Instance()->GetParticle(kElectron); 

            // Change 4-momentum of surviving particle
            track->px = gamma4m.X(); 
            track->py = gamma4m.Y();
            track->pz = gamma4m.Z(); 
               
            // All tracks survive - for now.  TODO: kill gamma below threshold
            if (trackout) trackout[nout] = itrack;
            nout++;

            // Create track for electron
            GeantTrack* eTrk= new GeantTrack();
            eTrk->px = gamma4m.X(); 
            eTrk->py = gamma4m.Y();
            eTrk->pz = gamma4m.Z(); 
            eTrk->e = gamma4m.E(); 
            eTrk->event= track->event; 
            eTrk->species= kLepton;     // electronPDG; 
            eTrk->pdg= 11; 
            eTrk->xpos= track->xpos;
            eTrk->ypos= track->ypos;
            eTrk->zpos= track->zpos;
            eTrk->safety= track->safety; 

            // Always add secondary electron -- TODO: do not add if E < threshold
            if (trackout) 
            {
               Int_t itracknew = gPropagator->AddTrack(eTrk);
               trackout[nout] = itracknew;
               nout++;
               // if (gPropagator->fWMgr->GetCurrentBasket(tid)) 
               //    gPropagator->fWMgr->GetCurrentBasket(tid)->AddPendingTrack(itracknew);
               std::cout << " WANT to add gamma of Energy= " << eTrk->e << " to stack. " << std::endl;
            }
         }
      }
   }   
   StepManager(0, ntracks, trackin, nout, trackout);
}



//  Adapted from  G4KleinNishinaModel::SampleSecondaries()

void GammaCompton::SampleSecondaries(
                     const TGeoMaterial&   tMaterial,
                     TLorentzVector&        gamma4Mom,     // In/Out: 4-mom of gamma
                     // Int_t&                electronOut,   // Out: true if secondary created
                     TLorentzVector&        electron4Mom,  // Out: 4-mom of outgoing e-
                     Double_t&             enDeposit)     // Out: Energy Deposit
{
  // primary gamma
  TVector3       direction= gamma4Mom.Vect(); 
  Double_t       energy=    gamma4Mom.E(); 
  TLorentzVector lv1, lv2; 

  const Int_t  MaxRandom= 12; 
  Double_t     UniformRandom[MaxRandom]; 
  Int_t randCount= 0; 

  // select atom
  const TGeoElement* elm = fComptonXS->SelectRandomAtom(tMaterial, energy);

  // select shell first
  Int_t Z= elm->Z(); 
  Int_t nShells = TAtomicShells::GetNumberOfShells(Z);
  if(nShells > (Int_t) fProbabilities.size() ) { fProbabilities.resize(nShells); }
  Double_t totprob = 0.0;
  Int_t i;

  for(i=0; i<nShells; ++i) {
     Double_t bindingEnergy =  TAtomicShells::GetBindingEnergy(Z, i);
     Double_t eth = sqrt(bindingEnergy*(bindingEnergy + electron_mass_c2)) -
        0.5*(sqrt(bindingEnergy*(bindingEnergy + 2*electron_mass_c2)) - bindingEnergy);
     Double_t prob = 1.0 - eth/energy;
     if(prob > 0.0) { totprob += prob*TAtomicShells::GetNumberOfShellElectrons(i); } 
     fProbabilities[i] = totprob; 
  }
  if(totprob == 0.0) { return; }

  // Loop on sampling
  Double_t eKinEnergy;
  const Int_t nlooplim = 100;
  Int_t nloop = 0;
  Int_t firstGo = true;

  do {
    ++nloop;

    //  Retry 
    const Int_t randsPerLoopIter= 5; 
    if( firstGo || (randCount + randsPerLoopIter >= MaxRand) ) {
       // Refill the Vector with Random numbers
       fRndEngine->RndmArray( MaxRandom, UniformRand );
       randCount=0; 
       firstGo= false; 
    }

    Double_t xprob = totprob*UniformRandom[randCount++];    // Rand #1

    // select shell
    for(i=0; i<nShells; ++i) { if(xprob <= fProbabilities[i]) {break;} }
   
    Double_t bindingEnergy = elm->GetAtomicShell(i);

    // shortcut if the loop is too long
    if(nloop >= nlooplim) {
      lv1.set(0.0,0.0,0.0,0.0);
      eKinEnergy = energy - bindingEnergy;
      if(eKinEnergy < 0.0) { eKinEnergy = 0.0; }
      Double_t eTotMomentum = sqrt(eKinEnergy*(eKinEnergy + 2*electron_mass_c2));
      Double_t phi = UniformRandom[randCount++]*twopi;      // Rand #2
      Double_t costet = 2*UniformRand[randCount++] - 1;     // Rand #3
      Double_t sintet = sqrt((1 - costet)*(1 + costet));
      lv2.set(eTotMomentum*sintet*cos(phi),eTotMomentum*sintet*sin(phi),
	      eTotMomentum*costet,eKinEnergy + electron_mass_c2);
      break;
    }

    Double_t limitEnergy = limitFactor*bindingEnergy;
    Double_t gamEnergy0 = energy;
    lv1.set(0.0,0.0,energy,energy);

    //std::cout << "nShells= " << nShells << " i= " << i 
    //   << " Egamma= " << energy << " Ebind= " << bindingEnergy
    //   << " Elim= " << limitEnergy 
    //   << std::endl;

    // for low energy rest frame of the electron
    if(energy < limitEnergy) { 
      Double_t eTotMomentum = sqrt(bindingEnergy*(bindingEnergy + 2*electron_mass_c2));
      Double_t phi = UniformRand[randCount++]*twopi;              // Rand #4
      Double_t costet = 2*UniformRand[randCount++] - 1;           // Rand #5
      Double_t sintet = sqrt((1 - costet)*(1 + costet));
      lv2.set(eTotMomentum*sintet*cos(phi),eTotMomentum*sintet*sin(phi),
	      eTotMomentum*costet,bindingEnergy + electron_mass_c2);
      bst = lv2.boostVector();
      lv1.boost(-bst);
      gamEnergy0 = lv1.e();
    }

    // In the rest frame of the electron
    // The scattered gamma energy is sampled according to Klein - Nishina formula.
    // The random number techniques of Butcher & Messel are used 
    // (Nuc Phys 20(1960),15).
 
    Double_t E0_m = gamEnergy0/electron_mass_c2;

    //
    // sample the energy rate of the scattered gamma 
    //

    Double_t epsilon, epsilonsq, onecost, sint2, greject ;

    Double_t epsilon0   = 1./(1 + 2*E0_m);
    Double_t epsilon0sq = epsilon0*epsilon0;
    Double_t alpha1     = - log(epsilon0);
    Double_t alpha2     = 0.5*(1 - epsilon0sq);

    do {
       const Int_t randsPerIteration= 3; 
       if( randCount + randsPerIteration >= MaxRand ) {
          // Refill the Vector with Random numbers
          fRndEngine->RndmArray( MaxRandom, UniformRand );
          randCount=0; 
       }

       if ( alpha1/(alpha1+alpha2) > UniformRand[randCount++] ) {   // Rand A1
          epsilon   = exp(-alpha1*UniformRand[randCount++]);        // Rand A2
                            // epsilon0**r
          epsilonsq = epsilon*epsilon; 

       } else {                                                  // Rand A1'
          epsilonsq = epsilon0sq + (1.- epsilon0sq)*UniformRand[randCount++];
          epsilon   = sqrt(epsilonsq);
       };

       onecost = (1.- epsilon)/(epsilon*E0_m);
       sint2   = onecost*(2.-onecost);
       greject = 1. - epsilon*sint2/(1.+ epsilonsq);
       
    } while (greject < UniformRand[randCount++]);            // Rand A3  (max)

    Double_t gamEnergy1 = epsilon*gamEnergy0;
 
    // before scattering total 4-momentum in e- system
    lv2.set(0.0, 0.0, 0.0, electron_mass_c2);
    lv2 += lv1;
 
    //
    // scattered gamma angles. ( Z - axis along the parent gamma)
    //
    if(sint2 < 0.0) { sint2 = 0.0; }
    Double_t cosTeta = 1. - onecost; 
    Double_t sinTeta = sqrt(sint2);
    Double_t Phi  = twopi * G4UniformRand();

    // e- recoil
    //
    // in  rest frame of the electron
    if(energy < limitEnergy) { 
      TVector3 gamDir = lv1.Vect().unit();
      TVector3 v = TVector3(sinTeta*cos(Phi),sinTeta*sin(Phi),cosTeta);
      v.rotateUz(gamDir);
      lv1.set(gamEnergy1*v.x(),gamEnergy1*v.y(),gamEnergy1*v.z(),gamEnergy1);
      lv2 -= lv1;
      //cout << "Egam= " << lv1.e() << "  Ee= " << lv2.e()-electron_mass_c2 << std::endl;
      lv2.boost(bst);
      lv1.boost(bst);
      eKinEnergy = lv2.e() - electron_mass_c2 - 2*bindingEnergy;
      
    } else {
      lv1.set(gamEnergy1*sinTeta*cos(Phi),gamEnergy1*sinTeta*sin(Phi),
	      gamEnergy1*cosTeta,gamEnergy1);
      lv2 -= lv1;
      eKinEnergy = lv2.e() - electron_mass_c2 - bindingEnergy;
    }
   
    // std::cout << "eKinEnergy= " << eKinEnergy << std::endl;

  } while ( eKinEnergy < 0.0 );

  //
  //  Prepare the scattered gamma
  //
  Double_t gamEnergy1 = lv1.e();
  TVector3 gamMomentum = lv1.Vect();
  gamMomentum.rotateUz(direction);
  gamma4Mom->SetPxPyPzE( gamMomentum.X(), gamMomentum.Y(), gamMomentum.Z(), gamEnergy1 );

  //
  // kinematic of the scattered electron
  //
  //  if(eKinEnergy > lowestGammaEnergy) {  .. leave the decision to calling method
  TVector3 eMomentum = lv2.Vect();
  eMomentum.rotateUz(direction);
  // G4DynamicParticle* dp = new G4DynamicParticle(theElectron, eDirection, eKinEnergy);
  electron4Mom->SetPxPyPzE( eMomentum.X(), eMomentum.Y(), eMomentum.Z(), eKinEnergy);

  enDeposit = energy - gamEnergy1 - eKinEnergy;
  
  // No sampling of deexcitation - at this stage   
  //

  // energy balance
  if(edep < 0.0) { edep = 0.0; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

