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

#include "WorkloadManager.h"   
#include "ComptonCrossSection.h"    // For ComptonCrossSection - not separate
#include "TAtomicShells.h" 

#include "physical_constants.h"  

#include <iostream>

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
   Double_t density = 0.0;  // No material = zero density

   // std::cout << " GammaCompton::ComputeIntLen called " 
   //           << " with " << ntracks << " tracks. " << std::endl;

   TGeoMaterial *mat = vol->GetMaterial();
   // std::cout << " GC: Material = " << mat->GetName() << std::endl;
   if (mat) density = mat->GetDensity();

   // Write in the thread space for the current basket  ... ???
   Double_t *rndArray = &gPropagator->fDblArray[2*tid*ntracks];
   Int_t irnd = 0;
   TRandom *rngEngine= gPropagator->fRndm[tid]; 
   rngEngine->RndmArray(ntracks, rndArray);

   GeantTrack** gTracks= gPropagator->fTracks; 
   static const Int_t PdgGamma= 22; 

   Int_t noGammas= 0; 
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
         noGammas++; 
      }
      if( lengthC < lengths[itrack] )
      {
         lengths[itrack] = lengthC; 
         // limited[itrack] = kComptonProcess;  // Flag that Compton limited this one
      }
   }
   if( noGammas > 0 ) { 
      std::cout << " GammaCompton::ComputeIntLen:  "; 
      std::cout << " GC: Material = " << mat->GetName() << std::endl;
      std::cout << " Found " << noGammas << "  gammas. " << std::endl;
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
   TRandom      *RngEngine= gPropagator->fRndm[tid]; 
   TGeoMaterial *tMaterial = vol->GetMaterial();

   RngEngine->RndmArray(2*ntracks, rndArray);
   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];
      track  = gPropagator->fTracks[itrack];
      if( 1 ) // should be ( limited[i] == kComptonProcess ) 
      {
         //  Only if Compton limited the step does it create a secondary
         TLorentzVector  gamma4m(  track->px, track->py, track->pz, track->e);
         TLorentzVector  electron4m( 0.0, 0.0, 0.0, 0.0);            //  
         // Int_t           electronOut;  CHANGED -> Always creates 'out' electron
         // Double_t        enDeposit; 

         SampleSecondaries( *tMaterial, gamma4m, electron4m, RngEngine ); 
         //****************

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
            
            GeantTrack* eTrk= gPropagator->AddTrack(track->evslot);
            eTrk->evslot = track->evslot;
            *eTrk->path = *track->path;
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
               Int_t itracknew = eTrk->particle;
               trackout[nout] = itracknew;
               nout++;
               std::cout << " WANT to add gamma of Energy= " << eTrk->e << " to stack. " << std::endl;

               GeantVolumeBasket *basket = gPropagator->fWMgr->GetCurrentBasket(tid);
               if (basket) gPropagator->fWMgr->AddPendingTrack(itracknew, basket, tid);
            }
         }
      }
   }   
   StepManager(3, ntracks, trackin, nout, trackout);
}




// Adapted from G4KleinNishinaCompton::SampleSecondaries 


void GammaCompton::SampleSecondaries(
			const TGeoMaterial&  tMaterial,
     		        TLorentzVector&      gamma4Mom,     // In/Out: 4-mom of gamma
                        TLorentzVector&      electron4Mom,  // Out: 4-mom of outgoing e-
                        // Double_t&            enDeposit,     // Out: Energy Deposit - None
                        TRandom*             RngEngine)
{
  // The scattered gamma energy is sampled according to Klein - Nishina formula.
  // The random number techniques of Butcher & Messel are used 
  // (Nuc Phys 20(1960),15).
  // Note : Effects due to binding of atomic electrons are negliged.
 
  Double_t gamEnergy0 =  gamma4Mom.E(); 

  // extra protection
  // if(gamEnergy0 < lowestGammaEnergy) {
  //  fParticleChange->ProposeTrackStatus(fStopAndKill);
  //  fParticleChange->ProposeLocalEnergyDeposit(gamEnergy0);
  //  fParticleChange->SetProposedKineticEnergy(0.0);
  //  return;
  // }

  Double_t E0_m = gamEnergy0 / electron_mass_c2 ;

  TVector3 gamDirection0 = gamma4Mom.P(); // GetMomentumDirection();

  //
  // sample the energy rate of the scattered gamma 
  //

  Double_t epsilon, epsilonsq, onecost, sint2, greject ;

  Double_t epsilon0   = 1./(1. + 2.*E0_m);
  Double_t epsilon0sq = epsilon0*epsilon0;
  Double_t alpha1     = - log(epsilon0);
  Double_t alpha2     = 0.5*(1.- epsilon0sq);

  const  Int_t RandsPerIter=3; 
  Double_t     dUniformRand[RandsPerIter]; 
  do {
     RngEngine->RndmArray( RandsPerIter, dUniformRand );
     if ( alpha1/(alpha1+alpha2) > dUniformRand[0] ) {
        epsilon   = exp(-alpha1*dUniformRand[1]);   // epsilon0**r
        epsilonsq = epsilon*epsilon; 
     } else {
        epsilonsq = epsilon0sq + (1.- epsilon0sq)*dUniformRand[1];
        epsilon   = sqrt(epsilonsq);
     };
     
     onecost = (1.- epsilon)/(epsilon*E0_m);
     sint2   = onecost*(2.-onecost);
     greject = 1. - epsilon*sint2/(1.+ epsilonsq);
     
  } while (greject < dUniformRand[2] );
 
  // Accept/Reject is NOT well suited to Vectorisation

  Double_t  aUniformRand;
  RngEngine->RndmArray( 1, &aUniformRand );
  //
  // scattered gamma angles. ( Z - axis along the parent gamma)
  //

  if(sint2 < 0.0) { sint2 = 0.0; }
  Double_t cosTeta = 1. - onecost; 
  Double_t sinTeta = sqrt (sint2);
  Double_t Phi     = 2.0 * PI * aUniformRand;
  //
  // update G4VParticleChange for the scattered gamma
  //
   
  TVector3 gamDirection1(sinTeta*cos(Phi), sinTeta*sin(Phi), cosTeta);
  gamDirection1.RotateUz(gamDirection0);
  Double_t gamEnergy1 = epsilon*gamEnergy0;   // Final Energy of Gamma

  TVector3 gamMomentum= gamEnergy1 * gamDirection1; 
  gamma4Mom.SetPxPyPzE( gamMomentum.X(), gamMomentum.Y(), gamMomentum.Z(), gamEnergy1 );

  // kinematic of the scattered electron
  //
  Double_t eKinEnergy = gamEnergy0 - gamEnergy1;
  // eKinEnergy = max( 0.0, eKinEnergy ); 

  // Double_t eMomentumMag = sqrt( eKinEnergy * ( electron_restMass_c2+ 2.0* eKinEnergy  ) ); 
  TVector3 eVector      = gamEnergy0*gamDirection0 - gamEnergy1*gamDirection1;
  // In case of problem with rounding or units, can normalise and fix
  //   eVector  = eVector.unit(); 
  //   eVector *= eMomentumMag;     
  
  electron4Mom.SetPxPyPzE( eVector.X(), eVector.Y(), eVector.Z(), eKinEnergy+electron_mass_c2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
