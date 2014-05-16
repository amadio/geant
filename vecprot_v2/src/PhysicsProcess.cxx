// Toy physics processes for our propagator prototype. Currently including:
// - single scattering as a discrete process
// - energy loss as continuous process
// - generic interaction as discrete process, producing secondaries

#include "PhysicsProcess.h"
#include "GeantThreadData.h"
#include "GeantVApplication.h"
#include "WorkloadManager.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "globals.h"
#include "GeantTrack.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/navigationstate.h"
#else
#include "TGeoBranchArray.h"
#endif
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "TGenPhaseSpace.h"

ClassImp(PhysicsProcess)
ClassImp(ScatteringProcess)

//______________________________________________________________________________
void ScatteringProcess::ComputeIntLen(TGeoMaterial *mat, 
                                      Int_t ntracks, 
                                      GeantTrack_v &tracks,
                                      Double_t *lengths, 
                                      Int_t tid)
{
// Generates an interaction length for the scattering process. Nothing physical,
// just generate something comparable with the size of the current volume.
// 
//
// trackin and lengths arrays should be be correctly set by the caller
   const Double_t kC1 = 500.;
   const Double_t xlen = TMath::Limits<double>::Max();
   Double_t density = 1.e-5;
   if (mat) density = mat->GetDensity();
   density = TMath::Max(density, 1.E-3);
   // Make sure we write in the thread space for the current basket
   Double_t *rndArray = gPropagator->fThreadData[tid]->fDblArray;
   Int_t irnd = 0;
   gPropagator->fThreadData[tid]->fRndm->RndmArray(ntracks, rndArray);
   for (Int_t i=0; i<ntracks; i++) {
      lengths[i] =  0.5*xlen;
      if (tracks.fStatusV[i] != kKilled)
         lengths[i] = kC1*tracks.fEV[i]*rndArray[irnd++]/density;
   }
}

//______________________________________________________________________________
void ScatteringProcess::PostStep(TGeoMaterial */*mat*/,
                                 Int_t ntracks,
                                 GeantTrack_v &tracks, 
                                 Int_t &nout, 
                                 Int_t tid)
{
// Do post-step actions on particle after scattering process. Surviving tracks
// copied in trackout
   // Compute the max theta angle opening after scattering.
   const Double_t ctmax = TMath::Cos(1.*TMath::DegToRad()); // 1 degree
   const Double_t de = gPropagator->fEmax-gPropagator->fEmin;
   Double_t theta, phi, scale,thetav,phiv; 
   Double_t dir[3];
   Double_t dirnew[3];
   Double_t *rndArray = gPropagator->fThreadData[tid]->fDblArray;
   Int_t irnd = 0;
   gPropagator->fThreadData[tid]->fRndm->RndmArray(2*ntracks, rndArray);
   for (Int_t i=0; i<ntracks; i++) {
      if (!tracks.fChargeV[i] || (tracks.fProcessV[i]!=0)) {
         nout++;
         continue;
      }   
      theta = TMath::ACos((1.-rndArray[irnd++]*(1.-ctmax)));
      // Re-scale from emin to emax
      scale = (tracks.fEV[i]-gPropagator->fEmin)/de;
      theta *= 1-scale;  // hi-energy don't scatter much
      phi = TMath::TwoPi()*rndArray[irnd++];
      // System along the (px,py,pz)
      thetav = TMath::ACos(tracks.fZdirV[i])*TMath::RadToDeg();
      phiv = TMath::ATan2(tracks.fYdirV[i],tracks.fXdirV[i])*TMath::RadToDeg();
      gPropagator->fThreadData[tid]->fRotation->SetAngles(phiv-90,-thetav,0);
      dir[0] = TMath::Sin(theta)*TMath::Cos(phi);
      dir[1] = TMath::Sin(theta)*TMath::Sin(phi);
      dir[2] = TMath::Cos(theta);
      gPropagator->fThreadData[tid]->fRotation->LocalToMaster(dir, dirnew);
      tracks.fXdirV[i] = dirnew[0];
      tracks.fYdirV[i] = dirnew[1];
      tracks.fZdirV[i] = dirnew[2];
      // All tracks survive
      nout++;
   }   
   gPropagator->fApplication->StepManager(tid, ntracks, tracks);
}

ClassImp(ElossProcess)

//______________________________________________________________________________
void ElossProcess::ComputeIntLen(TGeoMaterial *mat, 
                                 Int_t ntracks, 
                                 GeantTrack_v &tracks,
                                 Double_t *lengths, 
                                 Int_t /*tid*/)
{
// Energy loss process. Continuous process. Compute step limit for losing
// maximum dw per step.
   const Double_t dw = 1.E-3;  // 1 MEV
   Double_t mata = mat->GetA();
   Double_t matz = mat->GetZ();
   Double_t matr = mat->GetDensity();
   Bool_t invalid_material = kFALSE;
   if (matz<1 || mata<1 || matr<1.E-8) invalid_material = kTRUE;
   for (Int_t i=0; i<ntracks; i++) {
      if(tracks.fChargeV[i] && !invalid_material && tracks.fStatusV[i] != kKilled) {
         Double_t dedx = BetheBloch(tracks,i,matz,mata,matr);
         Double_t stepmax = (dedx>1.E-32)?dw/dedx:0.5*TMath::Limits<double>::Max();
         lengths[i] = stepmax;
      } else {
         lengths[i]=0.5*TMath::Limits<double>::Max();
      }      
   }
}

//______________________________________________________________________________
void ElossProcess::PostStep(TGeoMaterial *mat,
                                 Int_t ntracks,
                                 GeantTrack_v &tracks, 
                                 Int_t &nout, 
                                 Int_t tid)
{
// Do post-step actions after energy loss process. 
   Double_t eloss, dedx;
   Double_t mata = mat->GetA();
   Double_t matz = mat->GetZ();
   Double_t matr = mat->GetDensity();
   Bool_t invalid_material = kFALSE;
   if (matz<1 || mata<1 || matr<1.E-8) invalid_material = kTRUE;

   for (Int_t i=0; i<ntracks; i++) {
      if (tracks.fStatusV[i] == kKilled) continue;
      if (!tracks.fChargeV[i] || tracks.fStepV[i]==0 || invalid_material) {
         nout++;
         continue;
      }   
      if (tracks.fEV[i]-tracks.fMassV[i] < gPropagator->fEmin) {
         tracks.fStatusV[i] = kKilled;
//         gPropagator->StopTrack(track);
         continue;
      }   
      dedx = BetheBloch(tracks,i,matz,mata,matr);
      eloss = tracks.fStepV[i]*dedx;
      if (tracks.fEV[i]-tracks.fMassV[i]-eloss < gPropagator->fEmin) eloss = tracks.fEV[i]-tracks.fMassV[i];
      Double_t gammaold = tracks.Gamma(i);
      Double_t bgold = TMath::Sqrt((gammaold-1)*(gammaold+1));
      tracks.fEV[i] -= eloss;
      if (tracks.fEV[i]-tracks.fMassV[i] < gPropagator->fEmin) {
         tracks.fStatusV[i] = kKilled;
//         gPropagator->StopTrack(track);
         continue;
      }   
      nout++;

      Double_t gammanew = tracks.Gamma(i);
      Double_t bgnew = TMath::Sqrt((gammanew-1)*(gammanew+1));
      Double_t pnorm = bgnew/bgold;
      tracks.fPV[i] *= pnorm;
   }   
   gPropagator->fApplication->StepManager(tid, nout, tracks);
}

//______________________________________________________________________________
Double_t ElossProcess::BetheBloch(const GeantTrack_v &tracks, Int_t i, Double_t tz, Double_t ta, Double_t rho)
{
// Energy loss given by Bethe formula.
  if (tz<1. || ta<1.) return 0.;
  const Double_t konst = 0.1535; // MeV cm2/g
  const Double_t emass = 1000*TDatabasePDG::Instance()->GetParticle(kElectron)->Mass();
  const Double_t beta = tracks.Beta(i);
  const Double_t gamma = tracks.Gamma(i);
  const Double_t bg = beta*gamma;
  const Double_t wmax = 2*emass*bg*bg;
  Double_t ioniz;
  if(tz<13) ioniz = 12 + 7/tz;
  else ioniz = 9.76 + 58.8*TMath::Power(tz,-1.19);

  Double_t bethe = (konst * tz * rho * tracks.fChargeV[i] * tracks.fChargeV[i])/(ta * beta * beta);
//  Printf("ioniz %f",ioniz);
  bethe *= TMath::Log(2*emass*bg*bg*wmax*1e12/(ioniz*ioniz))-2*beta*beta;
//  Printf("bethe %f",bethe);
  return 1.e-3*bethe;
}

/*
//______________________________________________________________________________
Double_t ElossProcess::Bbf1(Double_t *x, Double_t *par)
{
  Double_t pimass = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
  GeantTrack t;
  Double_t bg = TMath::Power(10.,*x);
  Double_t gamma = TMath::Sqrt(bg*bg+1);

  t.px = bg*pimass;
  t.py = 0;
  t.pz = 0;
  t.e = gamma*pimass;
  t.charge = 1;
  t.mass = pimass;
  return 1000*BetheBloch(&t,par[0],par[1],par[2]);
}
//______________________________________________________________________________
void ElossProcess::PlotBB(Double_t z, Double_t a, Double_t rho, Double_t bgmin, Double_t bgmax)
{
  TF1 *f=new TF1("bb",ElossProcess::Bbf1,TMath::Log10(bgmin),TMath::Log10(bgmax),3);
  TH1F *h=new TH1F("hh","Bethe Bloch",100,TMath::Log10(bgmin),TMath::Log10(bgmax));
  h->SetMinimum(1.);
  h->SetMaximum(500.);
  f->SetParameter(0,z);
  f->SetParameter(1,a);
  f->SetParameter(2,rho);
  h->Draw();
  f->Draw("same");
}
*/

ClassImp(InteractionProcess)

//______________________________________________________________________________
InteractionProcess::InteractionProcess(const char *name)
                   :PhysicsProcess(name),
                    fNthreads(gPropagator->fNthreads),
                    fGen(0),
                    fMutex()
                   
{
// Ctor
   TObject::SetBit(kDiscrete);
   fGen = new TGenPhaseSpace[fNthreads];
}   

//______________________________________________________________________________
InteractionProcess::~InteractionProcess()
{
// dtor
   delete [] fGen;
}   

//______________________________________________________________________________
void InteractionProcess::ComputeIntLen(TGeoMaterial *mat, 
                                 Int_t ntracks, 
                                 GeantTrack_v &tracks,
                                 Double_t *lengths, 
                                 Int_t /*tid*/)
{
   Double_t fact = 1.0E-10;
   const Double_t nabarn = fact*TMath::Na()*1e-24;
   Double_t xlen = TMath::Limits<double>::Max();
   Double_t mata = mat->GetA();
   Double_t matz = mat->GetZ();
   Double_t matr = mat->GetDensity();
   Bool_t invalid_material = kFALSE;
   if (matz<1 || mata<1 || matr<1.E-8) invalid_material = kTRUE;
   if (!invalid_material) {
      Double_t density = TMath::Max(matr,1e-5);
      Double_t sigma = 28.5*TMath::Power(mata,0.75);
      xlen = mat->GetA()/(sigma*density*nabarn);
   } else {
      for (Int_t itr=0; itr<ntracks; itr++) lengths[itr] = 0.5*TMath::Limits<double>::Max();
      return;
   }   
 
   for (Int_t i=0; i<ntracks; i++) {
      if(tracks.fSpeciesV[i] == kHadron && tracks.fStatusV[i] != kKilled) {
         Double_t ek = tracks.fEV[i] - tracks.fMassV[i];
         lengths[i] = 10000.*xlen*(0.007+0.1*TMath::Log(ek)/ek+0.2/(ek*ek));
      } else {
         lengths[i] = 0.5*TMath::Limits<double>::Max();
      }
   }
}

//______________________________________________________________________________
void InteractionProcess::PostStep(TGeoMaterial *mat,
                                 Int_t ntracks,
                                 GeantTrack_v &tracks, 
                                 Int_t &nout, 
                                 Int_t tid)
{
// Do post-step actions on particle after interaction process. 
//   if (gUseDebug) Printf("PostStepInteraction %d tracks", ntracks);
// We calculate the CMS energy
// We suppose at first that the available energy is the Kin cms energy
// We produce equal number of pos and neg pions

   Double_t *rndArray = gPropagator->fThreadData[tid]->fDblArray;
   const Double_t pimass = TDatabasePDG::Instance()->GetParticle(kPiMinus)->Mass();
   const Double_t prodm[18] = {pimass, pimass, pimass, pimass, pimass, pimass,
			       pimass, pimass, pimass, pimass, pimass, pimass,
			       pimass, pimass, pimass, pimass, pimass, pimass};
   gPropagator->fThreadData[tid]->fRndm->RndmArray(ntracks, rndArray);

   Int_t nprod = 0;
   Int_t ngen  = 0;
   for (Int_t i=0; i<ntracks; i++) {
      if (tracks.fProcessV[i]!=2) {
         nout++;
         continue;
      }   
      Double_t en = tracks.fEV[i];
      Double_t m1 = tracks.fMassV[i];
      Double_t m2 = mat->GetA();
      Double_t cmsen = TMath::Sqrt(m1*m1+m2*m2+2*en*m2)-m1-m2;
      // Calculate the number of pions as a poisson distribution leaving half of the cms energy
      // for phase space momentum
      Int_t npi = 0.5*gPropagator->fThreadData[tid]->fRndm->Rndm()*cmsen/pimass+0.5;
      if(npi>1) {
         do { nprod = TMath::Min(gPropagator->fThreadData[tid]->fRndm->Poisson(npi),9); } 
         while(nprod*pimass*2>cmsen || nprod==0);
//         Printf("Inc en = %f, cms en = %f produced pis = %d",en,cmsen,nprod);
         TLorentzVector pcms(tracks.Px(i), tracks.Py(i), tracks.Pz(i), tracks.fEV[i] + m2);
         if(!fGen[tid].SetDecay(pcms,2*nprod,prodm)) Printf("Forbidden decay!");
         fMutex.Lock();
         fGen[tid].Generate();
         fMutex.UnLock();
         //Double_t pxtot=track->px;
         //Double_t pytot=track->py;
         //Double_t pztot=track->pz;
         // The mother particle dies
         tracks.fStatusV[i] = kKilled;
         GeantTrack &trackg = gPropagator->GetTempTrack(tid);
         for(Int_t j=0; j<2*nprod; ++j) {
            // Do not consider tracks below the production threshold. Normally the energy deposited should be taken into account
            TLorentzVector *lv = fGen[tid].GetDecay(j);
            if (lv->E()-pimass < gPropagator->fEmin) continue;
            *trackg.fPath = *tracks.fPathV[i];
            if(j%2) trackg.fPDG = kPiMinus;
            else trackg.fPDG = kPiPlus;
            trackg.fEvent = tracks.fEventV[i];
            trackg.fEvslot = tracks.fEvslotV[i];
            trackg.fSpecies = kHadron;
            trackg.fCharge = TDatabasePDG::Instance()->GetParticle(trackg.fPDG)->Charge()/3.;
            trackg.fMass = pimass;
//            trackg.fProcess = 0;
            trackg.fXpos = tracks.fXposV[i];
            trackg.fYpos = tracks.fYposV[i];
            trackg.fZpos = tracks.fZposV[i];
            Double_t oneoverp = 1./lv->P();
            trackg.fXdir = oneoverp*lv->Px();
            trackg.fYdir = oneoverp*lv->Py();
            trackg.fZdir = oneoverp*lv->Pz();
            trackg.fE = lv->E();
            trackg.fP = TMath::Sqrt((trackg.E()-trackg.Mass())*(trackg.E()+trackg.Mass()));
            if (TMath::IsNaN(trackg.fXdir)) {
               Printf("NaN");
            }   
            trackg.fStatus = kNew;
            ngen++;
            // Add track to the tracks vector
            gPropagator->AddTrack(trackg);
            tracks.AddTrack(trackg);
         }
      }
   }   
   gPropagator->fApplication->StepManager(tid, nout, tracks);
}
