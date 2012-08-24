// Toy physics processes for our propagator prototype. Currently including:
// - single scattering as a discrete process
// - energy loss as continuous process
// - generic interaction as discrete process, producing secondaries

#include "PhysicsProcess.h"
#include "GeantVolumeBasket.h"
#include "GeantThreadData.h"
#include "WorkloadManager.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "globals.h"
#include "GeantTrack.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoBranchArray.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "TGenPhaseSpace.h"

ClassImp(PhysicsProcess)

//______________________________________________________________________________
void PhysicsProcess::StepManager(Int_t iproc, Int_t npart, Int_t */*particles*/, Int_t nout, Int_t */*partnext*/)
{
// User stepping routine. <partnext> array can
// be null.
//   for (Int_t ipart=0; ipart<npart; ipart++) gTracks[particles[ipart]]->nsteps++;
   if (gPropagator->fUseDebug) {
      Printf("StepManager: process %s, npart=%d, nout=%d", gPropagator->Process(iproc)->GetName(), npart, nout);
   }   
}

ClassImp(ScatteringProcess)

//______________________________________________________________________________
void ScatteringProcess::ComputeIntLen(TGeoVolume *vol, 
                                      Int_t ntracks, 
                                      Int_t *trackin, 
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
   Int_t itrack;
   Double_t density = 1.e-5;
   TGeoMaterial *mat = vol->GetMaterial();
   if (mat) density = mat->GetDensity();
   density = TMath::Max(density, 1.E-3);
   // Make sure we write in the thread space for the current basket
   Double_t *rndArray = gPropagator->fThreadData[tid]->fDblArray;
   Int_t irnd = 0;
   gPropagator->fThreadData[tid]->fRndm->RndmArray(ntracks, rndArray);
   GeantTrack *track = 0;
   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];
      track = gPropagator->fTracks[itrack];
      if (!track->charge) lengths[i] =  0.5*xlen;
      else if (gPropagator->fTracks[itrack]->IsAlive())
        lengths[i] = kC1*gPropagator->fTracks[itrack]->e*rndArray[irnd++]/density;
      else
        lengths[i] =  0.5*xlen; 
   }
}

//______________________________________________________________________________
void ScatteringProcess::PostStep(TGeoVolume *,
                                 Int_t ntracks,
                                 Int_t *trackin, 
                                 Int_t &nout, 
                                 Int_t* trackout, 
                                 Int_t tid)
{
// Do post-step actions on particle after scattering process. Surviving tracks
// copied in trackout
   // Compute the max theta angle opening after scattering.
   const Double_t ctmax = TMath::Cos(1.*TMath::DegToRad()); // 1 degree
   Double_t theta, phi, scale,thetav,phiv; 
   Double_t dir[3];
   Double_t dirnew[3];
   GeantTrack *track = 0;
   Int_t itrack;
   Double_t p;
   Double_t *rndArray = gPropagator->fThreadData[tid]->fDblArray;
   Int_t irnd = 0;
   gPropagator->fThreadData[tid]->fRndm->RndmArray(2*ntracks, rndArray);
   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];
      track = gPropagator->fTracks[itrack];
      if (!track->charge) {
         if (trackout) trackout[nout] = itrack;
         nout++;
         continue;
      }   
      theta = TMath::ACos((1.-rndArray[irnd++]*(1.-ctmax)));
      // Re-scale from emin to emax
      scale = (track->e-gPropagator->fEmin)/(gPropagator->fEmax-gPropagator->fEmin);
      theta *= 1-scale;  // hi-energy don't scatter much
      phi = TMath::TwoPi()*rndArray[irnd++];
      // System along the (px,py,pz)
      p = track->P();
      thetav = TMath::ACos(track->pz/p)*TMath::RadToDeg();
      phiv = TMath::ATan2(track->py,track->px)*TMath::RadToDeg();
      gPropagator->fThreadData[tid]->fRotation->SetAngles(phiv-90,-thetav,0);
      dir[0] = TMath::Sin(theta)*TMath::Cos(phi);
      dir[1] = TMath::Sin(theta)*TMath::Sin(phi);
      dir[2] = TMath::Cos(theta);
      gPropagator->fThreadData[tid]->fRotation->LocalToMaster(dir, dirnew);
      track->px = p*dirnew[0];
      track->py = p*dirnew[1];
      track->pz = p*dirnew[2];
      // All tracks survive
      if (trackout) trackout[nout] = itrack;
      nout++;
   }   
   StepManager(0, ntracks, trackin, nout, trackout);
}

ClassImp(ElossProcess)

//______________________________________________________________________________
void ElossProcess::ComputeIntLen(TGeoVolume *vol, 
                                      Int_t ntracks, 
                                      Int_t *trackin, 
                                      Double_t *lengths, 
                                      Int_t /*tid*/)
{
// Energy loss process. Continuous process. Compute step limit for losing
// maximum dw per step.
   const Double_t dw = 1.E-3;  // 1 MEV
   TGeoMaterial *mat = vol->GetMaterial();
   Double_t mata = mat->GetA();
   Double_t matz = mat->GetZ();
   Double_t matr = mat->GetDensity();
   Bool_t invalid_material = kFALSE;
   if (matz<1 || mata<1 || matr<1.E-8) invalid_material = kTRUE;
   Int_t itrack;
   GeantTrack *track;
   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];  
      track = gPropagator->fTracks[itrack];
      if(track->charge && !invalid_material && track->IsAlive()) {
         Double_t dedx = BetheBloch(track,matz,mata,matr);
         Double_t stepmax = (dedx>1.E-32)?dw/dedx:0.5*TMath::Limits<double>::Max();
         lengths[i] = stepmax;
      } else {
         lengths[i]=0.5*TMath::Limits<double>::Max();
      }      
   }
}

//______________________________________________________________________________
void ElossProcess::PostStep(TGeoVolume *vol,
                                 Int_t ntracks,
                                 Int_t *trackin, 
                                 Int_t &nout, 
                                 Int_t* trackout, 
                                 Int_t tid)
{
// Do post-step actions after energy loss process. 
   Double_t eloss, dedx;
   GeantTrack *track;
   Int_t itrack;
   TGeoMaterial *mat = vol->GetMaterial();
   Double_t mata = mat->GetA();
   Double_t matz = mat->GetZ();
   Double_t matr = mat->GetDensity();
   Bool_t invalid_material = kFALSE;
   if (matz<1 || mata<1 || matr<1.E-8) invalid_material = kTRUE;

   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];   
      track = gPropagator->fTracks[itrack];
      if (!track->IsAlive()) continue;
      if (!track->charge || track->step==0 || invalid_material) {
         if (trackout) trackout[nout] = itrack;
         nout++;
         continue;
      }   
      if (track->e-track->mass < gPropagator->fEmin) {
         gPropagator->StopTrack(track);
         continue;
      }   
      dedx = BetheBloch(track,matz,mata,matr);
      eloss = track->step*dedx;
      if (track->e-track->mass-eloss < gPropagator->fEmin) eloss = track->e-track->mass;
      Double_t gammaold = track->Gamma();
      Double_t bgold = TMath::Sqrt((gammaold-1)*(gammaold+1));
      track->e -= eloss;
      if (track->e-track->mass < gPropagator->fEmin) {
         gPropagator->StopTrack(track);
         continue;
      }   
      if (trackout) trackout[nout] = itrack;
      nout++;

      Double_t gammanew = track->Gamma();
      Double_t bgnew = TMath::Sqrt((gammanew-1)*(gammanew+1));
      Double_t pnorm = bgnew/bgold;
      track->px *= pnorm;
      track->py *= pnorm;
      track->pz *= pnorm;
   }   
   StepManager(1, ntracks, trackin, nout, trackout);
}

//______________________________________________________________________________
Double_t ElossProcess::BetheBloch(GeantTrack* track, Double_t tz, Double_t ta, Double_t rho)
{
// Energy loss given by Bethe formula.
  if (tz<1. || ta<1.) return 0.;
  const Double_t konst = 0.1535; // MeV cm2/g
  const Double_t emass = 1000*TDatabasePDG::Instance()->GetParticle(kElectron)->Mass();
  const Double_t beta = track->Beta();
  const Double_t gamma = track->Gamma();
  const Double_t bg = beta*gamma;
  const Double_t wmax = 2*emass*bg*bg;
  Double_t ioniz;
  if(tz<13) ioniz = 12 + 7/tz;
  else ioniz = 9.76 + 58.8*TMath::Power(tz,-1.19);

  Double_t bethe = (konst * tz * rho * track->charge * track->charge)/(ta * beta * beta);
//  Printf("ioniz %f",ioniz);
  bethe *= TMath::Log(2*emass*bg*bg*wmax*1e12/(ioniz*ioniz))-2*beta*beta;
//  Printf("bethe %f",bethe);
  return 1.e-3*bethe;
}

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
void PlotBB(Double_t z, Double_t a, Double_t rho, Double_t bgmin=1e-2, Double_t bgmax=1e6)
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

ClassImp(InteractionProcess)

//______________________________________________________________________________
void InteractionProcess::ComputeIntLen(TGeoVolume *vol, 
                                      Int_t ntracks, 
                                      Int_t *trackin, 
                                      Double_t *lengths, 
                                      Int_t /*tid*/)
{
   Double_t fact = 1.;
   const Double_t nabarn = fact*TMath::Na()*1e-24;
   Int_t itrack;
   GeantTrack *track;
   Double_t xlen = TMath::Limits<double>::Max();
   TGeoMaterial *mat = vol->GetMaterial();
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
      for (itrack=0; itrack<ntracks; itrack++) lengths[itrack] = 0.5*TMath::Limits<double>::Max();
      return;
   }   
 
   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];
      track = gPropagator->fTracks[itrack];
      if(track->species == kHadron && track->IsAlive()) {
         Double_t ek = track->e - track->mass;
         lengths[i] = xlen*(0.007+0.1*TMath::Log(ek)/ek+0.2/(ek*ek));
      } else {
         lengths[i] = 0.5*TMath::Limits<double>::Max();
      }
   }
}

//______________________________________________________________________________
void InteractionProcess::PostStep(TGeoVolume *vol,
                                 Int_t ntracks,
                                 Int_t *trackin, 
                                 Int_t &nout, 
                                 Int_t* trackout, 
                                 Int_t tid)
{
// Do post-step actions on particle after interaction process. 
//   if (gUseDebug) Printf("PostStepInteraction %d tracks", ntracks);
// We calculate the CMS energy
// We suppose at first that the available energy is the Kin cms energy
// We produce equal number of pos and neg pions

   static TGenPhaseSpace gps;
   GeantTrack *track;
   Int_t itrack;
   Double_t *rndArray = gPropagator->fThreadData[tid]->fDblArray;
   const Double_t pimass = TDatabasePDG::Instance()->GetParticle(kPiMinus)->Mass();
   const Double_t prodm[18] = {pimass, pimass, pimass, pimass, pimass, pimass,
			       pimass, pimass, pimass, pimass, pimass, pimass,
			       pimass, pimass, pimass, pimass, pimass, pimass};
   gPropagator->fThreadData[tid]->fRndm->RndmArray(ntracks, rndArray);

   Int_t nprod = 0;
   Int_t ngen  = 0;
   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];
      track = gPropagator->fTracks[itrack];
      Double_t en = track->e;
      Double_t m1 = track->mass;
      Double_t m2 = gPropagator->fThreadData[tid]->fVolume->GetMaterial()->GetA();
      Double_t cmsen = TMath::Sqrt(m1*m1+m2*m2+2*en*m2)-m1-m2;
      // Calculate the number of pions as a poisson distribution leaving half of the cms energy
      // for phase space momentum
      Int_t npi = 0.5*gPropagator->fThreadData[tid]->fRndm->Rndm()*cmsen/pimass+0.5;
      if(npi>1) {
         do { nprod = TMath::Min(gPropagator->fThreadData[tid]->fRndm->Poisson(npi),9); } 
         while(nprod*pimass*2>cmsen || nprod==0);
//         Printf("Inc en = %f, cms en = %f produced pis = %d",en,cmsen,nprod);
         TLorentzVector pcms(track->px, track->py, track->pz, track->e + m2);
         if(!gps.SetDecay(pcms,2*nprod,prodm)) Printf("Forbidden decay!");
         gps.Generate();
         //Double_t pxtot=track->px;
         //Double_t pytot=track->py;
         //Double_t pztot=track->pz;
         TGeoBranchArray &a = *track->path;
         for(Int_t j=0; j<2*nprod; ++j) {
            GeantTrack *trackg = gPropagator->AddTrack(track->evslot);
            *trackg->path = a;
            TLorentzVector *lv = gps.GetDecay(j);
            if(j%2) trackg->pdg = kPiMinus;
            else trackg->pdg = kPiPlus;
            trackg->event = track->event;
            trackg->evslot = track->evslot;
            trackg->species = kHadron;
            trackg->charge = TDatabasePDG::Instance()->GetParticle(trackg->pdg)->Charge()/3.;
            trackg->mass = pimass;
            trackg->process = 0;
            trackg->xpos = track->xpos;
            trackg->ypos = track->ypos;
            trackg->zpos = track->zpos;
            trackg->px = lv->Px();
            trackg->py = lv->Py();
            trackg->pz = lv->Pz();
            trackg->e = lv->E();
//            Double_t mm2 = trackg->e*trackg->e-trackg->px*trackg->px-trackg->py*trackg->py-trackg->pz*trackg->pz;
            Int_t itracknew = trackg->particle;
            trackout[nout++] = itracknew;
            ngen++;
            GeantVolumeBasket *basket = gPropagator->fWMgr->GetCurrentBasket(tid);
            if (basket) gPropagator->fCollections[tid]->AddTrack(itracknew, basket);
           //check
           //pxtot -= trackg->px;
           //pytot -= trackg->py;
           //pztot -= trackg->pz;
         }
	//	Printf("pbal = %f %f %f",pxtot, pytot, pztot);
      }
   }   
   StepManager(2, ntracks, trackin, nout, trackout);
   if (ngen) {
      // Generated particles may be below threshold-> Call PostStepEloss
      Int_t nsurv = 0;
      Int_t *trackgen = new Int_t[ngen];
      gPropagator->Process(1)->PostStep(vol, ngen, &trackout[nout-ngen], nsurv, trackgen,tid);
      memcpy(&trackout[nout-ngen], trackgen, nsurv*sizeof(Int_t));
      nout += nsurv-ngen;
      delete [] trackgen;
   }
}
