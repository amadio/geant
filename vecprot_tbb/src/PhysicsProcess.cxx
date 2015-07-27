// Toy physics processes for our propagator prototype. Currently including:
// - single scattering as a discrete process
// - energy loss as continuous process
// - generic interaction as discrete process, producing secondaries

#include "PhysicsProcess.h"

#include "GeantVolumeBasket.h"
#include "GeantThreadData.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"

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

#include "base/Global.h"

using vecgeom::kTwoPi;
using vecgeom::kDegToRad;
using vecgeom::kRadToDeg;
using vecgeom::kAvogadro;
using std::numeric_limits;
using std::min;

ClassImp(PhysicsProcess)

    //______________________________________________________________________________
    void PhysicsProcess::StepManager(int iproc, int npart, int * /*particles*/, int nout,
                                     int * /*partnext*/) {
  GeantPropagator *gPropagator = GeantPropagator::Instance();
  // User stepping routine. <partnext> array can
  // be null.
  //   for (int ipart=0; ipart<npart; ipart++) gTracks[particles[ipart]]->nsteps++;
  if (gPropagator->fUseDebug) {
    Printf("StepManager: process %s, npart=%d, nout=%d", gPropagator->Process(iproc)->GetName(), npart, nout);
  }
}

ClassImp(ScatteringProcess)

    //______________________________________________________________________________
    void ScatteringProcess::ComputeIntLen(TGeoVolume *vol, int ntracks, int *trackin, double *lengths) {
  // Generates an interaction length for the scattering process. Nothing physical,
  // just generate something comparable with the size of the current volume.
  //
  //
  // trackin and lengths arrays should be be correctly set by the caller

  GeantPropagator *gPropagator = GeantPropagator::Instance();
  PerThread::reference TBBperThread = gPropagator->fTBBthreadData.local();

  const double kC1 = 500.;
  const double xlen = numeric_limits<double>.max();
  int itrack;
  double density = 1.e-5;
  TGeoMaterial *mat = vol->GetMaterial();
  if (mat)
    density = mat->GetDensity();
  density = max<double>(density, 1e-3);
  // Make sure we write in the thread space for the current basket
  double *rndArray = TBBperThread.fDblArray;
  int irnd = 0;
  TBBperThread.fRndm->RndmArray(ntracks, rndArray);
  GeantTrack *track = 0;
  for (int i = 0; i < ntracks; i++) {
    itrack = trackin[i];
    track = gPropagator->fTracks[itrack];
    if (!track->charge)
      lengths[i] = 0.5 * xlen;
    else if (gPropagator->fTracks[itrack]->IsAlive())
      lengths[i] = kC1 * gPropagator->fTracks[itrack]->e * rndArray[irnd++] / density;
    else
      lengths[i] = 0.5 * xlen;
  }
}

//______________________________________________________________________________
void ScatteringProcess::PostStep(TGeoVolume *, int ntracks, int *trackin, int &nout, int *trackout) {
  // Do post-step actions on particle after scattering process. Surviving tracks
  // copied in trackout

  GeantPropagator *gPropagator = GeantPropagator::Instance();
  PerThread::reference TBBperThread = gPropagator->fTBBthreadData.local();

  // Compute the max theta angle opening after scattering.
  const double ctmax = cos(1. * kDegToRad); // 1 degree
  double theta, phi, scale, thetav, phiv;
  double dir[3];
  double dirnew[3];
  GeantTrack *track = 0;
  int itrack;
  double p;
  double *rndArray = TBBperThread.fDblArray;
  int irnd = 0;
  TBBperThread.fRndm->RndmArray(2 * ntracks, rndArray);
  for (int i = 0; i < ntracks; i++) {
    itrack = trackin[i];
    track = gPropagator->fTracks[itrack];
    if (!track->charge) {
      if (trackout)
        trackout[nout] = itrack;
      nout++;
      continue;
    }
    theta = acos((1. - rndArray[irnd++] * (1. - ctmax)));
    // Re-scale from emin to emax
    scale = (track->e - gPropagator->fEmin) / (gPropagator->fEmax - gPropagator->fEmin);
    theta *= 1 - scale; // hi-energy don't scatter much
    phi = kTwoPi * rndArray[irnd++];
    // System along the (px,py,pz)
    p = track->P();
    thetav = acos(track->pz / p) * kRadToDeg;
    phiv = atan2(track->py, track->px) * kRadToDeg;
    TBBperThread.fRotation->SetAngles(phiv - 90, -thetav, 0);
    dir[0] = sin(theta) * cos(phi);
    dir[1] = sin(theta) * sin(phi);
    dir[2] = cos(theta);
    TBBperThread.fRotation->LocalToMaster(dir, dirnew);
    track->px = p * dirnew[0];
    track->py = p * dirnew[1];
    track->pz = p * dirnew[2];
    // All tracks survive
    if (trackout)
      trackout[nout] = itrack;
    nout++;
  }
  StepManager(0, ntracks, trackin, nout, trackout);
}

ClassImp(ElossProcess)

    //______________________________________________________________________________
    void ElossProcess::ComputeIntLen(TGeoVolume *vol, int ntracks, int *trackin, double *lengths) {
  // Energy loss process. Continuous process. Compute step limit for losing
  // maximum dw per step.
  GeantPropagator *gPropagator = GeantPropagator::Instance();
  const double dw = 1.E-3; // 1 MEV
  TGeoMaterial *mat = vol->GetMaterial();
  double mata = mat->GetA();
  double matz = mat->GetZ();
  double matr = mat->GetDensity();
  Bool_t invalid_material = kFALSE;
  if (matz < 1 || mata < 1 || matr < 1.E-8)
    invalid_material = kTRUE;
  int itrack;
  GeantTrack *track;
  for (int i = 0; i < ntracks; i++) {
    itrack = trackin[i];
    track = gPropagator->fTracks[itrack];
    if (track->charge && !invalid_material && track->IsAlive()) {
      double dedx = BetheBloch(track, matz, mata, matr);
      double stepmax = (dedx > 1.E-32) ? dw / dedx : 0.5 * numeric_limits<double>.max();
      lengths[i] = stepmax;
    } else {
      lengths[i] = 0.5 * numeric_limits<double>.max();
    }
  }
}

//______________________________________________________________________________
void ElossProcess::PostStep(TGeoVolume *vol, int ntracks, int *trackin, int &nout, int *trackout) {
  // Do post-step actions after energy loss process.
  GeantPropagator *gPropagator = GeantPropagator::Instance();
  double eloss, dedx;
  GeantTrack *track;
  int itrack;
  TGeoMaterial *mat = vol->GetMaterial();
  double mata = mat->GetA();
  double matz = mat->GetZ();
  double matr = mat->GetDensity();
  Bool_t invalid_material = kFALSE;
  if (matz < 1 || mata < 1 || matr < 1.E-8)
    invalid_material = kTRUE;

  for (int i = 0; i < ntracks; i++) {
    itrack = trackin[i];
    track = gPropagator->fTracks[itrack];
    if (!track->IsAlive())
      continue;
    if (!track->charge || track->step == 0 || invalid_material) {
      if (trackout)
        trackout[nout] = itrack;
      nout++;
      continue;
    }
    if (track->e - track->mass < gPropagator->fEmin) {
      gPropagator->StopTrack(track);
      continue;
    }
    dedx = BetheBloch(track, matz, mata, matr);
    eloss = track->step * dedx;
    if (track->e - track->mass - eloss < gPropagator->fEmin)
      eloss = track->e - track->mass;
    double gammaold = track->Gamma();
    double bgold = sqrt((gammaold - 1) * (gammaold + 1));
    track->e -= eloss;
    if (track->e - track->mass < gPropagator->fEmin) {
      gPropagator->StopTrack(track);
      continue;
    }
    if (trackout)
      trackout[nout] = itrack;
    nout++;

    double gammanew = track->Gamma();
    double bgnew = sqrt((gammanew - 1) * (gammanew + 1));
    double pnorm = bgnew / bgold;
    track->px *= pnorm;
    track->py *= pnorm;
    track->pz *= pnorm;
  }
  StepManager(1, ntracks, trackin, nout, trackout);
}

//______________________________________________________________________________
double ElossProcess::BetheBloch(GeantTrack *track, double tz, double ta, double rho) {
  // Energy loss given by Bethe formula.
  if (tz < 1. || ta < 1.)
    return 0.;
  const double konst = 0.1535; // MeV cm2/g
  const double emass = 1000 * TDatabasePDG::Instance()->GetParticle(kElectron)->Mass();
  const double beta = track->Beta();
  const double gamma = track->Gamma();
  const double bg = beta * gamma;
  const double wmax = 2 * emass * bg * bg;
  double ioniz;
  if (tz < 13)
    ioniz = 12 + 7 / tz;
  else
    ioniz = 9.76 + 58.8 * pow(tz, -1.19);

  double bethe = (konst * tz * rho * track->charge * track->charge) / (ta * beta * beta);
  //  Printf("ioniz %f",ioniz);
  bethe *= log(2 * emass * bg * bg * wmax * 1e12 / (ioniz * ioniz)) - 2 * beta * beta;
  //  Printf("bethe %f",bethe);
  return 1.e-3 * bethe;
}

//______________________________________________________________________________
double ElossProcess::Bbf1(double *x, double *par) {
  double pimass = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
  GeantTrack t;
  double bg = pow(10., *x);
  double gamma = sqrt(bg * bg + 1);

  t.px = bg * pimass;
  t.py = 0;
  t.pz = 0;
  t.e = gamma * pimass;
  t.charge = 1;
  t.mass = pimass;
  return 1000 * BetheBloch(&t, par[0], par[1], par[2]);
}

//______________________________________________________________________________
void ElossProcess::PlotBB(double z, double a, double rho, double bgmin, double bgmax) {
  TF1 *f = new TF1("bb", ElossProcess::Bbf1, log10(bgmin), log10(bgmax), 3);
  TH1F *h = new TH1F("hh", "Bethe Bloch", 100, log10(bgmin), log10(bgmax));
  h->SetMinimum(1.);
  h->SetMaximum(500.);
  f->SetParameter(0, z);
  f->SetParameter(1, a);
  f->SetParameter(2, rho);
  h->Draw();
  f->Draw("same");
}

ClassImp(InteractionProcess)

    //______________________________________________________________________________
    void InteractionProcess::ComputeIntLen(TGeoVolume *vol, int ntracks, int *trackin, double *lengths) {
  GeantPropagator *gPropagator = GeantPropagator::Instance();
  double fact = 1.;
  const double nabarn = fact * kAvogadro * 1e-24;
  int itrack;
  GeantTrack *track;
  double xlen = numeric_limits<double>.max();
  TGeoMaterial *mat = vol->GetMaterial();
  double mata = mat->GetA();
  double matz = mat->GetZ();
  double matr = mat->GetDensity();
  Bool_t invalid_material = kFALSE;
  if (matz < 1 || mata < 1 || matr < 1.E-8)
    invalid_material = kTRUE;
  if (!invalid_material) {
    double density = max<double>(matr, 1e-5);
    double sigma = 28.5 * pow(mata, 0.75);
    xlen = mat->GetA() / (sigma * density * nabarn);
  } else {
    for (itrack = 0; itrack < ntracks; itrack++)
      lengths[itrack] = 0.5 * numeric_limits<double>.max();
    return;
  }

  for (int i = 0; i < ntracks; i++) {
    itrack = trackin[i];
    track = gPropagator->fTracks[itrack];
    if (track->species == kHadron && track->IsAlive()) {
      double ek = track->e - track->mass;
      lengths[i] = xlen * (0.007 + 0.1 * log(ek) / ek + 0.2 / (ek * ek));
    } else {
      lengths[i] = 0.5 * numeric_limits<double>.max();
    }
  }
}

//______________________________________________________________________________
void InteractionProcess::PostStep(TGeoVolume *vol, int ntracks, int *trackin, int &nout, int *trackout) {
  // Do post-step actions on particle after interaction process.
  //   if (gUseDebug) Printf("PostStepInteraction %d tracks", ntracks);
  // We calculate the CMS energy
  // We suppose at first that the available energy is the Kin cms energy
  // We produce equal number of pos and neg pions

  GeantPropagator *gPropagator = GeantPropagator::Instance();
  PerThread::reference TBBperThread = gPropagator->fTBBthreadData.local();

  static TGenPhaseSpace gps;
  GeantTrack *track;
  int itrack;
  double *rndArray = TBBperThread.fDblArray;
  const double pimass = TDatabasePDG::Instance()->GetParticle(kPiMinus)->Mass();
  const double prodm[18] = {pimass, pimass, pimass, pimass, pimass, pimass, pimass, pimass, pimass,
                              pimass, pimass, pimass, pimass, pimass, pimass, pimass, pimass, pimass};
  TBBperThread.fRndm->RndmArray(ntracks, rndArray);

  int nprod = 0;
  int ngen = 0;
  for (int i = 0; i < ntracks; i++) {
    itrack = trackin[i];
    track = gPropagator->fTracks[itrack];
    double en = track->e;
    double m1 = track->mass;
    double m2 = TBBperThread.fVolume->GetMaterial()->GetA();
    double cmsen = sqrt(m1 * m1 + m2 * m2 + 2 * en * m2) - m1 - m2;
    // Calculate the number of pions as a poisson distribution leaving half of the cms energy
    // for phase space momentum
    int npi = 0.5 * TBBperThread.fRndm->Rndm() * cmsen / pimass + 0.5;
    if (npi > 1) {
      do {
        nprod = min<int>(TBBperThread.fRndm->Poisson(npi), 9);
      } while (nprod * pimass * 2 > cmsen || nprod == 0);
      //         Printf("Inc en = %f, cms en = %f produced pis = %d",en,cmsen,nprod);
      TLorentzVector pcms(track->px, track->py, track->pz, track->e + m2);
      if (!gps.SetDecay(pcms, 2 * nprod, prodm))
        Printf("Forbidden decay!");
      gps.Generate();
      // double pxtot=track->px;
      // double pytot=track->py;
      // double pztot=track->pz;
      TGeoBranchArray &a = *track->path;
      for (int j = 0; j < 2 * nprod; ++j) {
        GeantTrack *trackg = gPropagator->AddTrack(track->evslot);
        *trackg->path = a;
        TLorentzVector *lv = gps.GetDecay(j);
        if (j % 2)
          trackg->pdg = kPiMinus;
        else
          trackg->pdg = kPiPlus;
        trackg->event = track->event;
        trackg->evslot = track->evslot;
        trackg->species = kHadron;
        trackg->charge = TDatabasePDG::Instance()->GetParticle(trackg->pdg)->Charge() / 3.;
        trackg->mass = pimass;
        trackg->process = 0;
        trackg->xpos = track->xpos;
        trackg->ypos = track->ypos;
        trackg->zpos = track->zpos;
        trackg->px = lv->Px();
        trackg->py = lv->Py();
        trackg->pz = lv->Pz();
        trackg->e = lv->E();
        //            double mm2 =
        //            trackg->e*trackg->e-trackg->px*trackg->px-trackg->py*trackg->py-trackg->pz*trackg->pz;
        int itracknew = trackg->particle;
        trackout[nout++] = itracknew;
        ngen++;

        /*GeantVolumeBasket *basket = gPropagator->fWMgr->GetCurrentBasket(tid);*/
        GeantVolumeBasket *basket = (GeantVolumeBasket *)TBBperThread.fVolume->GetField();

        if (basket)
          TBBperThread.fCollection->AddTrack(itracknew, basket);
        // check
        // pxtot -= trackg->px;
        // pytot -= trackg->py;
        // pztot -= trackg->pz;
      }
      //	Printf("pbal = %f %f %f",pxtot, pytot, pztot);
    }
  }
  StepManager(2, ntracks, trackin, nout, trackout);
  if (ngen) {
    // Generated particles may be below threshold-> Call PostStepEloss
    int nsurv = 0;
    int *trackgen = new int[ngen];
    gPropagator->Process(1)->PostStep(vol, ngen, &trackout[nout - ngen], nsurv, trackgen);
    memcpy(&trackout[nout - ngen], trackgen, nsurv * sizeof(int));
    nout += nsurv - ngen;
    delete[] trackgen;
  }
}
