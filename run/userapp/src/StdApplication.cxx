#include "StdApplication.h"
#include "Propagator.h"
#include "GeantTaskData.h"
#include "globals.h"
#include "Geant/Error.h"
#ifdef USE_ROOT
#include "TProfile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#endif

//______________________________________________________________________________
StdApplication::StdApplication(RunManager *runmgr)
  : geant::GeantVApplication(runmgr), fInitialized(false),
#ifdef USE_ROOT
    fHeta(0), fHpt(0), fHStep(0), fStepSize(0), fStepCnt(0),
#endif
    fMHist(), fScore(kScore) {
  // Ctor.

#ifdef USE_ROOT
  double *array = 0;
  TH1::AddDirectory(false);
  fHeta = new TH1F("hEta", "Eta distribution per step", 50, -8., 8.);
  fHpt = new TH1F("hPt", "Pt distribution per step", 50, 0.1, 100.);
  array = MakeUniformLogArray(100, 1.e-7, 1.E3);
  fHStep = new TH1D("hSteps", "Distribution of small steps", 100, array);
  delete[] array;
  fStepSize = new TProfile("hStepEta", "Profile of step size with eta", 50, -8, 8);
  fStepCnt = new TProfile("hStepCnt", "Profile of step count with eta", 50, -8, 8);
  TH1::AddDirectory(true);
#endif
}

//______________________________________________________________________________
double *StdApplication::MakeUniformLogArray(int nbins, double lmin, double lmax) {
  // Create and fill a log scale bin limits array with nbins between lmin and lmax
  // To be passed to TH1D constructor. User responsability to delete.
  const double l10 = log(10.);
  if ((lmin <= 0) || (lmax <= 0))
    return 0;
  double *array = new double[nbins + 1];
  double lminlog = log10(lmin);
  double lmaxlog = log10(lmax);
  double dstep = (lmaxlog - lminlog) / nbins;
  for (auto i = 0; i <= nbins; ++i) {
    array[i] = exp(l10 * (lminlog + i * dstep));
  }
  // Use array as:
  // TH1D *hist = new TH1D("name", "title", nbins, array);
  // delete [] array;
  return array;
}

//______________________________________________________________________________
bool StdApplication::Initialize() {
  // Initialize application. Geometry must be loaded.
  if (fInitialized)
    return true;
  geant::Printf("=== StdApplication::Initialize done");
  fInitialized = true;
  return true;
}

//______________________________________________________________________________
void StdApplication::SteppingActions(Track &track, GeantTaskData * td) {
  // Application stepping actions.
#ifdef USE_ROOT
  Propagator *propagator = td->fPropagator;
  if ((!fInitialized) || (fScore == kNoScore))
    return;
  // Loop all tracks, check if they are in the right volume and collect the
  // energy deposit and step length
  double theta, eta;
  if (track.Dz() == 1)
    eta = 1.E30;
  else {
    theta = acos(track.Dz());
    eta = -log(tan(0.5 * theta));
  }
  if (propagator->fNthreads > 1)
    fMHist.lock();
  fHeta->Fill(eta);
  fHpt->Fill(track.Pt());
  fHStep->Fill(track.GetStep());
  fStepSize->Fill(eta, track.GetStep());
  if ((track.Status() == kKilled) || (track.Status() == kExitingSetup) ||
      (track.Path()->IsOutside()))
      fStepCnt->Fill(eta, track.GetNsteps());
  if (propagator->fNthreads > 1)
    fMHist.unlock();
  }
#else
  (void)track;
  (void)td;
#endif

//______________________________________________________________________________
void StdApplication::FinishRun() {
#ifdef USE_ROOT
  if (fScore == kNoScore)
    return;
  TVirtualPad *pad;
  TCanvas *c3 = new TCanvas("Step size profile", "Standard GeantV scoring", 800, 1600);
  c3->Divide(2, 3);
  pad = c3->cd(1);
  fHeta->Sumw2();
  fHeta->Draw("E");
  pad = c3->cd(2);
  pad->SetLogx();
  pad->SetLogy();
  fHpt->Sumw2();
  fHpt->Draw("E");
  pad = c3->cd(3);
  pad->SetLogx();
  pad->SetLogy();
  fHStep->Sumw2();
  fHStep->Draw("E");
  pad = c3->cd(4);
  pad->SetLogy();
  fStepSize->Draw();
  pad = c3->cd(5);
  pad->SetLogy();
  fStepCnt->Draw();
#endif
}
