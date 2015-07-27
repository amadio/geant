#include "StdApplication.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "GeantTaskData.h"
#include "globals.h"
#include "TProfile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"

ClassImp(StdApplication)

    //______________________________________________________________________________
    StdApplication::StdApplication()
    : GeantVApplication(), fInitialized(kFALSE), fHeta(0), fHpt(0), fHStep(0), fStepSize(0), fStepCnt(0), fMHist(),
      fScore(kScore) {
  // Ctor.
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
}

//______________________________________________________________________________
double *StdApplication::MakeUniformLogArray(Int_t nbins, double lmin, double lmax) {
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
Bool_t StdApplication::Initialize() {
  // Initialize application. Geometry must be loaded.
  if (fInitialized)
    return kTRUE;
  Printf("=== StdApplication::Initialize done");
  fInitialized = kTRUE;
  return kTRUE;
}

//______________________________________________________________________________
void StdApplication::StepManager(Int_t npart, const GeantTrack_v &tracks, GeantTaskData * /*td*/) {
  // Application stepping manager. The thread id has to be used to manage storage
  // of hits independently per thread.
  static GeantPropagator *propagator = GeantPropagator::Instance();
  if ((!fInitialized) || (fScore == kNoScore))
    return;
  // Loop all tracks, check if they are in the right volume and collect the
  // energy deposit and step length
  double theta, eta;
  for (Int_t itr = 0; itr < npart; itr++) {
    if (tracks.fZdirV[itr] == 1)
      eta = 1.E30;
    else {
      theta = acos(tracks.fZdirV[itr]);
      eta = -log(tan(0.5 * theta));
    }
    if (propagator->fNthreads > 1)
      fMHist.lock();
    fHeta->Fill(eta);
    fHpt->Fill(tracks.Pt(itr));
    fHStep->Fill(tracks.fStepV[itr]);
    fStepSize->Fill(eta, tracks.fStepV[itr]);
    if ((tracks.fStatusV[itr] == kKilled) || (tracks.fStatusV[itr] == kExitingSetup) ||
        (tracks.fPathV[itr]->IsOutside()))
      fStepCnt->Fill(eta, tracks.fNstepsV[itr]);
    if (propagator->fNthreads > 1)
      fMHist.unlock();
  }
}

//______________________________________________________________________________
void StdApplication::Digitize(Int_t /* event */) {
  // User method to digitize a full event, which is at this stage fully transported
  //   printf("======= Statistics for event %d:\n", event);
}

//______________________________________________________________________________
void StdApplication::FinishRun() {
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
}
