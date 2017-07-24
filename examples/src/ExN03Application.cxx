#include "ExN03Application.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "management/GeoManager.h"
using vecgeom::GeoManager;
#endif
#include "GeantEvent.h"
#include "GeantFactoryStore.h"
#include "GeantTrackVec.h"
#include "GeantRunManager.h"
#include "GeantTaskData.h"
#include "globals.h"
#ifdef USE_ROOT
#include "TGeoNode.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#endif
#include <cassert>

#include "Geant/Error.h"

using std::min;
using std::max;

//______________________________________________________________________________
ExN03Application::ExN03Application(GeantRunManager *runmgr)
  : GeantVApplication(runmgr), fInitialized(false), fIdGap(0), fIdAbs(0), fFactory(0) {
  // Ctor..
  GeantFactoryStore *store = GeantFactoryStore::Instance();
  fFactory = store->GetFactory<MyHit>(16, runmgr->GetNthreadsTotal());
  memset(fEdepGap, 0, kNlayers * kMaxThreads * sizeof(float));
  memset(fLengthGap, 0, kNlayers * kMaxThreads * sizeof(float));
  memset(fEdepAbs, 0, kNlayers * kMaxThreads * sizeof(float));
  memset(fLengthAbs, 0, kNlayers * kMaxThreads * sizeof(float));
}

//______________________________________________________________________________
void ExN03Application::AttachUserData(GeantTaskData *td)
{
  // Create task data for handling digits. Provide number of event slots.
  ExN03ScoringData *data = new ExN03ScoringData(fRunMgr->GetConfig()->fNbuff, kNlayers);
  fDigitsHandle->AttachUserData(data, td);
  printf("Attached user data %p for tid=%d\n", data, td->fTid);
}

//______________________________________________________________________________
bool ExN03Application::Initialize() {
  // Initialize application. Geometry must be loaded.
  if (fInitialized)
    return true;
#ifndef USE_VECGEOM_NAVIGATOR
  if (!gGeoManager) {
#else
  if (!GeoManager::Instance().GetWorld()) {
#endif
    Geant::Error("ExN03Application::Initialize", "Geometry not loaded");
    return false;
  }
#ifndef USE_VECGEOM_NAVIGATOR
  Volume_t *lvGap = gGeoManager->GetVolume("liquidArgon");
  Volume_t *lvAbs = gGeoManager->GetVolume("Lead");
#else
  Volume_t *lvGap = GeoManager::Instance().FindLogicalVolume("liquidArgon");
  Volume_t *lvAbs = GeoManager::Instance().FindLogicalVolume("Lead");
#endif

  if (!lvGap || !lvAbs) {
    Geant::Error("ExN03Application::Initialize", "Logical volumes for gap and absorber not found - do you use the right geometry");
    return false;
  }
#ifndef USE_VECGEOM_NAVIGATOR
  fIdGap = lvGap->GetNumber();
  fIdAbs = lvAbs->GetNumber();
#else
  fIdGap = lvGap->id();
  fIdAbs = lvAbs->id();
#endif
  // Register user data and get a handle for it with the task data manager
  fDigitsHandle = fRunMgr->GetTDManager()->RegisterUserData<ExN03ScoringData>("ExN03digits");
  fInitialized = true;
  return true;
}

//______________________________________________________________________________
void ExN03Application::StepManager(int npart, const GeantTrack_v &tracks, GeantTaskData *td) {
  // Application stepping manager. The thread id has to be used to manage storage
  // of hits independently per thread.
  if (!fInitialized)
    return; // FOR NOW
  // Loop all tracks, check if they are in the right volume and collect the
  // energy deposit and step length
  int tid = td->fTid;
  Node_t const *current;
  int idvol = -1;
  int idnode = -1;
  int ilev = -1;
  for (int i = 0; i < npart; i++) {
//      printf("%d=>\n", i);
//      tracks.PrintTrack(i);
#ifndef USE_VECGEOM_NAVIGATOR
    ilev = tracks.fPathV[i]->GetLevel();
#else
    ilev = tracks.fPathV[i]->GetCurrentLevel() - 1;
#endif
    if (ilev < 1)
      continue;
#ifndef USE_VECGEOM_NAVIGATOR
    current = tracks.fPathV[i]->GetCurrentNode();
#else
    current = tracks.fPathV[i]->Top();
#endif
    if (!current)
      continue;
#ifndef USE_VECGEOM_NAVIGATOR
    idnode = tracks.fPathV[i]->GetNode(ilev - 1)->GetNumber();
    idvol = current->GetVolume()->GetNumber();
#else
    idnode = tracks.fPathV[i]->At(ilev - 1)->id();
    idvol = current->GetLogicalVolume()->id();
#endif
    if (idvol == fIdGap) {
      //         tracks.PrintTrack(i);
      fEdepGap[idnode][tid] += tracks.fEdepV[i];
      fLengthGap[idnode][tid] += tracks.fStepV[i];
    } else if (idvol == fIdAbs) {
      //         tracks.PrintTrack(i);
      fEdepAbs[idnode][tid] += tracks.fEdepV[i];
      fLengthAbs[idnode][tid] += tracks.fStepV[i];
    }
  }
  if (fRunMgr->GetConfig()->fFillTree) {
    MyHit *hit;
    //    int nhits = 0;
    for (int i = 0; i < npart; i++) {
      // Deposit hits in scintillator
      if (idvol==fIdGap && tracks.fEdepV[i]>0.00002)
	{
	  hit = fFactory->NextFree(tracks.fEvslotV[i], tid);

	  hit->fX = tracks.fXposV[i];
	  hit->fY = tracks.fYposV[i];
	  hit->fZ = tracks.fZposV[i];
	  hit->fEdep = tracks.fEdepV[i];
          hit->fTime = tracks.fTimeV[i];
          hit->fEvent = tracks.fEventV[i];
          hit->fTrack = tracks.fParticleV[i];
	  hit->fVolId = idvol;
	  hit->fDetId = idnode;

	  //      if (track->path && track->path->GetCurrentNode()) {
	  //         hit->fVolId = track->path->GetCurrentNode()->GetVolume()->GetNumber();
	  //         hit->fDetId = track->path->GetCurrentNode()->GetNumber();
	  //      }
	  //	  nhits++;
	}
    }
  }

  //  Printf("Thread %d produced %d hits", tid, nhits);
  //  Printf("StepManager: size of queue %zu", fFactory->fOutputs.size());

  return;
}

//______________________________________________________________________________
void ExN03Application::SteppingActions(GeantTrack &track, GeantTaskData *td)
{
  if (!fInitialized) return;
  // Loop all tracks, check if they are in the right volume and collect the
  // energy deposit and step length
  Node_t const *current;
  int idvol = -1;
  int idnode = -1;
  int ilev = -1;
  ilev = track.fPath->GetLevel();
  if (ilev < 1) return;
#ifndef USE_VECGEOM_NAVIGATOR
  current = track.fPath->GetCurrentNode();
#else
  current = track.fPath->Top();
#endif
  if (!current) return;
#ifndef USE_VECGEOM_NAVIGATOR
    idnode = track.fPath->GetNode(ilev - 1)->GetNumber();
    idvol = current->GetVolume()->GetNumber();
#else
    idnode = track.fPath->At(ilev - 1)->id();
    idvol = current->GetLogicalVolume()->id();
#endif
  ExN03LayerDigit &digit = (*fDigitsHandle)(td).GetDigits(track.fEvslot).GetDigit(idnode);
  if (idvol == fIdGap)
    digit.ScoreInGap(track.fEdep, track.fStep);
  else if (idvol == fIdAbs)
    digit.ScoreInAbs(track.fEdep, track.fStep);
}

//______________________________________________________________________________
void ExN03Application::Digitize(GeantEvent *event) {
  // User method to digitize a full event, which is at this stage fully transported
  //   printf("======= Statistics for event %d:\n", event);
  // Merge the digits for the event
  ExN03ScoringData *digits = fRunMgr->GetTDManager()->MergeUserData(event->GetSlot(), *fDigitsHandle);
  if (digits) {
    printf("=== Merged digits for event %d\n", event->GetEvent());
    //digits->PrintDigits(event);
    digits->Clear(event->GetSlot());
  }
  return;
  printf("Energy deposit [MeV/primary] and cumulated track length [cm/primary] per layer");
  printf("================================================================================");
  double nprim = (double)event->GetNtracks(); // fRunMgr->GetNprimaries()
  for (int i = 0; i < kNlayers; ++i) {
    for (int tid = 1; tid < kMaxThreads; ++tid) {
      fEdepGap[i][0] += fEdepGap[i][tid];
      fLengthGap[i][0] += fLengthGap[i][tid];
      fEdepAbs[i][0] += fEdepAbs[i][tid];
      fLengthAbs[i][0] += fLengthAbs[i][tid];
    }
    printf("Layer %d: Egap=%f   Lgap=%f   Eabs=%f   Labs=%f", i + 3, fEdepGap[i][0] * 1000. / nprim,
           fLengthGap[i][0] / nprim, fEdepAbs[i][0] * 1000. / nprim, fLengthAbs[i][0] / nprim);
  }
  printf("================================================================================");
  //   TCanvas *c1 = new TCanvas("Edep", "Energy deposition for ExN03", 700, 800);
  #ifdef USE_ROOT
  TCanvas *c1 = (TCanvas *)gROOT->GetListOfCanvases()->FindObject("capp");
  if (!c1)
    return;
  c1->Divide(1, 2);
  TVirtualPad *pad = c1->cd(1);
  pad->SetGridx();
  pad->SetGridy();
  pad->SetLogy();
  TH1F *histeg = new TH1F("Edep_gap", "Primary track energy deposition per layer", kNlayers, 0.5, kNlayers+0.5);
  histeg->SetMarkerColor(kRed);
  histeg->SetMarkerStyle(2);
  histeg->SetStats(false);
  TH1F *histea = new TH1F("Edep_abs", "Primary track energy deposition per layer in absorber", kNlayers, 0.5, kNlayers+0.5);
  histea->SetMarkerColor(kBlue);
  histea->SetMarkerStyle(4);
  histea->SetStats(false);
  for (int i = 0; i < kNlayers; i++) {
    histeg->SetBinContent(i+1, fEdepGap[i][0] * 1000. / nprim);
    histea->SetBinContent(i+1, fEdepAbs[i][0] * 1000. / nprim);
  }
  double minval =
      min<double>(histeg->GetBinContent(histeg->GetMinimumBin()), histea->GetBinContent(histea->GetMinimumBin()));
  minval = max<double>(minval, 1.E-5);
  double maxval =
      max<double>(histeg->GetBinContent(histeg->GetMaximumBin()), histea->GetBinContent(histea->GetMaximumBin()));
  histeg->GetXaxis()->SetTitle("Layer");
  histeg->GetYaxis()->SetTitle("Edep per layer [MeV]");
  histeg->GetYaxis()->SetRangeUser(minval - 0.1 * minval, maxval + 0.1 * maxval);
  histeg->Draw("P");
  histea->Draw("SAMEP");
  //   TCanvas *c2 = new TCanvas("Length", "Length in layers for ExN03", 700, 800);
  pad = c1->cd(2);
  pad->SetGridx();
  pad->SetGridy();
  pad->SetLogy();
  TH1F *histlg = new TH1F("Len_gap", "Length per layer normalized per primary", kNlayers, 0.5, kNlayers+0.5);
  histlg->SetMarkerColor(kRed);
  histlg->SetMarkerStyle(2);
  histlg->SetStats(false);
  TH1F *histla = new TH1F("Len_abs", "Length per layer normalized per primary", kNlayers, 0.5, kNlayers+0.5);
  histla->SetMarkerColor(kBlue);
  histla->SetMarkerStyle(4);
  histla->SetStats(false);
  for (int i = 0; i < 10; i++) {
    histlg->SetBinContent(i+1, fLengthGap[i][0] / nprim);
    histla->SetBinContent(i+1, fLengthAbs[i][0] / nprim);
  }
  histlg->GetXaxis()->SetTitle("Layer");
  histlg->GetYaxis()->SetTitle("Length per layer");
  minval = min<double>(histlg->GetBinContent(histlg->GetMinimumBin()), histla->GetBinContent(histla->GetMinimumBin()));
  minval = max<double>(minval, 1.E-5);
  maxval = max<double>(histlg->GetBinContent(histlg->GetMaximumBin()), histla->GetBinContent(histla->GetMaximumBin()));
  histlg->GetYaxis()->SetRangeUser(minval - 0.1 * minval, maxval + 0.1 * maxval);
  histlg->Draw("P");
  histla->Draw("SAMEP");
  #endif
}
