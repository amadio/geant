#include "GeantVolumeBasket.h"
#include "globals.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "PhysicsProcess.h"

#include "TThread.h"
#include "TArrayI.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TGeoBranchArray.h"

#include <iostream>

using std::min;
using std::max;

const double gTolerance = TGeoShape::Tolerance();

//______________________________________________________________________________
GeantVolumeBasket::GeantVolumeBasket(TGeoVolume *vol, int number) : TObject(), fVolume(vol), fNumber(number) {
  // Constructor
}

//______________________________________________________________________________
GeantVolumeBasket::~GeantVolumeBasket() {
  // Clean up
}

//______________________________________________________________________________
void GeantVolumeBasket::ComputeTransportLength(int ntracks, int *trackin) {
  // Computes snext and safety for an array of tracks. This is the transportation
  // process. Tracks are assumed to be inside gVolume.

  GeantPropagator *gPropagator = GeantPropagator::Instance();

  static int icalls = 0;
  double pdir[3];
  int itr;
  bool isOnBoundary = kFALSE;
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();

  if (!nav) {
    nav = gGeoManager->AddNavigator();
    std::cerr << "[ComputeTransportLength] Added navigator" << std::endl;
  }

  nav->SetOutside(kFALSE);

  for (itr = 0; itr < ntracks; itr++) {
    GeantTrack *track = gPropagator->fTracks[trackin[itr]];
    track->Direction(pdir);
    track->path->UpdateNavigator(nav);
    nav->SetCurrentPoint(&track->xpos);
    nav->SetCurrentDirection(pdir);
    isOnBoundary = track->frombdr;
    double pstep = min<double>(1.E20, track->pstep);
    track->snext = 0;
    if (track->charge) {
      nav->FindNextBoundary(pstep, "", isOnBoundary);
      track->safety = nav->GetSafeDistance();
      track->snext = max<double>(2 * gTolerance, nav->GetStep());
    } else {
      // Propagate to next volume without computing safety for neutrals
      if (gPropagator->fUseDebug && (gPropagator->fDebugTrk == trackin[itr] || gPropagator->fDebugTrk < 0)) {
        Printf("track %d in %s - calling FNBAS", trackin[itr], nav->GetCurrentNode()->GetName());
      }
      nav->FindNextBoundaryAndStep(pstep, kFALSE);
      track->safety = 0.;
      track->snext = nav->GetStep();
      if (gPropagator->fUseDebug && (gPropagator->fDebugTrk == trackin[itr] || gPropagator->fDebugTrk < 0)) {
        Printf("track %d ending in %s after snext=%g", trackin[itr], nav->GetCurrentNode()->GetName(), track->snext);
      }
      if (nav->IsOutside() || track->snext > 1.E19) {
        gPropagator->StopTrack(track);
        continue;
      }
      // Check if it was a short step
      track->nextpath->InitFromNavigator(nav);
      track->frombdr = nav->IsOnBoundary();
      if (nav->IsOnBoundary() && track->snext < 2. * gTolerance) {
        // Make sure track crossed
        track->izero++;
        if (track->izero > 10) {
          gPropagator->StopTrack(track);
          continue;
        }
        nav->FindNextBoundaryAndStep(1.E30, kFALSE);
        track->nextpath->InitFromNavigator(nav);
        track->snext += nav->GetStep();
      }
      if (track->snext > 2. * gTolerance)
        track->izero = 0;
    }

    if (gPropagator->fUseDebug && (gPropagator->fDebugTrk == trackin[itr] || gPropagator->fDebugTrk < 0)) {
      Printf(" CTL(%d):   track %d: %s  snext=%19.15f safe=%19.15f pstep=%f", icalls, trackin[itr], nav->GetPath(),
             track->snext, track->safety, track->pstep);
      track->Print(trackin[itr]);
    }
  }
  icalls++;
}

//______________________________________________________________________________
void GeantVolumeBasket::PropagateTracks(int ntracks, int *trackin, int &nout, int *trackout, int &ntodo, int *tracktodo,
                                        int &ncross, int *trackcross) {
  // Propagate the ntracks with their selected steps. If a boundary is
  // found in the way, the track is stopped. Nout must be initialized from outside.
  //     trackin = array of <ntracks> input tracks
  //     trackout = array of <nout> tracks propagated to physics processes
  //     tracktodo = array of <ntodo> tracks propagated with a safe step or crossing
  //                 inside the same volume. These have to be propagated  again.
  //     trackcross = array of <ncross> tracks that crossed the boundary. For these tracks
  //                 the continuous processes have to be applied after propagation

  GeantPropagator *gPropagator = GeantPropagator::Instance();

  GeantTrack *track;
  double step, snext, safety, c;
  ntodo = 0;

  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();

  if (!nav) {
    nav = gGeoManager->AddNavigator();
    std::cerr << "[PropagateTracks] Added navigator" << std::endl;
  }

  GeantVolumeBasket *basket = 0;
  //   Printf("===== PropagateTracks: ntracks=%d nout=%d", ntracks, nout);
  for (int itr = 0; itr < ntracks; itr++) {
    track = gPropagator->fTracks[trackin[itr]];
    // Skip neutral tracks for the time being (!!!)
    if (!track->IsAlive())
      continue;
    if (!track->charge) {
      track->PropagateStraight(track->snext, trackin[itr]);
      if (track->frombdr)
        trackcross[ncross++] = trackin[itr];
      else
        trackout[nout++] = trackin[itr];
      continue;
    }
    track->nsteps++;
    if (track->nsteps > gPropagator->fMaxSteps) {
      gPropagator->StopTrack(track);
      continue;
    }
    step = track->pstep;
    snext = track->snext;
    safety = track->safety;
    // If proposed step less than safety, just propagate to physics process
    if (step < safety) {
      track->izero = 0;
      track->PropagateInField(step, kFALSE, trackin[itr]);
      gPropagator->fNsafeSteps++; // increment-only, thread safe
      // track transported to physics process
      if (gPropagator->fUseDebug && (gPropagator->fDebugTrk == trackin[itr] || gPropagator->fDebugTrk < 0))
        Printf("   track %d process: %s", trackin[itr], gPropagator->Process(track->process)->GetName());
      trackout[nout++] = trackin[itr]; // <- survives geometry, stopped due to physics
      continue;                        // -> to next track
    }
    // Check if we can propagate to boundary
    c = track->Curvature();
    if (0.25 * c * snext < 1E-6 && snext < 1E-3 && snext < step - 1E-6) {
      // Propagate with snext and check if we crossed
      //   backup track position and momentum
      if (track->izero > 10)
        snext = 1.E-3;
      basket = track->PropagateInField(snext + 10 * gTolerance, kTRUE, trackin[itr]);
      if (snext < 1.E-6)
        track->izero++;
      else
        track->izero = 0;
      gPropagator->fNsnextSteps++;
      if (!basket) {
        // track exiting
        if (nav->IsOutside()) {
          if (gPropagator->fUseDebug && (gPropagator->fDebugTrk == trackin[itr] || gPropagator->fDebugTrk < 0))
            Printf("   track %d exiting geometry", trackin[itr]);
          trackcross[ncross++] = trackin[itr];
          gPropagator->StopTrack(track);
          continue;
        }
        // these tracks do not cross
        //            if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0))
        //            Printf("   track %d propagated with snext=%19.15f", trackin[itr], snext);
        tracktodo[ntodo++] = trackin[itr]; // <- survives partial geometry step
        continue;                          // -> next track
      }
      // These tracks are reaching boundaries
      trackcross[ncross++] = trackin[itr];
      //         if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0))
      //         Printf("   track %d pushed to boundary of %s", trackin[itr], basket->GetName());
      //         basket = track->PropagateStraight(snext, trackin[itr]);
      continue; // -> to next track
    }
    // Track has safety<pstep but next boundary not close enough.
    // We propagate in field with the safety value.
    if (safety < gTolerance) {
      // Track getting away from boundary. Work to be done here
      // In principle we need a safety value for the content of the current volume only
      // This does not behave well on corners...
      // ... so we peek a small value and chech if this crosses, than recompute safety
      safety = 1.E-3;
      track->izero++;
      if (track->izero > 10)
        safety = 0.5 * snext;
      basket = track->PropagateInField(safety, kTRUE, trackin[itr]);
    } else {
      if (track->izero > 10) {
        // Propagate with snext
        basket = track->PropagateInField(snext + 10 * gTolerance, kTRUE, trackin[itr]);
        track->izero = 0;
        gPropagator->fNsnextSteps++;
        if (!basket) {
          // track exiting geometry
          if (nav->IsOutside()) {
            if (gPropagator->fUseDebug && (gPropagator->fDebugTrk == trackin[itr] || gPropagator->fDebugTrk < 0))
              Printf("   track %d exiting geometry", trackin[itr]);
            trackcross[ncross++] = trackin[itr];
            gPropagator->StopTrack(track);
            continue;
          }
          // these tracks do not cross
          //               if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] ||
          //               gPropagator->fDebugTrk<0)) Printf("   track %d propagated with snext=%19.15f", trackin[itr],
          //               snext);
          tracktodo[ntodo++] = trackin[itr]; // <- survives partial geometry step
          continue;                          // -> next track
        }
        // These tracks are reaching boundaries
        trackcross[ncross++] = trackin[itr];
        continue;
      }
      if (safety < 1.E-3)
        track->izero++;
      // Propagate with safety without checking crossing
      basket = track->PropagateInField(safety, kFALSE, trackin[itr]);
    }
    gPropagator->fNsafeSteps++;
    if (!basket) {
      // check if track exiting
      if (nav->IsOutside()) {
        if (gPropagator->fUseDebug && (gPropagator->fDebugTrk == trackin[itr] || gPropagator->fDebugTrk < 0))
          Printf("   track %d exiting geometry", trackin[itr]);
        trackcross[ncross++] = trackin[itr];
        gPropagator->StopTrack(track);
        continue;
      }
      // these tracks do not cross
      //         if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0))
      //         Printf("   track %d propagated with safety=%19.15f", trackin[itr], safety);
      tracktodo[ntodo++] = trackin[itr]; // <- survives partial geometry step
      continue;                          // -> next track
    }
    // These tracks are reaching boundaries
    trackcross[ncross++] = trackin[itr];
  }
  // Recompute snext and safety for todo tracks
  if (ntodo)
    ComputeTransportLength(ntodo, tracktodo);
}

//______________________________________________________________________________
void GeantVolumeBasket::Print(const char *) const {
  // Print info about the basket content.
}

void GeantVolumeBasket::ResetStep(int ntracks, int *array) {
  GeantPropagator *gPropagator = GeantPropagator::Instance();
  // Reset current step for a list of tracks.
  for (int i = 0; i < ntracks; i++) {
    GeantTrack *track = gPropagator->fTracks[array[i]];
    track->step = 0.;
  }
}
