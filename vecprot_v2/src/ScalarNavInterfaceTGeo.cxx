#include "ScalarNavInterfaceTGeo.h"

#include "TGeoNavigator.h"
#include "TGeoManager.h"
#include "TGeoNode.h"

#ifdef BUG_HUNT
#include "GeantPropagator.h"
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void ScalarNavInterfaceTGeo::NavFindNextBoundaryAndStep(int ntracks, const double *pstep, 
         const double *x, const double *y, const double *z,
         const double *dirx, const double *diry, const double *dirz, 
         const VolumePath_t **instate, VolumePath_t **outstate, 
         double *step, double *safe, bool *isonbdr) {
// Find the next boundary and state after propagating to the boundary. 
// Input:  ntracks - number of tracks
//         pstep - proposed step
//         x, y, z, dirx, diry, dirz - particle position and direction
//         instate - input particle navigation state
//         safe - estimated safety value for the input point
//         isonbdr - starting point is on a boundary
// Output: outstate - navigation state after propagation. If boundary further
//           than proposed step, outstate has to match instate
//         step - propagation step for which the state is sampled
//         safety - calculated safety value for the input point
//         isonbdr - propagated point is on a boundary

  const double epserr = 1.E-3; // push value in case of repeated geom error
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  int ismall;
  double snext;
  TGeoNode *nextnode, *lastnode;
  double pt[3];
#ifdef BUG_HUNT
  int index = (int)(x - trk->fXposV);
#endif // BUG_HUNT
  for (int itr = 0; itr < ntracks; itr++) {
#ifdef BUG_HUNT
    index += itr;
#endif // BUG_HUNT
    ismall = 0;
    step[itr] = 0;
    // Check if current safety allows for the proposed step
    if (safe[itr] > pstep[itr]) {
      step[itr] = pstep[itr];
      *outstate[itr] = *instate[itr];
      isonbdr[itr] = false;
      continue;
    }
    // Reset navigation state flags and safety to start fresh
    nav->ResetState();
    // Setup start state
    nav->SetCurrentPoint(x[itr], y[itr], z[itr]);
    nav->SetCurrentDirection(dirx[itr], diry[itr], dirz[itr]);
    instate[itr]->UpdateNavigator(nav);
    nextnode = nav->GetCurrentNode();
    while (nextnode) {
      lastnode = nextnode;
      // Compute distance to next boundary and propagate internally
      nextnode = nav->FindNextBoundaryAndStep(Math::Min<double>(1.E20, pstep[itr]), !isonbdr[itr]);
      snext = nav->GetStep();
      // Adjust step to be non-negative and cross the boundary
      step[itr] = Math::Max<double>(2 * gTolerance, snext + 2 * gTolerance);
      // Check for repeated small steps starting from boundary
      if (isonbdr[itr] && (snext < 1.E-8) && (pstep[itr] > 1.E-8)) {
        ismall++;
        if ((ismall < 3) && (nextnode != lastnode)) {
          // Make sure we don't have a thin layer artefact so repeat search
          nextnode = nav->FindNextBoundaryAndStep(Math::Min<double>(1.E20, pstep[itr] - snext), !isonbdr[itr]);
          snext = nav->GetStep();
          step[itr] += snext;
          // If step still small, repeat
          if (snext < 1.E-8)
            continue;
          // We managed to cross with macroscopic step: reset error counter and exit loop
          ismall = 0;
          break;
        } else {
          if (ismall > 3) {
            // Mark track to be killed
            step[itr] = -1;
            break;
          }
          // The block below can only happen if crossing into the same node on different geometry
          // branch with small step. Try to relocate the next point by making an epserr push
          memcpy(pt, nav->GetCurrentPoint(), 3 * sizeof(double));
          const double *dir = gGeoManager->GetCurrentDirection();
          for (int j = 0; j < 3; j++)
            pt[j] += epserr * dir[j];
          step[itr] += epserr;
          nav->CdTop();
          nextnode = nav->FindNode(pt[0], pt[1], pt[2]);
          if (nav->IsOutside())
            break;
          continue;
        }
      }
      // All OK here, reset error counter and exit loop
      ismall = 0;
      break;
    }
    // Update safety, boundary flag and next path
    safe[itr] = isonbdr[itr] ? 0. : nav->GetSafeDistance();
    isonbdr[itr] = nav->IsOnBoundary();
    outstate[itr]->InitFromNavigator(nav);
#ifdef VERBOSE
    double bruteforces = nav->Safety();
    Geant::Print("","##TGEOM  ## TRACK %d BOUND %d PSTEP %lg STEP %lg SAFETY %lg BRUTEFORCES %lg TOBOUND %d", itr, isonbdr[itr],
           pstep[itr], step[itr], safe[itr], bruteforces, nav->IsOnBoundary());
// assert( safe[itr]<=bruteforces );
#endif // VERBOSE
  }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void ScalarNavInterfaceTGeo::NavFindNextBoundary(GeantTrack &track) {
  // Check if current safety allows for the proposed step
  if (track.fSafety > track.fPstep) {
    track.fSnext = track.fPstep;
    track.UpdateSameNextPath();
    track.fBoundary = false;
    return;
  }
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  // Reset navigation state flags and safety to start fresh
  nav->ResetState();
  // Setup start state
  nav->SetCurrentPoint(track.fXpos, track.fYpos, track.fZpos);
  nav->SetCurrentDirection(track.fXdir, track.fYdir, track.fZdir);
  track.Path()->UpdateNavigator(nav);
  nav->FindNextBoundary(Math::Min<double>(1.E20, track.fPstep), "", track.fBoundary);
  track.fSnext = nav->GetStep();
  track.fSafety = track.fBoundary ? 0. : nav->GetSafeDistance();
  track.fBoundary = track.fSnext < track.fPstep;
  track.fSnext = Math::Max<double>(2 * gTolerance, track.fSnext + 2 * gTolerance);
}
  
//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void ScalarNavInterfaceTGeo::NavFindNextBoundaryAndStep(GeantTrack &track) {

  const double epserr = 1.E-3; // push value in case of repeated geom error
  int ismall;
  double snext;
  TGeoNode *nextnode, *lastnode;
  double pt[3];
  ismall = 0;
  track.fStep = 0;
  // Check if current safety allows for the proposed step
  if (track.fSafety > track.fPstep) {
    track.fStep = track.fPstep;
    track.UpdateSameNextPath();
    track.fBoundary = false;
    return;
  }
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  // Reset navigation state flags and safety to start fresh
  nav->ResetState();
  // Setup start state
  nav->SetCurrentPoint(track.fXpos, track.fYpos, track.fZpos);
  nav->SetCurrentDirection(track.fXdir, track.fYdir, track.fZdir);
  track.Path()->UpdateNavigator(nav);
  nextnode = nav->GetCurrentNode();
  while (nextnode) {
    lastnode = nextnode;
    // Compute distance to next boundary and propagate internally
    nextnode = nav->FindNextBoundaryAndStep(Math::Min<double>(1.E20, track.fPstep), !track.fBoundary);
    snext = nav->GetStep();
    // Adjust step to be non-negative and cross the boundary
    track.fStep = Math::Max<double>(2 * gTolerance, snext + 2 * gTolerance);
    // Check for repeated small steps starting from boundary
    if (track.fBoundary && (snext < 1.E-8) && (track.fPstep > 1.E-8)) {
      ismall++;
      if ((ismall < 3) && (nextnode != lastnode)) {
        // Make sure we don't have a thin layer artefact so repeat search
        nextnode = nav->FindNextBoundaryAndStep(Math::Min<double>(1.E20, track.fPstep - snext), !track.fBoundary);
        snext = nav->GetStep();
        track.fStep += snext;
        // If step still small, repeat
        if (snext < 1.E-8)
          continue;
        // We managed to cross with macroscopic step: reset error counter and exit loop
        ismall = 0;
        break;
      } else {
        if (ismall > 3) {
          // Mark track to be killed
          track.fStep = -1;
          break;
        }
        // The block below can only happen if crossing into the same node on different geometry
        // branch with small step. Try to relocate the next point by making an epserr push
        memcpy(pt, nav->GetCurrentPoint(), 3 * sizeof(double));
        const double *dir = gGeoManager->GetCurrentDirection();
        for (int j = 0; j < 3; j++)
          pt[j] += epserr * dir[j];
        track.fStep += epserr;
        nav->CdTop();
        nextnode = nav->FindNode(pt[0], pt[1], pt[2]);
        if (nav->IsOutside())
          break;
        continue;
      }
    }
    // All OK here, reset error counter and exit loop
    ismall = 0;
    break;
  }
  // Update safety, boundary flag and next path
  track.fSafety = track.fBoundary ? 0. : nav->GetSafeDistance();
  track.fBoundary = nav->IsOnBoundary();
  track.NextPath()->InitFromNavigator(nav);
#ifdef VERBOSE
  double bruteforces = nav->Safety();
  Geant::Print("","##TGEOM  BOUND %d PSTEP %lg STEP %lg SAFETY %lg BRUTEFORCES %lg TOBOUND %d", track.fBoundary,
         track.fPstep, track.fStep, track.fSafety, bruteforces, nav->IsOnBoundary());
// assert( track.fSafety<=bruteforces );
#endif // VERBOSE
}

//______________________________________________________________________________
void ScalarNavInterfaceTGeo::NavIsSameLocation(int ntracks,
       const double *x, const double *y, const double *z,
       const double *dirx, const double *diry, const double *dirz,
       const VolumePath_t **start, VolumePath_t **end, bool *same) {
// 
// Checks if the navigation states corresponding to positions (x,y,z) are the
// same as the ones pointed by start. Update new states in end.
// Input:  ntracks - number of tracks to be checked
//         x,y,z   - arrays of positions
//         start   - starting navigation paths to compare with
// Output: end     - navigation paths corresponding to given positions
//         same    - flags showing if the end and start positions are matching

  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  for (int itr = 0; itr < ntracks; ++itr) {
    nav->ResetState();
    nav->SetLastSafetyForPoint(0, 0, 0, 0);
    nav->SetCurrentPoint(x[itr], y[itr], z[itr]);
    nav->SetCurrentDirection(dirx[itr], diry[itr], dirz[itr]);
    start[itr]->UpdateNavigator(nav);
    if (!nav->IsSameLocation(x[itr], y[itr], z[itr], true)) {
      end[itr]->InitFromNavigator(nav);
#ifdef BUG_HUNT
      GeantPropagator *prop = GeantPropagator::Instance();
      BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep, "NavIsSameLoc:CROSSED", itr);
#endif
      same[itr] = false;      
    } else {
    // Track not crossing
      same[itr] = true;
    }  
#ifdef BUG_HUNT
    BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep, "NavIsSameLoc:SAME", itr);
#endif
  }
}

//______________________________________________________________________________
void ScalarNavInterfaceTGeo::NavIsSameLocation(GeantTrack &track, bool &same) {
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  nav->ResetState();
  nav->SetLastSafetyForPoint(0, 0, 0, 0);
  nav->SetCurrentPoint(track.fXpos, track.fYpos, track.fZpos);
  nav->SetCurrentDirection(track.fXdir, track.fYdir, track.fZdir);
  track.Path()->UpdateNavigator(nav);
  if (!nav->IsSameLocation(track.fXpos, track.fYpos, track.fZpos, true)) {
    track.NextPath()->InitFromNavigator(nav);
    same = false;     
  } else {
  // Track not crossing
    same = true;
  }
}

} // GEANT_IMPL_NAMESPACE
} // Geant
