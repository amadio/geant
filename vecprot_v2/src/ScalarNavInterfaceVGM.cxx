#include "ScalarNavInterfaceVGM.h"

#include "backend/Backend.h"
#include "navigation/VNavigator.h"
#ifdef VECCORE_CUDA
#include "navigation/SimpleNavigator.h"
#else
#include "navigation/ABBoxNavigator.h"
#endif
#include "volumes/PlacedVolume.h"
#include "base/Vector3D.h"
#include "base/Transformation3D.h"
#include "base/Global.h"
#include "management/GeoManager.h"

#ifdef CROSSCHECK
#include "TGeoNavigator.h"
#include "TGeoNode.h"
#endif

#ifdef BUG_HUNT
#include "GeantPropagator.h"
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

using namespace VECGEOM_NAMESPACE;

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void ScalarNavInterfaceVGM::NavFindNextBoundaryAndStep(int ntracks, const double *pstep,
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

  typedef Vector3D<Precision> Vector3D_t;
#ifdef VECCORE_CUDA
  SimpleNavigator nav;
#else
  ABBoxNavigator nav;
#endif // VECCORE_CUDA

  for (int itr = 0; itr < ntracks; ++itr) {
    // If the basket is mixed volumes/navigators may be different
    VNavigator const * newnav = instate[itr]->Top()->GetLogicalVolume()->GetNavigator();
    // Check if current safety allows for the proposed step
    if (safe[itr] > pstep[itr]) {
      step[itr] = pstep[itr];
      isonbdr[itr] = false;
      *outstate[itr] = *instate[itr];
      continue;
    }

    step[itr] = newnav->ComputeStepAndSafetyAndPropagatedState(Vector3D_t(x[itr], y[itr], z[itr]),
                               Vector3D_t(dirx[itr], diry[itr], dirz[itr]),
                               Math::Min<double>(1.E20, pstep[itr]), *instate[itr], *outstate[itr] /* the paths */, !isonbdr[itr], safe[itr]);
    step[itr] = Math::Max<double>(2. * gTolerance, step[itr] + 2. * gTolerance);
    safe[itr] = Math::Max<double>(safe[itr], 0);
    // onboundary with respect to new point
    isonbdr[itr] = outstate[itr]->IsOnBoundary();

    //#### To add small step detection and correction - see ScalarNavInterfaceTGeo ####//

#if defined(CROSSCHECK) && !defined(VECCORE_CUDA)
    //************
    // CROSS CHECK USING TGEO
    //************
    TGeoNavigator *rootnav = gGeoManager->GetCurrentNavigator();
    rootnav->ResetState();
    rootnav->SetCurrentPoint(x[itr], y[itr], z[itr]);
    rootnav->SetCurrentDirection(dirx[itr], diry[itr], dirz[itr]);
    TGeoBranchArray *tmp = instate[itr]->ToTGeoBranchArray();
    tmp->UpdateNavigator(rootnav);
    delete tmp;
    rootnav->FindNextBoundaryAndStep(Math::Min<double>(1.E20, pstep[itr]), !isonbdr[itr]);
    double stepcmp = Math::Max<double>(2 * gTolerance, rootnav->GetStep());
    double safecmp = rootnav->GetSafeDistance();
    if (Math::Abs(step[itr] - stepcmp) > 1E-6) {
      Geant::Print("","## PSTEP %lf VECGEOMSTEP %lf ROOTSTEP %lf", pstep[itr], step[itr], stepcmp);
      Geant::Print("","## PSTEP %lf ONBOUND %d VECGEOMSAFETY %lf ROOTSAFETY %lf BRUTEFORCEROOT %lf", pstep[itr], isonbdr[itr],
             safe[itr], rootnav->Safety());

      // check nextpath
      tmp = outstate[itr]->ToTGeoBranchArray();
      tmp->InitFromNavigator(rootnav);
      Geant::Print("","## VECGEOMNEXTNODE %p ROOTNEXTNODE %p", outstate[itr]->GetCurrentNode(), tmp->GetCurrentNode());
      Geant::Print("","## VECGEOMBOUNDARY %d ROOTBOUNDARY %d", outstate[itr]->IsOnBoundary(), rootnav->IsOnBoundary());

      Geant::Print("","INCONSISTENT STEP");
      nav.InspectEnvironmentForPointAndDirection(Vector3D_t(x[itr], y[itr], z[itr]) /*global pos*/,
                                                 Vector3D_t(dirx[itr], diry[itr], dirz[itr]), *instate[itr]);
    }
#endif // CROSSCHECK

#ifdef VERBOSE
    Geant::Print("","navfindbound on %p track %d with pstep %lf yields step %lf and safety %lf\n", this, itr, pstep[itr], step[itr],
           safe[itr]);
#endif // VERBOSE
  }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void ScalarNavInterfaceVGM::NavFindNextBoundaryAndStep(GeantTrack &track) {

  typedef Vector3D<Precision> Vector3D_t;
#ifdef VECCORE_CUDA
  SimpleNavigator nav;
#else
  ABBoxNavigator nav;
#endif // VECCORE_CUDA

  // Retrieve navigator for the track
  VNavigator const * newnav = track.GetVolume()->GetNavigator();
  // Check if current safety allows for the proposed step
  if (track.fSafety > track.fPstep) {
    track.fSnext = track.fPstep;
    track.fBoundary = false;
    track.UpdateSameNextPath();
    return;
  }

  track.fSnext = newnav->ComputeStepAndSafetyAndPropagatedState(Vector3D_t(track.fXpos, track.fYpos, track.fZpos),
                             Vector3D_t(track.fXdir, track.fYdir, track.fZdir),
                             Math::Min<double>(1.E20, track.fPstep),
                             *track.Path(), *track.NextPath(), !track.fBoundary, track.fSafety);
  track.fSnext = Math::Max<double>(2. * gTolerance, track.fSnext + 2. * gTolerance);
  track.fSafety = Math::Max<double>(track.fSafety, 0);
  // onboundary with respect to new point
  track.fBoundary = track.NextPath()->IsOnBoundary();

  //#### To add small step detection and correction - see ScalarNavInterfaceTGeo ####//

#ifdef CROSSCHECK
  //************
  // CROSS CHECK USING TGEO
  //************
  TGeoNavigator *rootnav = gGeoManager->GetCurrentNavigator();
  rootnav->ResetState();
  rootnav->SetCurrentPoint(track.fXpos, track.fYpos, track.fZpos);
  rootnav->SetCurrentDirection(track.fXdir, track.fYdir, track.fZdir);
  TGeoBranchArray *tmp = track.fPath->ToTGeoBranchArray();
  tmp->UpdateNavigator(rootnav);
  delete tmp;
  rootnav->FindNextBoundaryAndStep(Math::Min<double>(1.E20, track.fPstep), !track.fBoundary);
  double stepcmp = Math::Max<double>(2 * gTolerance, rootnav->GetStep());
  double safecmp = rootnav->GetSafeDistance();
  if (Math::Abs(track.fStep - stepcmp) > 1E-6) {
    Geant::Print("","## PSTEP %lf VECGEOMSTEP %lf ROOTSTEP %lf", track.fPstep, track.fStep, stepcmp);
    Geant::Print("","## PSTEP %lf ONBOUND %d VECGEOMSAFETY %lf ROOTSAFETY %lf BRUTEFORCEROOT %lf", track.fPstep, track.fBoundary,
           track.fSafety, rootnav->Safety());

    // check nextpath
    tmp = track.fNextpath->ToTGeoBranchArray();
    tmp->InitFromNavigator(rootnav);
    Geant::Print("","## VECGEOMNEXTNODE %p ROOTNEXTNODE %p", track.fNextpath->GetCurrentNode(), tmp->GetCurrentNode());
    Geant::Print("","## VECGEOMBOUNDARY %d ROOTBOUNDARY %d", track.fNextpath->IsOnBoundary(), rootnav->IsOnBoundary());

    Geant::Print("","INCONSISTENT STEP");
    nav.InspectEnvironmentForPointAndDirection(Vector3D_t(track.fXpos, track.fYpos, track.fZpos) /*global pos*/,
                                               Vector3D_t(track.fXdir, track.fYdir, track.fZdir),
                                               *track.fPath);
  }
#endif // CROSSCHECK

#ifdef VERBOSE
  Geant::Print("","navfindbound on %p with pstep %lf yields step %lf and safety %lf\n",
               this, track.fPtrack.fStep, track.fPstep, track.fSafety);
#endif // VERBOSE
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void ScalarNavInterfaceVGM::NavFindNextBoundary(GeantTrack &track) {
  // back-up the pre-step point boundary flag
  track.fIsOnBoundaryPreStp = track.fBoundary;
  // Find distance to next boundary, within proposed step.
  typedef Vector3D<Precision> Vector3D_t;
  // Retrieve navigator for the track
  VNavigator const * newnav = track.GetVolume()->GetNavigator();
  // Check if current safety allows for the proposed step
  if (track.fSafety > track.fPstep) {
    track.fSnext = track.fPstep;
    track.fBoundary = false;
//    *track.fNextpath = *track.fPath;
    return;
  }
  track.fSnext = newnav->ComputeStepAndSafety(Vector3D_t(track.fXpos, track.fYpos, track.fZpos),
                             Vector3D_t(track.fXdir, track.fYdir, track.fZdir),
                             Math::Min<double>(1.E20, track.fPstep),
                             *track.Path(), !track.fBoundary, track.fSafety);
  track.fBoundary = track.fSnext < track.fPstep;
  track.fSnext = Math::Max<double>(2. * gTolerance, track.fSnext + 2. * gTolerance);
  track.fSafety = Math::Max<double>(track.fSafety, 0);
}

//______________________________________________________________________________
void ScalarNavInterfaceVGM::NavIsSameLocation(int ntracks,
       const double *x, const double *y, const double *z,
       const double */*dirx*/, const double */*diry*/, const double */*dirz*/,
       const VolumePath_t **start, VolumePath_t **end, bool *same, VolumePath_t *tmpstate) {
//
// Checks if the navigation states corresponding to positions (x,y,z) are the
// same as the ones pointed by start. Update new states in end.
// Input:  ntracks - number of tracks to be checked
//         x,y,z   - arrays of positions
//         start   - starting navigation paths to compare with
// Output: end     - navigation paths corresponding to given positions
//         same    - flags showing if the end and start positions are matching


  //#### NOT USING YET THE NEW NAVIGATORS ####//

  typedef Vector3D<Precision> Vector3D_t;

  SimpleNavigator nav;
  for (int itr = 0; itr < ntracks; ++itr) {

// cross check with answer from ROOT

    // TODO: not using the direction yet here !!
    bool samepath = nav.HasSamePath(
      Vector3D_t(x[itr], y[itr], z[itr]), *start[itr], *tmpstate);
    if (!samepath) tmpstate->CopyTo(end[itr]);

#if defined(CROSSCHECK) && !defined(VECCORE_CUDA)
    TGeoBranchArray *sb = start[itr]->ToTGeoBranchArray();
    TGeoBranchArray *eb = end[itr]->ToTGeoBranchArray();
    TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
    nav->ResetState();
    nav->SetLastSafetyForPoint(0, 0, 0, 0);
    nav->SetCurrentPoint(x[itr], y[itr], z[itr]);
    sb->UpdateNavigator(nav);
    bool rootsame = nav->IsSameLocation(x[itr], y[itr], z[itr], true);
    if (rootsame != samepath) {
      Geant::Print("","INCONSISTENT ANSWER ROOT(%d) VECGEOM(%d)", rootsame, samepath);
      std::cout << Vector3D_t(x[itr], y[itr], z[itr]) << "\n";
      Geant::Print("","old state");
      sb->Print();
      nav->ResetState();
      nav->SetLastSafetyForPoint(0, 0, 0, 0);
      nav->SetCurrentPoint(x[itr], y[itr], z[itr]);
//      nav->SetCurrentDirection(xdir[itr], ydir[itr], zdir[itr]);
      sb->UpdateNavigator(nav);
      nav->InspectState();
      bool rootsame = nav->IsSameLocation(x[itr], y[itr], z[itr], true);
      nav->InspectState();
      eb->InitFromNavigator(nav);
      Geant::Print("","new state");
      eb->Print();
      Geant::Print("","VERSUS VECGEOM OLD AND NEW");
      start[itr]->printVolumePath();
      end[itr]->printVolumePath();
    }
    delete sb;
    delete eb;
#endif // CROSSCHECK && !VECCORE_CUDA
    same[itr] = samepath;
  }
}

//______________________________________________________________________________
void ScalarNavInterfaceVGM::NavIsSameLocation(GeantTrack &track, bool &same, VolumePath_t *tmpstate) {

  //#### NOT USING YET THE NEW NAVIGATORS ####//

  typedef Vector3D<Precision> Vector3D_t;

  SimpleNavigator nav;
// cross check with answer from ROOT
#ifdef CROSSCHECK
  TGeoBranchArray *sb = track.fPath->ToTGeoBranchArray();
  TGeoBranchArray *eb = track.fNextpath->ToTGeoBranchArray();
#endif

  // TODO: not using the direction yet here !!
  bool samepath = nav.HasSamePath(
    Vector3D_t(track.fXpos, track.fYpos, track.fZpos), *track.Path(), *tmpstate);
  if (!samepath) tmpstate->CopyTo(track.NextPath());

#ifdef CROSSCHECK
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  nav->ResetState();
  nav->SetLastSafetyForPoint(0, 0, 0, 0);
  nav->SetCurrentPoint(track.fXpos, track.fYpos, track.fZpos);
  sb->UpdateNavigator(nav);
  bool rootsame = nav->IsSameLocation(track.fXpos, track.fYpos, track.fZpos, true);
  if (rootsame != samepath) {
    Geant::Print("","INCONSISTENT ANSWER ROOT(%d) VECGEOM(%d)", rootsame, samepath);
    std::cout << Vector3D_t(track.fXpos, track.fYpos, track.fZpos) << "\n";
    Geant::Print("","old state");
    sb->Print();
    nav->ResetState();
    nav->SetLastSafetyForPoint(0, 0, 0, 0);
    nav->SetCurrentPoint(track.fXpos, track.fYpos, track.fZpos);
//      nav->SetCurrentDirection(xdir, ydir, zdir);
    sb->UpdateNavigator(nav);
    nav->InspectState();
    bool rootsame = nav->IsSameLocation(track.fXpos, track.fYpos, track.fZpos, true);
    nav->InspectState();
    eb->InitFromNavigator(nav);
    Geant::Print("","new state");
    eb->Print();
    Geant::Print("","VERSUS VECGEOM OLD AND NEW");
    track.fPath->printVolumePath();
    track.fNextpath->printVolumePath();
  }
  delete sb;
  delete eb;
#endif // CROSSCHECK
  same = samepath;
}

//______________________________________________________________________________
void ScalarNavInterfaceVGM::DisplaceTrack(GeantTrack &track, const double dir[3], double disp, double mindisp)
{
  // before calling displacement it must be checked that the point is not on boundary! It is done in the caller.
  typedef Vector3D<Precision> Vector3D_t;
  constexpr double reduceFactor = 0.99;
  double realDisp = disp;
  bool   isDisp   = true;
  // check if displacement is beyond the current (post-step point i.e. updated) safety
  double postSafety = track.fSafety*reduceFactor; // make sure that we never reach the boundary in displacement
  if (disp>postSafety) {
    // Retrieve navigator for the track
    VNavigator const * newnav = track.GetVolume()->GetNavigator();
    // NOTE: we should have a navigator method to compute only the safety and the following call should be change to that
    /*double step =*/ newnav->ComputeStepAndSafety(Vector3D_t(track.fXpos, track.fYpos, track.fZpos),
                                               Vector3D_t(dir[0], dir[1], dir[2]),
                                               disp, *track.Path(), true, postSafety);
    postSafety *= reduceFactor; // make sure that we never reach the boundary in displacement
    if (postSafety>mindisp) {
      realDisp = std::min(postSafety,disp);
    } else {
      isDisp = false;
    }
  }
  // apply displacement if any
  if (isDisp) {
    track.fXpos   += realDisp * dir[0];
    track.fYpos   += realDisp * dir[1];
    track.fZpos   += realDisp * dir[2];
    track.fSafety -= realDisp;
  }
}

/*
void ScalarNavInterfaceVGM::NavFindNextBoundaryMSC(GeantTrack &track, double dist) {
//track.fIsEverythingWasDone = track.fBoundary;
// Find distance to next boundary, within proposed step.
   typedef Vector3D<Precision> Vector3D_t;
  // Retrieve navigator for the track
  VNavigator const * newnav = track.fPath->Top()->GetLogicalVolume()->GetNavigator();
  // Check if current safety allows for the proposed step
  // doesn't work well: I always need the precise safety and dist. to b.
  //if (track.fSafety > dist) {
  //  track.fSnext = dist;
  //    // track.fBoundary = false;
  //    // *track.fNextpath = *track.fPath;
  //  return;
  //}

  track.fSnext = newnav->ComputeStepAndSafety(Vector3D_t(track.fXpos, track.fYpos, track.fZpos),
                             Vector3D_t(track.fXdir, track.fYdir, track.fZdir),
                             1.e+20, //Math::Min<double>(1.E20, dist),
                             *track.fPath, true, track.fSafety);
//  track.fBoundary = track.fPath->IsOnBoundary();
//  if (track.fSafety==0.) track.fBoundary=true;
}
*/
} // GEANT_IMPL_NAMESPACE
} // Geant
