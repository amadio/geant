#include "Geant/ScalarNavInterfaceVGM.h"

#include "navigation/VNavigator.h"
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
#include "Geant/Propagator.h"
#endif

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

using namespace VECGEOM_NAMESPACE;

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void ScalarNavInterfaceVGM::NavFindNextBoundaryAndStep(int ntracks, const double *pstep, const double *x,
                                                       const double *y, const double *z, const double *dirx,
                                                       const double *diry, const double *dirz,
                                                       const VolumePath_t **instate, VolumePath_t **outstate,
                                                       double *step, double *safe, bool *isonbdr)
{
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
  constexpr double gTolerance = 1.e-9;
  for (int itr = 0; itr < ntracks; ++itr) {
    // If the basket is mixed volumes/navigators may be different
    VNavigator const *nav = instate[itr]->Top()->GetLogicalVolume()->GetNavigator();
    // Check if current safety allows for the proposed step
    if (safe[itr] > pstep[itr]) {
      step[itr]      = pstep[itr];
      isonbdr[itr]   = false;
      *outstate[itr] = *instate[itr];
      continue;
    }

    step[itr] = nav->ComputeStepAndSafetyAndPropagatedState(
        Vector3D_t(x[itr], y[itr], z[itr]), Vector3D_t(dirx[itr], diry[itr], dirz[itr]),
        Math::Min<double>(1.E20, pstep[itr]), *instate[itr], *outstate[itr] /* the paths */, !isonbdr[itr], safe[itr]);
    step[itr] = Math::Max<double>(2. * gTolerance, step[itr] + 2. * gTolerance);
    safe[itr] = Math::Max<double>(safe[itr], 0);
    // onboundary with respect to new point
    isonbdr[itr] = outstate[itr]->IsOnBoundary();

    //#### To add small step detection and correction - see ScalarNavInterfaceTGeo ####//
  }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void ScalarNavInterfaceVGM::NavFindNextBoundaryAndStep(Track &track)
{

  typedef Vector3D<Precision> Vector3D_t;
  constexpr double gTolerance = 1.e-9;

  // Check if current safety allows for the proposed step
  if (track.GetSafety() > track.GetPstep()) {
    track.SetSnext(track.GetPstep());
    track.SetBoundary(false);
    track.UpdateSameNextPath();
    return;
  }

  // Retrieve navigator for the track
  VNavigator const *nav = track.GetVolume()->GetNavigator();
  double safety         = track.GetSafety();
  double snext          = nav->ComputeStepAndSafetyAndPropagatedState(
      Vector3D_t(track.X(), track.Y(), track.Z()), Vector3D_t(track.Dx(), track.Dy(), track.Dz()),
      Math::Min<double>(1.E20, track.GetPstep()), *track.Path(), *track.NextPath(), !track.Boundary(), safety);
  snext = Math::Max<double>(2. * gTolerance, snext + 2. * gTolerance);
  track.SetSnext(snext);
  track.SetSafety(Math::Max<double>(safety, 0));
  // onboundary with respect to new point
  track.SetBoundary(track.NextPath()->IsOnBoundary());

  //#### To add small step detection and correction - see ScalarNavInterfaceTGeo ####//
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void ScalarNavInterfaceVGM::NavFindNextBoundary(Track &track)
{
  constexpr double gTolerance = 1.e-9;
  // back-up the pre-step point boundary flag
  track.SetBoundaryPreStep(track.Boundary());
  // Find distance to next boundary, within proposed step.
  typedef Vector3D<Precision> Vector3D_t;
  // Retrieve navigator for the track
  VNavigator const *newnav = track.GetVolume()->GetNavigator();
  // Check if current safety allows for the proposed step
  if (track.GetSafety() > track.GetPstep()) {
    track.SetSnext(track.GetPstep());
    track.SetBoundary(false);
    return;
  }
  double safety = track.GetSafety();
  double snext  = newnav->ComputeStepAndSafety(
      Vector3D_t(track.X(), track.Y(), track.Z()), Vector3D_t(track.Dx(), track.Dy(), track.Dz()),
      Math::Min<double>(1.E20, track.GetPstep()), *track.Path(), !track.Boundary(), safety);
  track.SetBoundary(snext < track.GetPstep());
  track.SetSnext(Math::Max<double>(2. * gTolerance, snext + 2. * gTolerance));
  track.SetSafety(Math::Max<double>(safety, 0));
}

//______________________________________________________________________________
void ScalarNavInterfaceVGM::NavIsSameLocation(int ntracks, const double *x, const double *y, const double *z,
                                              const double * /*dirx*/, const double * /*diry*/, const double * /*dirz*/,
                                              const VolumePath_t **start, VolumePath_t **end, bool *same,
                                              VolumePath_t *tmpstate)
{
  //
  // Checks if the navigation states corresponding to positions (x,y,z) are the
  // same as the ones pointed by start. Update new states in end.
  // Input:  ntracks - number of tracks to be checked
  //         x,y,z   - arrays of positions
  //         start   - starting navigation paths to compare with
  // Output: end     - navigation paths corresponding to given positions
  //         same    - flags showing if the end and start positions are matching

  for (int itr = 0; itr < ntracks; ++itr) {

    // cross check with answer from ROOT

    // TODO: not using the direction yet here !!
    bool samepath = GlobalLocator::HasSamePath(Vector3D<double>(x[itr], y[itr], z[itr]), *start[itr], *tmpstate);
    if (!samepath) {
      tmpstate->CopyTo(end[itr]);
#ifdef VECGEOM_CACHED_TRANS
      end[itr]->UpdateTopMatrix();
#endif
    }
    same[itr] = samepath;
  }
}

//______________________________________________________________________________
void ScalarNavInterfaceVGM::NavIsSameLocation(Track &track, bool &same, VolumePath_t *tmpstate)
{

  //#### NOT USING YET THE NEW NAVIGATORS ####//

  typedef Vector3D<Precision> Vector3D_t;

  // TODO: not using the direction yet here !!
  bool samepath = GlobalLocator::HasSamePath(Vector3D_t(track.X(), track.Y(), track.Z()), *track.Path(), *tmpstate);
  if (!samepath) {
    tmpstate->CopyTo(track.NextPath());
#ifdef VECGEOM_CACHED_TRANS
    track.NextPath()->UpdateTopMatrix();
#endif
  }
  same = samepath;
}

//______________________________________________________________________________
void ScalarNavInterfaceVGM::DisplaceTrack(Track &track, const double dir[3], double disp, double mindisp)
{
  // before calling displacement it must be checked that the point is not on boundary! It is done in the caller.
  typedef Vector3D<Precision> Vector3D_t;
  constexpr double reduceFactor = 0.99;
  double realDisp               = disp;
  bool isDisp                   = true;
  // check if displacement is beyond the current (post-step point i.e. updated) safety
  double postSafety = track.GetSafety() * reduceFactor; // make sure that we never reach the boundary in displacement
  if (disp > postSafety) {
    // Retrieve navigator for the track
    VNavigator const *newnav = track.GetVolume()->GetNavigator();
    // NOTE: we should have a navigator method to compute only the safety and the following call should be change to
    // that
    /*double step =*/newnav->ComputeStepAndSafety(Vector3D_t(track.X(), track.Y(), track.Z()),
                                                  Vector3D_t(dir[0], dir[1], dir[2]), disp, *track.Path(), true,
                                                  postSafety);
    postSafety *= reduceFactor; // make sure that we never reach the boundary in displacement
    if (postSafety > mindisp) {
      realDisp = std::min(postSafety, disp);
    } else {
      isDisp = false;
    }
  }
  // apply displacement if any
  if (isDisp) {
    track.SetPosition(track.X() + realDisp * dir[0], track.Y() + realDisp * dir[1], track.Z() + realDisp * dir[2]);
    track.DecreaseSafety(realDisp);
  }
}
} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
