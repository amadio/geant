#include "Geant/ScalarNavInterfaceVG.h"

#include "volumes/PlacedVolume.h"
#include "base/Vector3D.h"
#include "base/Transformation3D.h"
#include "base/Global.h"
#include "management/GeoManager.h"
#include "navigation/VNavigator.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

using namespace VECGEOM_NAMESPACE;

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void ScalarNavInterfaceVG::NavFindNextBoundaryAndStep(int ntracks, const double *pstep, const double *x,
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
    // Check if current safety allows for the proposed step
    if (safe[itr] > pstep[itr]) {
      step[itr]      = pstep[itr];
      isonbdr[itr]   = false;
      *outstate[itr] = *instate[itr];
      continue;
    }
    VNavigator const *nav = instate[itr]->Top()->GetLogicalVolume()->GetNavigator();
    step[itr]             = nav->ComputeStepAndSafetyAndPropagatedState(
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
void ScalarNavInterfaceVG::NavFindNextBoundaryAndStep(Track &track)
{
  typedef Vector3D<Precision> Vector3D_t;

  constexpr double gTolerance = 1.e-9;
  // Check if current safety allows for the proposed step
  if (track.GetSafety() > track.GetPstep()) {
    track.SetStep(track.GetPstep());
    track.SetBoundary(false);
    track.UpdateSameNextPath();
    return;
  }

  VNavigator const *nav = track.GetVolume()->GetNavigator();
  double step           = track.GetStep();
  nav->FindNextBoundaryAndStep(Vector3D_t(track.X(), track.Y(), track.Z()),
                               Vector3D_t(track.Dx(), track.Dy(), track.Dz()), *track.Path(), *track.NextPath(),
                               Math::Min<double>(1.E20, track.GetPstep()), step);
  double safety = track.GetSafety();
  double snext  = nav->ComputeStepAndSafetyAndPropagatedState(
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
void ScalarNavInterfaceVG::NavIsSameLocation(int ntracks, const double *x, const double *y, const double *z,
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
void ScalarNavInterfaceVG::NavIsSameLocation(Track &track, bool &same, VolumePath_t *tmpstate)
{
  // TODO: not using the direction yet here !!
  bool samepath =
      GlobalLocator::HasSamePath(Vector3D<double>(track.X(), track.Y(), track.Z()), *track.Path(), *tmpstate);
  if (!samepath) {
    tmpstate->CopyTo(track.NextPath());
#ifdef VECGEOM_CACHED_TRANS
    track.NextPath()->UpdateTopMatrix();
#endif
  }
  same = samepath;
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
