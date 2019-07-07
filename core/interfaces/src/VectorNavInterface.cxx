#include "Geant/VectorNavInterface.h"

#include "navigation/VNavigator.h"
#include "volumes/PlacedVolume.h"
#include "base/Vector3D.h"
#include "base/Transformation3D.h"
#include "base/Global.h"
#include "management/GeoManager.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

using namespace VECGEOM_NAMESPACE;

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void VectorNavInterface::NavFindNextBoundary(int ntracks, const double *pstep, const double *x, const double *y,
                                             const double *z, const double *dirx, const double *diry,
                                             const double *dirz, const VolumePath_t **instate, double *snext,
                                             double *safe, bool *compsafety)
{
  // Find the next boundary and state after propagating to the boundary.
  // Input:  ntracks - number of tracks
  //         pstep - proposed step
  //         x, y, z, dirx, diry, dirz - particle position and direction
  //         instate - input particle navigation state
  //         safe - estimated safety value for the input point
  //         compsafety - starting point is on a boundary, so, compute safety
  // Output: outstate - navigation state after propagation. If boundary further
  //           than proposed step, outstate has to match instate
  //         step - propagation step for which the state is sampled
  //         safety - calculated safety value for the input point
  //         compsafety - propagated point is not on a boundary

  typedef SOA3D<double> SOA3D_t;
  constexpr double gTolerance = 1.e-9;
  VNavigator const *newnav    = instate[0]->Top()->GetLogicalVolume()->GetNavigator();
  newnav->ComputeStepsAndSafeties(
      SOA3D_t(const_cast<double *>(x), const_cast<double *>(y), const_cast<double *>(z), ntracks),
      SOA3D_t(const_cast<double *>(dirx), const_cast<double *>(diry), const_cast<double *>(dirz), ntracks), pstep,
      instate, snext, compsafety, safe);

  for (int itr = 0; itr < ntracks; ++itr) {
    // onboundary with respect to new point
    if (!compsafety[itr])
      safe[itr] = 0.;
    else
      safe[itr] = vecCore::math::Max<double>(safe[itr], 0.);
    compsafety[itr] = snext[itr] >= pstep[itr];
    snext[itr]      = vecCore::math::Max<double>(2. * gTolerance, snext[itr] + 2. * gTolerance);
  }
}

//______________________________________________________________________________
void VectorNavInterface::NavIsSameLocation(int ntracks, const double *x, const double *y, const double *z,
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

  //#### NOT USING YET THE NEW NAVIGATORS, NOR BEING VECTORIZED ####//

  typedef Vector3D<Precision> Vector3D_t;

  for (int itr = 0; itr < ntracks; ++itr) {
    // TODO: not using the direction yet here !!
    bool samepath = GlobalLocator::HasSamePath(Vector3D_t(x[itr], y[itr], z[itr]), *start[itr], *tmpstate);
    if (!samepath) {
      tmpstate->CopyTo(end[itr]);
#ifdef VECGEOM_CACHED_TRANS
      end[itr]->UpdateTopMatrix();
#endif
    }
    same[itr] = samepath;
  }
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
