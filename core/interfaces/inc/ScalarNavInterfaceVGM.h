//===--- ScalarNavInterfaceVGM.h - GeantV ---------------------------------*- C++ -*-===//
//
//                     GeantV Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file ScalarNavInterfaceVGM.h
 * @brief Scalar navigation interface layer using VecGeom with multiple navigators
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_SCALARNAVINTERFACEVGM
#define GEANT_SCALARNAVINTERFACEVGM

#include "Geant/Config.h"
#include "Geant/Typedefs.h"
#include "GeantTrack.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Class ScalarNavInterfaceVGM
 */
class ScalarNavInterfaceVGM {
public:
  /** @brief ScalarNavInterfaceVGM default constructor */
  ScalarNavInterfaceVGM() {}

  /** @brief ScalarNavInterfaceVGM destructor */
  ~ScalarNavInterfaceVGM() {}

  /**
   * @brief Function for navigation that find next boundary and step
   *
   * @param pstep proposed step
   * @param x X position
   * @param y Y position
   * @param z Z position
   * @param dirx X direction
   * @param diry Y direction
   * @param dirz Z direction
   * @param instate Initial navigation state
   * @param outstate Final navigation state if next boundary would be crossed
   * @param step Step to be made to next boundary
   * @param safe Safety distance
   * @param isonbdr Boundary flag set if next boundaries closer than proposed step
   */
  VECCORE_ATT_HOST_DEVICE
  static
  void NavFindNextBoundaryAndStep(int ntracks, const double *pstep,
         const double *x, const double *y, const double *z,
         const double *dirx, const double *diry, const double *dirz,
         const VolumePath_t **instate, VolumePath_t **outstate,
         double *step, double *safe, bool *isonbdr);

  /**
  * @brief Single track version of the function above */
  VECCORE_ATT_HOST_DEVICE
  static
  void NavFindNextBoundaryAndStep(GeantTrack &track);

  /** @brief Find distance to next boundary */
  VECCORE_ATT_HOST_DEVICE
  static
  void NavFindNextBoundary(GeantTrack &track);

//  VECCORE_ATT_HOST_DEVICE
//  static
//  void NavFindNextBoundaryMSC(GeantTrack &track,double dist);


  /**
   * @brief Function for navigation that checks if location is the same or not
   *
   * @param ntracks Number of tracks
   * @param x X positions
   * @param y Y positions
   * @param z Z positions
   * @param dirx X directions
   * @param diry Y directions
   * @param dirz Z directions
   * @param start Start volume paths
   * @param end End volume paths
   * @param same Boolean return flags specifying if the location is same
   & @param tmpstate Temporary navigation state to be used internally
   */
  VECCORE_ATT_HOST_DEVICE
  static
  void NavIsSameLocation(int ntracks,
         const double *x, const double *y, const double *z,
         const double *dirx, const double *diry, const double *dirz,
         const VolumePath_t **start, VolumePath_t **end, bool *same, VolumePath_t *tmpstate);

  /** @brief Single track version of the function above */
  VECCORE_ATT_HOST_DEVICE
  static
  void NavIsSameLocation(GeantTrack &track, bool &same, VolumePath_t *tmpstate);

  /** @brief Displace track along given normalized direction with given step. If boundary is crossed,
       displace to boundary and relocate, setting the boundary flag
   * @param track Track to be displaced
   * @param dir New direction
   * @param step Displacement step
   * @return The actual displacement.
   */
  VECCORE_ATT_HOST_DEVICE
  static
  void DisplaceTrack(GeantTrack &track, const double dir[3], double step, double mindisp);

};
} // GEANT_IMPL_NAMESPACE

#ifdef GEANT_CUDA_ENABLED
#ifdef VECCORE_CUDA
namespace cxx {
class ScalarNavInterfaceVGM;
}
#else
namespace cuda {
class ScalarNavInterfaceVGM;
}
#endif
#endif
} // Geant
#endif
