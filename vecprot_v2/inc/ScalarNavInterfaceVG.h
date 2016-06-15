//===--- ScalarNavInterfaceVG.h - GeantV ---------------------------------*- C++ -*-===//
//
//                     GeantV Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file ScalarNavInterfaceVG.h
 * @brief Scalar navigation interface layer using VecGeom
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_SCALARNAVINTERFACEVG
#define GEANT_SCALARNAVINTERFACEVG

#include "Geant/Config.h"
#include "GeantTrackVec.h"
#include "Geant/Typedefs.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Class ScalarNavInterfaceVG
 */
class ScalarNavInterfaceVG {
public:
  /** @brief ScalarNavInterfaceVG default constructor */
  ScalarNavInterfaceVG() {}

  /** @brief ScalarNavInterfaceVG destructor */
  ~ScalarNavInterfaceVG() {}

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
  GEANT_CUDA_BOTH_CODE
  static
  void NavFindNextBoundaryAndStep(int ntracks, const double *pstep, 
         const double *x, const double *y, const double *z,
         const double *dirx, const double *diry, const double *dirz, 
         const VolumePath_t **instate, VolumePath_t **outstate, 
         double *step, double *safe, bool *isonbdr);
    
  /**
  * @brief Single track version of the function above */
  GEANT_CUDA_BOTH_CODE
  static
  void NavFindNextBoundaryAndStep(GeantTrack &track);

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
  GEANT_CUDA_BOTH_CODE
  static 
  void NavIsSameLocation(int ntracks, 
         const double *x, const double *y, const double *z, 
         const double *dirx, const double *diry, const double *dirz, 
         const VolumePath_t **start, VolumePath_t **end, bool *same, VolumePath_t *tmpstate);
  
  /**
  * @brief Single track version of the function above */
  GEANT_CUDA_BOTH_CODE
  static
  void NavIsSameLocation(GeantTrack &track, bool &same, VolumePath_t *tmpstate);
  
};
} // GEANT_IMPL_NAMESPACE

#ifdef GEANT_CUDA
#ifdef GEANT_NVCC
namespace cxx {
class ScalarNavInterfaceVG;
}
#else
namespace cuda {
class ScalarNavInterfaceVG;
}
#endif
#endif
} // Geant
#endif
