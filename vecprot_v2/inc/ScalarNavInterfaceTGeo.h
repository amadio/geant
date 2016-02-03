//===--- ScalarNavInterfaceTGeo.h - GeantV ---------------------------------*- C++ -*-===//
//
//                     GeantV Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file ScalarNavInterfaceTGeo.h
 * @brief Scalar navigation interface layer using ROOT TGeo
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_SCALARNAVINTERFACETGEO
#define GEANT_VECTORNAVINTERFACETGEO

#include "Geant/Config.h"
#include "GeantTrack.h"
#include "Geant/Typedefs.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Class ScalarNavInterfaceTGeo
 */
class ScalarNavInterfaceTGeo {
public:
  /** @brief ScalarNavInterfaceTGeo default constructor */
  ScalarNavInterfaceTGeo() {}

  /** @brief ScalarNavInterfaceTGeo destructor */
  ~ScalarNavInterfaceTGeo() {}

  /**
   * @brief Function for navigation that find next boundary and step
   *
   * @param ntracks Number of tracks
   * @param pstep proposed steps (input)
   * @param x X positions (input)
   * @param y Y positions (input)
   * @param z Z positions (input)
   * @param dirx X directions (input)
   * @param diry Y directions (input)
   * @param dirz Z directions (input)
   * @param instate Initial navigation states (input)
   * @param outstate Final navigation states (output)
   * @param step Steps to be made to next boundaries or to proposed steps (output)
   * @param safe Safety distances, used as input and output (input/output)
   * @param isonbdr Boundary flags set if next boundaries closer than proposed step (input/output)
   */
  GEANT_CUDA_BOTH_CODE
  void NavFindNextBoundaryAndStep(int ntracks, const double *pstep, 
         const double *x, const double *y, const double *z,
         const double *dirx, const double *diry, const double *dirz, 
         const VolumePath_t **instate, VolumePath_t **outstate, 
         double *step, double *safe, bool *isonbdr);
    
  /**
   * @brief Function for navigation that checks if location is the same or not
   *
   * @param ntracks Number of tracks
   * @param x X positions
   * @param y Y positions
   * @param z Z positions   
   * @param start Start volume path
   * @param end End volume path
   * @param same Boolean return flag specifying if the location is same
   */
  void NavIsSameLocation(int ntracks, 
         const double *x, const double *y, const double *z, 
         const VolumePath_t **start, VolumePath_t **end, bool *same, VolumePath_t *tmpstate);
  
  
};
} // GEANT_IMPL_NAMESPACE

#ifdef GEANT_CUDA
#ifdef GEANT_NVCC
namespace cxx {
class ScalarNavInterfaceTGeo;
}
#else
namespace cuda {
class ScalarNavInterfaceTGeo;
}
#endif
#endif
} // Geant
#endif
