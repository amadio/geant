//===--- VectorNavInterface.h - GeantV ---------------------------------*- C++ -*-===//
//
//                     GeantV Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file VectorNavInterface.h
 * @brief Vectorized navigation interface layer
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_VECTORNAVINTERFACE
#define GEANT_VECTORNAVINTERFACE

#include "Geant/Config.h"
#include "GeantTrackVec.h"
#include "Geant/Typedefs.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Class VectorNavInterface
 */
class VectorNavInterface {
public:
  /** @brief VectorNavInterface default constructor */
  VectorNavInterface() {}

  /** @brief VectorNavInterface destructor */
  ~VectorNavInterface() {}

  /**
   * @brief Function for navigation that find next boundary and step
   *
   * @param ntracks Number of tracks
   * @param pstep proposed steps
   * @param x X positions
   * @param y Y positions
   * @param z Z positions
   * @param dirx X directions
   * @param diry Y directions
   * @param dirz Z directions
   * @param instate Initial navigation states
   * @param outstate Final navigation states if next boundary would be crossed
   * @param step Steps to be made to next boundaries
   * @param safe Safety distances
   * @param isonbdr Boundary flags set if next boundaries closer than proposed steps
   */
  VECCORE_ATT_HOST_DEVICE
  static
  void NavFindNextBoundaryAndStep(int ntracks, const double *pstep, 
         const double *x, const double *y, const double *z,
         const double *dirx, const double *diry, const double *dirz, 
         const VolumePath_t **instate, VolumePath_t **outstate, 
         double *step, double *safe, bool *isonbdr);

  VECCORE_ATT_HOST_DEVICE
  static
  void NavFindNextBoundary(int ntracks, const double *pstep, 
         const double *x, const double *y, const double *z,
         const double *dirx, const double *diry, const double *dirz, 
         const VolumePath_t **instate, 
         double *step, double *safe, bool *isonbdr);
    
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
   */
  VECCORE_ATT_HOST_DEVICE
  static
  void NavIsSameLocation(int ntracks, 
         const double *x, const double *y, const double *z, 
         const double *dirx, const double *diry, const double *dirz, 
         const VolumePath_t **start, VolumePath_t **end, bool *same, VolumePath_t *tmpstate);
  
  
};
} // GEANT_IMPL_NAMESPACE

#ifdef GEANT_CUDA_ENABLED
#ifdef VECCORE_CUDA
namespace cxx {
class VectorNavInterface;
}
#else
namespace cuda {
class VectorNavInterface;
}
#endif
#endif
} // Geant
#endif
