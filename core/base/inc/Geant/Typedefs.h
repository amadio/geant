#ifndef GEANT_TYPEDEFS_H
#define GEANT_TYPEDEFS_H

#include "Geant/Config.h"
#include "Geant/VectorTypes.h"

#ifdef VECCORE_CUDA
#include "base/Vector.h"
  template <class T>
  using vector_t = vecgeom::Vector<T>;
#else
#include <vector>
#ifdef GEANT_USE_NUMA
#include <GeantNuma.h>
  template <class T>
  using vector_t = std::vector<T, Geant::NumaAllocator<T>>;
#else
  template <class T>
  using vector_t = std::vector<T>;
#endif
#endif

namespace geantphysics {
   class Particle;
}
typedef geantphysics::Particle Particle_t;

#include "navigation/NavigationState.h"
typedef VECGEOM_NAMESPACE::NavigationState VolumePath_t;
#include "Material.h"
typedef geantphysics::Material Material_t;
#include "volumes/LogicalVolume.h"
typedef VECGEOM_NAMESPACE::LogicalVolume Volume_t;
#include "volumes/PlacedVolume.h"
typedef VECGEOM_NAMESPACE::VPlacedVolume Node_t;
#endif
