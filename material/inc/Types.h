
#ifndef TYPES_H
#define TYPES_H

/**
 * @brief   Generic vector and map containers.
 * @class   Types
 * @author  M Novak
 * @date    March 2017
 *
 * Makes possible to switch between vecgeom and std containers with a -DUSE_VECGEOM_CONTAINERS cmake option.
 *
 */
#include "Geant/Config.h"

#ifdef USE_VECGEOM_CONTAINERS
#include "base/Vector.h"
#include "base/Map.h"
#else
#include <vector>
#include <map>
#endif

namespace geantphysics {

#ifdef USE_VECGEOM_CONTAINERS
  template <class T>
  using Vector_t = vecgeom::Vector<T>;
  template <class KeyT, class ValueT>
  using Map_t    = vecgeom::map<KeyT,ValueT>;
#else
  template <class T>
  using Vector_t = std::vector<T>;
  template <class KeyT, class ValueT>
  using Map_t    = std::map<KeyT,ValueT>;
#endif

} // namespace geantphysics

#endif // TYPES_H
