
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

#ifdef USE_VECGEOM_CONTAINERS
// vecgeom::Vector
#include "base/Vector.h"
template <typename T>
struct VectorHelper {
  typedef vecgeom::Vector<T> Vector_t;
};
// vecgeom::map
#include "base/Map.h"
template <typename TK, typename TV>
struct MapHelper {
  typedef vecgeom::map<TK,TV> Map_t;
};
#else
// std::vector
#include <vector>
template <typename T>
struct VectorHelper {
  typedef std::vector<T> Vector_t;
};
// std::map
#include <map>
template <typename TK, typename TV>
struct MapHelper {
  typedef std::map<TK,TV> Map_t;
};
#endif

#endif
