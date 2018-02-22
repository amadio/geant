//===--- BasketCounters.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file BasketCounters.h
 * @brief A simple utility for handling track counting in baskets.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_BASKET_COUNTERS
#define GEANT_BASKET_COUNTERS

#include "Geant/Config.h"

/**
 * @brief Struct BasketCounters
 * @details Basket counters per simulation stage
 *
 */
namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

struct BasketCounters {
  size_t fNhandlers = 0; ///< number of handlers
  size_t fNscalar = 0;   ///< number of scalar DoIt calls per stage
  size_t fNvector = 0;   ///< number of basketized DoIt calls per stage
  size_t *fCounters; ///< counters

  BasketCounters(size_t nhandlers)
  {
    fNhandlers = nhandlers;
    fCounters  = new size_t[nhandlers];
    Reset();
  }

  ~BasketCounters() { delete[] fCounters; }

  GEANT_FORCE_INLINE
  void Reset()
  {
    for (size_t i = 0; i < fNhandlers; ++i)
      fCounters[i] = 0;
  }

  GEANT_FORCE_INLINE
  BasketCounters &operator+=(const BasketCounters &other)
  {
    fNscalar += other.fNscalar;
    fNvector += other.fNvector;
    for (size_t i = 0; i < fNhandlers; ++i)
      fCounters[i] += other.fCounters[i];
    return *this;
  }
};

} // GEANT_IMPL_NAMESPACE
} // Geant
#endif // GEANT_BASKET_COUNTERS
