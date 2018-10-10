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
  volatile size_t fNhandlers = 0; ///< number of handlers
  volatile size_t fNscalar   = 0; ///< number of scalar DoIt calls per stage
  volatile size_t fNvector   = 0; ///< number of basketized DoIt calls per stage
  volatile size_t *fCounters;     ///< counters
  volatile size_t *fFired;        ///< fired baskets per handler
  volatile size_t *fFlushed;      ///< flushed baskets per handler

  BasketCounters(size_t nhandlers)
  {
    fNhandlers = nhandlers;
    fCounters  = new size_t[nhandlers];
    fFired     = new size_t[nhandlers];
    fFlushed   = new size_t[nhandlers];
    Reset();
  }

  ~BasketCounters()
  {
    delete[] fCounters;
    delete[] fFired;
    delete[] fFlushed;
  }

  GEANT_FORCE_INLINE
  void Reset()
  {
    for (size_t i = 0; i < fNhandlers; ++i) {
      fCounters[i] = 0;
      fFired[i]    = 0;
      fFlushed[i]  = 0;
    }
  }

  GEANT_FORCE_INLINE
  BasketCounters &operator+=(const BasketCounters &other)
  {
    fNscalar += other.fNscalar;
    fNvector += other.fNvector;
    for (size_t i = 0; i < fNhandlers; ++i) {
      fCounters[i] += other.fCounters[i];
      fFired[i] += other.fFired[i];
      fFlushed[i] += other.fFlushed[i];
    }
    return *this;
  }

  GEANT_FORCE_INLINE
  size_t GetNcalls() const { return (fNscalar + fNvector); }

  GEANT_FORCE_INLINE
  void Increment(size_t ihandler, size_t threshold)
  {
    fCounters[ihandler]++;
    fFired[ihandler] += (size_t)(fCounters[ihandler] == threshold);
    fCounters[ihandler] &= threshold - 1;
  }
};

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
#endif // GEANT_BASKET_COUNTERS
