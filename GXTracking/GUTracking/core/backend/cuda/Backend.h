/// \file vc/backend.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)
//
// temporary clone for GUTracking based on VecGeom - syjun
//
#ifndef VECPHYS_CUDABACKEND_H
#define VECPHYS_CUDABACKEND_H

#include "base/Global.h"
#include "backend/scalar/Backend.h"
#include "backend/cuda/Interface.h"

#ifdef __CUDACC__
#include <cuda.h>
#include <curand_kernel.h>
#endif

namespace vecphys {

#ifdef VECPHYS_NVCC
inline
#endif
namespace cuda {

struct kCuda {

  typedef int        Int_t;
  typedef Precision  Double_t;
  typedef bool       Bool_t;
  typedef int        Index_t;

  const static       bool kTrue = true;
  const static       bool kFalse = false;
  const static       bool early_returns = false;

  const static       int kSize = 1;
  static constexpr   Precision kOne = 1.0;
  static constexpr   Precision kZero = 0.0;

};

} // End of cuda namespace
} // End global namespace

#endif // VECPHYS_CUDABACKEND_H
