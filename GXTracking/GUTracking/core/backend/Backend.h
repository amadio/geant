/// \file vector/Backend.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECPHYS_BACKEND_BACKEND_H_
#define VECPHYS_BACKEND_BACKEND_H_

#include "base/Global.h"

#ifdef VECPHYS_NVCC
#include "backend/cuda/Backend.h"
#elif defined(VECPHYS_VC)
#include "backend/vc/Backend.h"
#elif defined(VECPHYS_CILK)
#include "backend/cilk/Backend.h"
#else
#include "backend/scalar/Backend.h"
#endif

#endif // VECPHYS_BACKEND_BACKEND_H_
