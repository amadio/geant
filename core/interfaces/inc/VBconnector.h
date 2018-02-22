#ifndef GEANT_VBCONNECTOR
#define GEANT_VBCONNECTOR

#ifdef USE_ROOT
#include "TGeoExtension.h"
#endif

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

/** @brief Volume-basket manager connector structure attached to volumes as extension */
#if defined(USE_ROOT) && !defined(VECCORE_CUDA)
class VBconnector : public TGeoExtension {
#else
class VBconnector {
#endif
public:
  int index;                      /** Index of basket manager */
  VECCORE_ATT_HOST_DEVICE
  VBconnector(int i) : index(i) {}
#if defined(USE_ROOT) && !defined(VECCORE_CUDA)
  virtual TGeoExtension *Grab() { return this; }
  virtual void Release() const {}
#endif
};

}
}

#endif
