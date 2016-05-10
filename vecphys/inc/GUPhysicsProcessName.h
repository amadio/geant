#ifndef GUPhysicsProcessName_H
#define GUPhysicsProcessName_H 

#include "GUPhysicsProcessIndex.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

static const char* GUPhysicsProcessName[kNumberPhysicsProcess] = {
  "PhotonProcess",
  "ElectronProcess"
};

}
} // end namespace vecphys

#endif
