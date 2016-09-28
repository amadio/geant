#include "GeantVApplication.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GeantVApplication::GeantVApplication(GeantRunManager *runmgr):fRunMgr(runmgr)
{
  // Ctor..

}

//______________________________________________________________________________
void GeantVApplication::SetRunManager(GeantRunManager *runmgr) {
  fRunMgr = runmgr;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
