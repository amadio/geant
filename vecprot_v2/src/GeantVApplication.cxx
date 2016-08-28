#include "GeantVApplication.h"

//______________________________________________________________________________
GeantVApplication::GeantVApplication(GeantRunManager *runmgr):fRunMgr(runmgr)
{
  // Ctor..

}

//______________________________________________________________________________
void GeantVApplication::SetRunManager(GeantRunManager *runmgr) {
  fRunMgr = runmgr;
}
