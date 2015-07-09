#ifdef USE_VECGEOM_NAVIGATOR
#define RESTORE_USE_VECGEOM_NAVIGATOR
#undef USE_VECGEOM_NAVIGATOR
#endif

#include "TGeoManager.h"
#include "WorkloadManager.h"
using namespace Geant;

ClassImp(WorkloadManager)


int WorkloadManager::ThreadId() {
     gGeoManager->SetMultiThread();
     return TGeoManager::ThreadId();
  }
#ifdef RESTORE_USE_VEGEOM_NAVIGATOR
#define USE_VECGEOM_NAVIGATOR
#endif


