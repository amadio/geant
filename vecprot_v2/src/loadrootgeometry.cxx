#ifndef USE_VECGEOM_NAVIGATOR
#define UNDEF_VECGEOM_NAVIGATOR
#define USE_VECGEOM_NAVIGATOR
#endif
#include "management/RootGeoManager.h"
void loadrootgeometry() { vecgeom::RootGeoManager::Instance().LoadRootGeometry(); }
#ifdef UNDEF_VECGEOM_NAVIGATOR
#undef UNDEF_VECGEOM_NAVIGATOR
#undef USE_VECGEOM_NAVIGATOR
#endif
