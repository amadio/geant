#include "GeantPropagator.h"
void loadvecgeomgeometry(GeantPropagator *prop) {
#ifdef USE_VECGEOM_NAVIGATOR
  prop->LoadVecGeomGeometry();
#endif
}
