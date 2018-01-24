#ifndef GEANT_NAVIGATION_INTERFACE_H
#define GEANT_NAVIGATION_INTERFACE_H

#ifdef USE_VECGEOM_NAVIGATOR
#include "ScalarNavInterfaceVGM.h"
#include "VectorNavInterface.h"
typedef Geant::ScalarNavInterfaceVGM ScalarNavInterface;
#else
#include "ScalarNavInterfaceTGeo.h"
typedef Geant::ScalarNavInterfaceTGeo ScalarNavInterface;
#endif

#endif
