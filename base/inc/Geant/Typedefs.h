#ifndef GEANT_TYPEDEFS_H
#define GEANT_TYPEDEFS_H

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/NavigationState.h"
typedef VECGEOM_NAMESPACE::NavigationState VolumePath_t;
#include "materials/Material.h"
typedef VECGEOM_NAMESPACE::Material Material_t;
#include "materials/Medium.h"
typedef VECGEOM_NAMESPACE::Medium Medium_t;
#include "volumes/LogicalVolume.h"
typedef VECGEOM_NAMESPACE::LogicalVolume Volume_t;
#include "volumes/PlacedVolume.h"
typedef VECGEOM_NAMESPACE::VPlacedVolume Node_t;
#include "materials/Particle.h"
typedef VECGEOM_NAMESPACE::Particle Particle_t;
#else
#include "TGeoBranchArray.h"
typedef TGeoBranchArray VolumePath_t;
#include "TGeoMaterial.h"
typedef TGeoMaterial Material_t;
#include "TGeoMedium.h"
typedef TGeoMedium Medium_t;
#include "TGeoVolume.h"
typedef TGeoVolume Volume_t;
#include "TGeoNode.h"
typedef TGeoNode Node_t;
#include "TDatabasePDG.h"
typedef TParticlePDG Particle_t;
#endif
#endif
