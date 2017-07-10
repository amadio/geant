//===--- GeantVDetectorConstruction.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantVDetectorConstruction.h
 * @brief Implementation of user detector construction in Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_VDETECTORCONSTRUCTION
#define GEANT_VDETECTORCONSTRUCTION

#include "Geant/Typedefs.h"

#ifdef USE_ROOT
class TGeoMaterial;
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantRunManager;
class TaskBroker;

/** @brief GeantVDetectorConstruction class */
class GeantVDetectorConstruction {
private:
  GeantRunManager *fRunMgr = nullptr;
// Material conversion callback function
#ifdef USE_VECGEOM_NAVIGATOR
#ifdef USE_ROOT
  static
  std::function<void*(TGeoMaterial const *)> CreateMaterialConversion();
#endif
#endif

public:  
  /** @brief GeantVDetectorConstruction constructor */	
  GeantVDetectorConstruction(GeantRunManager *runmgr) { fRunMgr = runmgr; }

  /** @brief GeantVDetectorConstruction destructor */
  virtual ~GeantVDetectorConstruction() {}

  /** @brief Creation of materials (optional) */
  virtual void CreateMaterials() {}

  /** @brief Creation of geometry (mandatory) */
  virtual void CreateGeometry() = 0;
  
  /** Setup geometry connections to user indices */
  int SetupGeometry(vector_t<Volume_t const *> &volumes, TaskBroker *broker = nullptr);
  
  /** Load geometry from a file. Can be called from the user-defined CreateGeometry */
  bool LoadGeometry(const char *filename);

  // these methods will become private and non-static after migration to 
  // user detector construction is complete
  static bool LoadVecGeomGeometry(TaskBroker *broker = nullptr);
  static void InitNavigators();

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
