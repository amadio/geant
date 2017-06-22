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

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/** @brief GeantVDetectorConstruction class */
class GeantVDetectorConstruction {
private:
  GeantRunManager *fRunMgr = nullptr; /*taskData*/
  
public:  
  /** @brief GeantVDetectorConstruction constructor */	
  GeantVDetectorConstruction(GeantRunManager *runmgr) : fRunMgr(runmgr) {}

  /** @brief GeantVDetectorConstruction destructor */
  virtual ~GeantVDetectorConstruction() {}

  /** @brief Creation of materials (optional) */
  virtual void CreateMaterials() {}

  /** @brief Creation of geometry (mandatory) */
  virtual void CreateGeometry() = 0;

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
