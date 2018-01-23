//===--- ExN03DetectorConstruction.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file ExN03DetectorConstruction.h
 * @brief Implementation of user detector construction for ExN03
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef EXN03_DETECTORCONSTRUCTION
#define EXN03_DETECTORCONSTRUCTION

#include "GeantVDetectorConstruction.h"
#include "Geant/Typedefs.h"
#include "GeantFwd.h"


/** @brief DetectorConstruction class for ExN03*/
class ExN03DetectorConstruction : public Geant::GeantVDetectorConstruction {
  std::string fGeomFileName;
public:
  ExN03DetectorConstruction(const char *filename, Geant::GeantRunManager *mgr)
    : GeantVDetectorConstruction(mgr), fGeomFileName(filename) {}
  /** @brief GeantVDetectorConstruction constructor */	
  virtual ~ExN03DetectorConstruction() {}

  /** @brief Creation of geometry (mandatory) */
  virtual void CreateGeometry() { LoadGeometry(fGeomFileName.c_str()); }
};

#endif
