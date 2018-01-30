//===--- CMSDetectorConstruction.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file CMSDetectorConstruction.h
 * @brief Implementation of user detector construction for CMS
 * @author Guilherme Lima, Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef CMSDETECTORCONSTRUCTION
#define CMSDETECTORCONSTRUCTION

#include "GeantVDetectorConstruction.h"
//#include "Geant/Typedefs.h"
//#include "GeantFwd.h"

/** @brief DetectorConstruction class for CMS*/
class CMSDetectorConstruction : public Geant::GeantVDetectorConstruction {

  std::string fGeomFileName;

public:
  /** @brief CMSDetectorConstruction constructor */	
  CMSDetectorConstruction(const char *filename, Geant::GeantRunManager *mgr)
    : GeantVDetectorConstruction(mgr), fGeomFileName(filename) {}

  /** @brief CMSDetectorConstruction destructor */	
  virtual ~CMSDetectorConstruction() {}

  /** @brief Creation of geometry (mandatory) */
  //virtual void CreateGeometry() { LoadGeometry(fGeomFileName.c_str()); }
  void CreateGeometry() final { LoadGeometry(fGeomFileName.c_str()); }
};

#endif
