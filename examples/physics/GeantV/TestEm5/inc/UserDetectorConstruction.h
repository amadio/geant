
#ifndef USERDETECTORCONSTRUCTION_H
#define USERDETECTORCONSTRUCTION_H

#include "GeantVDetectorConstruction.h"
#include "Geant/Typedefs.h"
#include "Geant/Config.h"
#include "GeantFwd.h"

#include <string>

namespace Geant {
  inline namespace GEANT_IMPL_NAMESPACE {
    class GeantVDetectorConstruction;
    class GeantRunManager;
  }
}

namespace geantphysics {
  inline namespace GEANT_IMPL_NAMESPACE {
    class Material;
  }
}


namespace userapplication {

/**
 * @brief User detector construction for TeatEm5.
 *
 * The detector is a simple slab (with configurable thickness and material) in a world. The whole detector is in a
 * user defined region with configurable secondary particle production threads.
 *
 * @class   UserDetectorConstruction
 * @author  M Novak
 * @date    July 2017
 */



class UserDetectorConstruction : public Geant::GeantVDetectorConstruction {
public:
  // CTR
  UserDetectorConstruction(Geant::GeantRunManager *runmgr);
  // DTR
 ~UserDetectorConstruction();

  // interface mathod to define a set of custom materials for the application
  virtual void CreateMaterials();
  // interface method to define the geometry for the application
  virtual void CreateGeometry();
  //
  // the following methods need to be called before CreateGeometry() (i.e. initialization of detector)
  //
  // method to set the name of the target material (NIST_MAT_* or one of the custom material names defined in the
  // CreateMaterials) method
  void SetTargetMaterialName(std::string &matname) { fTargetMatName = matname; }
  // method to set the target thickness (in internal [length] units)
  void SetTargetThickness(double val) { fTargetX = val; }
  // methods to set the secondary production cuts in length (internal [length] unit)
  void SetProductionCuts(double val)        { fGammaCut    = fElectronCut = fPositronCut = val; }
  void SetGammaProductionCut(double val)    { fGammaCut    = val; }
  void SetElectronProductionCut(double val) { fElectronCut = val; }
  void SetPositronProductionCut(double val) { fPositronCut = val; }
  //
  // the following methods will give the proper response only after the detector is Initialized (after CreateGeometry())
  //
  // method to get the target material of the actual detector
  const geantphysics::Material* GetTargetMaterial()        const { return fTargetMaterial; }
  // method to get the target start x coordinate of the actual detector
  double                        GetWorldXStart()           const { return -fWorldX; }
  // method to get the target end x coordinate of the actual detector
  double                        GetTargetXStart()          const { return -fTargetX; }
  // method to get the target logical volume ID of the actual detector
  int                           GetTargetLogicalVolumeID() const { return fTargetLogicalVolumeID; }
  // method to get the index of the reagion in which the target is located in the actual detector
  int                           GetTargetRegionIndex()     const { return fTargetRegionIndx; }

private:
  // internal method to comupte/set detector parameters
  void ComputeSetup();

private:
  std::string             fTargetMatName;
  int                     fTargetLogicalVolumeID;
  int                     fTargetRegionIndx;
  double                  fTargetYZ;
  double                  fTargetX;
  double                  fWorldYZ;
  double                  fWorldX;
  double                  fGammaCut;
  double                  fElectronCut;
  double                  fPositronCut;
  geantphysics::Material *fTargetMaterial;
  geantphysics::Material *fWorldMaterial;
};

}        // namespace userapplication

#endif   // DETECTORCONSTRUCTION_H
