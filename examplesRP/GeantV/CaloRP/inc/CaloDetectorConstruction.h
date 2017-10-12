//===--- CaloDetectorConstruction.h - GeantV ---------------------------------*- C++ -*-===//
//
//                     GeantV Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file CaloDetectorConstruction.h
 * @brief Implementation of geometry for GeantV calorimeter prototype
 * @author Ryan Schmitz
 * @date August 1, 2017
 */
//===----------------------------------------------------------------------===//

#ifndef CALO_DETECTOR_CONSTRUCTION
#define CALO_DETECTOR_CONSTRUCTION
#define NEW_DETECTOR
#include "GeantRunManager.h"
#include "GeantVDetectorConstruction.h"
#include "Material.h"
#include "Region.h"
namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {
class GeantRunManager;
class Material;
}
}
namespace userapplication {
class CaloDetectorConstruction : public Geant::GeantVDetectorConstruction {

  static const int kMaxAbsorbers = 10;

private:
  std::string fWorldMaterialName;
  std::string fAbsMaterialName[kMaxAbsorbers];

  geantphysics::Material *fAbsMaterial[kMaxAbsorbers];
  geantphysics::Material *fWorldMaterial;

  bool fUserLayerNum    = false;
  bool fUserAbsorberNum = false;
  bool fUserCaloYZ      = false;
  bool fUserThickness[kMaxAbsorbers];
  bool fUserMaterial[kMaxAbsorbers];
  bool fProdCutByLength = true;

  int fNumAbsorbers;
  int fNumLayers;
  int fAbsLogicVolumeID[kMaxAbsorbers];
  int fDetectorRegionIndex;
  double fGammaCut    = 0.1;
  double fElectronCut = 0.1;
  double fPositronCut = 0.1;

  double fAbsThickness[kMaxAbsorbers];
  double fLayerThickness;

  double fCaloSizeYZ;
  double fWorldSizeX;
  double fWorldSizeYZ;

  void SetDetectorMaterials();

public:
  CaloDetectorConstruction(Geant::GeantRunManager *runmgr) : GeantVDetectorConstruction(runmgr) {}
  ~CaloDetectorConstruction();

public:
  void SetAbsorberMaterialName(int absNum, std::string matName) { fUserMaterial[absNum] = true; fAbsMaterialName[absNum] = matName; }
  std::string GetAbsorberMaterialName(int absNum) const { return fAbsMaterialName[absNum]; }
  geantphysics::Material *GetAbsorberMaterial(int absNum) const { return fAbsMaterial[absNum]; }

  void SetNumLayers(int nLayers) { fUserLayerNum = true; fNumLayers = nLayers; }
  void SetNumAbsorbers(int nAbsorbers) { fUserAbsorberNum = true; fNumAbsorbers = nAbsorbers; }
  int GetNumLayers() const { return fNumLayers; }
  int GetNumAbsorbers() const { return fNumAbsorbers; }

  void SetProductionCutsByEnergy(double);
  void SetProductionCutsByLength(double);
  void SetDetectorGammaProductionCut(double);
  void SetDetectorElectronProductionCut(double);
  void SetDetectorPositronProductionCut(double);

  void SetAbsorberThickness(int absNum, double thickness) { fUserThickness[absNum] = true; fAbsThickness[absNum] = thickness; }
  double GetAbsorberThickness(int absNum) const { return fAbsThickness[absNum]; }

  int GetAbsorberLogicalVolumeID(int absorber) const { return fAbsLogicVolumeID[absorber]; }
  int GetDetectorRegionIndex() const { return fDetectorRegionIndex; }

  void SetDetectorYZ(double yz) { fUserCaloYZ  = true; fCaloSizeYZ = yz; }
  double GetDetectorX() const { return fLayerThickness * fNumLayers; }
  double GetDetectorYZ() const { return fCaloSizeYZ; }
  double GetWorldX() const { return fWorldSizeX; }
  double GetWorldYZ() const { return fWorldSizeYZ; }

  void CreateMaterials();
  void CreateGeometry();
};
}
#endif
