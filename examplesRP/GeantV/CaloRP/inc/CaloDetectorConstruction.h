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

static const int maxAbsorbers = 10;
  private:

	std::string fWorldMaterialName;
        std::string fAbsMaterialName[maxAbsorbers];

	geantphysics::Material *fAbsMaterial[maxAbsorbers];
	geantphysics::Material *fWorldMaterial;

	bool userLayerNum=false;
        bool userAbsorberNum=false;
        bool userCaloYZ=false;
	bool userThickness[maxAbsorbers];
	bool userMaterial[maxAbsorbers];
	bool prodCutByLength=true;

        int numAbsorbers;
        int numLayers;
	int fAbsLogicVolumeID[maxAbsorbers];
	int fDetectorRegionIndex;
	double fGammaCut=0.1;
	double fElectronCut=0.1;
	double fPositronCut=0.1;

        double fAbsThickness[maxAbsorbers];
        double fLayerThickness;

        double fCaloSizeYZ;
        double fWorldSizeX;
	double fWorldSizeYZ;

	void SetDetectorMaterials();

  public:
	CaloDetectorConstruction(Geant::GeantRunManager *runmgr) : GeantVDetectorConstruction(runmgr) {}
	~CaloDetectorConstruction();

  public:
        void SetAbsorberMaterialName(int,std::string);
        std::string GetAbsorberMaterialName(int);
        geantphysics::Material* GetAbsorberMaterial(int);

        void SetNumLayers(int);
        void SetNumAbsorbers(int);
        int GetNumLayers();
        int GetNumAbsorbers();

	void SetProductionCutsByEnergy(double);
	void SetProductionCutsByLength(double);
	void SetDetectorGammaProductionCut(double);
	void SetDetectorElectronProductionCut(double);
	void SetDetectorPositronProductionCut(double);

        void SetAbsorberThickness(int,double);
        double GetAbsorberThickness(int);

	int GetAbsorberLogicalVolumeID(int);
	int GetDetectorRegionIndex();

        void SetDetectorYZ(double);
        double GetDetectorX();
        double GetDetectorYZ();
        double GetWorldX();
        double GetWorldYZ();

	void CreateMaterials();
	void CreateGeometry();
};
}
#endif
