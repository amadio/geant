#include "GeantVDetectorConstruction.h"

#include "Geant/Error.h"
#include "Region.h"
#include "TaskBroker.h"
#include "VBconnector.h"
#include "GeantRunManager.h"

#ifdef USE_VECGEOM_NAVIGATOR

#include "navigation/VNavigator.h"
#include "navigation/SimpleNavigator.h"
#include "navigation/NewSimpleNavigator.h"
#include "navigation/SimpleABBoxNavigator.h"
#include "navigation/SimpleABBoxLevelLocator.h"
#include "navigation/HybridNavigator2.h"
#include "Material.h"
#include "Element.h"
#ifdef USE_ROOT
#include "TGeoRegion.h"
#include "management/RootGeoManager.h"
#endif
#include "volumes/LogicalVolume.h"
#include "volumes/PlacedVolume.h"

#else

#ifdef USE_ROOT
#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TGeoVoxelFinder.h"
#include "TGeoNode.h"
#include "TGeoMaterial.h"
#include "TGeoRegion.h"
#include "TROOT.h"
#endif

#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
int GeantVDetectorConstruction::SetupGeometry(vector_t<Volume_t const *> &volumes, TaskBroker *broker)
{
  // Setup geometry after being constructed by the user (or loaded from file)
  int nvolumes = 0;
#ifdef USE_VECGEOM_NAVIGATOR
  LoadVecGeomGeometry(broker);
  vecgeom::GeoManager::Instance().GetAllLogicalVolumes(volumes);
  nvolumes = volumes.size();
  if (!nvolumes) {
    Fatal("GeantVDetectorConstruction::SetupGeometry", "Geometry is empty");
    return 0;
  }  
#else
  (void)broker;
  if (!gGeoManager) {
    Fatal("GeantVDetectorConstruction::SetupGeometry", "ROOT geometry not loaded");
    return 0;
  }
  gGeoManager->SetMaxThreads(fRunMgr->GetNthreadsTotal());
  TObjArray *lvolumes = gGeoManager->GetListOfVolumes();
  nvolumes = lvolumes->GetEntries();
  for (auto ivol = 0; ivol < nvolumes; ivol++)
    volumes.push_back((TGeoVolume *)lvolumes->At(ivol));
#endif

  for (auto i=0; i<nvolumes; ++i) {
    Volume_t *vol = (Volume_t*)volumes[i];
    VBconnector *connector = new VBconnector(i);
#ifdef USE_VECGEOM_NAVIGATOR
    vol->SetBasketManagerPtr(connector);
#else
    vol->SetFWExtension(connector);
#endif
  }
  return nvolumes;
}

//______________________________________________________________________________
bool GeantVDetectorConstruction::LoadGeometry(const char *filename) {
// Load geometry from file. This is callable from the user detector construction class.
#ifdef USE_ROOT
  if (gGeoManager || TGeoManager::Import(filename))
    return true;
#else
  vecgeom::GeoManager *geom = &vecgeom::GeoManager::Instance();
  if (geom) {
    geom->LoadGeometryFromSharedLib(filename);
    return true;
  }
#endif
  Error("GeantVDetectorConstruction::LoadGeometry", "Cannot load geometry from file %s", filename);
  return false;
}

//______________________________________________________________________________
bool GeantVDetectorConstruction::LoadVecGeomGeometry(TaskBroker *broker) {
#ifdef USE_VECGEOM_NAVIGATOR
  if (vecgeom::GeoManager::Instance().GetWorld() == NULL) {
#ifdef USE_ROOT
    vecgeom::RootGeoManager::Instance().SetMaterialConversionHook(CreateMaterialConversion());
    printf("Now loading VecGeom geometry...\n");
    vecgeom::RootGeoManager::Instance().LoadRootGeometry();
    printf("Loading VecGeom geometry done\n");
    std::vector<vecgeom::LogicalVolume *> v1;
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(v1);
    printf("--- Have logical volumes %ld\n", v1.size());
    std::vector<vecgeom::VPlacedVolume *> v2;
    vecgeom::GeoManager::Instance().getAllPlacedVolumes(v2);
    printf("--- Have placed volumes %ld\n", v2.size());
    //    vecgeom::RootGeoManager::Instance().world()->PrintContent();
    // Create regions if any available in the ROOT geometry
    int nregions = GeantVDetectorConstruction::ImportRegions();
    printf("--- Imported %d regions\n", nregions);
#endif
  }
  if (broker) {
    printf("Now uploading VecGeom geometry to Coprocessor(s)...\n");
    return broker->UploadGeometry();
  }
  InitNavigators();
  return true;
#else
  (void)broker;
  return false;
#endif
}

//______________________________________________________________________________
int GeantVDetectorConstruction::ImportRegions() {
// Import regions if available in TGeo
  int nregions = 0;
#if defined(USE_VECGEOM_NAVIGATOR) && defined(USE_ROOT)
  using Region = vecgeom::Region;
  using LogicalVolume = vecgeom::LogicalVolume;
  
  double electronCut, positronCut, gammaCut, protonCut;
  // Loop on ROOT regions if any
  nregions = gGeoManager->GetNregions();
  for (int i=0; i<nregions; ++i) {
    TGeoRegion *region_root = gGeoManager->GetRegion(i);
    electronCut = positronCut = gammaCut = protonCut = -1.;
    // loop cuts for this region
    for (int icut=0; icut < region_root->GetNcuts(); ++icut) {
      std::string cutname = region_root->GetCut(icut)->GetName();
      double cutvalue = region_root->GetCut(icut)->GetCut();
      if (cutname == "gamcut") gammaCut = cutvalue;
      if (cutname == "ecut") electronCut = cutvalue;
      if (cutname == "poscut") positronCut = cutvalue;
      if (cutname == "pcut") protonCut = cutvalue;
    }
    Region *region = new Region(std::string(region_root->GetName()),
                                gammaCut, electronCut, positronCut, protonCut);
    printf("Created region %s with: gammaCut = %g [cm], eleCut = %g [cm], posCut = %g [cm], protonCut = %g [cm]\n",
           region_root->GetName(), gammaCut, electronCut, positronCut, protonCut);
    // loop volumes in the region. Volumes should be already converted to LogicalVolume
    for (int ivol = 0; ivol < region_root->GetNvolumes(); ++ivol) {
      LogicalVolume *vol = vecgeom::RootGeoManager::Instance().Convert(region_root->GetVolume(ivol));
      vol->SetRegion(region);
      // printf("   added to volume %s\n", vol->GetName());
    }
  }
#endif
  return nregions;
}

//______________________________________________________________________________
void GeantVDetectorConstruction::InitNavigators() {
#ifdef USE_VECGEOM_NAVIGATOR
  for (auto &lvol : vecgeom::GeoManager::Instance().GetLogicalVolumesMap()) {
    if (lvol.second->GetDaughtersp()->size() < 4) {
      lvol.second->SetNavigator(vecgeom::NewSimpleNavigator<>::Instance());
    }
    if (lvol.second->GetDaughtersp()->size() >= 5) {
      lvol.second->SetNavigator(vecgeom::SimpleABBoxNavigator<>::Instance());
    }
    if (lvol.second->GetDaughtersp()->size() >= 10) {
      lvol.second->SetNavigator(vecgeom::HybridNavigator<>::Instance());
      vecgeom::HybridManager2::Instance().InitStructure((lvol.second));
    }
    lvol.second->SetLevelLocator(vecgeom::SimpleABBoxLevelLocator::GetInstance());
  }
#endif
}

//______________________________________________________________________________
#ifdef USE_VECGEOM_NAVIGATOR
#ifdef USE_ROOT
std::function<void*(TGeoMaterial const *)> GeantVDetectorConstruction::CreateMaterialConversion() {
  return [](TGeoMaterial const *rootmat) {
      //std::cout<<"     -->  Creating Material  "<<rootmat->GetName();
      int    numElem    = rootmat->GetNelements();
      double density    = rootmat->GetDensity()*geant::g/geant::cm3; // in g/cm3
      const std::string  name = rootmat->GetName();
      // check if it is a G4 NIST material
      std::string postName = "";
      bool isNistMaterial = false;
      if (name.substr(0,3)=="G4_") {
        postName = name.substr(3);
        isNistMaterial = true;
      }
      geantphysics::Material *gmat = nullptr;
      if (isNistMaterial) {
        std::string nistName = "NIST_MAT_"+postName;
        gmat = geantphysics::Material::NISTMaterial(nistName);
      } else {
        // find or create material
        gmat = geantphysics::Material::GetMaterial(name);
        if (gmat) {
          // std::cout<< " Material "<<name << " has already been created.!"<< std::endl;
          return gmat;
        }
        gmat = new geantphysics::Material(name, density, numElem);
        for (int j=0; j<numElem; ++j) {
          double va;
          double vz;
          double vw;
          const_cast<TGeoMaterial *>(rootmat)->GetElementProp(va, vz, vw, j);
          // create NIST element
          geantphysics::Element *elX = geantphysics::Element::NISTElement(vz);
          // add to the Material
          gmat->AddElement(elX, vw);
       }
     }
     // std::cout<< "  geantphysics::name = " << gmat->GetName() << std::endl;
     gmat->SetIsUsed(true);
     return gmat;
   };
}
#endif
#endif

} // GEANT_IMPL_NAMESPACE
} // Geant
