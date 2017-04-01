
#include "MaterialCuts.h"

#include "Material.h"

#include "CutConverter.h"
#include "CutConverterForGamma.h"
#include "CutConverterForElectron.h"
#include "CutConverterForPositron.h"

#include "PhysicsParameters.h"

#include "SystemOfUnits.h"

#include "Region.h"

#include <iostream>

// vecgeom GeoManager, LogicalVolume, VPlacedVolume, Material, Medium
#include "management/GeoManager.h"
#include "volumes/PlacedVolume.h"
#include "volumes/LogicalVolume.h"
#include "materials/Material.h"
#include "materials/Medium.h"

namespace geantphysics {

std::vector<MaterialCuts*> MaterialCuts::gTheMaterialCutsTable;
std::vector< std::vector<MaterialCuts*> > MaterialCuts::gTheMaterialCutsPerRegion;
std::vector< std::vector<MaterialCuts*> > MaterialCuts::gTheMaterialCutsPerRegionPerMaterial;

MaterialCuts::MaterialCuts(int regionindx, const Material *mat, bool iscutinlength, double gcut, double emcut, double epcut)
: fRegionIndex(regionindx), fIsProductionCutsGivenInLength(iscutinlength), fMaterial(mat) {
  if (iscutinlength) {
    fProductionCutsInLength[0] = gcut;  // for gamma in lenght
    fProductionCutsInLength[1] = emcut; // for e- in lenght
    fProductionCutsInLength[2] = epcut; // for e+ in lenght
    // they will be set by the converter
    fProductionCutsInEnergy[0] = -1.0;  // for gamma in internal energy units
    fProductionCutsInEnergy[1] = -1.0; // for e- in internal energy units
    fProductionCutsInEnergy[2] = -1.0; // for e+ in internal energy units
  } else {
    // they will be set by the converter
    fProductionCutsInLength[0] = -1.0;  // for gamma in lenght
    fProductionCutsInLength[1] = -1.0; // for e- in lenght
    fProductionCutsInLength[2] = -1.0; // for e+ in lenght

    fProductionCutsInEnergy[0] = gcut;  // for gamma in internal energy units
    fProductionCutsInEnergy[1] = emcut; // for e- in internal energy units
    fProductionCutsInEnergy[2] = epcut; // for e+ in internal energy units
  }

  fIndex = gTheMaterialCutsTable.size();
  // add the current MaterialCuts pointer to the global material cuts table
  gTheMaterialCutsTable.push_back(this);

  // it will be set only once to the number of regions
  unsigned long numRegions = vecgeom::Region::GetTheRegionTable().size();
  if (gTheMaterialCutsPerRegion.size()<numRegions) {
    gTheMaterialCutsPerRegion.resize(numRegions);
  }
  if (gTheMaterialCutsPerRegionPerMaterial.size()<numRegions) {
    gTheMaterialCutsPerRegionPerMaterial.resize(numRegions);
    int numMaterials = Material::GetTheMaterialTable().size();
    for (unsigned long i=0; i<numRegions; ++i) {
      gTheMaterialCutsPerRegionPerMaterial[i].resize(numMaterials,nullptr);
    }
  }
  // add the current MatrialCuts pointer to the list per region
  gTheMaterialCutsPerRegion[regionindx].push_back(this);
  // add the current MatrialCuts pointer to the list per region per material
  gTheMaterialCutsPerRegionPerMaterial[regionindx][mat->GetIndex()] = this;
}


void MaterialCuts::ClearAll() {
  for (unsigned long i=0; i<gTheMaterialCutsTable.size(); ++i) {
    delete gTheMaterialCutsTable[i];
  }
  gTheMaterialCutsTable.clear();
  for (unsigned long i=0; i<gTheMaterialCutsPerRegion.size(); ++i) {
    gTheMaterialCutsPerRegion[i].clear();
  }
  gTheMaterialCutsPerRegion.clear();
  for (unsigned long i=0; i<gTheMaterialCutsPerRegionPerMaterial.size(); ++i) {
    gTheMaterialCutsPerRegionPerMaterial[i].clear();
  }
  gTheMaterialCutsPerRegionPerMaterial.clear();
}

const MaterialCuts* MaterialCuts::GetMaterialCut(int indx) {
  if (indx>-1 && indx<gTheMaterialCutsTable.size()) {
    return gTheMaterialCutsTable[indx];
  } else {
    std::cerr << "  ***  ERROR:  MaterialCuts::GetMaterialCut() \n"
              << "        Requested MaterialCuts by index = " << indx << " cannot be find in the table! \n"
              << "        Index should be:   0 <= index < number of MaterialCuts = " << gTheMaterialCutsTable.size()
              << std::endl;
    exit(-1);
  }
}

// we should get rid of this
const MaterialCuts* MaterialCuts::GetMaterialCut(int regionindx, int materialindx) {
  if (regionindx>-1 && regionindx<gTheMaterialCutsPerRegionPerMaterial.size()
      && materialindx>-1 && materialindx<gTheMaterialCutsPerRegionPerMaterial[0].size()) {
    return gTheMaterialCutsPerRegionPerMaterial[regionindx][materialindx];
  } else {
    if (!(regionindx>-1 && regionindx<gTheMaterialCutsPerRegionPerMaterial.size())) {
      std::cerr << "  ***  ERROR:  MaterialCuts::GetMaterialCut() \n"
                << "        Requested MaterialCuts by Region and Material indices = "
                << regionindx << ", " << materialindx << " cannot be find in the table!"
                << "\n        Region index should be  :   0 <= index < region dimension of the table = "
                << gTheMaterialCutsPerRegionPerMaterial.size()
                << std::endl;
     } else {
       std::cerr << "  ***  ERROR:  MaterialCuts::GetMaterialCut() \n"
                 << "        Requested MaterialCuts by Region and Material indices = "
                 << regionindx << ", " << materialindx << " cannot be find in the table!"
                 << "\n        Region index should be  :   0 <= index < region dimension of the table   = "
                 << gTheMaterialCutsPerRegionPerMaterial.size()
                 << "\n        Material index should be:   0 <= index < material dimension of the table = "
                 << gTheMaterialCutsPerRegionPerMaterial[0].size()
                 << std::endl;
     }
    exit(-1);
  }
}

void MaterialCuts::ConvertAll() {
  CutConverter **converters = new CutConverter*[3];
  // get the min/max secondary production values; they are the same in each reagion
  converters[0] = new CutConverterForGamma(301, PhysicsParameters::GetMinAllowedGammaCutEnergy(), PhysicsParameters::GetMaxAllowedGammaCutEnergy());
  converters[1] = new CutConverterForElectron(301, PhysicsParameters::GetMinAllowedElectronCutEnergy(), PhysicsParameters::GetMaxAllowedElectronCutEnergy());
  converters[2] = new CutConverterForPositron(301, PhysicsParameters::GetMinAllowedPositronCutEnergy(), PhysicsParameters::GetMaxAllowedPositronCutEnergy());

  for (unsigned long i=0; i<gTheMaterialCutsTable.size(); ++i) {
    MaterialCuts   *matcut  = gTheMaterialCutsTable[i];
    const Material *mat     = matcut->GetMaterial();
    if (matcut->fIsProductionCutsGivenInLength) {
      for (int j=0; j<3; ++j) {
        const double *cuts = matcut->GetProductionCutsInLength();
        matcut->SetProductionCutEnergy(j,converters[j]->Convert(mat,cuts[j],true));
      }
    } else {
      for (int j=0; j<3; ++j) {
        const double *cuts      = matcut->GetProductionCutsInEnergy();
        matcut->SetProductionCutLength(j,converters[j]->Convert(mat,cuts[j],false));
      }
    }
  }
  delete converters[0];
  delete converters[1];
  delete converters[2];
  delete [] converters;
}

// create all MaterialCuts by using the Region table; this will be the standard way of automatically creating all
// MaterialCuts in the detetector
void MaterialCuts::CreateAll() {
  // clear all if there were any created before
  ClearAll();
  // get number of regions
  int   numRegions = vecgeom::Region::GetNumberOfRegions();
  // get the world LogicalVolume from vecgeom: this const_cast is very dirty !!
  vecgeom::LogicalVolume *logicWorld = const_cast<vecgeom::LogicalVolume *>(vecgeom::GeoManager::Instance().GetWorld()->GetLogicalVolume());
  // if there are no any regions created by the user create one, with default production cut values,
  // set to the world logical volume (that should set to all logical volumes)
  if (numRegions<1) {
    // create one region with default production cuts and assign all logical volumes to this
    vecgeom::Region *reg0 = new vecgeom::Region("DefaultWorld");
    logicWorld->SetRegion(reg0);
  }
  // it is still possible however, that the user defined at least one region but the world logical volume has no region;
  // we create a default world region and assign the world logical volume to that;
  // later all logical volumes, without region, will also assigned to this region
  if (!logicWorld->GetRegion()) {
    vecgeom::Region *reg0 = new vecgeom::Region("DefaultWorld");
    logicWorld->SetRegion(reg0,false); // set only to this logical volume
  }
  // keep the world region in hand
  vecgeom::Region *regionWorld = logicWorld->GetRegion();

  // loop over all logical volumes, take the associated region and add the material of the given
  // logical volume to the given region.
  // NOTE: we created (in PhysicsProcessHandler::BuildMaterials()) our geantphysics::Materials in the same order that
  // the vecgeom::Materials so their global index will be the same.
  std::vector<vecgeom::LogicalVolume*> theLogicVolumes;
  vecgeom::GeoManager::Instance().GetAllLogicalVolumes(theLogicVolumes);
  // get our Material table
  const VectorHelper<Material*>::Vector_t theMaterialTable = Material::GetTheMaterialTable();
  for (size_t i=0; i<theLogicVolumes.size(); ++i) {
    vecgeom::Material *vgMat = ((vecgeom::Medium*)theLogicVolumes[i]->GetTrackingMediumPtr())->GetMaterial();
    vecgeom::Region *reg     = theLogicVolumes[i]->GetRegion();
    // region must be set for each logical volumes
    // check it and assign all logical volumes without region to the world region and give warning
    if (!reg) {
      theLogicVolumes[i]->SetRegion(regionWorld,false); // set only for this logical volume
      reg = regionWorld;
      std::cerr << "  ***  WARNING:  MaterialCuts::CreateAll() \n"
                << "        LogicalVolume with name = " << theLogicVolumes[i]->GetLabel() << " was not associated to \n"
                << "        any Regions. It has been assigned to the world region ( " << regionWorld->GetName() << ")\n"
                << std::endl;
    }
    // get our material that corresponds to the vecgeom::Material vgMat
    // check if we have already created a MaterialCuts with this material in the current region
    // create a new if not yet
    Material *mat = theMaterialTable[vgMat->GetIndex()];
    // it returns with a pointer to the MaterialCuts that will be set into the logical volume later!
    //MaterialCuts *matCut = CheckMaterialForRegion(reg, mat);
    CheckMaterialForRegion(reg, mat);
  }
  // convert production cuts in lenght/energy to energy/lenght
  ConvertAll();
}


MaterialCuts* MaterialCuts::CheckMaterialForRegion(const vecgeom::Region *region, const Material *mat) {
  // get the region index
  int indxReg = region->GetIndex();
  int indxMat = mat->GetIndex();
  MaterialCuts *matCut = nullptr;
  if (indxReg>-1 && indxReg<gTheMaterialCutsPerRegionPerMaterial.size()
      && indxMat>-1 && indxMat<gTheMaterialCutsPerRegionPerMaterial[0].size()) {
    matCut = gTheMaterialCutsPerRegionPerMaterial[indxReg][indxMat];
  }
  // create it if it is still nullptr (i.e. has not been created yet)
  if (!matCut) {
    bool   iscutinlength = region->IsProductionCutInLength();
    double gcut  = region->GetGammaCut();
    double emcut = region->GetElectronCut();
    double epcut = region->GetPositronCut();
    // check if the cut was given by user or we need to take default values
    if (iscutinlength) {
      if (gcut<0.0)
        gcut  = PhysicsParameters::GetDefaultGammaCutInLength();
      if (emcut<0.0)
        emcut = PhysicsParameters::GetDefaultElectronCutInLength();
      if (epcut<0.0)
        epcut = PhysicsParameters::GetDefaultPositronCutInLength();
    } else {
      if (gcut<0.0)
        gcut  = PhysicsParameters::GetDefaultGammaCutInEnergy();
      if (emcut<0.0)
        emcut = PhysicsParameters::GetDefaultElectronCutInEnergy();
      if (epcut<0.0)
        epcut = PhysicsParameters::GetDefaultPositronCutInEnergy();
    }
    // create the MaterialCuts object
    matCut = new MaterialCuts(indxReg, mat, iscutinlength, gcut, emcut, epcut);
  }
  return matCut;
}

// printouts
std::ostream& operator<<(std::ostream& flux, const MaterialCuts* matcut) {
  //std::ios::fmtflags mode = flux.flags();
  flux.setf(std::ios::fixed,std::ios::floatfield);
  //long prec = flux.precision(6);
  const double *cutslenght = matcut->GetProductionCutsInLength();
  const double *cutsenergy = matcut->GetProductionCutsInEnergy();
  const std::vector<vecgeom::Region*> vregions = vecgeom::Region::GetTheRegionTable();
  int   regionIndx             = matcut->GetRegionIndex();
  const std::string regionName = vregions[regionIndx]->GetName();
  std::string str = " in length.";
  if (!matcut->IsProductionCutsGivenInLength())
    str = " in energy.";

  flux << "\n   Material name      : " << matcut->GetMaterial()->GetName();
  flux << "\n   Belongs to region  : " << "Region name = " << vregions[regionIndx]->GetName()
                                       << ", Region index = " << regionIndx;
  flux << "\n   Material-cut index : " << matcut->GetIndex() << ", Production cuts were given" << str;
  flux << "\n   Production cuts    : " << "     gamma         e-            e+"
       << "\n    In length         : "
       << cutslenght[0]/geant::mm << " [mm]  "
       << cutslenght[1]/geant::mm << " [mm]  "
       << cutslenght[2]/geant::mm << " [mm]  "
       << "\n    In energy         : "
       << cutsenergy[0]/geant::MeV << " [MeV] "
       << cutsenergy[1]/geant::MeV << " [MeV] "
       << cutsenergy[2]/geant::MeV << " [MeV] "
       << "\n";

  return flux;
}


std::ostream& operator<<(std::ostream& flux, const MaterialCuts& matcut)
{
  flux << &matcut;
  return flux;
}


std::ostream& operator<<(std::ostream& flux, std::vector<MaterialCuts*> matcuttable)
{
  //Dump info for all material-cut couples
  flux << "\n  ========= Table: table of registered material-cuts ================="
       << "\n   Number of material-cuts in the table = " << matcuttable.size()
       << "\n";

  for (unsigned long i = 0; i < matcuttable.size(); ++i) {
    flux << matcuttable[i];
    if (i<matcuttable.size()-1)
      flux << "  --------------------------------------------------------------------\n";
  }
  flux << "  ====================================================================\n\n";

  return flux;
}


}  // namespace geantphysics
