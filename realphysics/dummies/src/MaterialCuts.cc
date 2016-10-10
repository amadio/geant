
#include "MaterialCuts.h"

#include "Material.h"

#include "CutConverter.h"
#include "CutConverterForGamma.h"
#include "CutConverterForElectron.h"
#include "CutConverterForPositron.h"

#include "PhysicsParameters.h"

#include "SystemOfUnits.h"

// dummy
#include "Region.h"

#include <iostream>


namespace geantphysics {

std::vector<MaterialCuts*> MaterialCuts::gTheMaterialCutsTable;
std::vector< std::vector<MaterialCuts*> > MaterialCuts::gTheMaterialCutsPerRegion;
std::vector< std::vector<MaterialCuts*> > MaterialCuts::gTheMaterialCutsPerRegionPerMaterial;

/*
// NOTE: we might drop this Ctr and make the other private i.e. let to create MaterialCuts only by CreateAll method.
MaterialCuts::MaterialCuts(int regionindx, const Material *mat, bool iscutinlength)
: fRegionIndex(regionindx), fMaterial(mat), fIsProductionCutsGivenInLength(iscutinlength) {
  // set default values both for length and energy;
  // depending on the fIsProductionCutInLength flag, one set will be converted during ConvertAll
  fProductionCutsInLength[0] = PhysicsParameters::GetDefaultGammaCutInLength();    // for gamma in lenght
  fProductionCutsInLength[1] = PhysicsParameters::GetDefaultElectronCutInLength(); // for e- in lenght
  fProductionCutsInLength[2] = PhysicsParameters::GetDefaultPositronCutInLength(); // for e+ in lenght

  fProductionCutsInEnergy[0] = PhysicsParameters::GetDefaultGammaCutInEnergy();    // for gamma in internal energy units
  fProductionCutsInEnergy[1] = PhysicsParameters::GetDefaultElectronCutInEnergy(); // for e- in internal energy units
  fProductionCutsInEnergy[2] = PhysicsParameters::GetDefaultPositronCutInEnergy(); // for e+ in internal energy units

  fIndex = gTheMaterialCutsTable.size();
  // add the current MaterialCuts pointer to the global material cuts table
  gTheMaterialCutsTable.push_back(this);

  // it will be set only once to the number of regions
  int numRegions = Region::GetTheRegionTable().size();
  if (gTheMaterialCutsPerRegion.size()<numRegions) {
    gTheMaterialCutsPerRegion.resize(numRegions);
  }
  if (gTheMaterialCutsPerRegionPerMaterial.size()<numRegions) {
    gTheMaterialCutsPerRegionPerMaterial.resize(numRegions);
    int numMaterials = Material::GetTheMaterialTable().size();
    for (int i=0; i<numRegions; ++i) {
      gTheMaterialCutsPerRegionPerMaterial[i].resize(numMaterials,nullptr);
    }
  }
  // add the current MatrialCuts pointer to the list per region
  gTheMaterialCutsPerRegion[regionindx].push_back(this);
  // add the current MatrialCuts pointer to the list per region per material
  gTheMaterialCutsPerRegionPerMaterial[regionindx][mat->GetIndex()] = this;
}
*/

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
  unsigned long numRegions = Region::GetTheRegionTable().size();
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


void MaterialCuts::CleanUp() {
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

// create all MaterialCuts by using the Resgion table; this will be the standard way of automatically creating all
// MaterialCuts in the detetector
void MaterialCuts::CreateAll() {
  // get the global region table
  const std::vector<Region*> theRegions = Region::GetTheRegionTable();
  int   numRegions = theRegions.size();
  // loop over the regions
  for (int ir=0; ir<numRegions; ++ir) {
    // get the cuts, the flag if cut is given in length or energy;
    bool   iscutinlength = theRegions[ir]->IsProductionCutInLength();
    double gcut  = theRegions[ir]->GetGammaCut();
    double emcut = theRegions[ir]->GetElectronCut();
    double epcut = theRegions[ir]->GetPositronCut();
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

    // list of materials that are in this region; NOTE: assume that there might be duplications
    std::vector<const Material*> matList    = theRegions[ir]->GetMaterialList();
    std::vector<const Material*> difMatList;
    for (unsigned long i=0; i<matList.size(); ++i) {
      const Material *curMat = matList[i];
      bool isthere = false;
      for (unsigned long j=0; j<difMatList.size() && !isthere; ++j) {
        if (difMatList[j]==curMat) {
          isthere = true;
        }
      }
      if (!isthere) {
        difMatList.push_back(curMat);
      }
    }
    // create all different MatrialCuts
    for (unsigned long i=0; i<difMatList.size(); ++i) {
      new MaterialCuts(theRegions[ir]->GetIndex(), difMatList[i], iscutinlength, gcut, emcut, epcut);
    }
  }

  // convert production cuts in lenght/energy to energy/lenght
  ConvertAll();
}

// printouts
std::ostream& operator<<(std::ostream& flux, const MaterialCuts* matcut) {
  //std::ios::fmtflags mode = flux.flags();
  flux.setf(std::ios::fixed,std::ios::floatfield);
  //long prec = flux.precision(6);
  const double *cutslenght = matcut->GetProductionCutsInLength();
  const double *cutsenergy = matcut->GetProductionCutsInEnergy();
  const std::vector<Region*> vregions = Region::GetTheRegionTable();
  int   regionIndx             = matcut->GetRegionIndex();
  const std::string regionName = vregions[regionIndx]->GetName();
  std::string str = " in length.";
  if (!matcut->fIsProductionCutsGivenInLength)
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
