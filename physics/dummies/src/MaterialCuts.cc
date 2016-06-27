
#include "MaterialCuts.h"
#include "Material.h"

#include "SystemOfUnits.h"

namespace geant {

std::vector<MaterialCuts*> MaterialCuts::gTheMaterialCutsTable;

MaterialCuts::MaterialCuts(Material *mat) : fMaterial(mat) {
  fProductionCutsInLength[0] = 1.0*geant::mm;  // for gamma in lenght
  fProductionCutsInLength[1] = 1.0*geant::mm;  // for e- in lenght
  fProductionCutsInLength[2] = 1.0*geant::mm;  // for e+ in lenght

  fProductionCutsInEnergy[0] = 1.0*geant::keV;  // for gamma in internal energy units
  fProductionCutsInEnergy[1] = 1.0*geant::keV;  // for e- in internal energy units
  fProductionCutsInEnergy[2] = 1.0*geant::keV;  // for e+ in internal energy units

  fIndex = gTheMaterialCutsTable.size();
  gTheMaterialCutsTable.push_back(this);
}

MaterialCuts::MaterialCuts(Material *mat, double gcutlength, double emcutlength, double epcutlength,
                           double gcutenergy, double emcutenergy , double epcutenergy)
                           : fMaterial(mat) {
  fProductionCutsInLength[0] = gcutlength;  // for gamma in lenght
  fProductionCutsInLength[1] = emcutlength; // for e- in lenght
  fProductionCutsInLength[2] = epcutlength; // for e+ in lenght

  fProductionCutsInEnergy[0] = gcutenergy;  // for gamma in internal energy units
  fProductionCutsInEnergy[1] = emcutenergy; // for e- in internal energy units
  fProductionCutsInEnergy[2] = epcutenergy; // for e+ in internal energy units

  fIndex = gTheMaterialCutsTable.size();
  gTheMaterialCutsTable.push_back(this);
}

}
