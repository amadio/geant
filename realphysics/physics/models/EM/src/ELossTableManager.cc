
#include "ELossTableManager.h"

#include "PhysicsParameters.h"
#include "ELossTableRegister.h"
#include "ELossTable.h"
#include "Material.h"
#include "MaterialCuts.h"
#include "Particle.h"

#include <iostream>

namespace geantphysics {

ELossTableManager& ELossTableManager::Instance() {
  static ELossTableManager instance;
  return instance;
}


void ELossTableManager::BuildELossTables() {
  Clear(); // will delete all ELossTable and clear all vectors;
  const std::vector<MaterialCuts*> matCutTable = MaterialCuts::GetTheMaterialCutsTable();
  fElossTablePerMaterialCuts.resize(matCutTable.size(),nullptr);
  // set up and build each ELossTable
  for (unsigned long i=0; i<PhysicsParameters::GetThePhysicsParametersTable().size(); ++i) {
    // each PhysicsParameters will create an ELossTable object that will handle only those
    // MaterialCuts and Material-s that belong to regions where the PhysicsParameters object is active
    ELossTable  *lossTable = new ELossTable(PhysicsParameters::GetThePhysicsParametersTable()[i]);
    lossTable->BuildELossTable(fElossTablePerMaterialCuts);
    fTheElossTables.push_back(lossTable);
  }
}


void ELossTableManager::Clear() {
  for (unsigned long i=0; i<fTheElossTables.size(); ++i) {
    delete fTheElossTables[i];
  }
  fTheElossTables.clear();
  fElossTablePerMaterialCuts.clear();
}


double ELossTableManager::GetRestrictedDEDX(const MaterialCuts *matcut, const Particle *part, double kinenergy) {
  double dedx = 0.0;
  // should work properly without this check; we need to put it under verbose build
  if (fElossTablePerMaterialCuts[matcut->GetIndex()]) {
    dedx = fElossTablePerMaterialCuts[matcut->GetIndex()]->GetRestrictedDEDX(matcut->GetIndex(),part->GetInternalCode(), kinenergy);
  } else { // there is no any ELossTable active in the region where the MatrialCut belongs to; should never happen
    std::cerr<<"  ====  ELossTableManager:  No ELossTable for MaterialCuts: \n";
    std::cerr<<matcut<<std::endl;
  }
  return dedx;
}


double ELossTableManager::GetRestrictedRange(const MaterialCuts *matcut, const Particle *part, double kinenergy) {
  double range = 1.0e+20;
  // should work properly without this check; we need to put it under verbose build
  if (fElossTablePerMaterialCuts[matcut->GetIndex()]) {
    range = fElossTablePerMaterialCuts[matcut->GetIndex()]->GetRestrictedRange(matcut->GetIndex(),part->GetInternalCode(), kinenergy);
  } else { // there is no any ELossTable active in the region where the MatrialCut belongs to; should never happen
    std::cerr<<"  ====  ELossTableManager:  No ELossTable for MaterialCuts: \n";
    std::cerr<<matcut<<std::endl;
  }
  return range;
}


double ELossTableManager::GetEnergyForRestrictedRange(const MaterialCuts *matcut, const Particle *part, double range) {
  double energy = 0.0;
  // should work properly without this check; we need to put it under verbose build
  if (fElossTablePerMaterialCuts[matcut->GetIndex()]) {
    energy = fElossTablePerMaterialCuts[matcut->GetIndex()]->GetEnergyForRestrictedRange(matcut->GetIndex(),part->GetInternalCode(), range);
  } else { // there is no any ELossTable active in the region where the MatrialCut belongs to; should never happen
    std::cerr<<"  ====  ELossTableManager:  No ELossTable for MaterialCuts: \n";
    std::cerr<<matcut<<std::endl;
  }
  return energy;
}


double ELossTableManager::GetRange(const MaterialCuts *matcut, const Particle *part, double kinenergy) {
  double range = 1.0e+20;
  // should work properly without this check; we need to put it under verbose build
  if (fElossTablePerMaterialCuts[matcut->GetIndex()]) {
    range = fElossTablePerMaterialCuts[matcut->GetIndex()]->GetRange(matcut->GetMaterial()->GetIndex(),part->GetInternalCode(), kinenergy);
  } else { // there is no any ELossTable active in the region where the MatrialCut belongs to; should never happen
    std::cerr<<"  ====  ELossTableManager:  No ELossTable for MaterialCuts: \n";
    std::cerr<<matcut<<std::endl;
  }
  return range;
}


} // namespace geant
