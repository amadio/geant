
#include "LambdaTable.h"

#include "PhysicsProcess.h"
#include "PhysicsParameters.h"
#include "Particle.h"

#include "Types.h"
#include "Material.h"
#include "MaterialCuts.h"

#include "Spline.h"

#include <iostream>

namespace geantphysics {

LambdaTable::LambdaTable(const PhysicsProcess *process, bool ispermaterial)
: fProcess(process), fIsLambdaTablesPerMaterial(ispermaterial) {
  // these will be set properly when the tables are built
  fIsSpecialLambdaTableBinNum    = false;
  fNumLambdaTableBins            = 3;
  fNumSpecialLambdaTableBins     = 3;
  fMinLambdaTableEnergy          = 100.*geant::eV;
  fMaxLambdaTableEnergy          = 100.*geant::TeV;
  fLogMinLambdaTableEnergy       = 1.;
  fEnergyILDelta                 = 1.;
}

LambdaTable::~LambdaTable() {
  std::cerr<< "  deleting lambda tables for process = " << fProcess->GetName() << std::endl;
  ClearAllTables();
}


void LambdaTable::ClearAllTables() {
  if (fIsLambdaTablesPerMaterial) {
    for (size_t i=0; i<fLambdaTablesPerMaterial.size(); ++i) {
      if (fLambdaTablesPerMaterial[i]) {
        fLambdaTablesPerMaterial[i]->fOneLambdaTable.clear();
        delete fLambdaTablesPerMaterial[i]->fSpline;
        delete fLambdaTablesPerMaterial[i];
      }
    }
    fLambdaTablesPerMaterial.clear();
    fEnergyGrid.clear();
  } else {
    for (size_t i=0; i<fLambdaTablesPerMaterialCuts.size(); ++i) {
      if (fLambdaTablesPerMaterialCuts[i]) {
        fLambdaTablesPerMaterialCuts[i]->fEnergyGrid.clear();
        fLambdaTablesPerMaterialCuts[i]->fLambdaTable.clear();
        delete fLambdaTablesPerMaterialCuts[i]->fSpline;
        delete fLambdaTablesPerMaterialCuts[i];
      }
    }
    fLambdaTablesPerMaterialCuts.clear();
  }
}


double LambdaTable::GetMacroscopicXSection(const MaterialCuts *matcut, double ekin) {
  double macXsec = 0.;
  if (fIsLambdaTablesPerMaterial) {
    // return with zero if below or above the min/max lambda table energy
    if (ekin>=fEnergyGrid[0] && ekin<=fEnergyGrid[fNumLambdaTableBins-1]) {
      double logE     = std::log(ekin);
      int    lowEIndx = (int) ((logE-fLogMinLambdaTableEnergy)*fEnergyILDelta);
      if (lowEIndx>=fNumLambdaTableBins-1) --lowEIndx;
      // we might put it under verbose build since
      // protection against very small numerical uncertainties
  //      if (lowEIndx>0 && kinenergy<fEnergyGrid[lowEIndx]) {
  //        --lowEIndx;
  //      } else if (kinenergy>fEnergyGrid[lowEIndx+1]) {
  //        ++lowEIndx;
  //      }
      macXsec = fLambdaTablesPerMaterial[matcut->GetMaterial()->GetIndex()]->fSpline->GetValueAt(ekin,lowEIndx);
      if (macXsec<0.0) {
        macXsec = 0.0;
      }
    }
  } else {
    // return with zero if below or above the min/max lambda table energy
    struct LambdaTableForAMaterialCuts *data = fLambdaTablesPerMaterialCuts[matcut->GetIndex()];
    if (ekin>=data->fEnergyGrid[0] && ekin<=data->fEnergyGrid[data->fNumLambdaTableBins-1]) {
      double logE     = std::log(ekin);
      int    lowEIndx = (int) ((logE-data->fLogMinLambdaTableEnergy)*data->fEnergyILDelta);
      if (lowEIndx>=data->fNumLambdaTableBins-1) --lowEIndx;
      // we might put it under verbose build since
      // protection against very small numerical uncertainties
  //      if (lowEIndx>0 && kinenergy<fEnergyGrid[lowEIndx]) {
  //        --lowEIndx;
  //      } else if (kinenergy>fEnergyGrid[lowEIndx+1]) {
  //        ++lowEIndx;
  //      }
      macXsec = data->fSpline->GetValueAt(ekin,lowEIndx);
      if (macXsec<0.0) {
        macXsec = 0.0;
      }
    }
  }
  return macXsec;
}


double  LambdaTable::GetMacroscopicXSectionMaximumEnergy(const MaterialCuts *matcut) const {
  if (fIsLambdaTablesPerMaterial) {
    return fLambdaTablesPerMaterial[matcut->GetMaterial()->GetIndex()]->fLambdaMaxEnergy;
  } else {
    return fLambdaTablesPerMaterialCuts[matcut->GetIndex()]->fLambdaMaxEnergy;
  }
}


double  LambdaTable::GetMacroscopicXSectionMaximum(const MaterialCuts *matcut) const {
  if (fIsLambdaTablesPerMaterial) {
    return fLambdaTablesPerMaterial[matcut->GetMaterial()->GetIndex()]->fLambdaMax;
  } else {
    return fLambdaTablesPerMaterialCuts[matcut->GetIndex()]->fLambdaMax;
  }
}


void LambdaTable::BuildLambdaTables() {
  ClearAllTables();
  // tables per material
  if (fIsLambdaTablesPerMaterial) {
    std::cerr<< "      = Building lambda tables per material for process = " << fProcess->GetName() <<std::endl;
    const std::vector<MaterialCuts*> &matCutTable = MaterialCuts::GetTheMaterialCutsTable();
    const Vector_t<Material*>        &matTable    = Material::GetTheMaterialTable();
    fLambdaTablesPerMaterial.resize(matTable.size(),nullptr);
    for (size_t imatc=0; imatc<matCutTable.size(); ++imatc) {
      const MaterialCuts *matCut = matCutTable[imatc];
      if (!fProcess->IsActiveRegion(matCut->GetRegionIndex()) || !matCut->GetMaterial()->IsUsed()) {
        continue;
      }
      int matIndx = matCut->GetMaterial()->GetIndex();
      if (!fLambdaTablesPerMaterial[matIndx]) { // build lambda table for this material
        if (fEnergyGrid.size()<1) {  // build the energy grid: only once
          GenerateEnergyGrid(matCut);
        }
//        std::cerr << "        ===> Mat = " << matCut->GetMaterial()->GetName()
//                  << "   emin  = " << fMinLambdaTableEnergy/geant::MeV
//                  << "   emax  = " << fMaxLambdaTableEnergy/geant::MeV
//                  << "   bins  = " << fNumLambdaTableBins << std::endl;
        fLambdaTablesPerMaterial[matIndx]                   = new ALambdaTable();
        fLambdaTablesPerMaterial[matIndx]->fLambdaMax       = -1.0;
        fLambdaTablesPerMaterial[matIndx]->fLambdaMaxEnergy = -1.0;
        fLambdaTablesPerMaterial[matIndx]->fOneLambdaTable.resize(fNumLambdaTableBins,0.);
        for (int iener=0; iener<fNumLambdaTableBins; ++iener) {
          double ekin    = fEnergyGrid[iener];
//          if (iener==0) {
//            ekin += geant::eV;
//          }
          // dynamic mass of the particle is not considered !
          double macXsec = fProcess->ComputeMacroscopicXSection(matCut, ekin, fProcess->GetParticle(), 0.);
          if (macXsec<0.0) {
            macXsec = 0.0;
          }
          if (macXsec>fLambdaTablesPerMaterial[matIndx]->fLambdaMax) {
            fLambdaTablesPerMaterial[matIndx]->fLambdaMax       = macXsec;
            fLambdaTablesPerMaterial[matIndx]->fLambdaMaxEnergy = ekin;
          }
          fLambdaTablesPerMaterial[matIndx]->fOneLambdaTable[iener] = macXsec;
        }

//std::cerr<< " particle = " << fProcess->GetParticle()->GetName() << " proc = "<< fProcess->GetName()<< "  Sigma_max E (MeV) = " <<  fLambdaTablesPerMaterial[matIndx]->fLambdaMaxEnergy/geant::MeV
//<< " Sigma_Max (1/mm) = " <<  fLambdaTablesPerMaterial[matIndx]->fLambdaMax*geant::mm << std::endl;

        fLambdaTablesPerMaterial[matIndx]->fSpline =
            new Spline(&fEnergyGrid[0], &(fLambdaTablesPerMaterial[matIndx]->fOneLambdaTable[0]), fNumLambdaTableBins);
      }
    }
  } else { // tables per material-cuts
    const std::vector<MaterialCuts*> &matCutTable = MaterialCuts::GetTheMaterialCutsTable();
    std::cerr<< "      = Building lambda tables per material-cuts for process = " << fProcess->GetName() <<std::endl;
    fLambdaTablesPerMaterialCuts.resize(matCutTable.size(),nullptr);
    for (size_t imatc=0; imatc<matCutTable.size(); ++imatc) {
      const MaterialCuts *matCut = matCutTable[imatc];
      if (!fProcess->IsActiveRegion(matCut->GetRegionIndex())) {
        continue;
      }
      int matCutIndx = matCut->GetIndex();
      if (!fLambdaTablesPerMaterialCuts[matCutIndx]) { // build lambda table for this material-cuts
        fLambdaTablesPerMaterialCuts[matCutIndx]                   = new LambdaTableForAMaterialCuts();
        fLambdaTablesPerMaterialCuts[matCutIndx]->fLambdaMax       = -1.0;
        fLambdaTablesPerMaterialCuts[matCutIndx]->fLambdaMaxEnergy = -1.0;
        GenerateEnergyGrid(matCut, fLambdaTablesPerMaterialCuts[matCutIndx]);
//        std::cerr << "        ===> MatCut = "<< matCut->GetMaterial()->GetName()
//                  << "   emin  = " << fLambdaTablesPerMaterialCuts[matCutIndx]->fMinLambdaTableEnergy/geant::MeV
//                  << "   emax  = " << fLambdaTablesPerMaterialCuts[matCutIndx]->fMaxLambdaTableEnergy/geant::MeV
//                  << "   bins  = " << fLambdaTablesPerMaterialCuts[matCutIndx]->fNumLambdaTableBins << std::endl;
        for (int iener=0; iener<fLambdaTablesPerMaterialCuts[matCutIndx]->fNumLambdaTableBins; ++iener) {
          double ekin    = fLambdaTablesPerMaterialCuts[matCutIndx]->fEnergyGrid[iener];
//          if (iener==0) {
//            ekin += geant::eV;
//          }
          // dynamic mass of the particle is not considered
          double macXsec = fProcess->ComputeMacroscopicXSection(matCut, ekin, fProcess->GetParticle(), 0.);
          if (macXsec<0.0) {
            macXsec = 0.0;
          }
          if (macXsec>fLambdaTablesPerMaterialCuts[matCutIndx]->fLambdaMax) {
            fLambdaTablesPerMaterialCuts[matCutIndx]->fLambdaMax       = macXsec;
            fLambdaTablesPerMaterialCuts[matCutIndx]->fLambdaMaxEnergy = ekin;
          }
          fLambdaTablesPerMaterialCuts[matCutIndx]->fLambdaTable[iener] = macXsec;
        }
        fLambdaTablesPerMaterialCuts[matCutIndx]->fSpline =
            new Spline(&(fLambdaTablesPerMaterialCuts[matCutIndx]->fEnergyGrid[0]),
                       &(fLambdaTablesPerMaterialCuts[matCutIndx]->fLambdaTable[0]),
                       fLambdaTablesPerMaterialCuts[matCutIndx]->fNumLambdaTableBins);
//std::cerr<< " particle = " << fProcess->GetParticle()->GetName() << " proc = "<< fProcess->GetName()<< "  Sigma_max E (MeV) = " <<  fLambdaTablesPerMaterialCuts[matCutIndx]->fLambdaMaxEnergy/geant::MeV
//<< " Sigma_Max (1/mm) = " <<  fLambdaTablesPerMaterialCuts[matCutIndx]->fLambdaMax*geant::mm << std::endl;
      }
    }
  }
}


void LambdaTable::GenerateEnergyGrid(const MaterialCuts *matcut, struct LambdaTableForAMaterialCuts *data) {
  if (!data) {
    // per material
    // generate the kinetic energy grid common for each material
    double emin = fProcess->GetPhysicsParameters()->GetMinLambdaTableEnergy();
    double emax = fProcess->GetPhysicsParameters()->GetMaxLambdaTableEnergy();
    double nbin = fProcess->GetPhysicsParameters()->GetNumLambdaTableBins();
    if (fIsSpecialLambdaTableBinNum) {
      nbin = fNumSpecialLambdaTableBins;
    }
    double scl  = std::log(emax/emin);
//    std::cerr<< "   ++++ emin = " << emin/geant::MeV << "  emax = " <<emax/geant::MeV << " nbin = " << nbin << " scale = "<< scl << std::endl;
    fMinLambdaTableEnergy = std::max(fProcess->GetMinimumLambdaTableKineticEnergy(matcut, fProcess->GetParticle()), emin);
    fMaxLambdaTableEnergy = emax;
    if (fMinLambdaTableEnergy>=fMaxLambdaTableEnergy) {
      std::cerr<< "   **** ERROR:: LambdaTable::GenerateEnergyGrid() \n"
               << "   fMinLambdaTableEnergy = " << fMinLambdaTableEnergy
               << " >= fMaxLambdaTableEnergy = " <<fMaxLambdaTableEnergy
               << " for process = " << fProcess->GetName()
               << std::endl;
      exit(-1);
    }
    fNumLambdaTableBins = std::lrint(nbin*std::log(fMaxLambdaTableEnergy/fMinLambdaTableEnergy)/scl)+1;
                          //fProcess->GetPhysicsParameters()->GetNumLambdaTableBinsPerDecade()
                          //*std::lrint(std::log10(fMaxLambdaTableEnergy/fMinLambdaTableEnergy))+1;
    if (fNumLambdaTableBins<3) {
      fNumLambdaTableBins = 3;
    }
    fEnergyGrid.resize(fNumLambdaTableBins);
    fLogMinLambdaTableEnergy           = std::log(fMinLambdaTableEnergy);
    double delta                       = std::log(fMaxLambdaTableEnergy/fMinLambdaTableEnergy)/(fNumLambdaTableBins-1.0);
    fEnergyILDelta                     = 1.0/delta;
    fEnergyGrid[0]                     = fMinLambdaTableEnergy;
    fEnergyGrid[fNumLambdaTableBins-1] = fMaxLambdaTableEnergy;
    for (int i=1; i<fNumLambdaTableBins-1; ++i) {
      fEnergyGrid[i] = std::exp(fLogMinLambdaTableEnergy+i*delta);
    }
  } else {
    // generate the kinetic energy grid common for each material
    double emin = fProcess->GetPhysicsParameters()->GetMinLambdaTableEnergy();
    double emax = fProcess->GetPhysicsParameters()->GetMaxLambdaTableEnergy();
    double nbin = fProcess->GetPhysicsParameters()->GetNumLambdaTableBins();
    if (fIsSpecialLambdaTableBinNum) {
      nbin = fNumSpecialLambdaTableBins;
    }
    double scl  = std::log(emax/emin);
//    std::cerr<< "   ++++ emin = " << emin/geant::MeV << "  emax = " <<emax/geant::MeV << " nbin = " << nbin << " scale = "<< scl << std::endl;
    data->fMinLambdaTableEnergy = std::max(fProcess->GetMinimumLambdaTableKineticEnergy(matcut, fProcess->GetParticle()), emin);
    data->fMaxLambdaTableEnergy = emax;
    if (data->fMinLambdaTableEnergy>=data->fMaxLambdaTableEnergy) {
      std::cerr<< "   **** ERROR:: LambdaTable::GenerateEnergyGrid() \n"
               << "   fMinLambdaTableEnergy = " << data->fMinLambdaTableEnergy
               << " >= fMaxLambdaTableEnergy = " <<data->fMaxLambdaTableEnergy
               << " for process = " << fProcess->GetName()
               << std::endl;
      exit(-1);
    }
    data->fNumLambdaTableBins = std::lrint(nbin*std::log(data->fMaxLambdaTableEnergy/data->fMinLambdaTableEnergy)/scl)+1;
//                                fProcess->GetPhysicsParameters()->GetNumLambdaTableBinsPerDecade()
//                                *std::lrint(std::log10(data->fMaxLambdaTableEnergy/data->fMinLambdaTableEnergy))+1;
    if (data->fNumLambdaTableBins<3) {
      data->fNumLambdaTableBins = 3;
    }
    data->fEnergyGrid.resize(data->fNumLambdaTableBins);
    data->fLogMinLambdaTableEnergy                 = std::log(data->fMinLambdaTableEnergy);
    double delta                                   = std::log(data->fMaxLambdaTableEnergy/data->fMinLambdaTableEnergy)
                                                     /(data->fNumLambdaTableBins-1.0);
    data->fEnergyILDelta                           = 1.0/delta;
    data->fEnergyGrid[0]                           = data->fMinLambdaTableEnergy;
    data->fEnergyGrid[data->fNumLambdaTableBins-1] = data->fMaxLambdaTableEnergy;
    for (int i=1; i<data->fNumLambdaTableBins-1; ++i) {
      data->fEnergyGrid[i] = std::exp(data->fLogMinLambdaTableEnergy+i*delta);
    }
    data->fLambdaTable.resize(data->fNumLambdaTableBins,0.0);
  }
}


}  // namespace geantphysics
