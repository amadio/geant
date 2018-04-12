
#ifndef VECMOLLERBHABHAIONIZATIONMODELc1_H
#define VECMOLLERBHABHAIONIZATIONMODELc1_H

// from geantV
#include "Geant/Config.h"
#include "Geant/MollerBhabhaIonizationModel.h"
#include "Geant/AliasTableAlternative.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {
class TaskData;
}
}

namespace geantphysics {
inline namespace GEANT_IMPL_NAMESPACE {
class Material;
class Element;
}
}

#include <string>

namespace geantphysics {

// class Material;
class MaterialCuts;
// class Element;
class AliasTable;
class Particle;
class LightTrack;

class VecMollerBhabhaIonizationModel : public MollerBhabhaIonizationModel {
public:
  VecMollerBhabhaIonizationModel(bool iselectron, const std::string &modelname = "ieioniMollerBhabhaVec");

  virtual ~VecMollerBhabhaIonizationModel(){};

  void Initialize() override;

  void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td) override;

  bool IsModelUsable(const MaterialCuts *matCut, double ekin) override;
  ;

private:
  PhysDV SampleEnergyTransfer(PhysDV elProdCut, PhysDI mcLocalIdx, double *tableEmin, double *tableILDeta,
                              PhysDV primekin, PhysDV r1, PhysDV r2, PhysDV r3);

  void SampleEnergyTransfer(const double *elProdCut, const double *primekin, double *epsOut, int N,
                            const geant::TaskData *td);

  struct AliasDataForMatCut {
    AliasDataForMatCut(int ntables, double lemin, double ildel) : fNData(ntables), fLogEmin(lemin), fILDelta(ildel)
    {
      fAliasData.reserve(ntables);
    }
    int fNData;
    double fLogEmin;
    double fILDelta;
    std::vector<LinAliasCached> fAliasData;
  };

  struct AliasDataForAllMatCuts {
    std::vector<std::unique_ptr<AliasDataForMatCut>> fTablesPerMatCut;
    std::vector<int> fNData;
    std::vector<double> fLogEmin;
    std::vector<double> fILDelta;
  };

private:
  AliasDataForAllMatCuts fAliasData;
};

} // namespace geantphysics

#endif
