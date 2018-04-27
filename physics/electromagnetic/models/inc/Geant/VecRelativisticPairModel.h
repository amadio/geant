#ifndef GEANTV_VECRELATIVISTICPAIRMODEL_H
#define GEANTV_VECRELATIVISTICPAIRMODEL_H

#include "Geant/RelativisticPairModel.h"
#include "Geant/VectorTypes.h"
#include "Geant/AliasTableAlternative.h"

#include "Geant/Config.h"

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
#include <vector>

namespace geantphysics {

class MaterialCuts;
// class Element;
class AliasTable;
class Particle;
class LightTrack;

class VecRelativisticPairModel : public RelativisticPairModel {
public:
  VecRelativisticPairModel(const std::string &modelname = "PairRelativisticLPMVec") : RelativisticPairModel(modelname)
  {
  }

  virtual ~VecRelativisticPairModel(){};

  void Initialize() override;

  void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td) override;

  bool IsModelUsable(const MaterialCuts *, double ekin) override
  {
    return ekin < GetHighEnergyUsageLimit() && ekin > GetLowEnergyUsageLimit();
  };

protected:
  struct RatinAliasTablePerMaterial {
    std::vector<RatinAliasDataTrans> fTablePerEn;
    int fILowestZ;
    RatinAliasTablePerMaterial() : fTablePerEn(0), fILowestZ(200) {}
  };
  std::vector<RatinAliasTablePerMaterial> fAliasTablesPerMaterial;

  geant::Double_v SampleTotalEnergyTransferAliasOneShot(const geant::Double_v egamma, const int *matIDX,
                                                        const geant::Double_v r1, const geant::Double_v r2,
                                                        const geant::Double_v r3);

  void SampleTotalEnergyTransferRejVec(const double *egamma, const double *lpmEnergy, const int *izet, double *epsOut,
                                       int N, geant::TaskData *td);

  void ScreenFunction12(geant::Double_v &val1, geant::Double_v &val2, const geant::Double_v delta, const bool istsai);
  geant::Double_v ScreenFunction1(const geant::Double_v delta, const bool istsai);
  geant::Double_v ScreenFunction2(const geant::Double_v delta, const bool istsai);

  void ComputeScreeningFunctions(geant::Double_v &phi1, geant::Double_v &phi2, const geant::Double_v delta,
                                 const bool istsai);

  void ComputeLPMfunctions(geant::Double_v &funcXiS, geant::Double_v &funcGS, geant::Double_v &funcPhiS,
                           geant::Double_v lpmenergy, geant::Double_v eps, geant::Double_v egamma,
                           geant::Double_v varS1Cond, geant::Double_v ilVarS1Cond);
};
}

#endif
