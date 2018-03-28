#ifndef GEANTV_VECRELATIVISTICPAIRMODEL_H
#define GEANTV_VECRELATIVISTICPAIRMODEL_H

#include "Geant/RelativisticPairModel.h"
#include "Geant/VectorPhysicsTypes.h"
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

  virtual void Initialize();

  virtual void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td);

protected:
  struct RatinAliasTablePerMaterial {
    std::vector<RatinAliasDataTrans> fTablePerEn;
    int fILowestZ;
    RatinAliasTablePerMaterial() : fTablePerEn(0), fILowestZ(200) {}
  };
  std::vector<RatinAliasTablePerMaterial> fAliasTablesPerMaterial;

  PhysDV SampleTotalEnergyTransferAliasOneShot(const PhysDV egamma, const int *matIDX, const PhysDV r1, const PhysDV r2,
                                               const PhysDV r3);

  void SampleTotalEnergyTransferRejVec(const double *egamma, const double *lpmEnergy, const int *izet, double *epsOut,
                                       int N, geant::TaskData *td);

  void ScreenFunction12(PhysDV &val1, PhysDV &val2, const PhysDV delta, const bool istsai);
  PhysDV ScreenFunction1(const PhysDV delta, const bool istsai);
  PhysDV ScreenFunction2(const PhysDV delta, const bool istsai);

  void ComputeScreeningFunctions(PhysDV &phi1, PhysDV &phi2, const PhysDV delta, const bool istsai);

  void ComputeLPMfunctions(PhysDV &funcXiS, PhysDV &funcGS, PhysDV &funcPhiS, PhysDV lpmenergy, PhysDV eps,
                           PhysDV egamma, PhysDV varS1Cond, PhysDV ilVarS1Cond);
};
}

#endif
