#ifndef GEANTV_VECBETHEHEITLERPAIRMODEL_H
#define GEANTV_VECBETHEHEITLERPAIRMODEL_H

#include "Geant/BetheHeitlerPairModel.h"
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

class VecBetheHeitlerPairModel : public BetheHeitlerPairModel {
public:
  VecBetheHeitlerPairModel(const std::string &modelname = "PairBetheHeitlerVec") : BetheHeitlerPairModel(modelname) {}

  virtual ~VecBetheHeitlerPairModel(){};

  void Initialize() override;

  void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td) override;

protected:
  struct RatinAliasTablePerElem {
    std::vector<RatinAliasDataTrans> fTablePerEn;
    RatinAliasTablePerElem() : fTablePerEn(0) {}
  };
  std::vector<RatinAliasTablePerElem> fAliasTablesPerZ;

  PhysDV SampleTotalEnergyTransferAliasOneShot(const PhysDV egamma, const int *izet, const PhysDV r1, const PhysDV r2,
                                               const PhysDV r3);

  void SampleTotalEnergyTransferRejVec(const double *egamma, const int *izet, double *epsOut, int N,
                                       geant::TaskData *td);

  void ScreenFunction12(PhysDV &val1, PhysDV &val2, const PhysDV delta, const bool istsai);
  PhysDV ScreenFunction1(const PhysDV delta, const bool istsai);
  PhysDV ScreenFunction2(const PhysDV delta, const bool istsai);
};
}

#endif // GEANTV_VECBETHEHEITLERPAIRMODEL_H
