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

  virtual void Initialize();

  virtual void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td);

protected:
  struct RatinAliasTablePerElem {
    std::vector<RatinAliasDataTrans> fTablePerEn;
    RatinAliasTablePerElem() : fTablePerEn(0) {}
  };
  std::vector<RatinAliasTablePerElem> fAliasTablesPerZ;

  void SampleTotalEnergyTransferAliasVec(const double *egamma, const int *izet, const double *r1, const double *r2,
                                         const double *r3, int N, double *eps);

  PhysDV SampleTotalEnergyTransferAliasOneShot(const PhysDV egamma, const int *izet, const PhysDV r1, const PhysDV r2,
                                               const PhysDV r3);

  void SampleTotalEnergyTransferRejVec(const double *egamma, const int *izet, geant::TaskData *td);
};
}

#endif // GEANTV_VECBETHEHEITLERPAIRMODEL_H
