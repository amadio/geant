#ifndef GEANTV_VECBETHEHEITLERPAIRMODEL_H
#define GEANTV_VECBETHEHEITLERPAIRMODEL_H

#include "Geant/BetheHeitlerPairModel.h"
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

class VecBetheHeitlerPairModel : public BetheHeitlerPairModel {
public:
  VecBetheHeitlerPairModel(const std::string &modelname = "PairBetheHeitlerVec") : BetheHeitlerPairModel(modelname) {}

  virtual ~VecBetheHeitlerPairModel(){};

  void Initialize() override;

  void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td) override;

  bool IsModelUsable(const MaterialCuts *, double ekin) override;

protected:
  struct RatinAliasTablePerElem {
    std::vector<RatinAliasDataTrans> fTablePerEn;
    RatinAliasTablePerElem() : fTablePerEn(0) {}
  };
  std::vector<RatinAliasTablePerElem> fAliasTablesPerZ;

  geant::Double_v SampleTotalEnergyTransferAliasOneShot(const geant::Double_v egamma, const int *izet,
                                                        const geant::Double_v r1, const geant::Double_v r2,
                                                        const geant::Double_v r3);

  void SampleTotalEnergyTransferRejVec(const double *egamma, const int *izet, double *epsOut, int N,
                                       geant::TaskData *td);

  void ScreenFunction12(geant::Double_v &val1, geant::Double_v &val2, const geant::Double_v delta, const bool istsai);
  geant::Double_v ScreenFunction1(const geant::Double_v delta, const bool istsai);
  geant::Double_v ScreenFunction2(const geant::Double_v delta, const bool istsai);
};
}

#endif // GEANTV_VECBETHEHEITLERPAIRMODEL_H
