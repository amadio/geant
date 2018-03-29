#ifndef GEANTV_VECKLEINNISHINACOMPTONMODEL_H
#define GEANTV_VECKLEINNISHINACOMPTONMODEL_H

#include "Geant/KleinNishinaComptonModel.h"
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

class VecKleinNishinaComptonModel : public KleinNishinaComptonModel {
public:
  VecKleinNishinaComptonModel(const std::string &modelname = "ComptonKleinNishinaVec")
      : KleinNishinaComptonModel(modelname)
  {
  }

  ~VecKleinNishinaComptonModel() override{};

  void Initialize() override;

  void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td) override;

protected:
  std::vector<LinAliasCached> fAliasTablePerGammaEnergy;

  PhysDV SampleReducedPhotonEnergyVec(PhysDV egamma, PhysDV r1, PhysDV r2, PhysDV r3);

  void SampleReducedPhotonEnergyRej(const double *egamma, double *onemcost, double *sint2, double *eps, int N,
                                    const geant::TaskData *td);
};
}

#endif // GEANTV_VECKLEINNISHINACOMPTONMODEL_H
