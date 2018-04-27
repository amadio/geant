#ifndef GEANTV_VECKLEINNISHINACOMPTONMODEL_H
#define GEANTV_VECKLEINNISHINACOMPTONMODEL_H

#include "Geant/KleinNishinaComptonModel.h"
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

class VecKleinNishinaComptonModel : public KleinNishinaComptonModel {
public:
  VecKleinNishinaComptonModel(const std::string &modelname = "ComptonKleinNishinaVec")
      : KleinNishinaComptonModel(modelname)
  {
  }

  virtual ~VecKleinNishinaComptonModel(){};

  void Initialize() override;

  void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td) override;

  bool IsModelUsable(const MaterialCuts * /*cut*/, double ekin) override;

protected:
  std::vector<LinAliasCached> fAliasTablePerGammaEnergy;

  geant::Double_v SampleReducedPhotonEnergyVec(geant::Double_v egamma, geant::Double_v r1, geant::Double_v r2,
                                               geant::Double_v r3);

  void SampleReducedPhotonEnergyRej(const double *egamma, double *onemcost, double *sint2, double *eps, int N,
                                    const geant::TaskData *td);
};
}

#endif // GEANTV_VECKLEINNISHINACOMPTONMODEL_H
