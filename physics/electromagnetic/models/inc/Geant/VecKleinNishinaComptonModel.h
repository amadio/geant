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

  //  virtual ~VecKleinNishinaComptonModel();

  virtual void Initialize();

  virtual void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td);

protected:
  std::vector<LinAliasCached> fCachedAliasTable;
  void SampleSecondariesVectorAlias(LightTrack_v &tracks, geant::TaskData *td);

  void SampleReducedPhotonEnergyVec(const double *egamma, const double *r1, const double *r2, const double *r3,
                                    double *out, int N);
};
}

#endif // GEANTV_VECKLEINNISHINACOMPTONMODEL_H
