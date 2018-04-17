#ifndef GEANTV_VECSAUTERGAVRILAPHOTOELECTRICMODEL_H
#define GEANTV_VECSAUTERGAVRILAPHOTOELECTRICMODEL_H

#include "Geant/SauterGavrilaPhotoElectricModel.h"
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

class VecSauterGavrilaPhotoElectricModel : public SauterGavrilaPhotoElectricModel {
public:
  VecSauterGavrilaPhotoElectricModel(const std::string &modelname = "SauterGavrilaVec")
      : SauterGavrilaPhotoElectricModel(modelname)
  {
  }

    virtual ~VecSauterGavrilaPhotoElectricModel(){};

  void Initialize() override;

  void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td) override;

  virtual bool IsModelUsable(const MaterialCuts * /*cut*/, double ekin) override;

protected:
  std::vector<LinAliasCached> fAliasTablePerGammaEnergy;

  PhysDI SampleShellAliasVec(PhysDV egamma, PhysDI zed, PhysDV r1, PhysDV r2);
  PhysDV SamplePhotoElectronDirectionAliasVec(PhysDV egamma, PhysDV r1, PhysDV r2, PhysDV r3);
  void SamplePhotoElectronDirectionRejVec(const double *egamma, double *cosTheta, int N, const geant::TaskData *t);

};
}

#endif // GEANTV_VECSAUTERGAVRILAPHOTOELECTRICMODEL_H
