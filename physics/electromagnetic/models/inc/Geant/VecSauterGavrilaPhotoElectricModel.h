#ifndef GEANTV_VECSAUTERGAVRILAPHOTOELECTRICMODEL_H
#define GEANTV_VECSAUTERGAVRILAPHOTOELECTRICMODEL_H

#include "Geant/SauterGavrilaPhotoElectricModel.h"
#include "Geant/VectorTypes.h"
#include "Geant/AliasTableAlternative.h"

#include "Geant/Config.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {
class TaskData;
}
} // namespace geant

namespace geantphysics {
inline namespace GEANT_IMPL_NAMESPACE {
class Material;
class Element;
} // namespace GEANT_IMPL_NAMESPACE
} // namespace geantphysics

#include <string>
#include <vector>

namespace geantphysics {

using geant::Double_v;
using geant::IndexD_v;
using geant::kVecLenD;
using geant::MaskD_v;
using MaskDI_v = vecCore::Mask<IndexD_v>;
using vecCore::AssignMaskLane;
using vecCore::Gather;
using vecCore::Get;
using vecCore::MaskedAssign;
using vecCore::MaskEmpty;
using vecCore::MaskFull;
using vecCore::Set;

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

  void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td);

  virtual bool IsModelUsable(const MaterialCuts * /*cut*/, double ekin) override;

protected:
  std::vector<LinAliasCached> fAliasTablePerGammaEnergy;

  IndexD_v SampleShellAliasVec(Double_v egamma, IndexD_v zed, Double_v r1, Double_v r2);
  void SampleShellVec(double *egamma, int *zed, int *ss, int N, const geant::TaskData *td, double *randoms);
  Double_v SamplePhotoElectronDirectionAliasVec(Double_v egamma, Double_v r1, Double_v r2, Double_v r3);
  void SamplePhotoElectronDirectionRejVec(const double *egamma, double *cosTheta, int N, const geant::TaskData *t);
};
} // namespace geantphysics

#endif // GEANTV_VECSAUTERGAVRILAPHOTOELECTRICMODEL_H
