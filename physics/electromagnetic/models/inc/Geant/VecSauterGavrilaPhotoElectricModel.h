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

using geant::Double_v;
using geant::Int_v;
using geant::IndexD_v;
using geant::kVecLenD;
using geant::MaskD_v;
using vecCore::Get;
using vecCore::Set;
using vecCore::AssignMaskLane;
using vecCore::MaskFull;
using vecCore::MaskEmpty;
using vecCore::Gather;
using vecCore::MaskedAssign;
    
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

  void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td) ;

  virtual bool IsModelUsable(const MaterialCuts * /*cut*/, double ekin) override;

protected:
  std::vector<LinAliasCached> fAliasTablePerGammaEnergy;

  Int_v SampleShellAliasVec(Double_v egamma, Int_v zed, Double_v r1, Double_v r2);
  void SampleShellVec(double *egamma, int * zed, int* ss, int N, const geant::TaskData *td, double* randoms);
  Double_v SamplePhotoElectronDirectionAliasVec(Double_v egamma, Double_v r1, Double_v r2, Double_v r3);
  void SamplePhotoElectronDirectionRejVec(const double *egamma, double *cosTheta, int N, const geant::TaskData *t);

};
}

#endif // GEANTV_VECSAUTERGAVRILAPHOTOELECTRICMODEL_H
