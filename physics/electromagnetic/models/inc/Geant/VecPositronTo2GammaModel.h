#ifndef GEANTV_VECPOSITRONTWOGAMMAMODEL_H
#define GEANTV_VECPOSITRONTWOGAMMAMODEL_H

#include "Geant/PositronTo2GammaModel.h"
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

class VecPositronTo2GammaModel : public PositronTo2GammaModel {
public:
  VecPositronTo2GammaModel(const std::string &modelname = "e2GammaAnnihVec") : PositronTo2GammaModel(modelname) {}

  virtual ~VecPositronTo2GammaModel(){};

  virtual void Initialize();

  virtual void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td);

protected:
  std::vector<LinAliasCached> fCachedAliasTable;

  PhysDV SampleEnergyTransferAlias(PhysDV pekin, PhysDV r1, PhysDV r2, PhysDV r3, PhysDV gamma);

  void SampleEnergyTransferRej(const double *gamma, double *eps, int N, const geant::TaskData *td);
};
}

#endif
