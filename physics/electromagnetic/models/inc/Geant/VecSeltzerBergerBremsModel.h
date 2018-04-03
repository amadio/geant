#ifndef VECSELTZERBERGERBREMSMODELc1_H
#define VECSELTZERBERGERBREMSMODELc1_H

// from geantV
#include "Geant/Config.h"
#include "Geant/SeltzerBergerBremsModel.h"
#include "Geant/AliasTableAlternative.h"

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

namespace geantphysics {

// class Material;
class MaterialCuts;
// class Element;
class AliasTable;
class Particle;
class LightTrack;

class VecSeltzerBergerBremsModel : public SeltzerBergerBremsModel {
public:
  VecSeltzerBergerBremsModel(bool iselectron, const std::string &modelname = "eSeltzerBergerBremsVec");

  virtual ~VecSeltzerBergerBremsModel(){};

  virtual void Initialize();

  void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td) override;

private:
  PhysDV SampleEnergyTransfer(PhysDV gammaCut, PhysDV densityCor, PhysDI mcLocalIdx, double *tableEmin,
                              double *tableILDeta, PhysDV primekin, PhysDV r1, PhysDV r2, PhysDV r3);

  void SampleEnergyTransfer(const double *eEkin, const double *gammaCut, const int *IZet, const double *zet,
                            double *gammaEn, const double *densityCor, int N, const geant::TaskData *td);

  struct AliasDataForMatCut {
    AliasDataForMatCut(int ntables, double lemin, double ildel) : fNData(ntables), fLogEmin(lemin), fILDelta(ildel)
    {
      fAliasData.reserve(ntables);
    }
    int fNData;
    double fLogEmin;
    double fILDelta;
    std::vector<LinAliasCached> fAliasData;
  };

  struct AliasDataForAllMatCuts {
    std::vector<std::unique_ptr<AliasDataForMatCut>> fTablesPerMatCut;
    std::vector<int> fNData;
    std::vector<double> fLogEmin;
    std::vector<double> fILDelta;
  };

private:
  void SamplePhotonDirection(PhysDV elenergy, PhysDV &sinTheta, PhysDV &cosTheta, PhysDV rndm);
  PhysDV PositronCorrection1(PhysDV ekinelectron, PhysDV ephoton, PhysDV gcutener, PhysDV z);
  AliasDataForAllMatCuts fAliasData;
};

} // namespace geantphysics

#endif
