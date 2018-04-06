#ifndef VECRELATIVISTICBREMSMODELc1_H
#define VECRELATIVISTICBREMSMODELc1_H

// from geantV
#include "Geant/Config.h"
#include "Geant/RelativisticBremsModel.h"
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

class VecRelativisticBremsModel : public RelativisticBremsModel {
public:
  VecRelativisticBremsModel(const std::string &modelname = "eRelativisticBremsVec");

  virtual ~VecRelativisticBremsModel(){};

  void Initialize() override;

  void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td) override;

  virtual bool IsModelUsable(const MaterialCuts *matCut, double ekin);

private:
  PhysDV SampleEnergyTransfer(PhysDV gammaCut, PhysDV densityCor, PhysDI mcLocalIdx, double *tableEmin,
                              double *tableILDeta, PhysDV primekin, PhysDV r1, PhysDV r2, PhysDV r3);

  void SampleEnergyTransfer(const double *eEkin, const double *gammaCut, const double *zet, const double *densityCor,
                            const double *lpmEnergy, double *gammaEn, int N, const geant::TaskData *td);

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
  void GetLPMFunctions(PhysDV &lpmGs, PhysDV &lpmPhis, const PhysDV s);
  void ComputeLPMfunctions(PhysDV &funcXiS, PhysDV &funcGS, PhysDV &funcPhiS, const PhysDV lpmenergy,
                           const PhysDV egamma, const PhysDV etot, const PhysDV densitycor,
                           const std::array<int, kPhysDVWidth> izet);
  PhysDV ComputeURelDXSecPerAtom(PhysDV egamma, PhysDV etotal, PhysDV lpmenergy, PhysDV densitycor,
                                 std::array<int, kPhysDVWidth> izet);
  PhysDV ComputeDXSecPerAtom(PhysDV egamma, PhysDV etotal, PhysDV zet);
  void ComputeScreeningFunctions(PhysDV &phi1, PhysDV &phi1m2, PhysDV &xsi1, PhysDV &xsi1m2, const PhysDV gamma,
                                 const PhysDV epsilon);
  AliasDataForAllMatCuts fAliasData;
};

} // namespace geantphysics

#endif
