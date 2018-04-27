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

  bool IsModelUsable(const MaterialCuts *matCut, double ekin) override;

private:
  geant::Double_v SampleEnergyTransfer(geant::Double_v gammaCut, geant::Double_v densityCor, geant::IndexD_v mcLocalIdx,
                                       double *tableEmin, double *tableILDeta, geant::Double_v primekin,
                                       geant::Double_v r1, geant::Double_v r2, geant::Double_v r3);

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
  void SamplePhotonDirection(geant::Double_v elenergy, geant::Double_v &sinTheta, geant::Double_v &cosTheta,
                             geant::Double_v rndm);
  geant::Double_v PositronCorrection1(geant::Double_v ekinelectron, geant::Double_v ephoton, geant::Double_v gcutener,
                                      geant::Double_v z);
  void GetLPMFunctions(geant::Double_v &lpmGs, geant::Double_v &lpmPhis, const geant::Double_v s);
  void ComputeLPMfunctions(geant::Double_v &funcXiS, geant::Double_v &funcGS, geant::Double_v &funcPhiS,
                           const geant::Double_v lpmenergy, const geant::Double_v egamma, const geant::Double_v etot,
                           const geant::Double_v densitycor, const std::array<int, geant::kVecLenD> izet);
  geant::Double_v ComputeURelDXSecPerAtom(geant::Double_v egamma, geant::Double_v etotal, geant::Double_v lpmenergy,
                                          geant::Double_v densitycor, std::array<int, geant::kVecLenD> izet);
  geant::Double_v ComputeDXSecPerAtom(geant::Double_v egamma, geant::Double_v etotal, geant::Double_v zet);
  void ComputeScreeningFunctions(geant::Double_v &phi1, geant::Double_v &phi1m2, geant::Double_v &xsi1,
                                 geant::Double_v &xsi1m2, const geant::Double_v gamma, const geant::Double_v epsilon);
  AliasDataForAllMatCuts fAliasData;
};

} // namespace geantphysics

#endif
