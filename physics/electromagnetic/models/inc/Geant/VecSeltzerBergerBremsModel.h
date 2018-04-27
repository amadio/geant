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

  void Initialize() override;

  void SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td) override;

  bool IsModelUsable(const MaterialCuts *matCut, double ekin) override;

private:
  geant::Double_v SampleEnergyTransfer(geant::Double_v gammaCut, geant::Double_v densityCor, geant::IndexD_v mcLocalIdx,
                                       double *tableEmin, double *tableILDeta, geant::Double_v primekin,
                                       geant::Double_v r1, geant::Double_v r2, geant::Double_v r3);

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
  void SamplePhotonDirection(geant::Double_v elenergy, geant::Double_v &sinTheta, geant::Double_v &cosTheta,
                             geant::Double_v rndm);
  geant::Double_v PositronCorrection1(geant::Double_v ekinelectron, geant::Double_v ephoton, geant::Double_v gcutener,
                                      geant::Double_v z);
  AliasDataForAllMatCuts fAliasData;
};

} // namespace geantphysics

#endif
