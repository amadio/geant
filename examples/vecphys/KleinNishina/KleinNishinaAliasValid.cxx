#include <Geant/VecKleinNishinaComptonModel.h>
#include "KleinNishinaTestCommon.h"

#include "Geant/VecRngWrapper.h"

class KleinNishinaAliasTester : public KleinNishinaComptonModel {
public:
  KleinNishinaAliasTester() : KleinNishinaComptonModel(){};
  using KleinNishinaComptonModel::SampleReducedPhotonEnergy;
};

using geantphysics::VecKleinNishinaComptonModel;
class KleinNishinaAliasVec : public VecKleinNishinaComptonModel {
public:
  KleinNishinaAliasVec() : VecKleinNishinaComptonModel(){};
  using VecKleinNishinaComptonModel::SampleReducedPhotonEnergyVec;
};

const int kTestSize = 1024 * 10;
int main()
{
  KleinNishinaAliasVec *kln = new KleinNishinaAliasVec;
  kln->SetLowEnergyUsageLimit(minEn);
  kln->SetHighEnergyUsageLimit(maxEn);
  kln->SetUseSamplingTables(true);
  kln->Initialize();

  KleinNishinaAliasTester *klnStd = new KleinNishinaAliasTester;
  klnStd->SetLowEnergyUsageLimit(minEn);
  klnStd->SetHighEnergyUsageLimit(maxEn);
  klnStd->SetUseSamplingTables(true);
  klnStd->Initialize();

  geant::VecRngWrapper rng;
  std::vector<double> energy;
  std::vector<double> out1;
  std::vector<double> out2;
  std::vector<double> r1;
  std::vector<double> r2;
  std::vector<double> r3;

  for (int i = 0; i < kTestSize; ++i) {
    energy.push_back(rng.uniform(minEn, maxEn));
    out1.push_back(0.0);
    out2.push_back(0.0);
    r1.push_back(rng.uniform());
    r2.push_back(rng.uniform());
    r3.push_back(rng.uniform());
  }

  for (int i = 0; i < kTestSize; ++i) {
    out1[i] = klnStd->SampleReducedPhotonEnergy(energy[i], r1[i], r2[i], r3[i]);
  }
  kln->SampleReducedPhotonEnergyVec(energy.data(), r1.data(), r2.data(), r3.data(), out2.data(), kTestSize);

  double cumError = 0.0;
  for (int i = 0; i < kTestSize; ++i) {
    cumError += std::abs(out1[i] - out2[i]);
    Printf("TestSize: %d Cumulative error: %f", kTestSize, cumError);
  }

  delete kln;
  delete klnStd;
}
