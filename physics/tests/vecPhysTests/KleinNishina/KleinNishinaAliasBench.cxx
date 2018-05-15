#include <benchmark/benchmark.h>
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

static void SampleAliasScalar(benchmark::State &state)
{
  KleinNishinaAliasTester *kln = new KleinNishinaAliasTester;
  kln->SetLowEnergyUsageLimit(minEn);
  kln->SetHighEnergyUsageLimit(maxEn);
  kln->SetUseSamplingTables(true);
  kln->Initialize();

  geant::VecRngWrapper rng;
  std::vector<double> energy;
  std::vector<double> out;
  std::vector<double> r1;
  std::vector<double> r2;
  std::vector<double> r3;

  for (int i = 0; i < kMaxBasket; ++i) {
    energy.push_back(rng.uniform(minEn, maxEn));
    out.push_back(0.0);
    r1.push_back(rng.uniform());
    r2.push_back(rng.uniform());
    r3.push_back(rng.uniform());
  }

  for (auto _ : state) {
    for (int i = 0; i < state.range(0); ++i) {
      out[i] = kln->SampleReducedPhotonEnergy(energy[i], r1[i], r2[i], r3[i]);
    }
  }

  benchmark::DoNotOptimize(out.data());

  delete kln;
}
BENCHMARK(SampleAliasScalar)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void SampleAliasVector(benchmark::State &state)
{
  KleinNishinaAliasVec *kln = new KleinNishinaAliasVec;
  kln->SetLowEnergyUsageLimit(minEn);
  kln->SetHighEnergyUsageLimit(maxEn);
  kln->SetUseSamplingTables(true);
  kln->Initialize();

  geant::VecRngWrapper rng;
  std::vector<double> energy;
  std::vector<double> out;
  std::vector<double> r1;
  std::vector<double> r2;
  std::vector<double> r3;

  for (int i = 0; i < kMaxBasket; ++i) {
    energy.push_back(rng.uniform(minEn, maxEn));
    out.push_back(0.0);
    r1.push_back(rng.uniform());
    r2.push_back(rng.uniform());
    r3.push_back(rng.uniform());
  }

  for (auto _ : state) {
    for (int i = 0; i < state.range(0); i += kPhysDVWidth) {
      PhysDV en, r1v, r2v, r3v;
      vecCore::Load(en, energy.data() + i);
      vecCore::Load(r1v, r1.data() + i);
      vecCore::Load(r2v, r2.data() + i);
      vecCore::Load(r3v, r3.data() + i);
      PhysDV eps = kln->SampleReducedPhotonEnergyVec(en, r1v, r2v, r3v);
      vecCore::Store(eps, out.data() + i);
    }
  }

  benchmark::DoNotOptimize(out.data());

  delete kln;
}
BENCHMARK(SampleAliasVector)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

BENCHMARK_MAIN();
