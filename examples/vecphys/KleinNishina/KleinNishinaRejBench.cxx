#include <benchmark/benchmark.h>
#include <Geant/VecKleinNishinaComptonModel.h>
#include "KleinNishinaTestCommon.h"

#include "Geant/VecRngWrapper.h"

class KleinNishinaRejTesterScalar : public KleinNishinaComptonModel {
public:
  KleinNishinaRejTesterScalar() : KleinNishinaComptonModel(){};
  using KleinNishinaComptonModel::SampleReducedPhotonEnergy;
};

using geantphysics::VecKleinNishinaComptonModel;
class KleinNishinaRejTesterVector : public VecKleinNishinaComptonModel {
public:
  KleinNishinaRejTesterVector() : VecKleinNishinaComptonModel(){};
  using VecKleinNishinaComptonModel::SampleReducedPhotonEnergyRej;
};

static void SampleRejScalar(benchmark::State &state)
{
  KleinNishinaRejTesterScalar *kln = new KleinNishinaRejTesterScalar;
  kln->SetLowEnergyUsageLimit(minEn);
  kln->SetHighEnergyUsageLimit(maxEn);
  kln->SetUseSamplingTables(false);
  kln->Initialize();

  geant::VecRngWrapper rng;
  std::vector<double> energy;
  std::vector<double> out;
  std::vector<double> oneMinCost;
  std::vector<double> sin2t;

  TaskData *td = PrepareTaskData();

  for (int i = 0; i < kMaxBasket; ++i) {
    energy.push_back(rng.uniform(minEn, maxEn));
    out.push_back(0.0);
    oneMinCost.push_back(0.0);
    sin2t.push_back(0.0);
  }

  for (auto _ : state) {
    for (int i = 0; i < state.range(0); ++i) {
      double tmp1, tmp2;
      out[i]        = kln->SampleReducedPhotonEnergy(energy[i], tmp1, tmp2, td);
      oneMinCost[i] = tmp1;
      sin2t[i]      = tmp2;
    }
  }

  benchmark::DoNotOptimize(out.data());

  CleanTaskData(td);

  delete kln;
}
BENCHMARK(SampleRejScalar)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void SampleRejVector(benchmark::State &state)
{
  KleinNishinaRejTesterVector *kln = new KleinNishinaRejTesterVector;
  kln->SetLowEnergyUsageLimit(minEn);
  kln->SetHighEnergyUsageLimit(maxEn);
  kln->SetUseSamplingTables(false);
  kln->Initialize();

  geant::VecRngWrapper rng;
  std::vector<double> energy;
  std::vector<double> out;
  std::vector<double> oneMinCost;
  std::vector<double> sin2t;

  TaskData *td = PrepareTaskData();

  for (int i = 0; i < kMaxBasket; ++i) {
    energy.push_back(rng.uniform(minEn, maxEn));
    out.push_back(0.0);
    oneMinCost.push_back(0.0);
    sin2t.push_back(0.0);
  }

  for (auto _ : state) {
    kln->SampleReducedPhotonEnergyRej(energy.data(), oneMinCost.data(), sin2t.data(), out.data(), state.range(0), td);
  }

  benchmark::DoNotOptimize(out.data());
  CleanTaskData(td);

  delete kln;
}
BENCHMARK(SampleRejVector)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

BENCHMARK_MAIN();
