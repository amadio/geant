#include <benchmark/benchmark.h>
#include <Geant/VecSauterGavrilaPhotoElectricModel.h>
#include "SauterGavrilaTestCommon.h"

#include "Geant/RngWrapper.h"
#include <random>

class SauterGavrilaRejShellTester : public SauterGavrilaPhotoElectricModel {
public:
  SauterGavrilaRejShellTester() : SauterGavrilaPhotoElectricModel(){};
  using SauterGavrilaPhotoElectricModel::SampleShell;
};

using geantphysics::VecSauterGavrilaPhotoElectricModel;
class SauterGavrilaRejShellVec : public VecSauterGavrilaPhotoElectricModel {
public:
  SauterGavrilaRejShellVec() : VecSauterGavrilaPhotoElectricModel(){};
  using VecSauterGavrilaPhotoElectricModel::SampleShellVec;
};

static void SampleRejShellScalar(benchmark::State &state)
{
  SauterGavrilaRejShellTester *sgt = new SauterGavrilaRejShellTester;
  sgt->SetLowEnergyUsageLimit(minEn);
  sgt->SetHighEnergyUsageLimit(maxEn);
  sgt->SetUseSamplingTables(true);
  sgt->Initialize();

  geant::RngWrapper rngGV;
  std::vector<double> energy;
  std::vector<int> zed;
  std::vector<size_t> out;
  std::vector<double> r1;
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rngMT(rd());    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uni(3,98); // guaranteed unbiased


  for (int i = 0; i < kMaxBasket; ++i) {
    energy.push_back(rngGV.uniform(minEn, maxEn));
    auto random_integer = uni(rngMT);
    zed.push_back(random_integer);
    out.push_back(0);
    r1.push_back(rngGV.uniform());
  }

  for (auto _ : state) {
    for (int i = 0; i < state.range(0); ++i) {
       sgt->SampleShell(energy[i], zed[i], r1[i], out[i]);
    }
  }

  benchmark::DoNotOptimize(out.data());

  delete sgt;
}
BENCHMARK(SampleRejShellScalar)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void SampleRejShellVector(benchmark::State &state)
{
  SauterGavrilaRejShellVec *sgt = new SauterGavrilaRejShellVec;
  sgt->SetLowEnergyUsageLimit(minEn);
  sgt->SetHighEnergyUsageLimit(maxEn);
  sgt->SetUseSamplingTables(true);
  sgt->Initialize();

  geant::RngWrapper rng;
  std::vector<double> energy;
  std::vector<int> zed;
  std::vector<int> out;
  std::vector<double> r1;
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rngMT(rd());    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uni(3,98); // guaranteed unbiased

  for (int i = 0; i < kMaxBasket; ++i) {
    energy.push_back(rng.uniform(minEn, maxEn));
    auto random_integer = uni(rngMT);
    zed.push_back(random_integer);
    out.push_back(0);
    r1.push_back(rng.uniform());
  }
  auto td = PrepareTaskData();
  for (auto _ : state) {
      sgt->SampleShellVec(energy.data(), zed.data(), out.data(), state.range(0), td.get(), r1.data());
    }
  
  benchmark::DoNotOptimize(out.data());

  delete sgt;
}
BENCHMARK(SampleRejShellVector)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

BENCHMARK_MAIN();
