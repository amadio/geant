#include <benchmark/benchmark.h>
#include <Geant/VecSauterGavrilaPhotoElectricModel.h>
#include "SauterGavrilaTestCommon.h"

#include "Geant/RngWrapper.h"
#include <random>

class SauterGavrilaAliasTester : public SauterGavrilaPhotoElectricModel {
public:
  SauterGavrilaAliasTester() : SauterGavrilaPhotoElectricModel(){};
  using SauterGavrilaPhotoElectricModel::SampleShellAlias;
};

using geantphysics::VecSauterGavrilaPhotoElectricModel;
class SauterGavrilaAliasVec : public VecSauterGavrilaPhotoElectricModel {
public:
  SauterGavrilaAliasVec() : VecSauterGavrilaPhotoElectricModel(){};
  using VecSauterGavrilaPhotoElectricModel::SampleShellAliasVec;
};

static void SampleShellAliasScalar(benchmark::State &state)
{
  SauterGavrilaAliasTester *sgt = new SauterGavrilaAliasTester;
  sgt->SetLowEnergyUsageLimit(minEn);
  sgt->SetHighEnergyUsageLimit(maxEn);
  sgt->SetUseSamplingTables(true);
  sgt->Initialize();

  geant::RngWrapper rngGV;
  std::vector<double> energy;
  std::vector<int> zed;
  std::vector<size_t> out;
  std::vector<double> r1;
  std::vector<double> r2;
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rngMT(rd());    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uni(3,98); // guaranteed unbiased


  for (int i = 0; i < kMaxBasket; ++i) {
    energy.push_back(rngGV.uniform(minEn, maxEn));
    auto random_integer = uni(rngMT);
    zed.push_back(random_integer);
    out.push_back(0);
    r1.push_back(rngGV.uniform());
    r2.push_back(rngGV.uniform());
  }

  for (auto _ : state) {
    for (int i = 0; i < state.range(0); ++i) {
       size_t tmp=(size_t) zed[i];
       sgt->SampleShellAlias(energy[i], tmp, r1[i], r2[i], out[i]);
    }
  }

  benchmark::DoNotOptimize(out.data());

  delete sgt;
}
BENCHMARK(SampleShellAliasScalar)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void SampleShellAliasVector(benchmark::State &state)
{
  SauterGavrilaAliasVec *sgt = new SauterGavrilaAliasVec;
  sgt->SetLowEnergyUsageLimit(minEn);
  sgt->SetHighEnergyUsageLimit(maxEn);
  sgt->SetUseSamplingTables(true);
  sgt->Initialize();

  geant::RngWrapper rng;
  std::vector<double> energy;
  std::vector<int> zed;
  std::vector<int> out;
  std::vector<double> r1;
  std::vector<double> r2;
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rngMT(rd());    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uni(3,98); // guaranteed unbiased

  for (int i = 0; i < kMaxBasket; ++i) {
    energy.push_back(rng.uniform(minEn, maxEn));
    auto random_integer = uni(rngMT);
    zed.push_back(random_integer);
    out.push_back(0);
    r1.push_back(rng.uniform());
    r2.push_back(rng.uniform());
  }

  for (auto _ : state) {
      for (int i = 0; i < state.range(0); i += geant::kVecLenD) {
      geant::Double_v en, r1v, r2v;
      geant::IndexD_v z;
      vecCore::Load(en, energy.data() + i);
      vecCore::Load(z, zed.data() + i);
      vecCore::Load(r1v, r1.data() + i);
      vecCore::Load(r2v, r2.data() + i);
      geant::IndexD_v eps = sgt->SampleShellAliasVec(en, z, r1v, r2v);
      vecCore::Store(eps, out.data() + i);
    }
  }

  benchmark::DoNotOptimize(out.data());

  delete sgt;
}
BENCHMARK(SampleShellAliasVector)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

BENCHMARK_MAIN();
