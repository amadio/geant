#include <benchmark/benchmark.h>
#include <Geant/VecSauterGavrilaPhotoElectricModel.h>
#include "SauterGavrilaTestCommon.h"

#include "Geant/RngWrapper.h"

class SauterGavrilaRejTesterScalar : public SauterGavrilaPhotoElectricModel {
public:
  SauterGavrilaRejTesterScalar() : SauterGavrilaPhotoElectricModel(){};
  using SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirection_Rejection;
};

using geantphysics::VecSauterGavrilaPhotoElectricModel;
class SauterGavrilaRejTesterVector : public VecSauterGavrilaPhotoElectricModel {
public:
  SauterGavrilaRejTesterVector() : VecSauterGavrilaPhotoElectricModel(){};
  using VecSauterGavrilaPhotoElectricModel::SamplePhotoElectronDirectionRejVec;
};

static void SampleAngleRejScalar(benchmark::State &state)
{
  SauterGavrilaRejTesterScalar *sgt = new SauterGavrilaRejTesterScalar;
  sgt->SetLowEnergyUsageLimit(minEn);
  sgt->SetHighEnergyUsageLimit(maxEn);
  sgt->SetUseSamplingTables(false);
  sgt->Initialize();

  geant::RngWrapper rng;
  std::vector<double> energy;
  std::vector<double> cosTheta;

  auto td = PrepareTaskData();

  for (int i = 0; i < kMaxBasket; ++i) {
    energy.push_back(rng.uniform(minEn, maxEn));
    cosTheta.push_back(0.0);
  }

  for (auto _ : state) {
    for (int i = 0; i < state.range(0); ++i) {
      double ct;
      sgt->SamplePhotoElectronDirection_Rejection(energy[i], ct, td.get());
      cosTheta[i] = ct;
    }
  }

  //benchmark::DoNotOptimize(out.data());

    CleanTaskData(td.get());

  delete sgt;
}
BENCHMARK(SampleAngleRejScalar)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void SampleAngleRejVector(benchmark::State &state)
{
  SauterGavrilaRejTesterVector *sgt = new SauterGavrilaRejTesterVector;
  sgt->SetLowEnergyUsageLimit(minEn);
  sgt->SetHighEnergyUsageLimit(maxEn);
  sgt->SetUseSamplingTables(false);
  sgt->Initialize();

  geant::RngWrapper rng;
  std::vector<double> energy;
  std::vector<double> costheta;
  //double costheta[kMaxBasket];

  auto td = PrepareTaskData();

  for (int i = 0; i < kMaxBasket; ++i) {
    energy.push_back(rng.uniform(minEn, maxEn));
    costheta.push_back(0.0);
  }

  for (auto _ : state) {
    sgt->SamplePhotoElectronDirectionRejVec(energy.data(), costheta.data(), state.range(0), td.get());
  }

  //benchmark::DoNotOptimize(out.data());
  CleanTaskData(td.get());
  delete sgt;
}
BENCHMARK(SampleAngleRejVector)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

BENCHMARK_MAIN();
