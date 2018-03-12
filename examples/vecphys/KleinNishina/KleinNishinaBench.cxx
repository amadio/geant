#include "KleinNishinaBench.h"
#include <benchmark/benchmark.h>

#include "Geant/VecRngWrapper.h"

static void KleinNishinaScalarAlias(benchmark::State &state)
{
  KleinNishinaComptonModel *kNish = PrepareKnishinaModel(true);

  auto Td = PrepareTaskData();
  std::vector<LightTrack> primaries;
  PreparePrimaries(primaries);

  for (auto _ : state) {
    int basketSize = state.range(0);
    Td->fPhysicsData->SetNumUsedSecondaries(0);
    for (int i = 0; i < basketSize; ++i) {
      kNish->SampleSecondaries(primaries[i], Td);
    }
  }

  delete kNish;
  CleanTaskData(Td);
}
BENCHMARK(KleinNishinaScalarAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void KleinNishinaScalarRej(benchmark::State &state)
{
  KleinNishinaComptonModel *kNish = PrepareKnishinaModel(false);

  auto Td = PrepareTaskData();
  std::vector<LightTrack> primaries;
  PreparePrimaries(primaries);

  for (auto _ : state) {
    int basketSize = state.range(0);
    Td->fPhysicsData->SetNumUsedSecondaries(0);
    for (int i = 0; i < basketSize; ++i) {
      kNish->SampleSecondaries(primaries[i], Td);
    }
  }

  delete kNish;
  CleanTaskData(Td);
}
BENCHMARK(KleinNishinaScalarRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

BENCHMARK_MAIN();

void PreparePrimaries(std::vector<LightTrack> &output)
{
  geant::VecRngWrapper rng;
  output.clear();
  for (int i = 0; i < kMaxBasket; ++i) {
    LightTrack gamma;
    double phi = geant::units::kTwoPi * rng.uniform(); // NOT uniform on shpere
    double th  = geant::units::kPi * rng.uniform();
    gamma.SetDirX(sin(phi) * cos(th));
    gamma.SetDirY(cos(phi) * cos(th));
    gamma.SetDirZ(sin(th));
    double eKin = minEn + (maxEn - minEn) * rng.uniform();
    gamma.SetKinE(eKin);
    output.push_back(gamma);
  }
}

TaskData *PrepareTaskData()
{
  auto PhysData    = new geantphysics::PhysicsData();
  auto Td          = new TaskData(1, kMaxBasket);
  Td->fPhysicsData = PhysData;
  return Td;
}

void CleanTaskData(TaskData *Td)
{
  delete Td->fPhysicsData;
  delete Td;
}
