#include "KleinNishinaBench.h"
#include <benchmark/benchmark.h>

#include "Geant/VecRngWrapper.h"

static void KleinNishinaScalarAlias(benchmark::State &state)
{
  KleinNishinaComptonModel *kNish = PrepareKnishinaModel(true);

  auto Td = PrepareTaskData();
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries,basketSize);
    state.ResumeTiming();
    Td->fPhysicsData->ClearSecondaries();
    for (int i = 0; i < basketSize; ++i) {
      kNish->SampleSecondaries(primaries[i], Td);
    }
  }

  delete kNish;
  CleanTaskData(Td);
}
BENCHMARK(KleinNishinaScalarAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void KleinNishinaVectorAlias(benchmark::State &state)
{
  VecKleinNishinaComptonModel *kNish = PrepareVecKnishinaModel(true);

  auto Td = PrepareTaskData();

  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries,basketSize);
    primaries.SetNtracks(basketSize);
    state.ResumeTiming();

    Td->fPhysicsData->GetSecondarySOA().ClearTracks();
    kNish->SampleSecondariesVector(primaries, Td);
  }

  delete kNish;
  CleanTaskData(Td);
}
BENCHMARK(KleinNishinaVectorAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void KleinNishinaScalarRej(benchmark::State &state)
{
  KleinNishinaComptonModel *kNish = PrepareKnishinaModel(false);

  auto Td = PrepareTaskData();
  std::vector<LightTrack> primaries;


  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries,basketSize);
    state.ResumeTiming();
    Td->fPhysicsData->ClearSecondaries();
    for (int i = 0; i < basketSize; ++i) {
      kNish->SampleSecondaries(primaries[i], Td);
    }
  }

  delete kNish;
  CleanTaskData(Td);
}
BENCHMARK(KleinNishinaScalarRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

BENCHMARK_MAIN();

void PreparePrimaries(std::vector<LightTrack> &output,int N)
{
    geant::VecRngWrapper rng;
    output.clear();
    for (int i = 0; i < N; ++i) {
      LightTrack gamma;
      double phi = geant::units::kTwoPi * rng.uniform(); // NOT uniform on shpere
      double th = geant::units::kPi * rng.uniform();
      gamma.SetDirX(sin(phi) * cos(th));
      gamma.SetDirY(cos(phi) * cos(th));
      gamma.SetDirZ(sin(th));
      double eKin = minEn + (maxEn - minEn) * rng.uniform();
      gamma.SetKinE(eKin);
      output.push_back(gamma);
    }
}


void PreparePrimaries(LightTrack_v &output,int N)
{
  geant::VecRngWrapper rng;
  output.SetNtracks(N);
  for (int i = 0; i < N; ++i) {
    double phi = geant::units::kTwoPi * rng.uniform(); // NOT uniform on shpere
    double th  = geant::units::kPi * rng.uniform();
    output.SetDirX(sin(phi) * cos(th),i);
    output.SetDirY(cos(phi) * cos(th),i);
    output.SetDirZ(sin(th),i);
    double eKin = minEn + (maxEn - minEn) * rng.uniform();
    output.SetKinE(eKin,i);
    output.SetTrackIndex(i,i);
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
