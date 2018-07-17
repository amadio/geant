#include <benchmark/benchmark.h>
#include "SauterGavrilaTestCommon.h"

#include "Geant/RngWrapper.h"

SauterGavrilaPhotoElectricModel *sg;
SauterGavrilaPhotoElectricModel *vsg;
SauterGavrilaPhotoElectricModel *sgrej;
SauterGavrilaPhotoElectricModel *vsgrej;
auto td  = PrepareTaskData();

static void SauterGavrilaSampleSecondariesAliasScalar(benchmark::State &state)
{
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries, basketSize);
    td->fPhysicsData->ClearSecondaries();
    state.ResumeTiming();
    for (int i = 0; i < basketSize; ++i) {
      sg->SampleSecondaries(primaries[i], td.get());
    }
  }
}
BENCHMARK(SauterGavrilaSampleSecondariesAliasScalar)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void SauterGavrilaSampleSecondariesAliasVector(benchmark::State &state)
{
  LightTrack_v primaries;
  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    td->fPhysicsData->GetSecondarySOA().ClearTracks();
    state.ResumeTiming();
    vsg->SampleSecondaries(primaries, td.get());
  }
}
BENCHMARK(SauterGavrilaSampleSecondariesAliasVector)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void SauterGavrilaSampleSecondariesRejScalar(benchmark::State &state)
{
  std::vector<LightTrack> primaries;
  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries, basketSize);
    td->fPhysicsData->ClearSecondaries();
    state.ResumeTiming();
    for (int i = 0; i < basketSize; ++i) {
      sgrej->SampleSecondaries(primaries[i], td.get());
    }
  }

}
BENCHMARK(SauterGavrilaSampleSecondariesRejScalar)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void SauterGavrilaSampleSecondariesRejVector(benchmark::State &state)
{
  LightTrack_v primaries;
  int basketSize = state.range(0);
  
    for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    td->fPhysicsData->GetSecondarySOA().ClearTracks();
    state.ResumeTiming();
    vsgrej->SampleSecondaries(primaries, td.get());
  }

}
BENCHMARK(SauterGavrilaSampleSecondariesRejVector)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

int main(int argc, char **argv)
{
    SetUpSimulation();
    sg  = PrepareSauterGavrilaModel(true);
    vsg = PrepareVecSauterGavrilaModel(true);
    sgrej = PrepareSauterGavrilaModel(false);
    vsgrej= PrepareVecSauterGavrilaModel(false);
    
    
    ::benchmark::Initialize(&argc, argv);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
    ::benchmark::RunSpecifiedBenchmarks();
    
    CleanTaskData(td.get());
}
