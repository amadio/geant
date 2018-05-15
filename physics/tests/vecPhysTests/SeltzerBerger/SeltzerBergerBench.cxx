#include <benchmark/benchmark.h>
#include "SeltzerBergerTestCommon.h"

#include "Geant/VecRngWrapper.h"

SeltzerBergerBremsModel *aliasEl;
SeltzerBergerBremsModel *rejEl;
VecSeltzerBergerBremsModel *vecAliasEl;
VecSeltzerBergerBremsModel *vecRejEl;
TaskData *td;

static void SeltzerBergerScalarAlias(benchmark::State &state)
{
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    td->fPhysicsData->ClearSecondaries();

    state.ResumeTiming();

    for (int i = 0; i < basketSize; ++i) {
      aliasEl->SampleSecondaries(primaries[i], td);
    }
  }
}
BENCHMARK(SeltzerBergerScalarAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void SeltzerBergerVectorAlias(benchmark::State &state)
{
  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    td->fPhysicsData->GetSecondarySOA().ClearTracks();

    state.ResumeTiming();

    vecAliasEl->SampleSecondariesVector(primaries, td);
  }
}
BENCHMARK(SeltzerBergerVectorAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void SeltzerBergerScalarRej(benchmark::State &state)
{
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    td->fPhysicsData->ClearSecondaries();

    state.ResumeTiming();

    for (int i = 0; i < basketSize; ++i) {
      rejEl->SampleSecondaries(primaries[i], td);
    }
  }
}
BENCHMARK(SeltzerBergerScalarRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void SeltzerBergerVectorRej(benchmark::State &state)
{
  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    td->fPhysicsData->GetSecondarySOA().ClearTracks();

    state.ResumeTiming();

    vecRejEl->SampleSecondariesVector(primaries, td);
  }
}
BENCHMARK(SeltzerBergerVectorRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

int main(int argc, char **argv)
{
  PrepareWorld();
  rejEl      = PrepareSBModel(false, true);
  aliasEl    = PrepareSBModel(true, true);
  vecRejEl   = PrepareVecSBModel(false, true);
  vecAliasEl = PrepareVecSBModel(true, true);
  td         = PrepareTaskData();

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td);
}
