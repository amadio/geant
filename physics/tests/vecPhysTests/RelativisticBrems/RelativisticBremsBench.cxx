#include <benchmark/benchmark.h>
#include "RelativisticBremsTestCommon.h"

#include "Geant/VecRngWrapper.h"

RelativisticBremsModel *alias;
RelativisticBremsModel *rej;
VecRelativisticBremsModel *vecAlias;
VecRelativisticBremsModel *vecRej;
TaskData *td;

static void RelativisticScalarAlias(benchmark::State &state)
{
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    td->fPhysicsData->ClearSecondaries();

    state.ResumeTiming();

    for (int i = 0; i < basketSize; ++i) {
      alias->SampleSecondaries(primaries[i], td);
    }
  }
}
BENCHMARK(RelativisticScalarAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void RelativisticVectorAlias(benchmark::State &state)
{
  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    td->fPhysicsData->GetSecondarySOA().ClearTracks();

    state.ResumeTiming();

    vecAlias->SampleSecondariesVector(primaries, td);
  }
}
BENCHMARK(RelativisticVectorAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void RelativisticScalarRej(benchmark::State &state)
{
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    td->fPhysicsData->ClearSecondaries();

    state.ResumeTiming();

    for (int i = 0; i < basketSize; ++i) {
      rej->SampleSecondaries(primaries[i], td);
    }
  }
}
BENCHMARK(RelativisticScalarRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void RelativisticVectorRej(benchmark::State &state)
{
  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    td->fPhysicsData->GetSecondarySOA().ClearTracks();

    state.ResumeTiming();

    vecRej->SampleSecondariesVector(primaries, td);
  }
}
BENCHMARK(RelativisticVectorRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

int main(int argc, char **argv)
{
  PrepareWorld();
  rej      = PrepareRelBremsModel(false);
  alias    = PrepareRelBremsModel(true);
  vecRej   = PrepareVecRelBremsModel(false);
  vecAlias = PrepareVecRelBremsModel(true);
  td       = PrepareTaskData();

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td);
}
