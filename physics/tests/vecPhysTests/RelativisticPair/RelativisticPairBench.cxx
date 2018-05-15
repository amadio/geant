#include <benchmark/benchmark.h>
#include "RelativisticPairTestCommon.h"

#include "Geant/VecRngWrapper.h"

RelativisticPairModel *rpModelRej;
VecRelativisticPairModel *rpModelRejVec;
RelativisticPairModel *rpModelAlias;
VecRelativisticPairModel *rpModelAliasVec;
TaskData *td;

static void RelativisticPairScalarAlias(benchmark::State &state)
{
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    td->fPhysicsData->ClearSecondaries();

    state.ResumeTiming();

    for (int i = 0; i < basketSize; ++i) {
      rpModelAlias->SampleSecondaries(primaries[i], td);
    }
  }
}
BENCHMARK(RelativisticPairScalarAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void RelativisticPairVectorAlias(benchmark::State &state)
{
  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    td->fPhysicsData->GetSecondarySOA().ClearTracks();

    state.ResumeTiming();

    rpModelAliasVec->SampleSecondariesVector(primaries, td);
  }
}
BENCHMARK(RelativisticPairVectorAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void RelativisticPairScalarRej(benchmark::State &state)
{
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    td->fPhysicsData->ClearSecondaries();

    state.ResumeTiming();

    for (int i = 0; i < basketSize; ++i) {
      rpModelRej->SampleSecondaries(primaries[i], td);
    }
  }
}
BENCHMARK(RelativisticPairScalarRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void RelativisticPairVectorRej(benchmark::State &state)
{
  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    td->fPhysicsData->GetSecondarySOA().ClearTracks();

    state.ResumeTiming();

    rpModelRejVec->SampleSecondariesVector(primaries, td);
  }
}
BENCHMARK(RelativisticPairVectorRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

int main(int argc, char **argv)
{
  PrepareWorld();
  rpModelAlias    = PrepareRPModel(true);
  rpModelRej      = PrepareRPModel(false);
  rpModelAliasVec = PrepareVecRPModel(true);
  rpModelRejVec   = PrepareVecRPModel(false);
  td              = PrepareTaskData();

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td);
}
