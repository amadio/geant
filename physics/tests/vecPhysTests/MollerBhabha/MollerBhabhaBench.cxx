#include <benchmark/benchmark.h>
#include "MollerBhabhaTestCommon.h"

#include "Geant/VecRngWrapper.h"

MollerBhabhaIonizationModel *aliasEl;
MollerBhabhaIonizationModel *rejEl;
VecMollerBhabhaIonizationModel *vecAliasEl;
VecMollerBhabhaIonizationModel *vecRejEl;
TaskData *td;

static void MollerBhabhaScalarAlias(benchmark::State &state)
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
BENCHMARK(MollerBhabhaScalarAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void MollerBhabhaVectorAlias(benchmark::State &state)
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
BENCHMARK(MollerBhabhaVectorAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void MollerBhabhaScalarRej(benchmark::State &state)
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
BENCHMARK(MollerBhabhaScalarRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void MollerBhabhaVectorRej(benchmark::State &state)
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
BENCHMARK(MollerBhabhaVectorRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

int main(int argc, char **argv)
{
  PrepareWorld();
  rejEl      = PrepareMBModel(false, true);
  aliasEl    = PrepareMBModel(true, true);
  vecRejEl   = PrepareVecMBModel(false, true);
  vecAliasEl = PrepareVecMBModel(true, true);
  td         = PrepareTaskData();

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td);
}
