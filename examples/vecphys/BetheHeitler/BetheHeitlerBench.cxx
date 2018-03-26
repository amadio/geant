#include <benchmark/benchmark.h>
#include "BetheHeitlerTestCommon.h"

#include "Geant/VecRngWrapper.h"

BetheHeitlerPairModel *bHModelRej;
VecBetheHeitlerPairModel *bHVecModelRej;
BetheHeitlerPairModel *bHModelAlias;
VecBetheHeitlerPairModel *bHVecModelAlias;
TaskData *td;

static void BetheHeitlerScalarAlias(benchmark::State &state)
{
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    td->fPhysicsData->ClearSecondaries();

    state.ResumeTiming();

    for (int i = 0; i < basketSize; ++i) {
      bHModelAlias->SampleSecondaries(primaries[i], td);
    }
  }
}
BENCHMARK(BetheHeitlerScalarAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void BetheHeitlerVectorAlias(benchmark::State &state)
{
  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    td->fPhysicsData->GetSecondarySOA().ClearTracks();

    state.ResumeTiming();

    bHVecModelAlias->SampleSecondariesVector(primaries, td);
  }
}
BENCHMARK(BetheHeitlerVectorAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

// static void KleinNishinaScalarRejBaseline(benchmark::State &state)
//{
//  KleinNishinaComptonModel *kNish = PrepareKnishinaModel(false);
//
//  auto Td = PrepareTaskData();
//  std::vector<LightTrack> primaries;
//
//  int basketSize = state.range(0);
//
//  for (auto _ : state) {
//    PreparePrimaries(primaries, basketSize);
//    Td->fPhysicsData->ClearSecondaries();
//  }
//
//  benchmark::DoNotOptimize(primaries.data());
//
//  delete kNish;
//  CleanTaskData(Td);
//}
// BENCHMARK(KleinNishinaScalarRejBaseline)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);
//
// static void KleinNishinaScalarRej(benchmark::State &state)
//{
//  KleinNishinaComptonModel *kNish = PrepareKnishinaModel(false);
//
//  auto Td = PrepareTaskData();
//  std::vector<LightTrack> primaries;
//
//  int basketSize = state.range(0);
//
//  for (auto _ : state) {
//    PreparePrimaries(primaries, basketSize);
//    Td->fPhysicsData->ClearSecondaries();
//    for (int i = 0; i < basketSize; ++i) {
//      kNish->SampleSecondaries(primaries[i], Td);
//    }
//  }
//
//  delete kNish;
//  CleanTaskData(Td);
//}
// BENCHMARK(KleinNishinaScalarRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);
//
// static void KleinNishinaVectorRejBaseline(benchmark::State &state)
//{
//  VecKleinNishinaComptonModel *kNish = PrepareVecKnishinaModel(false);
//
//  auto Td = PrepareTaskData();
//
//  LightTrack_v primaries;
//
//  int basketSize = state.range(0);
//
//  for (auto _ : state) {
//    PreparePrimaries(primaries, basketSize);
//    primaries.SetNtracks(basketSize);
//
//    Td->fPhysicsData->GetSecondarySOA().ClearTracks();
//  }
//  benchmark::DoNotOptimize(&primaries);
//
//  delete kNish;
//  CleanTaskData(Td);
//}
// BENCHMARK(KleinNishinaVectorRejBaseline)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);
//
// static void KleinNishinaVectorRej(benchmark::State &state)
//{
//  VecKleinNishinaComptonModel *kNish = PrepareVecKnishinaModel(false);
//
//  auto Td = PrepareTaskData();
//
//  LightTrack_v primaries;
//
//  int basketSize = state.range(0);
//
//  for (auto _ : state) {
//    PreparePrimaries(primaries, basketSize);
//    primaries.SetNtracks(basketSize);
//
//    Td->fPhysicsData->GetSecondarySOA().ClearTracks();
//    kNish->SampleSecondariesVector(primaries, Td);
//  }
//
//  delete kNish;
//  CleanTaskData(Td);
//}
// BENCHMARK(KleinNishinaVectorRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

int main(int argc, char **argv)
{
  PrepareWorld();
  bHModelAlias    = PrepareBHModel(true);
  bHModelRej      = PrepareBHModel(false);
  bHVecModelAlias = PrepareVecBHModel(true);
  bHVecModelRej   = PrepareVecBHModel(false);
  td              = PrepareTaskData();

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td);
}
