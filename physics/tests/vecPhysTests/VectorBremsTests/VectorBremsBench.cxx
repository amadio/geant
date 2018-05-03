#include "VectorBremsTestCommon.h"
#include "CommonBenchmark.h"

void PrepareSBScalarPrims(std::vector<LightTrack> &out, int N)
{
  CreateParticles(kSBminEn, kSBmaxEn, true, TestParticleType::Em, out, N);
}

void PrepareSBVectorPrims(LightTrack_v &out, int N)
{
  CreateParticles(kSBminEn, kSBmaxEn, true, TestParticleType::Em, out, N);
}

void PrepareRBScalarPrims(std::vector<LightTrack> &out, int N)
{
  CreateParticles(kRelBremsMinEn, kRelBremsMaxEn, true, TestParticleType::Em, out, N);
}

void PrepareRBVectorPrims(LightTrack_v &out, int N)
{
  CreateParticles(kRelBremsMinEn, kRelBremsMaxEn, true, TestParticleType::Em, out, N);
}

int main(int argc, char **argv)
{
  PrepareWorld();

  auto td = PrepareTaskData();

  std::unique_ptr<EMModel> sbRej   = InitEMModel(new SeltzerBergerBremsModel(true), kSBminEn, kSBmaxEn, false);
  std::unique_ptr<EMModel> sbTable = InitEMModel(new SeltzerBergerBremsModel(true), kSBminEn, kSBmaxEn, true);

  std::unique_ptr<EMModel> rbRej   = InitEMModel(new RelativisticBremsModel(), kRelBremsMinEn, kRelBremsMaxEn, false);
  std::unique_ptr<EMModel> rbTable = InitEMModel(new RelativisticBremsModel(), kRelBremsMinEn, kRelBremsMaxEn, true);

  benchmark::RegisterBenchmark("SeltzerBergerAliasScal", ScalarModelBenchmark, sbTable.get(), PrepareSBScalarPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("SeltzerBergerAliasVec", VectorModelBenchmark, sbTable.get(), PrepareSBVectorPrims,
                               td.get(), kBasketSize);

  benchmark::RegisterBenchmark("SeltzerBergerRejScal", ScalarModelBenchmark, sbRej.get(), PrepareSBScalarPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("SeltzerBergerRejVec", VectorModelBenchmark, sbRej.get(), PrepareSBVectorPrims, td.get(),
                               kBasketSize);

  benchmark::RegisterBenchmark("RelBremsAliasScal", ScalarModelBenchmark, rbTable.get(), PrepareRBScalarPrims, td.get(),
                               kBasketSize);
  benchmark::RegisterBenchmark("RelBremsAliasVec", VectorModelBenchmark, rbTable.get(), PrepareRBVectorPrims, td.get(),
                               kBasketSize);

  benchmark::RegisterBenchmark("RelBremsRejScal", ScalarModelBenchmark, rbRej.get(), PrepareRBScalarPrims, td.get(),
                               kBasketSize);
  benchmark::RegisterBenchmark("RelBremsRejVec", VectorModelBenchmark, rbRej.get(), PrepareRBVectorPrims, td.get(),
                               kBasketSize);

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td.get());
  return 0;
}
