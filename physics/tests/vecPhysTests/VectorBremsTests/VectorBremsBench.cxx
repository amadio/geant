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

  std::unique_ptr<EMModel> sbScalarRej   = InitEMModel(new SeltzerBergerBremsModel(true), kSBminEn, kSBmaxEn, false);
  std::unique_ptr<EMModel> sbVectorRej   = InitEMModel(new SeltzerBergerBremsModel(true), kSBminEn, kSBmaxEn, false);
  std::unique_ptr<EMModel> sbScalarTable = InitEMModel(new SeltzerBergerBremsModel(true), kSBminEn, kSBmaxEn, true);
  std::unique_ptr<EMModel> sbVectorTable = InitEMModel(new SeltzerBergerBremsModel(true), kSBminEn, kSBmaxEn, true);

  std::unique_ptr<EMModel> rbScalarRej =
      InitEMModel(new RelativisticBremsModel(), kRelBremsMinEn, kRelBremsMaxEn, false);
  std::unique_ptr<EMModel> rbVectorRej =
      InitEMModel(new RelativisticBremsModel(), kRelBremsMinEn, kRelBremsMaxEn, false);
  std::unique_ptr<EMModel> rbScalarTable =
      InitEMModel(new RelativisticBremsModel(), kRelBremsMinEn, kRelBremsMaxEn, true);
  std::unique_ptr<EMModel> rbVectorTable =
      InitEMModel(new RelativisticBremsModel(), kRelBremsMinEn, kRelBremsMaxEn, true);

  benchmark::RegisterBenchmark("SeltzerBergerAliasScal", ScalarModelBenchmark, sbScalarTable.get(),
                               PrepareSBScalarPrims, td.get(), kBasketSize);
  benchmark::RegisterBenchmark("SeltzerBergerAliasVec", VectorModelBenchmark, sbVectorTable.get(), PrepareSBVectorPrims,
                               td.get(), kBasketSize);

  benchmark::RegisterBenchmark("SeltzerBergerRejScal", ScalarModelBenchmark, sbScalarRej.get(), PrepareSBScalarPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("SeltzerBergerRejVec", VectorModelBenchmark, sbVectorRej.get(), PrepareSBVectorPrims,
                               td.get(), kBasketSize);

  benchmark::RegisterBenchmark("RelBremsAliasScal", ScalarModelBenchmark, rbScalarTable.get(), PrepareRBScalarPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("RelBremsAliasVec", VectorModelBenchmark, rbVectorTable.get(), PrepareRBVectorPrims,
                               td.get(), kBasketSize);

  benchmark::RegisterBenchmark("RelBremsRejScal", ScalarModelBenchmark, rbScalarRej.get(), PrepareRBScalarPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("RelBremsRejVec", VectorModelBenchmark, rbVectorRej.get(), PrepareRBVectorPrims,
                               td.get(), kBasketSize);

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td.get());
  return 0;
}
