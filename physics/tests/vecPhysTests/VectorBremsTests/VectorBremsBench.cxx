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

  EMMODEL_REGISTER_SCALAR_BENCHMARK("SeltzerBergerAliasScal", sbTable.get(), PrepareSBScalarPrims, td.get(),
                                    kBasketSize);
  EMMODEL_REGISTER_VECTOR_BENCHMARK("SeltzerBergerAliasVec", sbTable.get(), PrepareSBVectorPrims, td.get(),
                                    kBasketSize);

  EMMODEL_REGISTER_SCALAR_BENCHMARK("SeltzerBergerRejScal", sbRej.get(), PrepareSBScalarPrims, td.get(), kBasketSize);
  EMMODEL_REGISTER_VECTOR_BENCHMARK("SeltzerBergerRejVec", sbRej.get(), PrepareSBVectorPrims, td.get(), kBasketSize);

  EMMODEL_REGISTER_SCALAR_BENCHMARK("RelBremsAliasScal", rbTable.get(), PrepareRBScalarPrims, td.get(), kBasketSize);
  EMMODEL_REGISTER_VECTOR_BENCHMARK("RelBremsAliasVec", rbTable.get(), PrepareRBVectorPrims, td.get(), kBasketSize);

  EMMODEL_REGISTER_SCALAR_BENCHMARK("RelBremsRejScal", rbRej.get(), PrepareRBScalarPrims, td.get(), kBasketSize);
  EMMODEL_REGISTER_VECTOR_BENCHMARK("RelBremsRejVec", rbRej.get(), PrepareRBVectorPrims, td.get(), kBasketSize);

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td.get());
  return 0;
}
