#include "VectorPairTestCommon.h"
#include "CommonBenchmark.h"

void PrepareBHScalarPrims(std::vector<LightTrack> &out, int N)
{
  CreateParticles(kBHminEn, kBHmaxEn, true, TestParticleType::Gamma, out, N);
}

void PrepareBHVectorPrims(LightTrack_v &out, int N)
{
  CreateParticles(kBHminEn, kBHmaxEn, true, TestParticleType::Gamma, out, N);
}

void PrepareRPScalarPrims(std::vector<LightTrack> &out, int N)
{
  CreateParticles(kRelPairMinEn, kRelPairMaxEn, true, TestParticleType::Gamma, out, N);
}

void PrepareRPVectorPrims(LightTrack_v &out, int N)
{
  CreateParticles(kRelPairMinEn, kRelPairMaxEn, true, TestParticleType::Gamma, out, N);
}

int main(int argc, char **argv)
{
  PrepareWorld();

  auto td = PrepareTaskData();

  std::unique_ptr<EMModel> bhRej   = InitEMModel(new BetheHeitlerPairModel, kBHminEn, kBHmaxEn, false);
  std::unique_ptr<EMModel> bhTable = InitEMModel(new BetheHeitlerPairModel, kBHminEn, kBHmaxEn, true);

  std::unique_ptr<EMModel> relPairRej   = InitEMModel(new RelativisticPairModel, kRelPairMinEn, kRelPairMaxEn, false);
  std::unique_ptr<EMModel> relPairTable = InitEMModel(new RelativisticPairModel, kRelPairMinEn, kRelPairMaxEn, true);

  benchmark::RegisterBenchmark("BetheHeitlerAliasScal", ScalarModelBenchmark, bhTable.get(), PrepareBHScalarPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("BetheHeitlerAliasVec", VectorModelBenchmark, bhTable.get(), PrepareBHVectorPrims,
                               td.get(), kBasketSize);

  benchmark::RegisterBenchmark("BetheHeitlerRejScal", ScalarModelBenchmark, bhRej.get(), PrepareBHScalarPrims, td.get(),
                               kBasketSize);
  benchmark::RegisterBenchmark("BetheHeitlerRejVec", VectorModelBenchmark, bhRej.get(), PrepareBHVectorPrims, td.get(),
                               kBasketSize);

  benchmark::RegisterBenchmark("RelativisticPairRejScal", ScalarModelBenchmark, relPairRej.get(), PrepareRPScalarPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("RelativisticPairRejVec", VectorModelBenchmark, relPairRej.get(), PrepareRPVectorPrims,
                               td.get(), kBasketSize);

  benchmark::RegisterBenchmark("RelativisticPairAliasScal", ScalarModelBenchmark, relPairTable.get(),
                               PrepareRPScalarPrims, td.get(), kBasketSize);
  benchmark::RegisterBenchmark("RelativisticPairAliasVec", VectorModelBenchmark, relPairTable.get(),
                               PrepareRPVectorPrims, td.get(), kBasketSize);

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td.get());
  return 0;
}
