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

  std::unique_ptr<EMModel> bhScalarRej   = InitEMModel(new BetheHeitlerPairModel, kBHminEn, kBHmaxEn, false);
  std::unique_ptr<EMModel> bhVectorRej   = InitEMModel(new VecBetheHeitlerPairModel, kBHminEn, kBHmaxEn, false);
  std::unique_ptr<EMModel> bhScalarTable = InitEMModel(new BetheHeitlerPairModel, kBHminEn, kBHmaxEn, true);
  std::unique_ptr<EMModel> bhVectorTable = InitEMModel(new VecBetheHeitlerPairModel, kBHminEn, kBHmaxEn, true);

  std::unique_ptr<EMModel> relPairScalarRej =
      InitEMModel(new RelativisticPairModel, kRelPairMinEn, kRelPairMaxEn, false);
  std::unique_ptr<EMModel> relPairVectorRej =
      InitEMModel(new VecRelativisticPairModel, kRelPairMinEn, kRelPairMaxEn, false);
  std::unique_ptr<EMModel> relPairScalarTable =
      InitEMModel(new RelativisticPairModel, kRelPairMinEn, kRelPairMaxEn, true);
  std::unique_ptr<EMModel> relPairVectorTable =
      InitEMModel(new VecRelativisticPairModel, kRelPairMinEn, kRelPairMaxEn, true);

  benchmark::RegisterBenchmark("BetheHeitlerAliasScal", ScalarModelBenchmark, bhScalarTable.get(), PrepareBHScalarPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("BetheHeitlerAliasVec", VectorModelBenchmark, bhVectorTable.get(), PrepareBHVectorPrims,
                               td.get(), kBasketSize);

  benchmark::RegisterBenchmark("BetheHeitlerRejScal", ScalarModelBenchmark, bhScalarRej.get(), PrepareBHScalarPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("BetheHeitlerRejVec", VectorModelBenchmark, bhVectorRej.get(), PrepareBHVectorPrims,
                               td.get(), kBasketSize);

  benchmark::RegisterBenchmark("RelativisticPairRejScal", ScalarModelBenchmark, relPairScalarRej.get(),
                               PrepareRPScalarPrims, td.get(), kBasketSize);
  benchmark::RegisterBenchmark("RelativisticPairRejVec", VectorModelBenchmark, relPairVectorRej.get(),
                               PrepareRPVectorPrims, td.get(), kBasketSize);

  benchmark::RegisterBenchmark("RelativisticPairAliasScal", ScalarModelBenchmark, relPairScalarTable.get(),
                               PrepareRPScalarPrims, td.get(), kBasketSize);
  benchmark::RegisterBenchmark("RelativisticPairAliasVec", VectorModelBenchmark, relPairVectorTable.get(),
                               PrepareRPVectorPrims, td.get(), kBasketSize);

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();
}
