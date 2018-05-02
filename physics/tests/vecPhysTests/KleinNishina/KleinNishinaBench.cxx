#include "CommonBenchmark.h"
#include "KleinNishinaTestCommon.h"

void PrepareKNScalarPrims(std::vector<LightTrack> &out, int N)
{
  CreateParticles(kKNminEn, kKNmaxEn, true, TestParticleType::Em, out, N);
}

void PrepareKNVectorPrims(LightTrack_v &out, int N)
{
  CreateParticles(kKNminEn, kKNmaxEn, true, TestParticleType::Em, out, N);
}

int main(int argc, char **argv)
{
  PrepareWorld();

  auto td = PrepareTaskData();

  std::unique_ptr<EMModel> knScalarRej   = InitEMModel(new KleinNishinaComptonModel, kKNminEn, kKNmaxEn, false);
  std::unique_ptr<EMModel> knVectorRej   = InitEMModel(new KleinNishinaComptonModel, kKNminEn, kKNmaxEn, false);
  std::unique_ptr<EMModel> knScalarTable = InitEMModel(new KleinNishinaComptonModel, kKNminEn, kKNmaxEn, true);
  std::unique_ptr<EMModel> knVectorTable = InitEMModel(new KleinNishinaComptonModel, kKNminEn, kKNmaxEn, true);

  benchmark::RegisterBenchmark("KleinNishinaAliasScal", ScalarModelBenchmark, knScalarTable.get(), PrepareKNScalarPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("KleinNishinaAliasVec", VectorModelBenchmark, knVectorTable.get(), PrepareKNVectorPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("KleinNishinaRejScal", ScalarModelBenchmark, knScalarRej.get(), PrepareKNScalarPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("KleinNishinaRejVec", VectorModelBenchmark, knVectorRej.get(), PrepareKNVectorPrims,
                               td.get(), kBasketSize);

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td.get());
}
