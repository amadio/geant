#include "PositronTo2GammaTestCommon.h"
#include "CommonBenchmark.h"

void PreparePosAnihScalarPrims(std::vector<LightTrack> &out, int N)
{
  CreateParticles(kPos2GammaMinEn, kPos2GammaMaxEn, true, TestParticleType::Ep, out, N);
}

void PreparePosAnihVectorPrims(LightTrack_v &out, int N)
{
  CreateParticles(kPos2GammaMinEn, kPos2GammaMaxEn, true, TestParticleType::Ep, out, N);
}

int main(int argc, char **argv)
{
  PrepareWorld();

  auto td = PrepareTaskData();

  std::unique_ptr<EMModel> pos2gRej   = InitEMModel(new PositronTo2GammaModel, kPos2GammaMinEn, kPos2GammaMaxEn, false);
  std::unique_ptr<EMModel> pos2gTable = InitEMModel(new PositronTo2GammaModel, kPos2GammaMinEn, kPos2GammaMaxEn, true);

  benchmark::RegisterBenchmark("PosTo2GammaAliasScal", ScalarModelBenchmark, pos2gTable.get(),
                               PreparePosAnihScalarPrims, td.get(), kBasketSize);
  benchmark::RegisterBenchmark("PosTo2GammaAliasVec", VectorModelBenchmark, pos2gTable.get(), PreparePosAnihVectorPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("PosTo2GammaRejScal", ScalarModelBenchmark, pos2gRej.get(), PreparePosAnihScalarPrims,
                               td.get(), kBasketSize);
  benchmark::RegisterBenchmark("PosTo2GammaRejVec", VectorModelBenchmark, pos2gRej.get(), PreparePosAnihVectorPrims,
                               td.get(), kBasketSize);

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td.get());
}
