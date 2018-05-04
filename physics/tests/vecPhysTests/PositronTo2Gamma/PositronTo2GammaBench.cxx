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

  EMMODEL_REGISTER_SCALAR_BENCHMARK("PosTo2GammaAliasScal", pos2gTable.get(), PreparePosAnihScalarPrims, td.get(),
                                    kBasketSize);
  EMMODEL_REGISTER_VECTOR_BENCHMARK("PosTo2GammaAliasVec", pos2gTable.get(), PreparePosAnihVectorPrims, td.get(),
                                    kBasketSize);
  EMMODEL_REGISTER_SCALAR_BENCHMARK("PosTo2GammaRejScal", pos2gRej.get(), PreparePosAnihScalarPrims, td.get(),
                                    kBasketSize);
  EMMODEL_REGISTER_VECTOR_BENCHMARK("PosTo2GammaRejVec", pos2gRej.get(), PreparePosAnihVectorPrims, td.get(),
                                    kBasketSize);

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td.get());
}
