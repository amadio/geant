#include "CommonBenchmark.h"
#include "SGTestCommon.h"

void PrepareMBScalarPrims(std::vector<LightTrack> &out, int N)
{
  CreateParticles(kSauGavminEn, kSauGavmaxEn, true, TestParticleType::Gamma, out, N);
}

void PrepareMBVectorPrims(LightTrack_v &out, int N)
{
  CreateParticles(kSauGavminEn, kSauGavmaxEn, true, TestParticleType::Gamma, out, N);
}

int main(int argc, char **argv)
{
  PrepareWorld();

  auto td = PrepareTaskData();

  std::unique_ptr<EMModel> mbRej =
      InitEMModel(new SauterGavrilaPhotoElectricModel(), kSauGavminEn, kSauGavmaxEn, false);
  std::unique_ptr<EMModel> mbTable =
      InitEMModel(new SauterGavrilaPhotoElectricModel(), kSauGavminEn, kSauGavmaxEn, true);

  EMMODEL_REGISTER_SCALAR_BENCHMARK("SauterGavrilaAliasScal", mbTable.get(), PrepareMBScalarPrims, td.get(),
                                    kBasketSize);
  EMMODEL_REGISTER_VECTOR_BENCHMARK("SauterGavrilaAliasVec", mbTable.get(), PrepareMBVectorPrims, td.get(),
                                    kBasketSize);
  EMMODEL_REGISTER_SCALAR_BENCHMARK("SauterGavrilaRejScal", mbRej.get(), PrepareMBScalarPrims, td.get(), kBasketSize);
  EMMODEL_REGISTER_VECTOR_BENCHMARK("SauterGavrilaRejVec", mbRej.get(), PrepareMBVectorPrims, td.get(), kBasketSize);

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td.get());
}
