#include "CommonBenchmark.h"
#include "MollerBhabhaTestCommon.h"

void PrepareMBScalarPrims(std::vector<LightTrack> &out, int N)
{
  CreateParticles(kMollBhminEn, kMollBhmaxEn, true, TestParticleType::Em, out, N);
}

void PrepareMBVectorPrims(LightTrack_v &out, int N)
{
  CreateParticles(kMollBhminEn, kMollBhmaxEn, true, TestParticleType::Em, out, N);
}

int main(int argc, char **argv)
{
  PrepareWorld();

  auto td = PrepareTaskData();

  std::unique_ptr<EMModel> mbRej =
      InitEMModel(new MollerBhabhaIonizationModel(true), kMollBhminEn, kMollBhmaxEn, false);
  std::unique_ptr<EMModel> mbTable =
      InitEMModel(new MollerBhabhaIonizationModel(true), kMollBhminEn, kMollBhmaxEn, true);

  EMMODEL_REGISTER_SCALAR_BENCHMARK("MollerBhabhaAliasScal", mbTable.get(), PrepareMBScalarPrims, td.get(),
                                    kBasketSize);
  EMMODEL_REGISTER_VECTOR_BENCHMARK("MollerBhabhaAliasVec", mbTable.get(), PrepareMBVectorPrims, td.get(), kBasketSize);
  EMMODEL_REGISTER_SCALAR_BENCHMARK("MollerBhabhaRejScal", mbRej.get(), PrepareMBScalarPrims, td.get(), kBasketSize);
  EMMODEL_REGISTER_VECTOR_BENCHMARK("MollerBhabhaRejVec", mbRej.get(), PrepareMBVectorPrims, td.get(), kBasketSize);

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td.get());
}
