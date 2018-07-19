#include "SauterGavrilaTestCommon.h"

#include "Geant/RngWrapper.h"

class SauterGavrilaAliasTester : public SauterGavrilaPhotoElectricModel {
public:
  SauterGavrilaAliasTester() : SauterGavrilaPhotoElectricModel(){};
  using SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirection_Alias;
};

class SauterGavrilaAliasVec : public SauterGavrilaPhotoElectricModel {
public:
  SauterGavrilaAliasVec() : SauterGavrilaPhotoElectricModel(){};
  using SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirectionAliasVec;
};

const int kTestSize = 1024 * 10;
int main()
{
  SauterGavrilaAliasVec *sgv = new SauterGavrilaAliasVec;
  sgv->SetLowEnergyUsageLimit(minEn);
  sgv->SetHighEnergyUsageLimit(maxEn);
  sgv->SetUseSamplingTables(true);
  sgv->Initialize();

  SauterGavrilaAliasTester *slnStd = new SauterGavrilaAliasTester;
  slnStd->SetLowEnergyUsageLimit(minEn);
  slnStd->SetHighEnergyUsageLimit(maxEn);
  slnStd->SetUseSamplingTables(true);
  slnStd->Initialize();

  geant::RngWrapper rng;
  std::vector<double> energy;
  std::vector<double> out1;
  std::vector<double> out2;
  std::vector<double> r1;
  std::vector<double> r2;
  std::vector<double> r3;

  for (int i = 0; i < kTestSize; ++i) {
    energy.push_back(rng.uniform(minEn, maxEn));
    out1.push_back(0.0);
    out2.push_back(0.0);
    r1.push_back(rng.uniform());
    r2.push_back(rng.uniform());
    r3.push_back(rng.uniform());
  }

  for (int i = 0; i < kTestSize; ++i) {
    out1[i] = slnStd->SamplePhotoElectronDirection_Alias(energy[i], r1[i], r2[i], r3[i]);
    // std::cout<<out1[i] <<std::endl;
  }
  for (int i = 0; i < kTestSize; i += geant::kVecLenD) {
    geant::Double_v en, r1v, r2v, r3v;
    vecCore::Load(en, energy.data() + i);
    vecCore::Load(r1v, r1.data() + i);
    vecCore::Load(r2v, r2.data() + i);
    vecCore::Load(r3v, r3.data() + i);
    geant::Double_v eps = sgv->SamplePhotoElectronDirectionAliasVec(en, r1v, r2v, r3v);
    // std::cout<<eps<<std::endl;
    vecCore::Store(eps, out2.data() + i);
  }

  double cumError = 0.0;
  for (int i = 0; i < kTestSize; ++i) {
    cumError += std::abs(out1[i] - out2[i]);
  }
  Printf("TestSize: %d Cumulative error: %f", kTestSize, cumError);

  delete sgv;
  delete slnStd;
}
