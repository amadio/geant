#include <Geant/VecSauterGavrilaPhotoElectricModel.h>
#include "SauterGavrilaTestCommon.h"

#include "Geant/RngWrapper.h"
#include <random>
#include <cmath>

class SauterGavrilaAliasTester : public SauterGavrilaPhotoElectricModel {
public:
  SauterGavrilaAliasTester() : SauterGavrilaPhotoElectricModel(){};
  using SauterGavrilaPhotoElectricModel::SampleShell;
};

using geantphysics::SauterGavrilaPhotoElectricModel;
class SauterGavrilaAliasVec : public SauterGavrilaPhotoElectricModel {
public:
  SauterGavrilaAliasVec() : SauterGavrilaPhotoElectricModel(){};
  using SauterGavrilaPhotoElectricModel::SampleShellAliasVec;
};

const int kTestSize = 1024 * 10;
int main()
{
  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<int> uni(3,98); // guaranteed unbiased
  
  SauterGavrilaAliasTester *slnStd = new SauterGavrilaAliasTester;
  slnStd->SetLowEnergyUsageLimit(minEn);
  slnStd->SetHighEnergyUsageLimit(maxEn);
  slnStd->SetUseSamplingTables(true);
  slnStd->Initialize();
    
  SauterGavrilaAliasVec *sgv = new SauterGavrilaAliasVec;
  sgv->SetLowEnergyUsageLimit(minEn);
  sgv->SetHighEnergyUsageLimit(maxEn);
  sgv->SetUseSamplingTables(true);
  sgv->Initialize();

  geant::RngWrapper rngGV;
  std::vector<double> energy;
  std::vector<int> zed;
  std::vector<size_t> out1;
  std::vector<int> out2;
  std::vector<double> r1;
  std::vector<double> r2;

  for (int i = 0; i < kTestSize; ++i) {
    energy.push_back(rngGV.uniform(minEn, maxEn));
    auto random_integer = uni(rng);
    zed.push_back(random_integer);
    out1.push_back(0);
    out2.push_back(0);
    r1.push_back(rngGV.uniform());
    r2.push_back(rngGV.uniform());
  }

  for (int i = 0; i < kTestSize; ++i) {
      size_t tmp=(size_t) zed[i];
      slnStd->SampleShellAlias(energy[i], tmp, r1[i], r2[i],out1[i]);
      //std::cout<<out1[i] <<std::endl;
  }
    for (int i = 0; i < kTestSize; i += geant::kVecLenD) {
    geant::Double_v en, r1v, r2v;
    geant::IndexD_v  z;
    vecCore::Load(en, energy.data() + i);
    vecCore::Load(z,  zed.data() + i);
    vecCore::Load(r1v, r1.data() + i);
    vecCore::Load(r2v, r2.data() + i);
    geant::IndexD_v eps = sgv->SampleShellAliasVec(en, z, r1v, r2v);
    //std::cout<<eps<<std::endl;
    vecCore::Store(eps, out2.data() + i);
  }

  double cumError = 0.0;
  for (int i = 0; i < kTestSize; ++i) {
    cumError += std::abs((int)(out1[i] - out2[i]));
  }
  Printf("TestSize: %d Cumulative error: %f", kTestSize, cumError);

  delete sgv;
  delete slnStd;
}
