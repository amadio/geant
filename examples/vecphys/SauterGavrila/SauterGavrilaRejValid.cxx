/*
 This sanity check doesn't make real sense until we will guarantee the reproducibility for random numbers
 (the ones 'wasted' in the vectorized rejection have to be put back and used again)
 */
#include <Geant/VecSauterGavrilaPhotoElectricModel.h>
#include "SauterGavrilaTestCommon.h"

#include "Geant/VecRngWrapper.h"

class SauterGavrilaRejTesterScalar : public SauterGavrilaPhotoElectricModel {
public:
  SauterGavrilaRejTesterScalar() : SauterGavrilaPhotoElectricModel(){};
  using SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirection_Rejection;
};

using geantphysics::VecSauterGavrilaPhotoElectricModel;
class SauterGavrilaRejTesterVec : public VecSauterGavrilaPhotoElectricModel {
public:
  SauterGavrilaRejTesterVec() : VecSauterGavrilaPhotoElectricModel(){};
  using VecSauterGavrilaPhotoElectricModel::SamplePhotoElectronDirectionRejVec;
};

const int kTestSize = 1024 * 10;
int main()
{
  SauterGavrilaRejTesterVec *sgv = new SauterGavrilaRejTesterVec;
  sgv->SetLowEnergyUsageLimit(minEn);
  sgv->SetHighEnergyUsageLimit(maxEn);
  sgv->SetUseSamplingTables(true);
  sgv->Initialize();

  SauterGavrilaRejTesterScalar *slnStd = new SauterGavrilaRejTesterScalar;
  slnStd->SetLowEnergyUsageLimit(minEn);
  slnStd->SetHighEnergyUsageLimit(maxEn);
  slnStd->SetUseSamplingTables(true);
  slnStd->Initialize();
  
  geant::VecRngWrapper rng;
  std::vector<double> energy;
  std::vector<double> cosTheta1;
  std::vector<double> cosTheta2;

  for (int i = 0; i < kTestSize; ++i) {
    energy.push_back(rng.uniform(minEn, maxEn));
    cosTheta1.push_back(0.0);
    cosTheta2.push_back(0.0);
  }

  TaskData *td = PrepareTaskData();
  for (int i = 0; i < kTestSize; ++i) {
    slnStd->SamplePhotoElectronDirection_Rejection(energy[i], cosTheta1[i] , td);
  }
  //for (int i = 0; i < kTestSize; i += kPhysDVWidth) {
  sgv->SamplePhotoElectronDirectionRejVec(energy.data(), cosTheta2.data(), kTestSize, td);
  //}

  double cumError = 0.0;
  for (int i = 0; i < kTestSize; ++i) {
    cumError += std::abs(cosTheta1[i] - cosTheta2[i]);
  }
  Printf("TestSize: %d Cumulative error: %f", kTestSize, cumError);

  CleanTaskData(td);
  delete sgv;
  delete slnStd;
}
