#include "SGTestCommon.h"
#include "Hist.h"

const int kNumBins = 100;
struct SGValidData {
  Hist primEn;
  Hist primAngle;
  Hist primAzimuth;
  Hist secEn;
  Hist secAngle;
  Hist secAzimuth;
  SGValidData()
      : primEn(0.0, 1.0, kNumBins), primAngle(-1.0, 1.0, kNumBins), primAzimuth(-1.0, 1.0, kNumBins),
        secEn(0.0, 1.0, kNumBins), secAngle(-1.0, 1.0, kNumBins), secAzimuth(-1.0, 1.0, kNumBins)
  {
  }
};

void TestSGModel(geantphysics::EMModel *model, TestParticleType primary, geant::TaskData *td)
{
  SGValidData validDataVector;
  SGValidData validDataScalar;

  for (int btry = 0; btry < kBasketSamples; ++btry) {
    LightTrack_v primariesVec;
    CreateParticles(model->GetLowEnergyUsageLimit(), model->GetHighEnergyUsageLimit(), true, primary, primariesVec,
                    kBasketSize);
    std::vector<double> energyBeforeInteraction = GetE(primariesVec);

    double E0 = GetTotalE(primariesVec);
    auto P0   = GetTotalP(primariesVec);

    SampleSecondariesVector(model, primariesVec, td);

    //CheckEnergyMomentumConservation(E0, P0, primariesVec, td->fPhysicsData->GetSecondarySOA());

    // Fill histogram
    for (int i = 0; i < primariesVec.GetNtracks(); ++i) {
      double enNormed = primariesVec.GetKinE(i) / energyBeforeInteraction[primariesVec.GetTrackIndex(i)];
      double dirZ     = primariesVec.GetDirZ(i);
      double azimuth  = XYDirToAzumuth(primariesVec.GetDirX(i), primariesVec.GetDirY(i));

      validDataVector.primEn.Fill(enNormed);
      validDataVector.primAngle.Fill(dirZ);
      if (azimuth != kWrongVal) validDataVector.primAzimuth.Fill(azimuth);
    }
    auto &sec = td->fPhysicsData->GetSecondarySOA();
    for (int i = 0; i < sec.GetNtracks(); ++i) {
      double enNormed = sec.GetKinE(i) / energyBeforeInteraction[sec.GetTrackIndex(i)];
      double dirZ     = sec.GetDirZ(i);
      double azimuth  = XYDirToAzumuth(sec.GetDirX(i), sec.GetDirY(i));

      validDataVector.secEn.Fill(enNormed);
      validDataVector.secAngle.Fill(dirZ);
      if (azimuth != kWrongVal) validDataVector.secAzimuth.Fill(azimuth);
    }
  }

  for (int btry = 0; btry < kBasketSamples; ++btry) {
    std::vector<LightTrack> primariesVec;
    CreateParticles(model->GetLowEnergyUsageLimit(), model->GetHighEnergyUsageLimit(), true, primary, primariesVec,
                    kBasketSize);
    std::vector<double> energyBeforeInteraction = GetE(primariesVec);

    double E0 = GetTotalE(primariesVec.data(), primariesVec.size());
    auto P0   = GetTotalP(primariesVec.data(), primariesVec.size());

    SampleSecondariesScalar(model, primariesVec, td);

    //CheckEnergyMomentumConservation(E0, P0, primariesVec.data(), primariesVec.size(),
                                    //td->fPhysicsData->GetListOfSecondaries(), td->fPhysicsData->GetNumOfSecondaries());
    //exit(-1);

    for (size_t i = 0; i < primariesVec.size(); ++i) {
      double enNormed = primariesVec[i].GetKinE() / energyBeforeInteraction[primariesVec[i].GetTrackIndex()];
      double dirZ     = primariesVec[i].GetDirZ();
      double azimuth  = XYDirToAzumuth(primariesVec[i].GetDirX(), primariesVec[i].GetDirY());

      validDataScalar.primEn.Fill(enNormed);
      validDataScalar.primAngle.Fill(dirZ);
      if (azimuth != kWrongVal) validDataScalar.primAzimuth.Fill(azimuth);
    }
    LightTrack *sec = td->fPhysicsData->GetListOfSecondaries();
    for (int i = 0; i < td->fPhysicsData->GetNumOfSecondaries(); ++i) {
      double enNormed = sec[i].GetKinE() / energyBeforeInteraction[sec[i].GetTrackIndex()];
      double dirZ     = sec[i].GetDirZ();
      double azimuth  = XYDirToAzumuth(sec[i].GetDirX(), sec[i].GetDirY());

      validDataScalar.secEn.Fill(enNormed);
      validDataScalar.secAngle.Fill(dirZ);
      if (azimuth != kWrongVal) validDataScalar.secAzimuth.Fill(azimuth);
    }
  }

  Printf("Comparing histograms");
  Printf("Primary reduced energy histogram");
  validDataVector.primEn.Compare(validDataScalar.primEn);
  Printf("Secondary electron energy histogram");
  validDataVector.secEn.Compare(validDataScalar.secEn);

  Printf("Primary transformed Theta angle histogram");
  validDataVector.primAngle.Compare(validDataScalar.primAngle);
  Printf("Secondary electron Theta angle histogram");
  validDataVector.secAngle.Compare(validDataScalar.secAngle);

  Printf("Primary transformed Phi angle histogram");
  validDataVector.primAzimuth.Compare(validDataScalar.primAzimuth);
  Printf("Secondary transformed Phi angle histogram");
  validDataVector.secAzimuth.Compare(validDataScalar.secAzimuth);
}

int main()
{

  PrepareWorld();
  auto td = PrepareTaskData();

  std::unique_ptr<EMModel> sgRej =
      InitEMModel(new SauterGavrilaPhotoElectricModel(), kSauGavminEn, kSauGavmaxEn, false);
  std::unique_ptr<EMModel> sgTable =
      InitEMModel(new SauterGavrilaPhotoElectricModel(), kSauGavminEn, kSauGavmaxEn, true);

  Printf("Testing SauterGavrila rejection model");
  TestSGModel(sgRej.get(), TestParticleType::Gamma, td.get());
  Printf("Testing SauterGavrila alias model");
  TestSGModel(sgTable.get(), TestParticleType::Gamma, td.get());

  CleanTaskData(td.get());
  return 0;
}
