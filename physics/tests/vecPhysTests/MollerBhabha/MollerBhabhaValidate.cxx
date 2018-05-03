#include "MollerBhabhaTestCommon.h"
#include "Hist.h"

const int kNumBins = 100;
struct MBValidData {
  Hist primEn;
  Hist primAngle;
  Hist primAzimuth;
  Hist secEn;
  Hist secAngle;
  Hist secAzimuth;
  MBValidData()
      : primEn(0.0, 1.0, kNumBins), primAngle(-1.0, 1.0, kNumBins), primAzimuth(-1.0, 1.0, kNumBins),
        secEn(0.0, 1.0, kNumBins), secAngle(-1.0, 1.0, kNumBins), secAzimuth(-1.0, 1.0, kNumBins)
  {
  }
};

void TestMBModel(geantphysics::EMModel *model, TestParticleType primary, geant::TaskData *td)
{
  MBValidData validDataVector;
  MBValidData validDataScalar;

  for (int btry = 0; btry < kBasketSamples; ++btry) {
    LightTrack_v primariesVec;
    CreateParticles(model->GetLowEnergyUsageLimit(), model->GetHighEnergyUsageLimit(), true, primary, primariesVec,
                    kBasketSize);
    std::vector<double> energyBeforeInteraction = GetE(primariesVec);

    double E0 = GetTotalE(primariesVec);
    auto P0   = GetTotalP(primariesVec);

    SampleSecondariesVector(model, primariesVec, td);

    CheckEnergyMomentumConservation(E0, P0, primariesVec, td->fPhysicsData->GetSecondarySOA());

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

    CheckEnergyMomentumConservation(E0, P0, primariesVec.data(), primariesVec.size(),
                                    td->fPhysicsData->GetListOfSecondaries(), td->fPhysicsData->GetNumOfSecondaries());

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
  Printf("Secondary reduced energy histogram");
  validDataVector.secEn.Compare(validDataScalar.secEn);

  Printf("Primary transformed Theta angle histogram");
  validDataVector.primAngle.Compare(validDataScalar.primAngle);
  Printf("Secondary transformed Theta angle histogram");
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

  std::unique_ptr<EMModel> mbRej_em =
      InitEMModel(new MollerBhabhaIonizationModel(true), kMollBhminEn, kMollBhmaxEn, false);
  std::unique_ptr<EMModel> mbTable_em =
      InitEMModel(new MollerBhabhaIonizationModel(true), kMollBhminEn, kMollBhmaxEn, true);

  std::unique_ptr<EMModel> mbRej_ep =
      InitEMModel(new MollerBhabhaIonizationModel(false), kMollBhminEn, kMollBhmaxEn, false);
  std::unique_ptr<EMModel> mbTable_ep =
      InitEMModel(new MollerBhabhaIonizationModel(false), kMollBhminEn, kMollBhmaxEn, true);

  Printf("Testing MollerBhabha rejection model for electron");
  TestMBModel(mbRej_em.get(), TestParticleType::Em, td.get());
  Printf("Testing MollerBhabha alias model for electron");
  TestMBModel(mbTable_em.get(), TestParticleType::Em, td.get());

  Printf("Testing MollerBhabha rejection model for positron");
  TestMBModel(mbRej_ep.get(), TestParticleType::Ep, td.get());
  Printf("Testing MollerBhabha alias model for positron");
  TestMBModel(mbTable_ep.get(), TestParticleType::Ep, td.get());

  CleanTaskData(td.get());
  return 0;
}
