#include "PositronTo2GammaTestCommon.h"
#include "Hist.h"

const int kNumBins = 100;
struct PositronAnihilValidData {
  Hist gammaE;
  Hist gammaTheta;
  Hist gammaAzimuth;
  PositronAnihilValidData()
      : gammaE(0.0, 1.0, kNumBins), gammaTheta(-1.0, 1.0, kNumBins), gammaAzimuth(-1.0, 1.0, kNumBins)
  {
  }
};

void TestPos2GammaModel(geantphysics::EMModel *model, geant::TaskData *td)
{
  PositronAnihilValidData validDataVector;
  PositronAnihilValidData validDataScalar;

  for (int btry = 0; btry < kBasketSamples; ++btry) {
    LightTrack_v primariesVec;
    CreateParticles(model->GetLowEnergyUsageLimit(), model->GetHighEnergyUsageLimit(), true, TestParticleType::Ep,
                    primariesVec, kBasketSize);
    std::vector<double> energyBeforeInteraction = GetE(primariesVec);

    double E0 = GetTotalE(primariesVec);
    auto P0   = GetTotalP(primariesVec);

    SampleSecondariesVector(model, primariesVec, td);

    CheckEnergyMomentumConservation(E0, P0, primariesVec, td->fPhysicsData->GetSecondarySOA());

    // Fill histogram
    auto &sec = td->fPhysicsData->GetSecondarySOA();
    for (int i = 0; i < sec.GetNtracks(); ++i) {
      double enNormed = sec.GetKinE(i) / energyBeforeInteraction[sec.GetTrackIndex(i)];
      double dirZ     = sec.GetDirZ(i);
      double azimuth  = XYDirToAzumuth(sec.GetDirX(i), sec.GetDirY(i));

      validDataVector.gammaE.Fill(enNormed);
      validDataVector.gammaTheta.Fill(dirZ);
      if (azimuth != kWrongVal) validDataVector.gammaAzimuth.Fill(azimuth);
    }
  }

  for (int btry = 0; btry < kBasketSamples; ++btry) {
    std::vector<LightTrack> primariesVec;
    CreateParticles(model->GetLowEnergyUsageLimit(), model->GetHighEnergyUsageLimit(), true, TestParticleType::Gamma,
                    primariesVec, kBasketSize);
    std::vector<double> energyBeforeInteraction = GetE(primariesVec);

    double E0 = GetTotalE(primariesVec.data(), primariesVec.size());
    auto P0   = GetTotalP(primariesVec.data(), primariesVec.size());

    SampleSecondariesScalar(model, primariesVec, td);

    CheckEnergyMomentumConservation(E0, P0, primariesVec.data(), primariesVec.size(),
                                    td->fPhysicsData->GetListOfSecondaries(), td->fPhysicsData->GetNumOfSecondaries());

    LightTrack *sec = td->fPhysicsData->GetListOfSecondaries();
    for (int i = 0; i < td->fPhysicsData->GetNumOfSecondaries(); ++i) {
      double enNormed = sec[i].GetKinE() / energyBeforeInteraction[sec[i].GetTrackIndex()];
      double dirZ     = sec[i].GetDirZ();
      double azimuth  = XYDirToAzumuth(sec[i].GetDirX(), sec[i].GetDirY());

      validDataScalar.gammaE.Fill(enNormed);
      validDataScalar.gammaTheta.Fill(dirZ);
      if (azimuth != kWrongVal) validDataScalar.gammaAzimuth.Fill(azimuth);
    }
  }

  Printf("Comparing histograms");
  Printf("Gamma reduced energy histogram");
  validDataVector.gammaE.Compare(validDataScalar.gammaE);
  Printf("Gamma transformed Theta angle histogram");
  validDataVector.gammaTheta.Compare(validDataScalar.gammaTheta);

  Printf("Lepton transformed Phi angle histogram");
  validDataVector.gammaAzimuth.Compare(validDataScalar.gammaAzimuth);
}

int main()
{
  PrepareWorld();
  auto td = PrepareTaskData();

  std::unique_ptr<EMModel> pos2gRej   = InitEMModel(new PositronTo2GammaModel, kPos2GammaMinEn, kPos2GammaMaxEn, false);
  std::unique_ptr<EMModel> pos2gTable = InitEMModel(new PositronTo2GammaModel, kPos2GammaMinEn, kPos2GammaMaxEn, true);

  Printf("Testing PositronTo2Gamma alias model");
  TestPos2GammaModel(pos2gTable.get(), td.get());
  Printf("Testing PositronTo2Gamma rej model");
  TestPos2GammaModel(pos2gRej.get(), td.get());

  CleanTaskData(td.get());
  return 0;
}
