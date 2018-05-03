#include <Geant/Electron.h>
#include <Geant/Positron.h>
#include "VectorPairTestCommon.h"
#include "Hist.h"

const int kNumBins = 100;
struct PairValidData {
  Hist primEn;
  Hist primAngle;
  Hist primAzimuth;
  Hist secEn;
  Hist secAngle;
  Hist secAzimuth;
  PairValidData()
      : primEn(0.0, 1.0, kNumBins), primAngle(-12.0, 0.0, kNumBins), primAzimuth(-1.0, 1.0, kNumBins),
        secEn(0.0, 1.0, kNumBins), secAngle(-12.0, 0.0, kNumBins), secAzimuth(-1.0, 1.0, kNumBins)
  {
  }
};

void TestPairModel(geantphysics::EMModel *model, geant::TaskData *td)
{
  PairValidData validDataVector;
  PairValidData validDataScalar;

  for (int btry = 0; btry < kBasketSamples; ++btry) {
    LightTrack_v primariesVec;
    CreateParticles(model->GetLowEnergyUsageLimit(), model->GetHighEnergyUsageLimit(), true, TestParticleType::Gamma,
                    primariesVec, kBasketSize);
    std::vector<double> energyBeforeInteraction = GetE(primariesVec);

    double E0 = GetTotalE(primariesVec);
    auto P0   = GetTotalP(primariesVec);

    SampleSecondariesVector(model, primariesVec, td);

    CheckEnergyMomentumConservation(E0, P0, primariesVec, td->fPhysicsData->GetSecondarySOA());

    // Fill histogram
    auto &sec = td->fPhysicsData->GetSecondarySOA();
    for (int i = 0; i < sec.GetNtracks(); ++i) {
      double enNormed = (sec.GetKinE(i) + sec.GetMass(i)) / energyBeforeInteraction[sec.GetTrackIndex(i)];
      double logElz   = ZDirLog(sec.GetDirZ(i));
      double azimuth  = XYDirToAzumuth(sec.GetDirX(i), sec.GetDirY(i));

      if (sec.GetGVcode(i) == geantphysics::Electron::Definition()->GetInternalCode()) {
        validDataVector.primEn.Fill(enNormed);
        if (logElz != kWrongVal) validDataVector.primAngle.Fill(logElz);
        if (azimuth != kWrongVal) validDataVector.primAzimuth.Fill(azimuth);
      }
      if (sec.GetGVcode(i) == geantphysics::Electron::Definition()->GetInternalCode()) {
        validDataVector.secEn.Fill(enNormed);
        if (logElz != kWrongVal) validDataVector.secAngle.Fill(logElz);
        if (azimuth != kWrongVal) validDataVector.secAzimuth.Fill(azimuth);
      }
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

    // Fill histogram
    LightTrack *sec = td->fPhysicsData->GetListOfSecondaries();
    for (int i = 0; i < td->fPhysicsData->GetNumOfSecondaries(); ++i) {
      double enNormed = (sec[i].GetKinE() + sec[i].GetMass()) / energyBeforeInteraction[sec[i].GetTrackIndex()];
      double logElz   = ZDirLog(sec[i].GetDirZ());
      double azimuth  = XYDirToAzumuth(sec[i].GetDirX(), sec[i].GetDirY());

      if (sec[i].GetGVcode() == geantphysics::Electron::Definition()->GetInternalCode()) {
        validDataScalar.primEn.Fill(enNormed);
        if (logElz != kWrongVal) validDataScalar.primAngle.Fill(logElz);
        if (azimuth != kWrongVal) validDataScalar.primAzimuth.Fill(azimuth);
      }
      if (sec[i].GetGVcode() == geantphysics::Electron::Definition()->GetInternalCode()) {
        validDataScalar.secEn.Fill(enNormed);
        if (logElz != kWrongVal) validDataScalar.secAngle.Fill(logElz);
        if (azimuth != kWrongVal) validDataScalar.secAzimuth.Fill(azimuth);
      }
    }
  }

  Printf("Comparing histograms");
  Printf("Electron reduced energy histogram");
  validDataVector.primEn.Compare(validDataScalar.primEn);
  Printf("Positron reduced energy histogram");
  validDataVector.secEn.Compare(validDataScalar.secEn);

  Printf("Electron transformed Theta angle histogram");
  validDataVector.primAngle.Compare(validDataScalar.primAngle);
  Printf("Positron transformed Theta angle histogram");
  validDataVector.secAngle.Compare(validDataScalar.secAngle);

  Printf("Electron transformed Phi angle histogram");
  validDataVector.primAzimuth.Compare(validDataScalar.primAzimuth);
  Printf("Positron transformed Phi angle histogram");
  validDataVector.secAzimuth.Compare(validDataScalar.secAzimuth);
}

int main()
{
  PrepareWorld();
  auto td = PrepareTaskData();

  std::unique_ptr<EMModel> bhRej   = InitEMModel(new BetheHeitlerPairModel, kBHminEn, kBHmaxEn, false);
  std::unique_ptr<EMModel> bhTable = InitEMModel(new BetheHeitlerPairModel, kBHminEn, kBHmaxEn, true);

  std::unique_ptr<EMModel> relPairRej   = InitEMModel(new RelativisticPairModel, kRelPairMinEn, kRelPairMaxEn, false);
  std::unique_ptr<EMModel> relPairTable = InitEMModel(new RelativisticPairModel, kRelPairMinEn, kRelPairMaxEn, true);

  Printf("Testing RelativisticPair alias model");
  TestPairModel(relPairTable.get(), td.get());
  Printf("Testing RelativisticPair rej model");
  TestPairModel(relPairRej.get(), td.get());

  Printf("Testing BetheHeitler alias model");
  TestPairModel(bhTable.get(), td.get());
  Printf("Testing BetheHeitler rej model");
  TestPairModel(bhRej.get(), td.get());

  CleanTaskData(td.get());
  return 0;
}
