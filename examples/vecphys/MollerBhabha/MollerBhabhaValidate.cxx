#include <Geant/Electron.h>
#include <Geant/Positron.h>
#include "MollerBhabhaTestCommon.h"
#include "Hist.h"

const int kBasketTries = 100000;

const int kNumBins = 100;
struct MBValidData {
  Hist primEn;
  Hist primAngle;
  Hist secEn;
  Hist secAngle;
  MBValidData()
      : primEn(0.0, 1.0, kNumBins), primAngle(-1.0, 1.0, kNumBins), secEn(0.0, 1.0, kNumBins),
        secAngle(-1.0, 1.0, kNumBins)
  {
  }
};

MollerBhabhaIonizationModel *aliasEl;
MollerBhabhaIonizationModel *aliasPos;
MollerBhabhaIonizationModel *rejEl;
MollerBhabhaIonizationModel *rejPos;
VecMollerBhabhaIonizationModel *vecAliasEl;
VecMollerBhabhaIonizationModel *vecAliasPos;
VecMollerBhabhaIonizationModel *vecRejEl;
VecMollerBhabhaIonizationModel *vecRejPos;

TaskData *td;

void FillDataVector(MBValidData &data, bool useAlias)
{
  LightTrack_v primaries;
  VecMollerBhabhaIonizationModel *model = useAlias ? vecAliasEl : vecRejEl;

  for (int bask = 0; bask < kBasketTries; ++bask) {
    PreparePrimaries(primaries, kMaxBasket);

    std::vector<double> enBeforeInteraction;
    enBeforeInteraction.insert(enBeforeInteraction.begin(), primaries.GetKinEArr(),
                               primaries.GetKinEArr() + kMaxBasket);

    primaries.SetNtracks(kMaxBasket);

    td->fPhysicsData->GetSecondarySOA().ClearTracks();
    model->SampleSecondariesVector(primaries, td);

    for (int i = 0; i < primaries.GetNtracks(); ++i) {
      int primIdx     = primaries.GetTrackIndex(i);
      double enNormed = primaries.GetKinE(i) / enBeforeInteraction[primIdx];
      data.primEn.Fill(enNormed);
      data.primAngle.Fill(primaries.GetDirZ(i));
    }

    auto &secondaries = td->fPhysicsData->GetSecondarySOA();
    for (int i = 0; i < secondaries.GetNtracks(); ++i) {
      int primIdx     = secondaries.GetTrackIndex(i);
      double enNormed = secondaries.GetKinE(i) / enBeforeInteraction[primIdx];
      data.secEn.Fill(enNormed);
      data.secAngle.Fill(secondaries.GetDirZ(i));
    }
  }
}

void FillDataScalar(MBValidData &data, bool useAlias)
{
  std::vector<LightTrack> primaries;
  MollerBhabhaIonizationModel *model = useAlias ? aliasEl : rejEl;

  for (int bask = 0; bask < kBasketTries; ++bask) {
    PreparePrimaries(primaries, kMaxBasket);

    std::vector<double> enBeforeInteraction;
    for (auto &lep : primaries) {
      enBeforeInteraction.push_back(lep.GetKinE());
    }

    for (int i = 0; i < kMaxBasket; ++i) {
      td->fPhysicsData->ClearSecondaries();
      model->SampleSecondaries(primaries[i], td);
      int primIdx = primaries[i].GetTrackIndex();
      data.primEn.Fill(primaries[i].GetKinE() / enBeforeInteraction[primIdx]);
      data.primAngle.Fill(primaries[i].GetDirZ());
      LightTrack *secondaries = td->fPhysicsData->GetListOfSecondaries();
      int numSec              = td->fPhysicsData->GetNumOfSecondaries();
      for (int sec = 0; sec < numSec; ++sec) {
        int primIdx     = secondaries[sec].GetTrackIndex();
        double enNormed = secondaries[sec].GetKinE() / enBeforeInteraction[primIdx];
        data.secEn.Fill(enNormed);
        data.secAngle.Fill(secondaries[sec].GetDirZ());
      }
    }
  }
}

int main()
{

  PrepareWorld();
  rejEl       = PrepareMBModel(false, true);
  rejPos      = PrepareMBModel(false, false);
  aliasEl     = PrepareMBModel(true, true);
  aliasPos    = PrepareMBModel(true, false);
  vecRejEl    = PrepareVecMBModel(false, true);
  vecRejPos   = PrepareVecMBModel(false, false);
  vecAliasEl  = PrepareVecMBModel(true, true);
  vecAliasPos = PrepareVecMBModel(true, false);
  td          = PrepareTaskData();

  Printf("Number of leptons for each test %d", kMaxBasket * kBasketTries);
  Printf("Relative histograms of kinematics (difference in percents)");
  {
    Printf("Test for alias method");

    MBValidData scalar;
    FillDataScalar(scalar, true);
    MBValidData vector;
    FillDataVector(vector, true);

    Printf("====Prim EN====");
    vector.primEn.Compare(scalar.primEn);
    Printf("====Sec EN====");
    vector.secEn.Compare(scalar.secEn);
    Printf("====Prim Dir====");
    vector.primAngle.Compare(scalar.primAngle);
    Printf("====Sec Dir====");
    vector.secAngle.Compare(scalar.secAngle);
  }
  {
    Printf("Test for rej method");

    MBValidData scalar;
    FillDataScalar(scalar, false);
    MBValidData vector;
    FillDataVector(vector, false);

    Printf("====Prim EN====");
    vector.primEn.Compare(scalar.primEn);
    Printf("====Sec EN====");
    vector.secEn.Compare(scalar.secEn);
    Printf("====Prim Dir====");
    vector.primAngle.Compare(scalar.primAngle);
    Printf("====Sec Dir====");
    vector.secAngle.Compare(scalar.secAngle);
  }

  return 0;
}
