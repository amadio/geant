#include <Geant/Electron.h>
#include <Geant/Positron.h>
#include "SeltzerBergerTestCommon.h"
#include "Hist.h"

const int kBasketTries = 100000;

const int kNumBins = 100;
struct SBValidData {
  Hist primEn;
  Hist primAngle;
  Hist secEn;
  Hist secAngle;
  SBValidData()
      : primEn(0.0, 1.0, kNumBins), primAngle(-1.0, 1.0, kNumBins), secEn(0.0, 1.0, kNumBins),
        secAngle(-1.0, 1.0, kNumBins)
  {
  }
};

SeltzerBergerBremsModel *aliasEl;
SeltzerBergerBremsModel *aliasPos;
SeltzerBergerBremsModel *rejEl;
SeltzerBergerBremsModel *rejPos;
VecSeltzerBergerBremsModel *vecAliasEl;
VecSeltzerBergerBremsModel *vecAliasPos;
VecSeltzerBergerBremsModel *vecRejEl;
VecSeltzerBergerBremsModel *vecRejPos;

TaskData *td;

void FillDataVector(SBValidData &data, bool useAlias, bool forElectron)
{
  LightTrack_v primaries;
  VecSeltzerBergerBremsModel *model;
  if (forElectron) {
    model = useAlias ? vecAliasEl : vecRejEl;
  } else {
    model = useAlias ? vecAliasPos : vecRejPos;
  }

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

void FillDataScalar(SBValidData &data, bool useAlias, bool forElectron)
{
  std::vector<LightTrack> primaries;
  SeltzerBergerBremsModel *model;
  if (forElectron) {
    model = useAlias ? aliasEl : rejEl;
  } else {
    model = useAlias ? aliasPos : rejPos;
  }

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
  rejEl       = PrepareSBModel(false, true);
  rejPos      = PrepareSBModel(false, false);
  aliasEl     = PrepareSBModel(true, true);
  aliasPos    = PrepareSBModel(true, false);
  vecRejEl    = PrepareVecSBModel(false, true);
  vecRejPos   = PrepareVecSBModel(false, false);
  vecAliasEl  = PrepareVecSBModel(true, true);
  vecAliasPos = PrepareVecSBModel(true, false);
  td          = PrepareTaskData();

  Printf("Number of leptons for each test %d", kMaxBasket * kBasketTries);
  {
    Printf("Test for alias method electron");

    SBValidData scalar;
    FillDataScalar(scalar, true, true);
    SBValidData vector;
    FillDataVector(vector, true, true);

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
    Printf("Test for alias method positron");

    SBValidData scalar;
    FillDataScalar(scalar, true, false);
    SBValidData vector;
    FillDataVector(vector, true, false);

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
    Printf("Test for rej method electron");

    SBValidData scalar;
    FillDataScalar(scalar, false, true);
    SBValidData vector;
    FillDataVector(vector, false, true);

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
    Printf("Test for rej method positron");

    SBValidData scalar;
    FillDataScalar(scalar, false, false);
    SBValidData vector;
    FillDataVector(vector, false, false);

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
