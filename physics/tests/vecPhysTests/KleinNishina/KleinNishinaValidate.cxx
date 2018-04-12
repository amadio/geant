/*
 * Basic sanity check for KleinNishina physics.
 * Checks first two moments of gamma and electron energy distributions
 */
#include "KleinNishinaTestCommon.h"
#include "Hist.h"

const int kBasketTries = 100000;

const int kNumBins = 100;
struct ComptonValidData {
  Hist gammaE;
  Hist gammaTheta;
  Hist emE;
  Hist emTheta;
  ComptonValidData()
      : gammaE(0.0, 1.0, kNumBins), gammaTheta(-1.0, 1.0, kNumBins), emE(0.0, 1.0, kNumBins),
        emTheta(-1.0, 1.0, kNumBins)
  {
  }
};

void FillDataVector(ComptonValidData &data, bool useAlias)
{
  auto Td = PrepareTaskData();
  LightTrack_v primaries;
  VecKleinNishinaComptonModel *kNish = PrepareVecKnishinaModel(useAlias);

  for (int bask = 0; bask < kBasketTries; ++bask) {
    PreparePrimaries(primaries, kMaxBasket);

    std::vector<double> enBeforeInteraction;
    enBeforeInteraction.insert(enBeforeInteraction.begin(), primaries.GetKinEArr(),
                               primaries.GetKinEArr() + kMaxBasket);

    primaries.SetNtracks(kMaxBasket);

    Td->fPhysicsData->GetSecondarySOA().ClearTracks();
    kNish->SampleSecondariesVector(primaries, Td);

    for (int i = 0; i < kMaxBasket; ++i) {
      double enNormed = primaries.GetKinE(i) / enBeforeInteraction[i];
      if (enNormed > 0.0) {
        data.gammaE.Fill(enNormed);
        double gammaCost = primaries.GetDirZ(i);
        data.gammaTheta.Fill(gammaCost);
      }
    }
    auto &secondaries = Td->fPhysicsData->GetSecondarySOA();
    for (int i = 0; i < secondaries.GetNtracks(); ++i) {
      double enNormed = secondaries.GetKinE(i) / enBeforeInteraction[i];
      data.emE.Fill(enNormed);
      double emCost = secondaries.GetDirZ(i);
      data.emTheta.Fill(emCost);
    }
  }

  //  delete kNish;
  CleanTaskData(Td);
}

void FillDataScalar(ComptonValidData &data, bool useAlias)
{
  auto Td = PrepareTaskData();
  std::vector<LightTrack> primaries;
  KleinNishinaComptonModel *kNish = PrepareKnishinaModel(useAlias);

  for (int bask = 0; bask < kBasketTries; ++bask) {
    PreparePrimaries(primaries, kMaxBasket);

    std::vector<double> enBeforeInteraction;
    for (auto &gamma : primaries) {
      enBeforeInteraction.push_back(gamma.GetKinE());
    }

    for (int i = 0; i < kMaxBasket; ++i) {
      Td->fPhysicsData->ClearSecondaries();
      kNish->SampleSecondaries(primaries[i], Td);

      double enNormed = primaries[i].GetKinE() / enBeforeInteraction[i];
      if (enNormed > 0.0) {
        data.gammaE.Fill(enNormed);
        double gammaCost = primaries[i].GetDirZ();
        data.gammaTheta.Fill(gammaCost);
      }

      LightTrack *secondaries = Td->fPhysicsData->GetListOfSecondaries();
      int numSec              = Td->fPhysicsData->GetNumOfSecondaries();
      for (int sec = 0; sec < numSec; ++sec) {
        double enNormed = secondaries[sec].GetKinE() / enBeforeInteraction[i];
        data.emE.Fill(enNormed);
        double emCost = secondaries[sec].GetDirZ();
        data.emTheta.Fill(emCost);
      }
    }
  }

  delete kNish;
  CleanTaskData(Td);
}

int main()
{
  Printf("Number of gamma for each test %d", kMaxBasket * kBasketTries);
  {
    Printf("Test for alias method");

    ComptonValidData scalar;
    FillDataScalar(scalar, true);
    ComptonValidData vector;
    FillDataVector(vector, true);

    Printf("Normilized gamma energy");
    vector.gammaE.Compare(scalar.gammaE);

    Printf("Normilized electron energy");
    vector.emE.Compare(scalar.emE);

    Printf("Gamma z direction");
    vector.gammaTheta.Compare(scalar.gammaTheta);

    Printf("Electron z direction");
    vector.emTheta.Compare(scalar.emTheta);
  }
  {
    Printf("Test for rej method");

    ComptonValidData scalar;
    FillDataScalar(scalar, false);
    ComptonValidData vector;
    FillDataVector(vector, false);

    Printf("Normilized gamma energy");
    vector.gammaE.Compare(scalar.gammaE);

    Printf("Normilized electron energy");
    vector.emE.Compare(scalar.emE);

    Printf("Gamma z direction");
    vector.gammaTheta.Compare(scalar.gammaTheta);

    Printf("Electron z direction");
    vector.emTheta.Compare(scalar.emTheta);
  }

  return 0;
}
