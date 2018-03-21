/*
 * Basic sanity check for KleinNishina physics.
 * Checks first two moments of gamma and electron energy distributions
 */
#include "KleinNishinaTestCommon.h"

const double kGammaEn0 = geant::units::GeV;
const int kBasketTries = 10000;

void FixedGammaTestVector(double &gammaEnMean, double &gammaEn2mean, double &emEnMean, double &emEn2Mean, bool useAlias)
{
  auto Td = PrepareTaskData();
  LightTrack_v primaries;

  int Ngamma   = 0;
  int Nem      = 0;
  gammaEnMean  = 0.0;
  gammaEn2mean = 0.0;
  emEnMean     = 0.0;
  emEn2Mean    = 0.0;

  VecKleinNishinaComptonModel *kNish = PrepareVecKnishinaModel(useAlias);

  for (int t = 0; t < kBasketTries; t++) {
    PreparePrimaries(primaries, kMaxBasket);
    primaries.SetNtracks(kMaxBasket);
    for (int i = 0; i < kMaxBasket; ++i) {
      primaries.SetKinE(kGammaEn0, i);
    }
    Td->fPhysicsData->GetSecondarySOA().ClearTracks();
    kNish->SampleSecondariesVector(primaries, Td);

    for (int i = 0; i < kMaxBasket; ++i) {
      double en = primaries.GetKinE(i);
      gammaEnMean += en;
      gammaEn2mean += en * en;
      Ngamma++;
    }
    auto &secondaries = Td->fPhysicsData->GetSecondarySOA();
    for (int i = 0; i < secondaries.GetNtracks(); ++i) {
      double en = secondaries.GetKinE(i);
      emEnMean += en;
      emEn2Mean += en * en;
      Nem++;
    }
  }

  gammaEnMean /= Ngamma;
  gammaEn2mean /= Ngamma;

  emEn2Mean /= Nem;
  emEnMean /= Nem;

  delete kNish;
  CleanTaskData(Td);
}

void FixedGammaTestScalar(double &gammaEnMean, double &gammaEn2mean, double &emEnMean, double &emEn2Mean, bool useAlias)
{
  auto Td = PrepareTaskData();
  std::vector<LightTrack> primaries;

  int Ngamma   = 0;
  int Nem      = 0;
  gammaEnMean  = 0.0;
  gammaEn2mean = 0.0;
  emEnMean     = 0.0;
  emEn2Mean    = 0.0;

  KleinNishinaComptonModel *kNish = PrepareVecKnishinaModel(useAlias);

  for (int t = 0; t < kBasketTries; t++) {
    PreparePrimaries(primaries, kMaxBasket);
    for (int i = 0; i < kMaxBasket; ++i) {
      primaries[i].SetKinE(kGammaEn0);
    }
    Td->fPhysicsData->ClearSecondaries();

    for (int i = 0; i < kMaxBasket; ++i) {
      kNish->SampleSecondaries(primaries[i], Td);
    }

    for (int i = 0; i < kMaxBasket; ++i) {
      double en = primaries[i].GetKinE();
      gammaEnMean += en;
      gammaEn2mean += en * en;
      Ngamma++;
    }
    LightTrack *secondaries = Td->fPhysicsData->GetListOfSecondaries();
    for (int i = 0; i < Td->fPhysicsData->GetNumOfSecondaries(); ++i) {
      double en = secondaries[i].GetKinE();
      emEnMean += en;
      emEn2Mean += en * en;
      Nem++;
    }
  }

  gammaEnMean /= Ngamma;
  gammaEn2mean /= Ngamma;

  emEn2Mean /= Nem;
  emEnMean /= Nem;

  delete kNish;
  CleanTaskData(Td);
}

int main()
{
  {
    Printf("NTracks: %d", kMaxBasket * kBasketTries);
    double vectorE    = 0.0;
    double vectorE2   = 0.0;
    double vectorEmE  = 0.0;
    double vectorEmE2 = 0.0;

    FixedGammaTestVector(vectorE, vectorE2, vectorEmE, vectorEmE2, true);
    Printf("VectorKNishAlias En after gamma <En>: %f <En^2>: %f em <En>: %f <En^2>: %f", vectorE, vectorE2, vectorEmE,
           vectorEmE2);

    double scalarE    = 0.0;
    double scalarE2   = 0.0;
    double scalarEmE  = 0.0;
    double scalarEmE2 = 0.0;
    FixedGammaTestScalar(scalarE, scalarE2, scalarEmE, scalarEmE2, true);
    Printf("ScalarKNishAlias En after gamma <En>: %f <En^2>: %f em <En>: %f <En^2>: %f", scalarE, scalarE2, scalarEmE,
           scalarEmE2);
  }
  {
    Printf("NTracks: %d", kMaxBasket * kBasketTries);
    double vectorE    = 0.0;
    double vectorE2   = 0.0;
    double vectorEmE  = 0.0;
    double vectorEmE2 = 0.0;

    FixedGammaTestVector(vectorE, vectorE2, vectorEmE, vectorEmE2, false);
    Printf("VectorKNishRej En after gamma <En>: %f <En^2>: %f em <En>: %f <En^2>: %f", vectorE, vectorE2, vectorEmE,
           vectorEmE2);

    double scalarE    = 0.0;
    double scalarE2   = 0.0;
    double scalarEmE  = 0.0;
    double scalarEmE2 = 0.0;
    FixedGammaTestScalar(scalarE, scalarE2, scalarEmE, scalarEmE2, false);
    Printf("ScalarKNishRej En after gamma <En>: %f <En^2>: %f em <En>: %f <En^2>: %f", scalarE, scalarE2, scalarEmE,
           scalarEmE2);
  }
  return 0;
}
