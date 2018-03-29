#include "PositronTo2GammaTestCommon.h"
#include "Hist.h"

const int kBasketTries = 1000;

const int kNumBins = 100;
struct PositronAnihilValidData {
  Hist gammaE;
  Hist gammaTheta;
  PositronAnihilValidData() : gammaE(0.0, 1.0, kNumBins), gammaTheta(-1.0, 1.0, kNumBins) {}
};

PositronTo2GammaModel *baseAlias;
PositronTo2GammaModel *baseRej;
VecPositronTo2GammaModel *vectorAlias;
VecPositronTo2GammaModel *vectorRej;

void FillDataVector(PositronAnihilValidData &data, bool useAlias)
{
  auto Td = PrepareTaskData();
  LightTrack_v primaries;
  VecPositronTo2GammaModel *model = useAlias ? vectorAlias : vectorRej;

  for (int bask = 0; bask < kBasketTries; ++bask) {
    PreparePrimaries(primaries, kMaxBasket);

    std::vector<double> enBeforeInteraction;
    for (int i = 0; i < primaries.GetNtracks(); ++i) {
      enBeforeInteraction.push_back(primaries.GetKinE(i) + 2. * geant::units::kElectronMassC2);
    }

    primaries.SetNtracks(kMaxBasket);

    Td->fPhysicsData->GetSecondarySOA().ClearTracks();
    model->SampleSecondariesVector(primaries, Td);

    auto &secondaries = Td->fPhysicsData->GetSecondarySOA();
    for (int i = 0; i < secondaries.GetNtracks(); ++i) {
      int primIndx    = secondaries.GetTrackIndex(i);
      double enNormed = secondaries.GetKinE(i) / enBeforeInteraction[primIndx];
      data.gammaE.Fill(enNormed);
      double gammaCost = secondaries.GetDirZ(i);
      data.gammaTheta.Fill(gammaCost);
    }
  }

  //  delete model;
  CleanTaskData(Td);
}

void FillDataScalar(PositronAnihilValidData &data, bool useAlias)
{
  auto Td = PrepareTaskData();
  std::vector<LightTrack> primaries;
  PositronTo2GammaModel *model = useAlias ? baseAlias : baseRej;

  for (int bask = 0; bask < kBasketTries; ++bask) {
    PreparePrimaries(primaries, kMaxBasket);

    std::vector<double> enBeforeInteraction;
    for (auto &positron : primaries) {
      enBeforeInteraction.push_back(positron.GetKinE() + 2. * geant::units::kElectronMassC2);
    }

    for (int i = 0; i < kMaxBasket; ++i) {
      Td->fPhysicsData->ClearSecondaries();
      model->SampleSecondaries(primaries[i], Td);

      LightTrack *secondaries = Td->fPhysicsData->GetListOfSecondaries();
      int numSec              = Td->fPhysicsData->GetNumOfSecondaries();
      int primIndex           = primaries[i].GetTrackIndex();
      for (int sec = 0; sec < numSec; ++sec) {
        double enNormed = secondaries[sec].GetKinE() / enBeforeInteraction[primIndex];
        data.gammaE.Fill(enNormed);
        double gammaCost = secondaries[sec].GetDirZ();
        data.gammaTheta.Fill(gammaCost);
      }
    }
  }

  delete model;
  CleanTaskData(Td);
}

int main()
{
  baseAlias   = PrepareAnihilModel(true);
  baseRej     = PrepareAnihilModel(false);
  vectorAlias = PrepareVecAnihilModel(true);
  vectorRej   = PrepareVecAnihilModel(false);

  Printf("Number of gamma for each test %d", kMaxBasket * kBasketTries);
  Printf("Relative histograms of kinematics (difference in percents)");
  {
    Printf("Test for alias method");

    PositronAnihilValidData scalar;
    FillDataScalar(scalar, true);
    PositronAnihilValidData vector;
    FillDataVector(vector, true);

    Printf("Normilized gamma energy");
    vector.gammaE.Compare(scalar.gammaE);

    Printf("Gamma z direction");
    vector.gammaTheta.Compare(scalar.gammaTheta);
  }
  {
    Printf("Test for rej method");

    PositronAnihilValidData scalar;
    FillDataScalar(scalar, false);
    PositronAnihilValidData vector;
    FillDataVector(vector, false);

    Printf("Normilized gamma energy");
    vector.gammaE.Compare(scalar.gammaE);

    Printf("Gamma z direction");
    vector.gammaTheta.Compare(scalar.gammaTheta);
  }

  return 0;
}
