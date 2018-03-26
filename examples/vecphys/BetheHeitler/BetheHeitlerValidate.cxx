#include <Geant/Electron.h>
#include <Geant/Positron.h>
#include "BetheHeitlerTestCommon.h"
#include "Hist.h"

const int kBasketTries = 1000;

const int kNumBins = 100;
struct BHValidData {
  Hist elEn;
  Hist elAngle;
  Hist posEn;
  Hist posAngle;
  BHValidData()
      : elEn(0.0, 1.0, kNumBins), elAngle(-12.0, 0.0, kNumBins), posEn(0.0, 1.0, kNumBins),
        posAngle(-12.0, 0.0, kNumBins)
  {
  }
};

BetheHeitlerPairModel *bHModelRej;
VecBetheHeitlerPairModel *bHVecModelRej;
BetheHeitlerPairModel *bHModelAlias;
VecBetheHeitlerPairModel *bHVecModelAlias;
TaskData *td;

void FillDataVector(BHValidData &data, bool useAlias)
{
  LightTrack_v primaries;
  BetheHeitlerPairModel *model = useAlias ? bHVecModelAlias : bHVecModelRej;

  for (int bask = 0; bask < kBasketTries; ++bask) {
    PreparePrimaries(primaries, kMaxBasket);

    std::vector<double> enBeforeInteraction;
    enBeforeInteraction.insert(enBeforeInteraction.begin(), primaries.GetKinEVec(),
                               primaries.GetKinEVec() + kMaxBasket);

    primaries.SetNtracks(kMaxBasket);

    td->fPhysicsData->GetSecondarySOA().ClearTracks();
    model->SampleSecondariesVector(primaries, td);

    auto &secondaries = td->fPhysicsData->GetSecondarySOA();
    for (int i = 0; i < secondaries.GetNtracks(); ++i) {
      int primIdx     = secondaries.GetTrackIndex(i);
      double enNormed = (secondaries.GetKinE(i) + secondaries.GetMass(i)) / enBeforeInteraction[primIdx];

      if (secondaries.GetGVcode(i) == geantphysics::Electron::Definition()->GetInternalCode()) {

        data.elEn.Fill(enNormed);
        double elZ = 0.5 * (1.0 - secondaries.GetDirZ(i));
        if (elZ > 0 && log10(elZ) > -12.0) {
          data.elAngle.Fill(log10(elZ));
        }

      } else if (secondaries.GetGVcode(i) == geantphysics::Positron::Definition()->GetInternalCode()) {

        data.posEn.Fill(enNormed);

        double posZ = 0.5 * (1.0 - secondaries.GetDirZ(i));
        if (posZ > 0 && log10(posZ) > -12.0) {
          data.posAngle.Fill(log10(posZ));
        }
      }
    }
  }
}

void FillDataScalar(BHValidData &data, bool useAlias)
{
  std::vector<LightTrack> primaries;
  BetheHeitlerPairModel *model = useAlias ? bHModelAlias : bHModelRej;

  for (int bask = 0; bask < kBasketTries; ++bask) {
    PreparePrimaries(primaries, kMaxBasket);

    std::vector<double> enBeforeInteraction;
    for (auto &gamma : primaries) {
      enBeforeInteraction.push_back(gamma.GetKinE());
    }

    for (int i = 0; i < kMaxBasket; ++i) {
      td->fPhysicsData->ClearSecondaries();
      model->SampleSecondaries(primaries[i], td);

      LightTrack *secondaries = td->fPhysicsData->GetListOfSecondaries();
      int numSec              = td->fPhysicsData->GetNumOfSecondaries();
      for (int sec = 0; sec < numSec; ++sec) {
        int primIdx     = secondaries[sec].GetTrackIndex();
        double enNormed = (secondaries[sec].GetKinE() + secondaries[sec].GetMass()) / enBeforeInteraction[primIdx];

        if (secondaries[sec].GetGVcode() == geantphysics::Electron::Definition()->GetInternalCode()) {

          data.elEn.Fill(enNormed);
          double elZ = 0.5 * (1.0 - secondaries[sec].GetDirZ());
          if (elZ > 0 && log10(elZ) > -12.0) {
            data.elAngle.Fill(log10(elZ));
          }

        } else if (secondaries[sec].GetGVcode() == geantphysics::Positron::Definition()->GetInternalCode()) {

          data.posEn.Fill(enNormed);

          double posZ = 0.5 * (1.0 - secondaries[sec].GetDirZ());
          if (posZ > 0 && log10(posZ) > -12.0) {
            data.posAngle.Fill(log10(posZ));
          }
        }
      }
    }
  }
}

int main()
{

  PrepareWorld();
  bHModelAlias    = PrepareBHModel(true);
  bHModelRej      = PrepareBHModel(false);
  bHVecModelAlias = PrepareVecBHModel(true);
  bHVecModelRej   = PrepareVecBHModel(false);
  td              = PrepareTaskData();

  Printf("Number of gamma for each test %d", kMaxBasket * kBasketTries);
  Printf("Relative histograms of kinematics (difference in percents)");
  {
    Printf("Test for alias method");

    BHValidData scalar;
    FillDataScalar(scalar, true);
    BHValidData vector;
    FillDataVector(vector, true);

    Printf("====EL EN====");
    vector.elEn.Compare(scalar.elEn);
    Printf("====POS EN====");
    vector.posEn.Compare(scalar.posEn);
    Printf("====EL Dir====");
    vector.elAngle.Compare(scalar.elAngle);
    Printf("====POS Dir====");
    vector.posAngle.Compare(scalar.posAngle);
  }
  {
    Printf("Test for rej method");

    BHValidData scalar;
    FillDataScalar(scalar, false);
    BHValidData vector;
    FillDataVector(vector, false);

    Printf("====EL EN====");
    vector.elEn.Compare(scalar.elEn);
    Printf("====POS EN====");
    vector.posEn.Compare(scalar.posEn);
    Printf("====EL Dir====");
    vector.elAngle.Compare(scalar.elAngle);
    Printf("====POS Dir====");
    vector.posAngle.Compare(scalar.posAngle);
  }

  return 0;
}
