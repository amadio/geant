#include <Geant/Electron.h>
#include <Geant/Positron.h>
#include "BetheHeitlerTestCommon.h"
#include "Hist.h"

const int kBasketTries = 100000;

const int kNumBins = 100;
struct BHValidData {
  Hist primEn;
  Hist primAngle;
  Hist secEn;
  Hist secAngle;
  BHValidData()
      : primEn(0.0, 1.0, kNumBins), primAngle(-12.0, 0.0, kNumBins), secEn(0.0, 1.0, kNumBins),
        secAngle(-12.0, 0.0, kNumBins)
  {
  }
};

BetheHeitlerPairModel *rpModelRej;
VecBetheHeitlerPairModel *rpModelRejVec;
BetheHeitlerPairModel *rpModelAlias;
VecBetheHeitlerPairModel *rpModelAliasVec;
TaskData *td;

void FillDataVector(BHValidData &data, bool useAlias)
{
  LightTrack_v primaries;
  BetheHeitlerPairModel *model = useAlias ? rpModelAliasVec : rpModelRejVec;

  for (int bask = 0; bask < kBasketTries; ++bask) {
    PreparePrimaries(primaries, kMaxBasket);

    std::vector<double> enBeforeInteraction;
    enBeforeInteraction.insert(enBeforeInteraction.begin(), primaries.GetKinEArr(),
                               primaries.GetKinEArr() + kMaxBasket);

    primaries.SetNtracks(kMaxBasket);

    td->fPhysicsData->GetSecondarySOA().ClearTracks();
    model->SampleSecondariesVector(primaries, td);

    auto &secondaries = td->fPhysicsData->GetSecondarySOA();
    for (int i = 0; i < secondaries.GetNtracks(); ++i) {
      int primIdx     = secondaries.GetTrackIndex(i);
      double enNormed = (secondaries.GetKinE(i) + secondaries.GetMass(i)) / enBeforeInteraction[primIdx];

      if (secondaries.GetGVcode(i) == geantphysics::Electron::Definition()->GetInternalCode()) {

        data.primEn.Fill(enNormed);
        double elZ = 0.5 * (1.0 - secondaries.GetDirZ(i));
        if (elZ > 0 && log10(elZ) > -12.0) {
          data.primAngle.Fill(log10(elZ));
        }

      } else if (secondaries.GetGVcode(i) == geantphysics::Positron::Definition()->GetInternalCode()) {

        data.secEn.Fill(enNormed);

        double posZ = 0.5 * (1.0 - secondaries.GetDirZ(i));
        if (posZ > 0 && log10(posZ) > -12.0) {
          data.secAngle.Fill(log10(posZ));
        }
      }
    }
  }
}

void FillDataScalar(BHValidData &data, bool useAlias)
{
  std::vector<LightTrack> primaries;
  BetheHeitlerPairModel *model = useAlias ? rpModelAlias : rpModelRej;

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

          data.primEn.Fill(enNormed);
          double elZ = 0.5 * (1.0 - secondaries[sec].GetDirZ());
          if (elZ > 0 && log10(elZ) > -12.0) {
            data.primAngle.Fill(log10(elZ));
          }

        } else if (secondaries[sec].GetGVcode() == geantphysics::Positron::Definition()->GetInternalCode()) {

          data.secEn.Fill(enNormed);

          double posZ = 0.5 * (1.0 - secondaries[sec].GetDirZ());
          if (posZ > 0 && log10(posZ) > -12.0) {
            data.secAngle.Fill(log10(posZ));
          }
        }
      }
    }
  }
}

int main()
{

  PrepareWorld();
  rpModelAlias    = PrepareBHModel(true);
  rpModelRej      = PrepareBHModel(false);
  rpModelAliasVec = PrepareVecBHModel(true);
  rpModelRejVec   = PrepareVecBHModel(false);
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
    vector.primEn.Compare(scalar.primEn);
    Printf("====POS EN====");
    vector.secEn.Compare(scalar.secEn);
    Printf("====EL Dir====");
    vector.primAngle.Compare(scalar.primAngle);
    Printf("====POS Dir====");
    vector.secAngle.Compare(scalar.secAngle);
  }
  {
    Printf("Test for rej method");

    BHValidData scalar;
    FillDataScalar(scalar, false);
    BHValidData vector;
    FillDataVector(vector, false);

    Printf("====EL EN====");
    vector.primEn.Compare(scalar.primEn);
    Printf("====POS EN====");
    vector.secEn.Compare(scalar.secEn);
    Printf("====EL Dir====");
    vector.primAngle.Compare(scalar.primAngle);
    Printf("====POS Dir====");
    vector.secAngle.Compare(scalar.secAngle);
  }

  return 0;
}
