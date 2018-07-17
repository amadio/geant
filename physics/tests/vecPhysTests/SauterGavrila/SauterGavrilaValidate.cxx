/*
 * Basic sanity check for SauterGavrila physics.
 * Checks secondary energy and dirZ
 */
#include "SauterGavrilaTestCommon.h"
#include "Hist.h"

const int kBasketTries = 100000;

const int kNumBins = 100;
struct PhotoElectricValidData {
  Hist peE;
  Hist peTheta;
  PhotoElectricValidData()
      : peE(0.0, 1.0, kNumBins),
        peTheta(-1.0, 1.0, kNumBins)
  {
  }
};

SauterGavrilaPhotoElectricModel *sg;
SauterGavrilaPhotoElectricModel *vsg;
//TaskData *td;

void FillDataVector(PhotoElectricValidData &data,  geant::TaskData *td)
{

  LightTrack_v primaries;
  //VecSauterGavrilaPhotoElectricModel *sauterG = PrepareVecSauterGavrilaModel(useAlias);

  for (int bask = 0; bask < kBasketTries; ++bask) {
    PreparePrimaries(primaries, kMaxBasket); //to add the zed

    std::vector<double> enBeforeInteraction;
    enBeforeInteraction.insert(enBeforeInteraction.begin(), primaries.GetKinEArr(),
                               primaries.GetKinEArr() + kMaxBasket);

    primaries.SetNtracks(kMaxBasket);

    td->fPhysicsData->GetSecondarySOA().ClearTracks();
    vsg->SampleSecondaries(primaries, td);
    //sauterG->SampleSecondariesVector(primaries, td);
      
    auto &secondaries = td->fPhysicsData->GetSecondarySOA();
    for (int i = 0; i < secondaries.GetNtracks(); ++i) {
      double enNormed = secondaries.GetKinE(i) / enBeforeInteraction[i];
//        std::cout<<"Vector enNormed: "<<enNormed<<std::endl;
//        std::cout<<"secondaries.GetKinE("<<i<<"<<): "<<secondaries.GetKinE(i)<<std::endl;
      data.peE.Fill(enNormed);
      double emCost = secondaries.GetDirZ(i);
      data.peTheta.Fill(emCost);
      //std::cout<<"Vector emCost: "<<emCost<<std::endl;
    }
  }
}

void FillDataScalar(PhotoElectricValidData &data,  geant::TaskData *td)
{

  std::vector<LightTrack> primaries;
  for (int bask = 0; bask < kBasketTries; ++bask) {
    PreparePrimaries(primaries, kMaxBasket);

    std::vector<double> enBeforeInteraction;
    for (auto &gamma : primaries) {
      enBeforeInteraction.push_back(gamma.GetKinE());
    }

    for (int i = 0; i < kMaxBasket; ++i) {
      td->fPhysicsData->ClearSecondaries();

      //sauterG->SampleSecondaries(primaries[i], Td);
      sg->SampleSecondaries(primaries[i], td);
      LightTrack *secondaries = td->fPhysicsData->GetListOfSecondaries();
      int numSec              = td->fPhysicsData->GetNumOfSecondaries();
      for (int sec = 0; sec < numSec; ++sec) {
        double enNormed = secondaries[sec].GetKinE() / enBeforeInteraction[i];
        data.peE.Fill(enNormed);
        double emCost = secondaries[sec].GetDirZ();
        data.peTheta.Fill(emCost);
      }
    }
  }
}

int main()
{
  SetUpSimulation();
  sg  = PrepareSauterGavrilaModel(true);
  vsg = PrepareVecSauterGavrilaModel(true);
  auto td = PrepareTaskData();

  Printf("Number of gamma for each test %d", kMaxBasket * kBasketTries);
  Printf("Relative histograms of kinematics (difference in percents)");
  {
    Printf("*** Test for alias method *** ");

    PhotoElectricValidData scalar;
    FillDataScalar(scalar, td.get());
    PhotoElectricValidData vector;
    FillDataVector(vector, td.get());
    Printf("Normalized photoelectron energy");
    vector.peE.Compare(scalar.peE);

    Printf("PhotoElectron z direction");
    vector.peTheta.Compare(scalar.peTheta);
  }
//  {
//    Printf("*** Test for rej method ***");
//
//    PhotoElectricValidData scalar;
//    FillDataScalar(scalar, false);
//    PhotoElectricValidData vector;
//    FillDataVector(vector, false);
//
//    Printf("Normalized photoelectron energy");
//    vector.peE.Compare(scalar.peE);
//    Printf("Photoelectron z direction");
//    vector.peTheta.Compare(scalar.peTheta);
//  }

  return 0;
}
