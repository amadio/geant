
#include "Geant/VecMollerBhabhaIonizationModel.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/Material.h"
#include "Geant/Element.h"
#include "Geant/MaterialProperties.h"

#include "Geant/MaterialCuts.h"
#include "Geant/AliasTable.h"

#include "Geant/Electron.h"

#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"

// from geantV
#include "Geant/TaskData.h"

#include <iostream>
#include <cmath>
#include <Geant/VecMollerBhabhaIonizationModel.h>
#include "Geant/AliasTable.h"

namespace geantphysics {

using geant::Double_v;
using geant::IndexD_v;
using geant::kVecLenD;
using geant::MaskD_v;

VecMollerBhabhaIonizationModel::VecMollerBhabhaIonizationModel(bool iselectron, const std::string &modelname)
    : MollerBhabhaIonizationModel(iselectron, modelname)
{
}

void VecMollerBhabhaIonizationModel::Initialize()
{
  MollerBhabhaIonizationModel::Initialize();

  int NCutTables = (int)fSamplingTables.size();
  for (int i = 0; i < NCutTables; ++i) {
    // push empty data for equal spacing
    fAliasData.fILDelta.emplace_back();
    fAliasData.fLogEmin.emplace_back();
    fAliasData.fNData.emplace_back();
    fAliasData.fTablesPerMatCut.push_back(std::unique_ptr<AliasDataForMatCut>(nullptr));
    if (fSamplingTables[i] == nullptr) continue;

    AliasDataMaterialCuts *oldTable = fSamplingTables[i];
    fAliasData.fILDelta[i]          = oldTable->fILDelta;
    fAliasData.fLogEmin[i]          = oldTable->fLogEmin;
    fAliasData.fNData[i]            = oldTable->fNData;
    fAliasData.fTablesPerMatCut[i]  = std::unique_ptr<AliasDataForMatCut>(
        new AliasDataForMatCut((int)oldTable->fAliasData.size(), fAliasData.fLogEmin[i], fAliasData.fILDelta[i]));
    auto &newTable = *fAliasData.fTablesPerMatCut[i];

    for (int j = 0; j < (int)oldTable->fAliasData.size(); ++j) {
      auto &oldETable = oldTable->fAliasData[j];
      LinAliasCached newETable((int)oldETable->fAliasIndx.size());

      for (int k = 0; k < (int)oldETable->fAliasIndx.size() - 1; ++k) {
        auto &transAlias       = newETable.fLinAliasData[k];
        transAlias.fAliasIndx  = oldETable->fAliasIndx[k];
        transAlias.fAliasW     = oldETable->fAliasW[k];
        transAlias.fX          = oldETable->fXdata[k];
        transAlias.fYdata      = oldETable->fYdata[k];
        transAlias.fXdelta     = oldETable->fXdata[k + 1] - oldETable->fXdata[k];
        transAlias.fYdataDelta = (oldETable->fYdata[k + 1] - oldETable->fYdata[k]) / oldETable->fYdata[k];
        transAlias.fXdivYdelta = transAlias.fXdelta / transAlias.fYdataDelta;
      }

      newTable.fAliasData.push_back(newETable);
    }
  }
}

void VecMollerBhabhaIonizationModel::SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td)
{
  int N            = tracks.GetNtracks();
  double *epsArr   = td->fPhysicsData->fPhysicsScratchpad.fEps;
  double *elCutArr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr;

  for (int i = 0; i < N; i += kVecLenD) {

    Double_v electronCut;
    IndexD_v mcIndxLocal; // Used by alias method
    for (int l = 0; l < kVecLenD; ++l) {
      const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(tracks.GetMaterialCutCoupleIndex(i + l));
      electronCut[l]             = matCut->GetProductionCutsInEnergy()[1];
      if (GetUseSamplingTables()) {
        mcIndxLocal[l] = fGlobalMatECutIndxToLocal[matCut->GetIndex()];
      }
    }
    Double_v primEkin = tracks.GetKinEVec(i);
    assert((primEkin >= electronCut).isFull()); // Cut filtering should be applied up the call chain.
    if (GetUseSamplingTables()) {
      Double_v r1  = td->fRndm->uniformV();
      Double_v r2  = td->fRndm->uniformV();
      Double_v r3  = td->fRndm->uniformV();
      Double_v eps = SampleEnergyTransfer(electronCut, mcIndxLocal, fAliasData.fLogEmin.data(),
                                          fAliasData.fILDelta.data(), primEkin, r1, r2, r3);
      vecCore::Store(eps, epsArr + i);
    } else {
      vecCore::Store(electronCut, elCutArr + i);
    }
  }

  if (!GetUseSamplingTables()) {
    elCutArr[N]            = elCutArr[N - 1];
    tracks.GetKinEArr()[N] = tracks.GetKinEArr()[N - 1];
    SampleEnergyTransfer(elCutArr, tracks.GetKinEArr(), epsArr, N, td);
  }

  for (int i = 0; i < N; i += kVecLenD) {
    Double_v ekin = tracks.GetKinEVec(i);

    Double_v elInitTotalEnergy   = ekin + geant::units::kElectronMassC2; // initial total energy of the e-/e+
    Double_v elInitTotalMomentum = Math::Sqrt(ekin * (elInitTotalEnergy + geant::units::kElectronMassC2));

    Double_v deltaKinEnergy;
    vecCore::Load(deltaKinEnergy, epsArr + i);
    Double_v deltaTotalMomentum = Math::Sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0 * geant::units::kElectronMassC2));
    Double_v cost               = deltaKinEnergy * (elInitTotalEnergy + geant::units::kElectronMassC2) /
                    (deltaTotalMomentum * elInitTotalMomentum);
    Double_v cosTheta = Math::Min(cost, (Double_v)1.0);
    Double_v sinTheta = Math::Sqrt((1.0 - cosTheta) * (1.0 + cosTheta));

    Double_v phi = geant::units::kTwoPi * td->fRndm->uniformV();
    Double_v sinPhi, cosPhi;
    Math::SinCos(phi, sinPhi, cosPhi);

    Double_v deltaDirX = sinTheta * cosPhi;
    Double_v deltaDirY = sinTheta * sinPhi;
    Double_v deltaDirZ = cosTheta;
    RotateToLabFrame(deltaDirX, deltaDirY, deltaDirZ, tracks.GetDirXVec(i), tracks.GetDirYVec(i), tracks.GetDirZVec(i));

    LightTrack_v &secondaries = td->fPhysicsData->GetSecondarySOA();
    for (int l = 0; l < kVecLenD; ++l) {
      int idx = secondaries.InsertTrack();
      secondaries.SetKinE(deltaKinEnergy[l], idx);
      secondaries.SetDirX(deltaDirX[l], idx);
      secondaries.SetDirY(deltaDirY[l], idx);
      secondaries.SetDirZ(deltaDirZ[l], idx);
      secondaries.SetGVcode(fSecondaryInternalCode, idx);
      secondaries.SetMass(geant::units::kElectronMassC2, idx);
      secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx);
    }

    // compute the primary e-/e+ post interaction kinetic energy and direction: from momentum vector conservation
    // final momentum of the primary e-/e+ in the lab frame
    Double_v elDirX = elInitTotalMomentum * tracks.GetDirXVec(i) - deltaTotalMomentum * deltaDirX;
    Double_v elDirY = elInitTotalMomentum * tracks.GetDirYVec(i) - deltaTotalMomentum * deltaDirY;
    Double_v elDirZ = elInitTotalMomentum * tracks.GetDirZVec(i) - deltaTotalMomentum * deltaDirZ;
    // normalisation
    Double_v norm = 1.0 / Math::Sqrt(elDirX * elDirX + elDirY * elDirY + elDirZ * elDirZ);

    // update primary track direction
    for (int l = 0; l < kVecLenD; ++l) {
      tracks.SetDirX((elDirX * norm)[l], i + l);
      tracks.SetDirY((elDirY * norm)[l], i + l);
      tracks.SetDirZ((elDirZ * norm)[l], i + l);
      tracks.SetKinE((ekin - deltaKinEnergy)[l], i + l);
    }
  }
}

Double_v VecMollerBhabhaIonizationModel::SampleEnergyTransfer(Double_v elProdCut, IndexD_v mcLocalIdx,
                                                              double *tableEmin, double *tableILDeta, Double_v primekin,
                                                              Double_v r1, Double_v r2, Double_v r3)
{
  Double_v lPrimEkin    = Math::Log(primekin);
  Double_v logEmin      = vecCore::Gather<Double_v>(tableEmin, mcLocalIdx);
  Double_v ilDelta      = vecCore::Gather<Double_v>(tableILDeta, mcLocalIdx);
  Double_v val          = (lPrimEkin - logEmin) * ilDelta;
  IndexD_v indxPrimEkin = (IndexD_v)val; // lower electron energy bin index
  Double_v pIndxHigh    = val - indxPrimEkin;
  MaskD_v mask          = r1 < pIndxHigh;
  if (!mask.isEmpty()) {
    vecCore::MaskedAssign(indxPrimEkin, mask, indxPrimEkin + 1);
  }

  Double_v xiV;
  for (int l = 0; l < kVecLenD; ++l) {
    //    LinAliasCached& als = fAliasData.fTablesPerMatCut[mcLocalIdx[l]]->fAliasData[indxPrimEkin[l]];
    //    double xi = AliasTableAlternative::SampleLinear(als,fSTNumSamplingElecEnergies,r2[l],r3[l]);

    const LinAlias *als = fSamplingTables[mcLocalIdx[l]]->fAliasData[indxPrimEkin[l]];
    const double xi     = fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]),
                                                  &(als->fAliasIndx[0]), fSTNumSamplingElecEnergies, r2[l], r3[l]);
    xiV[l] = xi;
  }

  // sample the transformed variable
  Double_v dum1 = Math::Log(primekin / elProdCut);
  if (fIsElectron) {
    dum1 -= 0.693147180559945; // dum1 = dum1 + log(0.5)
  }
  return Math::Exp(xiV * dum1) * elProdCut;
}

void VecMollerBhabhaIonizationModel::SampleEnergyTransfer(const double *elProdCut, const double *primeKinArr,
                                                          double *epsOut, int N, const geant::TaskData *td)
{

  // assert(N>=kVecLenD)
  int currN         = 0;
  MaskD_v lanesDone = MaskD_v::Zero();
  IndexD_v idx;
  for (int l = 0; l < kVecLenD; ++l) {
    idx[l] = currN++;
  }

  while (currN < N || !lanesDone.isFull()) {
    Double_v primekin = vecCore::Gather<Double_v>(primeKinArr, idx);
    Double_v tmin     = vecCore::Gather<Double_v>(elProdCut, idx);
    Double_v tmax     = (fIsElectron) ? (0.5 * primekin) : (primekin);
    Double_v xmin     = tmin / primekin;
    Double_v xmax     = tmax / primekin;
    Double_v gamma    = primekin / geant::units::kElectronMassC2 + 1.0;
    Double_v gamma2   = gamma * gamma;
    Double_v beta2    = 1. - 1. / gamma2;
    Double_v xminmax  = xmin * xmax;

    MaskD_v accepted;
    Double_v deltaEkin;
    if (fIsElectron) {
      Double_v gg   = (2.0 * gamma - 1.0) / gamma2;
      Double_v y    = 1. - xmax;
      Double_v gf   = 1.0 - gg * xmax + xmax * xmax * (1.0 - gg + (1.0 - gg * y) / (y * y));
      Double_v rnd0 = td->fRndm->uniformV();
      Double_v rnd1 = td->fRndm->uniformV();
      deltaEkin     = xminmax / (xmin * (1.0 - rnd0) + xmax * rnd0);
      Double_v xx   = 1.0 - deltaEkin;
      Double_v dum  = 1.0 - gg * deltaEkin + deltaEkin * deltaEkin * (1.0 - gg + (1.0 - gg * xx) / (xx * xx));

      accepted = gf * rnd1 < dum;
    } else {
      Double_v y     = 1.0 / (1.0 + gamma);
      Double_v y2    = y * y;
      Double_v y12   = 1.0 - 2.0 * y;
      Double_v b1    = 2.0 - y2;
      Double_v b2    = y12 * (3.0 + y2);
      Double_v y122  = y12 * y12;
      Double_v b4    = y122 * y12;
      Double_v b3    = b4 + y122;
      Double_v xmax2 = xmax * xmax;
      Double_v gf    = 1.0 + (xmax2 * b4 - xmin * xmin * xmin * b3 + xmax2 * b2 - xmin * b1) * beta2;
      Double_v rnd0  = td->fRndm->uniformV();
      Double_v rnd1  = td->fRndm->uniformV();
      deltaEkin      = xminmax / (xmin * (1.0 - rnd0) + xmax * rnd0);
      Double_v xx    = deltaEkin * deltaEkin;
      Double_v dum   = 1.0 + (xx * xx * b4 - deltaEkin * xx * b3 + xx * b2 - deltaEkin * b1) * beta2;
      accepted       = gf * rnd1 < dum;
    }
    deltaEkin *= primekin;

    if (accepted.isNotEmpty()) {
      vecCore::Scatter(deltaEkin, epsOut, idx);
    }

    lanesDone = lanesDone || accepted;
    for (int l = 0; l < kVecLenD; ++l) {
      auto laneDone = accepted[l];
      if (laneDone) {
        if (currN < N) {
          idx[l]       = currN++;
          lanesDone[l] = false;
        } else {
          idx[l] = N;
        }
      }
    }
  }
}

bool VecMollerBhabhaIonizationModel::IsModelUsable(const MaterialCuts *matCut, double ekin)
{
  const double electronCut = matCut->GetProductionCutsInEnergy()[1];
  double maxETransfer      = ekin;
  if (fIsElectron) {
    maxETransfer *= 0.5;
  }
  return ekin < GetHighEnergyUsageLimit() && ekin > GetLowEnergyUsageLimit() && maxETransfer > electronCut;
}

} // namespace geantphysics
