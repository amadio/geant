
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

  for (int i = 0; i < N; i += kPhysDVWidth) {

    PhysDV electronCut;
    PhysDI mcIndxLocal; // Used by alias method
    for (int l = 0; l < kPhysDVWidth; ++l) {
      const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(tracks.GetMaterialCutCoupleIndex(i + l));
      electronCut[l]             = matCut->GetProductionCutsInEnergy()[1];
      if (GetUseSamplingTables()) {
        mcIndxLocal[l] = fGlobalMatECutIndxToLocal[matCut->GetIndex()];
      }
    }
    if (GetUseSamplingTables()) {
      PhysDV r1       = td->fRndm->uniformV();
      PhysDV r2       = td->fRndm->uniformV();
      PhysDV r3       = td->fRndm->uniformV();
      PhysDV primEkin = tracks.GetKinEVec(i);
      PhysDV eps      = SampleEnergyTransfer(electronCut, mcIndxLocal, fAliasData.fLogEmin.data(),
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

  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV ekin = tracks.GetKinEVec(i);

    PhysDV elInitTotalEnergy   = ekin + geant::units::kElectronMassC2; // initial total energy of the e-/e+
    PhysDV elInitTotalMomentum = Sqrt(ekin * (elInitTotalEnergy + geant::units::kElectronMassC2));

    PhysDV deltaKinEnergy;
    vecCore::Load(deltaKinEnergy, epsArr + i);
    PhysDV deltaTotalMomentum = Sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0 * geant::units::kElectronMassC2));
    PhysDV cost               = deltaKinEnergy * (elInitTotalEnergy + geant::units::kElectronMassC2) /
                  (deltaTotalMomentum * elInitTotalMomentum);
    PhysDV cosTheta = Min(cost, (PhysDV)1.0);
    PhysDV sinTheta = Sqrt((1.0 - cosTheta) * (1.0 + cosTheta));

    PhysDV phi = geant::units::kTwoPi * td->fRndm->uniformV();
    PhysDV sinPhi, cosPhi;
    SinCos(phi, &sinPhi, &cosPhi);

    PhysDV deltaDirX = sinTheta * cosPhi;
    PhysDV deltaDirY = sinTheta * sinPhi;
    PhysDV deltaDirZ = cosTheta;
    RotateToLabFrame(deltaDirX, deltaDirY, deltaDirZ, tracks.GetDirXVec(i), tracks.GetDirYVec(i), tracks.GetDirZVec(i));

    LightTrack_v &secondaries = td->fPhysicsData->GetSecondarySOA();
    for (int l = 0; l < kPhysDVWidth; ++l) {
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
    PhysDV elDirX = elInitTotalMomentum * tracks.GetDirXVec(i) - deltaTotalMomentum * deltaDirX;
    PhysDV elDirY = elInitTotalMomentum * tracks.GetDirYVec(i) - deltaTotalMomentum * deltaDirY;
    PhysDV elDirZ = elInitTotalMomentum * tracks.GetDirZVec(i) - deltaTotalMomentum * deltaDirZ;
    // normalisation
    PhysDV norm = 1.0 / Sqrt(elDirX * elDirX + elDirY * elDirY + elDirZ * elDirZ);

    // update primary track direction
    for (int l = 0; l < kPhysDVWidth; ++l) {
      tracks.SetDirX((elDirX * norm)[l], i + l);
      tracks.SetDirY((elDirY * norm)[l], i + l);
      tracks.SetDirZ((elDirZ * norm)[l], i + l);
      tracks.SetKinE((ekin - deltaKinEnergy)[l], i + l);
    }
  }
}

PhysDV VecMollerBhabhaIonizationModel::SampleEnergyTransfer(PhysDV elProdCut, PhysDI mcLocalIdx, double *tableEmin,
                                                            double *tableILDeta, PhysDV primekin, PhysDV r1, PhysDV r2,
                                                            PhysDV r3)
{
  PhysDV lPrimEkin    = Log(primekin);
  PhysDV logEmin      = vecCore::Gather<PhysDV>(tableEmin, mcLocalIdx);
  PhysDV ilDelta      = vecCore::Gather<PhysDV>(tableILDeta, mcLocalIdx);
  PhysDV val          = (lPrimEkin - logEmin) * ilDelta;
  PhysDI indxPrimEkin = (PhysDI)val; // lower electron energy bin index
  PhysDV pIndxHigh    = val - indxPrimEkin;
  PhysDM mask         = r1 < pIndxHigh;
  if (!mask.isEmpty()) {
    vecCore::MaskedAssign(indxPrimEkin, mask, indxPrimEkin + 1);
  }

  PhysDV xiV;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    //    LinAliasCached& als = fAliasData.fTablesPerMatCut[mcLocalIdx[l]]->fAliasData[indxPrimEkin[l]];
    //    double xi = AliasTableAlternative::SampleLinear(als,fSTNumSamplingElecEnergies,r2[l],r3[l]);

    const LinAlias *als = fSamplingTables[mcLocalIdx[l]]->fAliasData[indxPrimEkin[l]];
    const double xi     = fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]),
                                                  &(als->fAliasIndx[0]), fSTNumSamplingElecEnergies, r2[l], r3[l]);
    xiV[l] = xi;
  }

  // sample the transformed variable
  PhysDV dum1 = Log(primekin / elProdCut);
  if (fIsElectron) {
    dum1 -= 0.693147180559945; // dum1 = dum1 + log(0.5)
  }
  return Exp(xiV * dum1) * elProdCut;
}

void VecMollerBhabhaIonizationModel::SampleEnergyTransfer(const double *elProdCut, const double *primeKinArr,
                                                          double *epsOut, int N, const geant::TaskData *td)
{

  // assert(N>=kPhysDVWidth)
  int currN        = 0;
  PhysDM lanesDone = PhysDM::Zero();
  PhysDI idx;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    idx[l] = currN++;
  }

  while (currN < N || !lanesDone.isFull()) {
    PhysDV primekin = vecCore::Gather<PhysDV>(primeKinArr, idx);
    PhysDV tmin     = vecCore::Gather<PhysDV>(elProdCut, idx);
    PhysDV tmax     = (fIsElectron) ? (0.5 * primekin) : (primekin);
    PhysDV xmin     = tmin / primekin;
    PhysDV xmax     = tmax / primekin;
    PhysDV gamma    = primekin / geant::units::kElectronMassC2 + 1.0;
    PhysDV gamma2   = gamma * gamma;
    PhysDV beta2    = 1. - 1. / gamma2;
    PhysDV xminmax  = xmin * xmax;

    PhysDM accepted;
    PhysDV deltaEkin;
    if (fIsElectron) {
      PhysDV gg   = (2.0 * gamma - 1.0) / gamma2;
      PhysDV y    = 1. - xmax;
      PhysDV gf   = 1.0 - gg * xmax + xmax * xmax * (1.0 - gg + (1.0 - gg * y) / (y * y));
      PhysDV rnd0 = td->fRndm->uniformV();
      PhysDV rnd1 = td->fRndm->uniformV();
      deltaEkin   = xminmax / (xmin * (1.0 - rnd0) + xmax * rnd0);
      PhysDV xx   = 1.0 - deltaEkin;
      PhysDV dum  = 1.0 - gg * deltaEkin + deltaEkin * deltaEkin * (1.0 - gg + (1.0 - gg * xx) / (xx * xx));

      accepted = gf * rnd1 < dum;
    } else {
      PhysDV y     = 1.0 / (1.0 + gamma);
      PhysDV y2    = y * y;
      PhysDV y12   = 1.0 - 2.0 * y;
      PhysDV b1    = 2.0 - y2;
      PhysDV b2    = y12 * (3.0 + y2);
      PhysDV y122  = y12 * y12;
      PhysDV b4    = y122 * y12;
      PhysDV b3    = b4 + y122;
      PhysDV xmax2 = xmax * xmax;
      PhysDV gf    = 1.0 + (xmax2 * b4 - xmin * xmin * xmin * b3 + xmax2 * b2 - xmin * b1) * beta2;
      PhysDV rnd0  = td->fRndm->uniformV();
      PhysDV rnd1  = td->fRndm->uniformV();
      deltaEkin    = xminmax / (xmin * (1.0 - rnd0) + xmax * rnd0);
      PhysDV xx    = deltaEkin * deltaEkin;
      PhysDV dum   = 1.0 + (xx * xx * b4 - deltaEkin * xx * b3 + xx * b2 - deltaEkin * b1) * beta2;
      accepted     = gf * rnd1 < dum;
    }
    deltaEkin *= primekin;

    if (accepted.isNotEmpty()) {
      vecCore::Scatter(deltaEkin, epsOut, idx);
    }

    lanesDone = lanesDone || accepted;
    for (int l = 0; l < kPhysDVWidth; ++l) {
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

} // namespace geantphysics
