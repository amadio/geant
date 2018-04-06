#include <Geant/SystemOfUnits.h>
#include <Geant/PhysicalConstants.h>
#include <Geant/TaskData.h>
#include <Geant/PhysicsData.h>
#include <Geant/MaterialCuts.h>
#include <TError.h>
#include "Geant/VecBetheHeitlerPairModel.h"
#include "Geant/AliasTable.h"

namespace geantphysics {

void VecBetheHeitlerPairModel::Initialize()
{
  BetheHeitlerPairModel::Initialize();
  fAliasTablesPerZ.resize(fSamplingTables.size());

  // Steal sampling tables to new format
  for (size_t i = 0; i < fSamplingTables.size(); ++i) { // Loop over Z
    if (fSamplingTables[i] == nullptr) continue;

    auto &defaultGammaETables    = fSamplingTables[i]->fRatinAliasData;
    auto &transposedGammaETables = fAliasTablesPerZ[i].fTablePerEn;
    for (size_t j = 0; j < defaultGammaETables.size(); ++j) { // Loop over incoming gamma Es in element table
      transposedGammaETables.emplace_back(fSTNumDiscreteEnergyTransferVals);
      for (int k = 0; k < fSTNumDiscreteEnergyTransferVals;
           ++k) { // Loop over transferred E in corresponding gamma E tbl.
        auto &aliasTableData       = transposedGammaETables[j].fData[k];
        aliasTableData.fAliasIndx  = defaultGammaETables[j]->fAliasIndx[k];
        aliasTableData.fAliasW     = defaultGammaETables[j]->fAliasW[k];
        aliasTableData.fCumulative = defaultGammaETables[j]->fCumulative[k];
        aliasTableData.fParaA      = defaultGammaETables[j]->fParaA[k];
        aliasTableData.fParaB      = defaultGammaETables[j]->fParaB[k];
        aliasTableData.fXdata      = defaultGammaETables[j]->fXdata[k];
      }
    }
  }

  // fIsUseTsaisScreening = false;
}

void VecBetheHeitlerPairModel::SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td)
{
  // EMModel::SampleSecondariesVector(tracks,td);
  // return;
  int N     = tracks.GetNtracks();
  int *IZet = td->fPhysicsData->fPhysicsScratchpad.fIzet;

  short *sortKey = tracks.GetSortKeyV();
  for (int i = 0; i < N; ++i) {
    double ekin = tracks.GetKinE(i);
    sortKey[i]  = ekin < fGammaEneregyLimit ? (short)0 : (short)1;
  }
  int smallEtracksInd = tracks.SortTracks();
  for (int k = 0; k < smallEtracksInd; ++k) {
    double eps0                                  = geant::units::kElectronMassC2 / tracks.GetKinE(k);
    td->fPhysicsData->fPhysicsScratchpad.fEps[k] = eps0 + (0.5 - eps0) * td->fRndm->uniform();
  }

  // Chose Z
  for (int i = 0; i < N; ++i) {
    MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[tracks.GetMaterialCutCoupleIndex(i)];
    const Vector_t<Element *> &theElements = matCut->GetMaterial()->GetElementVector();
    double targetElemIndx                  = 0;
    if (theElements.size() > 1) {
      targetElemIndx = SampleTargetElementIndex(matCut, tracks.GetKinE(i), td->fRndm->uniform());
    }
    const double zet = theElements[targetElemIndx]->GetZ();
    const int izet   = std::min(std::lrint(zet), gMaxZet - 1);
    IZet[i]          = izet;
  }

  if (GetUseSamplingTables()) {
    {
      // First iteration of loop below manually extracted since in first SIMD pack we could have mixed
      // particles that were sampled with simplified procedure
      int i = (smallEtracksInd / kPhysDVWidth) * kPhysDVWidth;
      if (i < N) {
        PhysDV ekin    = tracks.GetKinEVec(i);
        PhysDM smallE  = ekin < fGammaEneregyLimit;
        PhysDV r1      = td->fRndm->uniformV();
        PhysDV r2      = td->fRndm->uniformV();
        PhysDV r3      = td->fRndm->uniformV();
        PhysDV tmpEkin = ekin;
        // To be sure that energy that is too small for alias table will not slip in alias sampling
        vecCore::MaskedAssign(tmpEkin, smallE, (PhysDV)fGammaEneregyLimit * (1.1));

        PhysDV sampledEps = SampleTotalEnergyTransferAliasOneShot(tmpEkin, &IZet[i], r1, r2, r3);
        PhysDV eps;
        vecCore::Load(eps, &td->fPhysicsData->fPhysicsScratchpad.fEps[i]);
        vecCore::MaskedAssign(eps, !smallE, sampledEps);
        vecCore::Store(eps, &td->fPhysicsData->fPhysicsScratchpad.fEps[i]);
      }
    }
    for (int i = (smallEtracksInd / kPhysDVWidth) * kPhysDVWidth + kPhysDVWidth; i < N; i += kPhysDVWidth) {
      PhysDV ekin;
      vecCore::Load(ekin, tracks.GetKinEArr(i));
      PhysDV r1 = td->fRndm->uniformV();
      PhysDV r2 = td->fRndm->uniformV();
      PhysDV r3 = td->fRndm->uniformV();

      PhysDV sampledEps = SampleTotalEnergyTransferAliasOneShot(ekin, &IZet[i], r1, r2, r3);
      vecCore::Store(sampledEps, &td->fPhysicsData->fPhysicsScratchpad.fEps[i]);
    }
  } else {

    int i = (smallEtracksInd / kPhysDVWidth) * kPhysDVWidth;
    if (i < N) {
      // Always create fake particle int the end of real particle array for rejection sampling
      tracks.GetKinEArr()[N] = tracks.GetKinEArr()[N - 1];
      IZet[N]                = IZet[N - 1];
      SampleTotalEnergyTransferRejVec(tracks.GetKinEArr(i), &IZet[i], &td->fPhysicsData->fPhysicsScratchpad.fEps[i],
                                      N - i, td);
    }
  }

  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV ekin = tracks.GetKinEVec(i);

    PhysDV eps;
    vecCore::Load(eps, &td->fPhysicsData->fPhysicsScratchpad.fEps[i]);

    PhysDV electronTotE, positronTotE;
    PhysDM tmpM = td->fRndm->uniformV() > 0.5;
    vecCore::MaskedAssign(electronTotE, tmpM, (1. - eps) * ekin);
    vecCore::MaskedAssign(positronTotE, tmpM, eps * ekin);

    vecCore::MaskedAssign(positronTotE, !tmpM, (1. - eps) * ekin);
    vecCore::MaskedAssign(electronTotE, !tmpM, eps * ekin);

    PhysDV r0 = td->fRndm->uniformV();
    PhysDV r1 = td->fRndm->uniformV();
    PhysDV r2 = td->fRndm->uniformV();

    PhysDV uvar  = -Log(r0 * r1);
    PhysDM tmpM2 = 9. > 36. * r2;
    vecCore::MaskedAssign(uvar, tmpM2, uvar * 1.6);
    vecCore::MaskedAssign(uvar, !tmpM2, uvar * 0.53333);

    const PhysDV thetaElectron = uvar * geant::units::kElectronMassC2 / electronTotE;
    PhysDV sintEle, costEle;
    SinCos(thetaElectron, &sintEle, &costEle);
    const PhysDV thetaPositron = uvar * geant::units::kElectronMassC2 / positronTotE;
    PhysDV sintPos, costPos;
    SinCos(thetaPositron, &sintPos, &costPos);
    sintPos          = -sintPos;
    const PhysDV phi = geant::units::kTwoPi * td->fRndm->uniformV();
    PhysDV sinphi, cosphi;
    SinCos(phi, &sinphi, &cosphi);

    PhysDV eleDirX = sintEle * cosphi;
    PhysDV eleDirY = sintEle * sinphi;
    PhysDV eleDirZ = costEle;

    PhysDV posDirX = sintPos * cosphi;
    PhysDV posDirY = sintPos * sinphi;
    PhysDV posDirZ = costPos;

    for (int l = 0; l < kPhysDVWidth; ++l) {
      tracks.SetTrackStatus(LTrackStatus::kKill, i + l);
      tracks.SetKinE(0.0, i + l);
    }

    const PhysDV ekinElectron = Max((electronTotE - geant::units::kElectronMassC2), (PhysDV)0.);
    const PhysDV ekinPositron = Max((positronTotE - geant::units::kElectronMassC2), (PhysDV)0.);
    // 5. rotate direction back to the lab frame: current directions are relative to the photon dir as z-dir
    PhysDV gammaX = tracks.GetDirXVec(i);
    PhysDV gammaY = tracks.GetDirYVec(i);
    PhysDV gammaZ = tracks.GetDirZVec(i);
    RotateToLabFrame(eleDirX, eleDirY, eleDirZ, gammaX, gammaY, gammaZ);
    RotateToLabFrame(posDirX, posDirY, posDirZ, gammaX, gammaY, gammaZ);

    for (int l = 0; l < kPhysDVWidth; ++l) {
      auto &secondarySoA = td->fPhysicsData->GetSecondarySOA();

      // Add e-
      int idx = secondarySoA.InsertTrack();
      secondarySoA.SetDirX(eleDirX[l], idx);
      secondarySoA.SetDirY(eleDirY[l], idx);
      secondarySoA.SetDirZ(eleDirZ[l], idx);
      secondarySoA.SetKinE(ekinElectron[l], idx);
      secondarySoA.SetGVcode(fElectronInternalCode, idx);
      secondarySoA.SetMass(geant::units::kElectronMassC2, idx);
      secondarySoA.SetTrackIndex(tracks.GetTrackIndex(i + l), idx); // parent Track index

      // Add e+
      idx = secondarySoA.InsertTrack();
      secondarySoA.SetDirX(posDirX[l], idx);
      secondarySoA.SetDirY(posDirY[l], idx);
      secondarySoA.SetDirZ(posDirZ[l], idx);
      secondarySoA.SetKinE(ekinPositron[l], idx);
      secondarySoA.SetGVcode(fPositronInternalCode, idx);
      secondarySoA.SetMass(geant::units::kElectronMassC2, idx);
      secondarySoA.SetTrackIndex(tracks.GetTrackIndex(i + l), idx); // parent Track index
    }
  }
}

PhysDV VecBetheHeitlerPairModel::SampleTotalEnergyTransferAliasOneShot(PhysDV egamma, const int *izet, PhysDV r1,
                                                                       PhysDV r2, PhysDV r3)
{
  const PhysDV legamma = Log(egamma);

  PhysDV val        = (legamma - fSTLogMinPhotonEnergy) * fSTILDeltaPhotonEnergy;
  PhysDI indxEgamma = (PhysDI)val; // lower electron energy bin index
  PhysDV pIndxHigh  = val - indxEgamma;
  PhysDM mask       = r1 < pIndxHigh;
  if (!mask.isEmpty()) {
    vecCore::MaskedAssign(indxEgamma, mask, indxEgamma + 1);
  }

  PhysDV xiV;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    int idx                  = (int)indxEgamma[l];
    RatinAliasDataTrans &als = fAliasTablesPerZ[izet[l]].fTablePerEn[idx];
    double xi = AliasTableAlternative::SampleRatin(als, fSTNumDiscreteEnergyTransferVals, r2[l], r3[l], 0);

    xiV[l] = xi;
  }
  PhysDV deltaMax;
  PhysDV deltaFactor;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    deltaMax[l] = egamma[l] > 50. * geant::units::MeV ? gElementData[izet[l]]->fDeltaMaxHighTsai
                                                      : gElementData[izet[l]]->fDeltaMaxLowTsai;
    deltaFactor[l] = gElementData[izet[l]]->fDeltaFactor;
  }
  const PhysDV eps0   = geant::units::kElectronMassC2 / egamma;
  const PhysDV epsp   = 0.5 - 0.5 * Sqrt(1. - 4. * eps0 * deltaFactor / deltaMax);
  const PhysDV epsMin = Max(eps0, epsp);
  const PhysDV epsV   = epsMin * Exp(xiV * Log(0.5 / epsMin));
  return epsV;
}

void VecBetheHeitlerPairModel::SampleTotalEnergyTransferRejVec(const double *egamma, const int *izet, double *epsOut,
                                                               int N, geant::TaskData *td)
{
  // assert(N>=kPhysDVWidth)
  int currN        = 0;
  PhysDM lanesDone = PhysDM::Zero();
  PhysDI idx;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    idx[l] = currN++;
  }

  while (currN < N || !lanesDone.isFull()) {
    PhysDV fz;
    PhysDV deltaMax;
    PhysDV deltaFac;

    for (int l = 0; l < kPhysDVWidth; ++l) {
      int lIdx        = idx[l];
      int lZet        = izet[lIdx];
      bool egammaHigh = egamma[lIdx] > 50. * geant::units::MeV;

      fz[l] = egammaHigh ? gElementData[lZet]->fFzHigh : gElementData[lZet]->fFzLow;
      if (fIsUseTsaisScreening) {
        deltaMax[l] = egammaHigh ? gElementData[lZet]->fDeltaMaxHighTsai : gElementData[lZet]->fDeltaMaxLowTsai;
      } else {
        deltaMax[l] = egammaHigh ? gElementData[lZet]->fDeltaMaxHigh : gElementData[lZet]->fDeltaMaxLow;
      }

      deltaFac[l] = gElementData[lZet]->fDeltaFactor;
    }

    PhysDV egammaV        = vecCore::Gather<PhysDV>(egamma, idx);
    const PhysDV eps0     = geant::units::kElectronMassC2 / egammaV;
    const PhysDV deltaMin = 4. * eps0 * deltaFac;
    const PhysDV eps1     = 0.5 - 0.5 * Sqrt(1. - deltaMin / deltaMax);
    const PhysDV epsMin   = Max(eps0, eps1);
    const PhysDV epsRange = 0.5 - epsMin;

    PhysDV F10, F20;
    ScreenFunction12(F10, F20, deltaMin, fIsUseTsaisScreening);
    F10 -= fz;
    F20 -= fz;

    const PhysDV NormF1   = Max(F10 * epsRange * epsRange, (PhysDV)0.);
    const PhysDV NormF2   = Max(1.5 * F20, (PhysDV)0.);
    const PhysDV NormCond = NormF1 / (NormF1 + NormF2);

    PhysDV rnd0 = td->fRndm->uniformV();
    PhysDV rnd1 = td->fRndm->uniformV();
    PhysDV rnd2 = td->fRndm->uniformV();

    PhysDV eps     = 0.0;
    PhysDV greject = 0.0;
    PhysDM cond1   = NormCond > rnd0;
    vecCore::MaskedAssign(eps, cond1, 0.5 - epsRange * Pow(rnd1, (PhysDV)1. / 3.));
    vecCore::MaskedAssign(eps, !cond1, epsMin + epsRange * rnd1);
    const PhysDV delta = deltaFac * eps0 / (eps * (1. - eps));
    if (cond1.isNotEmpty()) {
      vecCore::MaskedAssign(greject, cond1, (ScreenFunction1(delta, fIsUseTsaisScreening) - fz) / F10);
    }
    if ((!cond1).isNotEmpty()) {
      vecCore::MaskedAssign(greject, !cond1, (ScreenFunction2(delta, fIsUseTsaisScreening) - fz) / F20);
    }

    PhysDM accepted = greject > rnd2;

    if (accepted.isNotEmpty()) {
      vecCore::Scatter(eps, epsOut, idx);
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

void VecBetheHeitlerPairModel::ScreenFunction12(PhysDV &val1, PhysDV &val2, const PhysDV delta, const bool istsai)
{
  if (istsai) {
    const PhysDV gamma  = delta * 0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const PhysDV gamma2 = gamma * gamma;
    const PhysDV dum1 = 33.726 - 4. * Log(1.0 + 0.311877 * gamma2) + 4.8 * Exp(-0.9 * gamma) + 3.2 * Exp(-1.5 * gamma);
    const PhysDV dum2 = 2. / (3. + 19.5 * gamma + 18. * gamma2);
    val1              = dum1 + dum2;
    val2              = dum1 - 0.5 * dum2;
  } else {
    PhysDM tmp = delta > 1.;
    vecCore::MaskedAssign(val1, tmp, 42.24 - 8.368 * Log(delta + 0.952));
    vecCore::MaskedAssign(val2, tmp, val1);
    vecCore::MaskedAssign(val1, !tmp, 42.392 - delta * (7.796 - 1.961 * delta));
    vecCore::MaskedAssign(val2, !tmp, 41.405 - delta * (5.828 - 0.8945 * delta));
  }
}

// 3xPhi_1 - Phi_2: used in case of rejection (either istsai or not)
PhysDV VecBetheHeitlerPairModel::ScreenFunction1(const PhysDV delta, const bool istsai)
{
  PhysDV val;
  if (istsai) {
    const PhysDV gamma  = delta * 0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const PhysDV gamma2 = gamma * gamma;
    val = 33.726 - 4. * Log(1.0 + 0.311877 * gamma2) + 4.8 * Exp(-0.9 * gamma) + 3.2 * Exp(-1.5 * gamma) +
          2. / (3. + 19.5 * gamma + 18. * gamma2);
  } else {
    PhysDM tmp = delta > 1.;
    vecCore::MaskedAssign(val, tmp, 42.24 - 8.368 * Log(delta + 0.952));
    vecCore::MaskedAssign(val, !tmp, 42.392 - delta * (7.796 - 1.961 * delta));
  }
  return val;
}

// 1.5*Phi_1 + 0.5*Phi_2: used in case of rejection (either istsai or not)
PhysDV VecBetheHeitlerPairModel::ScreenFunction2(const PhysDV delta, const bool istsai)
{
  PhysDV val;
  if (istsai) {
    const PhysDV gamma  = delta * 0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const PhysDV gamma2 = gamma * gamma;
    val = 33.726 - 4. * Log(1.0 + 0.311877 * gamma2) + 4.8 * Exp(-0.9 * gamma) + 3.2 * Exp(-1.5 * gamma) -
          1. / (3. + 19.5 * gamma + 18. * gamma2);
  } else {
    PhysDM tmp = delta > 1.;
    vecCore::MaskedAssign(val, tmp, 42.24 - 8.368 * Log(delta + 0.952));
    vecCore::MaskedAssign(val, !tmp, 41.405 - delta * (5.828 - 0.8945 * delta));
  }
  return val;
}

bool VecBetheHeitlerPairModel::IsModelUsable(const MaterialCuts *, double ekin)
{
  double invEps = ekin * geant::units::kInvElectronMassC2;
  return ekin > GetLowEnergyUsageLimit() && ekin < GetHighEnergyUsageLimit() && invEps > 2.0;
}
}
