#include <Geant/SystemOfUnits.h>
#include <Geant/PhysicalConstants.h>
#include <Geant/TaskData.h>
#include <Geant/PhysicsData.h>
#include <Geant/MaterialCuts.h>
#include <TError.h>
#include "Geant/VecBetheHeitlerPairModel.h"
#include "Geant/AliasTable.h"

namespace geantphysics {

using geant::Double_v;
using geant::IndexD_v;
using geant::kVecLenD;
using geant::MaskD_v;
using vecCore::Get;
using vecCore::Set;
using vecCore::AssignMaskLane;
using vecCore::MaskFull;
using vecCore::MaskEmpty;

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
      int i = (smallEtracksInd / kVecLenD) * kVecLenD;
      if (i < N) {
        Double_v ekin    = tracks.GetKinEVec(i);
        MaskD_v smallE   = ekin < fGammaEneregyLimit;
        Double_v r1      = td->fRndm->uniformV();
        Double_v r2      = td->fRndm->uniformV();
        Double_v r3      = td->fRndm->uniformV();
        Double_v tmpEkin = ekin;
        // To be sure that energy that is too small for alias table will not slip in alias sampling
        vecCore::MaskedAssign(tmpEkin, smallE, (Double_v)fGammaEneregyLimit * (1.1));

        Double_v sampledEps = SampleTotalEnergyTransferAliasOneShot(tmpEkin, &IZet[i], r1, r2, r3);
        Double_v eps;
        vecCore::Load(eps, &td->fPhysicsData->fPhysicsScratchpad.fEps[i]);
        vecCore::MaskedAssign(eps, !smallE, sampledEps);
        vecCore::Store(eps, &td->fPhysicsData->fPhysicsScratchpad.fEps[i]);
      }
    }
    for (int i = (smallEtracksInd / kVecLenD) * kVecLenD + kVecLenD; i < N; i += kVecLenD) {
      Double_v ekin;
      vecCore::Load(ekin, tracks.GetKinEArr(i));
      Double_v r1 = td->fRndm->uniformV();
      Double_v r2 = td->fRndm->uniformV();
      Double_v r3 = td->fRndm->uniformV();

      Double_v sampledEps = SampleTotalEnergyTransferAliasOneShot(ekin, &IZet[i], r1, r2, r3);
      vecCore::Store(sampledEps, &td->fPhysicsData->fPhysicsScratchpad.fEps[i]);
    }
  } else {

    int i = (smallEtracksInd / kVecLenD) * kVecLenD;
    if (i < N) {
      // Always create fake particle int the end of real particle array for rejection sampling
      tracks.GetKinEArr()[N] = tracks.GetKinEArr()[N - 1];
      IZet[N]                = IZet[N - 1];
      SampleTotalEnergyTransferRejVec(tracks.GetKinEArr(i), &IZet[i], &td->fPhysicsData->fPhysicsScratchpad.fEps[i],
                                      N - i, td);
    }
  }

  for (int i = 0; i < N; i += kVecLenD) {
    Double_v ekin = tracks.GetKinEVec(i);

    Double_v eps;
    vecCore::Load(eps, &td->fPhysicsData->fPhysicsScratchpad.fEps[i]);

    Double_v electronTotE, positronTotE;
    MaskD_v tmpM = td->fRndm->uniformV() > 0.5;
    vecCore::MaskedAssign(electronTotE, tmpM, (1. - eps) * ekin);
    vecCore::MaskedAssign(positronTotE, tmpM, eps * ekin);

    vecCore::MaskedAssign(positronTotE, !tmpM, (1. - eps) * ekin);
    vecCore::MaskedAssign(electronTotE, !tmpM, eps * ekin);

    Double_v r0 = td->fRndm->uniformV();
    Double_v r1 = td->fRndm->uniformV();
    Double_v r2 = td->fRndm->uniformV();

    Double_v uvar = -Math::Log(r0 * r1);
    MaskD_v tmpM2 = 9. > 36. * r2;
    vecCore::MaskedAssign(uvar, tmpM2, uvar * 1.6);
    vecCore::MaskedAssign(uvar, !tmpM2, uvar * 0.53333);

    const Double_v thetaElectron = uvar * geant::units::kElectronMassC2 / electronTotE;
    Double_v sintEle, costEle;
    Math::SinCos(thetaElectron, sintEle, costEle);
    const Double_v thetaPositron = uvar * geant::units::kElectronMassC2 / positronTotE;
    Double_v sintPos, costPos;
    Math::SinCos(thetaPositron, sintPos, costPos);
    sintPos            = -sintPos;
    const Double_v phi = geant::units::kTwoPi * td->fRndm->uniformV();
    Double_v sinphi, cosphi;
    Math::SinCos(phi, sinphi, cosphi);

    Double_v eleDirX = sintEle * cosphi;
    Double_v eleDirY = sintEle * sinphi;
    Double_v eleDirZ = costEle;

    Double_v posDirX = sintPos * cosphi;
    Double_v posDirY = sintPos * sinphi;
    Double_v posDirZ = costPos;

    for (int l = 0; l < kVecLenD; ++l) {
      tracks.SetTrackStatus(LTrackStatus::kKill, i + l);
      tracks.SetKinE(0.0, i + l);
    }

    const Double_v ekinElectron = Math::Max((electronTotE - geant::units::kElectronMassC2), (Double_v)0.);
    const Double_v ekinPositron = Math::Max((positronTotE - geant::units::kElectronMassC2), (Double_v)0.);
    // 5. rotate direction back to the lab frame: current directions are relative to the photon dir as z-dir
    Double_v gammaX = tracks.GetDirXVec(i);
    Double_v gammaY = tracks.GetDirYVec(i);
    Double_v gammaZ = tracks.GetDirZVec(i);
    RotateToLabFrame(eleDirX, eleDirY, eleDirZ, gammaX, gammaY, gammaZ);
    RotateToLabFrame(posDirX, posDirY, posDirZ, gammaX, gammaY, gammaZ);

    for (int l = 0; l < kVecLenD; ++l) {
      auto &secondarySoA = td->fPhysicsData->GetSecondarySOA();

      // Add e-
      int idx = secondarySoA.InsertTrack();
      secondarySoA.SetDirX(Get(eleDirX, l), idx);
      secondarySoA.SetDirY(Get(eleDirY, l), idx);
      secondarySoA.SetDirZ(Get(eleDirZ, l), idx);
      secondarySoA.SetKinE(Get(ekinElectron, l), idx);
      secondarySoA.SetGVcode(fElectronInternalCode, idx);
      secondarySoA.SetMass(geant::units::kElectronMassC2, idx);
      secondarySoA.SetTrackIndex(tracks.GetTrackIndex(i + l), idx); // parent Track index

      // Add e+
      idx = secondarySoA.InsertTrack();
      secondarySoA.SetDirX(Get(posDirX, l), idx);
      secondarySoA.SetDirY(Get(posDirY, l), idx);
      secondarySoA.SetDirZ(Get(posDirZ, l), idx);
      secondarySoA.SetKinE(Get(ekinPositron, l), idx);
      secondarySoA.SetGVcode(fPositronInternalCode, idx);
      secondarySoA.SetMass(geant::units::kElectronMassC2, idx);
      secondarySoA.SetTrackIndex(tracks.GetTrackIndex(i + l), idx); // parent Track index
    }
  }
}

Double_v VecBetheHeitlerPairModel::SampleTotalEnergyTransferAliasOneShot(Double_v egamma, const int *izet, Double_v r1,
                                                                         Double_v r2, Double_v r3)
{
  const Double_v legamma = Math::Log(egamma);

  Double_v val        = (legamma - fSTLogMinPhotonEnergy) * fSTILDeltaPhotonEnergy;
  IndexD_v indxEgamma = (IndexD_v)val; // lower electron energy bin index
  Double_v pIndxHigh  = val - indxEgamma;
  MaskD_v mask        = r1 < pIndxHigh;
  if (!MaskEmpty(mask)) {
    vecCore::MaskedAssign(indxEgamma, mask, indxEgamma + 1);
  }

  Double_v xiV;
  for (int l = 0; l < kVecLenD; ++l) {
    int idx                  = (int)Get(indxEgamma, l);
    RatinAliasDataTrans &als = fAliasTablesPerZ[izet[l]].fTablePerEn[idx];
    double xi = AliasTableAlternative::SampleRatin(als, fSTNumDiscreteEnergyTransferVals, Get(r2, l), Get(r3, l), 0);

    Set(xiV, l, xi);
  }
  Double_v deltaMax;
  Double_v deltaFactor;
  for (int l = 0; l < kVecLenD; ++l) {
    Set(deltaMax, l, Get(egamma, l) > 50. * geant::units::MeV ? gElementData[izet[l]]->fDeltaMaxHighTsai
                                                              : gElementData[izet[l]]->fDeltaMaxLowTsai);
    Set(deltaFactor, l, gElementData[izet[l]]->fDeltaFactor);
  }
  const Double_v eps0   = geant::units::kElectronMassC2 / egamma;
  const Double_v epsp   = 0.5 - 0.5 * Math::Sqrt(1. - 4. * eps0 * deltaFactor / deltaMax);
  const Double_v epsMin = Math::Max(eps0, epsp);
  const Double_v epsV   = epsMin * Math::Exp(xiV * Math::Log(0.5 / epsMin));
  return epsV;
}

void VecBetheHeitlerPairModel::SampleTotalEnergyTransferRejVec(const double *egamma, const int *izet, double *epsOut,
                                                               int N, geant::TaskData *td)
{
  // assert(N>=kVecLenD)
  int currN = 0;
  MaskD_v lanesDone;
  IndexD_v idx;
  for (int l = 0; l < kVecLenD; ++l) {
    AssignMaskLane(lanesDone, l, false);
    Set(idx, l, currN++);
  }

  while (currN < N || !MaskFull(lanesDone)) {
    Double_v fz;
    Double_v deltaMax;
    Double_v deltaFac;

    for (int l = 0; l < kVecLenD; ++l) {
      int lIdx        = Get(idx, l);
      int lZet        = izet[lIdx];
      bool egammaHigh = egamma[lIdx] > 50. * geant::units::MeV;

      Set(fz, l, egammaHigh ? gElementData[lZet]->fFzHigh : gElementData[lZet]->fFzLow);
      if (fIsUseTsaisScreening) {
        Set(deltaMax, l, egammaHigh ? gElementData[lZet]->fDeltaMaxHighTsai : gElementData[lZet]->fDeltaMaxLowTsai);
      } else {
        Set(deltaMax, l, egammaHigh ? gElementData[lZet]->fDeltaMaxHigh : gElementData[lZet]->fDeltaMaxLow);
      }

      Set(deltaFac, l, gElementData[lZet]->fDeltaFactor);
    }

    Double_v egammaV        = vecCore::Gather<Double_v>(egamma, idx);
    const Double_v eps0     = geant::units::kElectronMassC2 / egammaV;
    const Double_v deltaMin = 4. * eps0 * deltaFac;
    const Double_v eps1     = 0.5 - 0.5 * Math::Sqrt(1. - deltaMin / deltaMax);
    const Double_v epsMin   = Math::Max(eps0, eps1);
    const Double_v epsRange = 0.5 - epsMin;

    Double_v F10, F20;
    ScreenFunction12(F10, F20, deltaMin, fIsUseTsaisScreening);
    F10 -= fz;
    F20 -= fz;

    const Double_v NormF1   = Math::Max(F10 * epsRange * epsRange, (Double_v)0.);
    const Double_v NormF2   = Math::Max(1.5 * F20, (Double_v)0.);
    const Double_v NormCond = NormF1 / (NormF1 + NormF2);

    Double_v rnd0 = td->fRndm->uniformV();
    Double_v rnd1 = td->fRndm->uniformV();
    Double_v rnd2 = td->fRndm->uniformV();

    Double_v eps     = 0.0;
    Double_v greject = 0.0;
    MaskD_v cond1    = NormCond > rnd0;
    vecCore::MaskedAssign(eps, cond1, 0.5 - epsRange * Math::Pow(rnd1, (Double_v)1. / 3.));
    vecCore::MaskedAssign(eps, !cond1, epsMin + epsRange * rnd1);
    const Double_v delta = deltaFac * eps0 / (eps * (1. - eps));
    if (!MaskEmpty(cond1)) {
      vecCore::MaskedAssign(greject, cond1, (ScreenFunction1(delta, fIsUseTsaisScreening) - fz) / F10);
    }
    if (!MaskEmpty(!cond1)) {
      vecCore::MaskedAssign(greject, !cond1, (ScreenFunction2(delta, fIsUseTsaisScreening) - fz) / F20);
    }

    MaskD_v accepted = greject > rnd2;

    if (!MaskEmpty(accepted)) {
      vecCore::Scatter(eps, epsOut, idx);
    }

    lanesDone = lanesDone || accepted;
    for (int l = 0; l < kVecLenD; ++l) {
      auto laneDone = Get(accepted, l);
      if (laneDone) {
        if (currN < N) {
          Set(idx, l, currN++);
          AssignMaskLane(lanesDone, l, false);
        } else {
          Set(idx, l, N);
        }
      }
    }
  }
}

void VecBetheHeitlerPairModel::ScreenFunction12(Double_v &val1, Double_v &val2, const Double_v delta, const bool istsai)
{
  if (istsai) {
    const Double_v gamma  = delta * 0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const Double_v gamma2 = gamma * gamma;
    const Double_v dum1   = 33.726 - 4. * Math::Log(1.0 + 0.311877 * gamma2) + 4.8 * Math::Exp(-0.9 * gamma) +
                          3.2 * Math::Exp(-1.5 * gamma);
    const Double_v dum2 = 2. / (3. + 19.5 * gamma + 18. * gamma2);
    val1                = dum1 + dum2;
    val2                = dum1 - 0.5 * dum2;
  } else {
    MaskD_v tmp = delta > 1.;
    vecCore::MaskedAssign(val1, tmp, 42.24 - 8.368 * Math::Log(delta + 0.952));
    vecCore::MaskedAssign(val2, tmp, val1);
    vecCore::MaskedAssign(val1, !tmp, 42.392 - delta * (7.796 - 1.961 * delta));
    vecCore::MaskedAssign(val2, !tmp, 41.405 - delta * (5.828 - 0.8945 * delta));
  }
}

// 3xPhi_1 - Phi_2: used in case of rejection (either istsai or not)
Double_v VecBetheHeitlerPairModel::ScreenFunction1(const Double_v delta, const bool istsai)
{
  Double_v val;
  if (istsai) {
    const Double_v gamma  = delta * 0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const Double_v gamma2 = gamma * gamma;
    val                   = 33.726 - 4. * Math::Log(1.0 + 0.311877 * gamma2) + 4.8 * Math::Exp(-0.9 * gamma) +
          3.2 * Math::Exp(-1.5 * gamma) + 2. / (3. + 19.5 * gamma + 18. * gamma2);
  } else {
    MaskD_v tmp = delta > 1.;
    vecCore::MaskedAssign(val, tmp, 42.24 - 8.368 * Math::Log(delta + 0.952));
    vecCore::MaskedAssign(val, !tmp, 42.392 - delta * (7.796 - 1.961 * delta));
  }
  return val;
}

// 1.5*Phi_1 + 0.5*Phi_2: used in case of rejection (either istsai or not)
Double_v VecBetheHeitlerPairModel::ScreenFunction2(const Double_v delta, const bool istsai)
{
  Double_v val;
  if (istsai) {
    const Double_v gamma  = delta * 0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const Double_v gamma2 = gamma * gamma;
    val                   = 33.726 - 4. * Math::Log(1.0 + 0.311877 * gamma2) + 4.8 * Math::Exp(-0.9 * gamma) +
          3.2 * Math::Exp(-1.5 * gamma) - 1. / (3. + 19.5 * gamma + 18. * gamma2);
  } else {
    MaskD_v tmp = delta > 1.;
    vecCore::MaskedAssign(val, tmp, 42.24 - 8.368 * Math::Log(delta + 0.952));
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
