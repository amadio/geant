#include <Geant/SystemOfUnits.h>
#include <Geant/PhysicalConstants.h>
#include <Geant/TaskData.h>
#include <Geant/PhysicsData.h>
#include <Geant/MaterialCuts.h>
#include <TError.h>
#include "Geant/VecRelativisticPairModel.h"
#include "Geant/AliasTable.h"
#include "Geant/MaterialProperties.h"

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

void VecRelativisticPairModel::Initialize()
{
  RelativisticPairModel::Initialize();

  fAliasTablesPerMaterial.resize(fSamplingTables.size());

  // Steal sampling tables to new format
  for (size_t i = 0; i < fSamplingTables.size(); ++i) { // Loop over material
    if (fSamplingTables[i] == nullptr) continue;

    fAliasTablesPerMaterial[i].fILowestZ = fSamplingTables[i]->fILowestZ;
    auto &defaultGammaETables            = fSamplingTables[i]->fRatinAliasData;
    auto &transposedGammaETables         = fAliasTablesPerMaterial[i].fTablePerEn;
    for (size_t j = 0; j < defaultGammaETables.size(); ++j) { // Loop over incoming gamma Es in material table
      transposedGammaETables.emplace_back(fSTNumDiscreteEnergyTransferVals);
      // Loop over transferred E in corresponding gamma E tbl.
      for (int k = 0; k < fSTNumDiscreteEnergyTransferVals; ++k) {
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

void VecRelativisticPairModel::SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td)
{
  int N             = tracks.GetNtracks();
  int *IZet         = td->fPhysicsData->fPhysicsScratchpad.fIzet;   // used by rej method
  int *MatIDX       = td->fPhysicsData->fPhysicsScratchpad.fMatIdx; // used by alias method
  double *LPMEnergy = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr;

  // Sort by LPM energies to simplify branching in rej. accept. sampling method
  if (!GetUseSamplingTables() && fIsUseLPM) {
    short *sortKey = tracks.GetSortKeyV();
    for (int i = 0; i < N; ++i) {
      double ekin = tracks.GetKinE(i);
      sortKey[i]  = ekin < fLPMEnergyLimit ? (short)0 : (short)1;
    }
    tracks.SortTracks();
  }

  // Populate material property arrays
  for (int i = 0; i < N; ++i) {
    MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[tracks.GetMaterialCutCoupleIndex(i)];
    const Material *mat  = matCut->GetMaterial();

    if (GetUseSamplingTables()) {
      MatIDX[i] = mat->GetIndex();
    } else {
      const double lpmEnergy = mat->GetMaterialProperties()->GetRadiationLength() * gLPMFactor;
      LPMEnergy[i]           = lpmEnergy;

      const Vector_t<Element *> &theElements = mat->GetElementVector();
      double targetElemIndx                  = 0;
      if (theElements.size() > 1) {
        targetElemIndx = SampleTargetElementIndex(matCut, tracks.GetKinE(i), td->fRndm->uniform());
      }
      const double zet = theElements[targetElemIndx]->GetZ();
      const int izet   = std::min(std::lrint(zet), gMaxZet - 1);
      IZet[i]          = izet;
    }
  }

  if (GetUseSamplingTables()) {
    for (int i = 0; i < N; i += kVecLenD) {
      Double_v ekin = tracks.GetKinEVec(i);
      Double_v r1   = td->fRndm->uniformV();
      Double_v r2   = td->fRndm->uniformV();
      Double_v r3   = td->fRndm->uniformV();

      Double_v sampledEps = SampleTotalEnergyTransferAliasOneShot(ekin, &MatIDX[i], r1, r2, r3);
      vecCore::Store(sampledEps, &td->fPhysicsData->fPhysicsScratchpad.fEps[i]);
    }
  } else {

    tracks.GetKinEArr()[N] = tracks.GetKinEArr()[N - 1];
    LPMEnergy[N]           = LPMEnergy[N - 1];
    IZet[N]                = IZet[N - 1];
    SampleTotalEnergyTransferRejVec(tracks.GetKinEArr(), LPMEnergy, IZet, td->fPhysicsData->fPhysicsScratchpad.fEps, N,
                                    td);
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
      tracks.SetKinE(0.0, i + l);
      tracks.SetTrackStatus(LTrackStatus::kKill, i + l);
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

      int idx = secondarySoA.InsertTrack();
      secondarySoA.SetDirX(Get(eleDirX, l), idx);
      secondarySoA.SetDirY(Get(eleDirY, l), idx);
      secondarySoA.SetDirZ(Get(eleDirZ, l), idx);
      secondarySoA.SetKinE(Get(ekinElectron, l), idx);
      secondarySoA.SetGVcode(fElectronInternalCode, idx);
      secondarySoA.SetMass(geant::units::kElectronMassC2, idx);
      secondarySoA.SetTrackIndex(tracks.GetTrackIndex(i + l), idx); // parent Track index
      // then set the e+
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

Double_v VecRelativisticPairModel::SampleTotalEnergyTransferAliasOneShot(const Double_v egamma, const int *matIDX,
                                                                         const Double_v r1, const Double_v r2,
                                                                         const Double_v r3)
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
    RatinAliasDataTrans &als = fAliasTablesPerMaterial[matIDX[l]].fTablePerEn[idx];
    double xi = AliasTableAlternative::SampleRatin(als, fSTNumDiscreteEnergyTransferVals, Get(r2, l), Get(r3, l), 0);
    Set(xiV, l, xi);
  }

  Double_v deltaMax;
  Double_v deltaFactor;
  for (int l = 0; l < kVecLenD; ++l) {
    int izet = fAliasTablesPerMaterial[matIDX[l]].fILowestZ;
    Set(deltaMax, l, gElementData[izet]->fDeltaMaxTsai);
    Set(deltaFactor, l, gElementData[izet]->fDeltaFactor);
  }
  const Double_v eps0   = geant::units::kElectronMassC2 / egamma;
  const Double_v epsp   = 0.5 - 0.5 * Math::Sqrt(1. - 4. * eps0 * deltaFactor / deltaMax);
  const Double_v epsMin = Math::Max(eps0, epsp);
  const Double_v epsV   = epsMin * Math::Exp(xiV * Math::Log(0.5 / epsMin));
  return epsV;
}

void VecRelativisticPairModel::SampleTotalEnergyTransferRejVec(const double *egamma, const double *lpmEnergy,
                                                               const int *izet, double *epsOut, int N,
                                                               geant::TaskData *td)
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
    Double_v varS1Cond, ilVarS1Cond;

    for (int l = 0; l < kVecLenD; ++l) {
      int lZet = izet[Get(idx, l)];
      Set(fz, l, gElementData[lZet]->fFz);
      Set(deltaMax, l, fIsUseTsaisScreening ? gElementData[lZet]->fDeltaMaxTsai : gElementData[lZet]->fDeltaMax);
      Set(deltaFac, l, gElementData[lZet]->fDeltaFactor);
      Set(varS1Cond, l, gElementData[lZet]->fVarS1Cond);
      Set(ilVarS1Cond, l, gElementData[lZet]->fILVarS1Cond);
    }

    Double_v egammaV        = vecCore::Gather<Double_v>(egamma, idx);
    const Double_v eps0     = geant::units::kElectronMassC2 / egammaV;
    const Double_v deltaMin = 4. * eps0 * deltaFac;
    const Double_v epsp     = 0.5 - 0.5 * Math::Sqrt(1. - deltaMin / deltaMax);
    const Double_v epsMin   = Math::Max(eps0, epsp);
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

    MaskD_v lpmMask     = MaskD_v(fIsUseLPM) && egammaV > fLPMEnergyLimit;
    Double_v lpmEnergyV = vecCore::Gather<Double_v>(lpmEnergy, idx);
    if (!MaskEmpty(cond1)) {
      if (!MaskEmpty(lpmMask)) {
        Double_v lpmPhiS, lpmGS, lpmXiS, phi1, phi2;
        ComputeScreeningFunctions(phi1, phi2, delta, fIsUseTsaisScreening);
        ComputeLPMfunctions(lpmXiS, lpmGS, lpmPhiS, lpmEnergyV, eps, egammaV, varS1Cond, ilVarS1Cond);
        Double_v tmpRej = lpmXiS * ((lpmGS + 2. * lpmPhiS) * phi1 - lpmGS * phi2 - lpmPhiS * fz) / F10;
        vecCore::MaskedAssign(greject, cond1 && lpmMask, tmpRej);
      }
      if (!MaskEmpty(!lpmMask)) {
        vecCore::MaskedAssign(greject, cond1 && !lpmMask, (ScreenFunction1(delta, fIsUseTsaisScreening) - fz) / F10);
      }
    }
    if (!MaskEmpty(!cond1)) {
      if (!MaskEmpty(lpmMask)) {
        Double_v lpmPhiS, lpmGS, lpmXiS, phi1, phi2;
        ComputeScreeningFunctions(phi1, phi2, delta, fIsUseTsaisScreening);
        ComputeLPMfunctions(lpmXiS, lpmGS, lpmPhiS, lpmEnergyV, eps, egammaV, varS1Cond, ilVarS1Cond);
        Double_v tmpRej =
            lpmXiS * ((0.5 * lpmGS + lpmPhiS) * phi1 + 0.5 * lpmGS * phi2 - 0.5 * (lpmGS + lpmPhiS) * fz) / F20;
        vecCore::MaskedAssign(greject, !cond1 && lpmMask, tmpRej);
      }
      if (!MaskEmpty(!lpmMask)) {
        vecCore::MaskedAssign(greject, !cond1 && !lpmMask, (ScreenFunction2(delta, fIsUseTsaisScreening) - fz) / F20);
      }
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

void VecRelativisticPairModel::ScreenFunction12(Double_v &val1, Double_v &val2, const Double_v delta, const bool istsai)
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
Double_v VecRelativisticPairModel::ScreenFunction1(const Double_v delta, const bool istsai)
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
Double_v VecRelativisticPairModel::ScreenFunction2(const Double_v delta, const bool istsai)
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

void VecRelativisticPairModel::ComputeScreeningFunctions(Double_v &phi1, Double_v &phi2, const Double_v delta,
                                                         const bool istsai)
{
  if (istsai) {
    const Double_v gamma  = delta * 0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const Double_v gamma2 = gamma * gamma;
    phi1                  = 16.863 - 2.0 * Math::Log(1.0 + 0.311877 * gamma2) + 2.4 * Math::Exp(-0.9 * gamma) +
           1.6 * Math::Exp(-1.5 * gamma);
    phi2 = phi1 - 2.0 / (3.0 + 19.5 * gamma + 18.0 * gamma2); // phi1-phi2
  } else {
    MaskD_v tmp = delta > 1.;
    vecCore::MaskedAssign(phi1, tmp, 21.12 - 4.184 * Math::Log(delta + 0.952));
    vecCore::MaskedAssign(phi2, tmp, phi1);

    Double_v delta2 = delta * delta;
    vecCore::MaskedAssign(phi1, !tmp, 20.867 - 3.242 * delta + 0.625 * delta2);
    vecCore::MaskedAssign(phi2, !tmp, 20.209 - 1.93 * delta - 0.086 * delta2);
  }
}

void VecRelativisticPairModel::ComputeLPMfunctions(Double_v &funcXiS, Double_v &funcGS, Double_v &funcPhiS,
                                                   Double_v lpmenergy, Double_v eps, Double_v egamma,
                                                   Double_v varS1Cond, Double_v ilVarS1Cond)
{

  //  1. y = E_+/E_{\gamma} with E_+ being the total energy transfered to one of the e-/e+ pair
  //     s'  = \sqrt{ \frac{1}{8} \frac{1}{y(1-y)}   \frac{E^{KL}_{LPM}}{E_{\gamma}}  }
  const Double_v varSprime = Math::Sqrt(0.125 * lpmenergy / (eps * egamma * (1.0 - eps)));
  funcXiS                  = 2.0;

  MaskD_v tmpM = varSprime > 1.0;
  vecCore::MaskedAssign(funcXiS, tmpM, (Double_v)1.0);
  tmpM = varSprime > varS1Cond;
  if (!MaskEmpty(tmpM)) {
    const Double_v funcHSprime = Math::Log(varSprime) * ilVarS1Cond;
    Double_v tmpFuncXiS =
        1.0 + funcHSprime - 0.08 * (1.0 - funcHSprime) * funcHSprime * (2.0 - funcHSprime) * ilVarS1Cond;
    vecCore::MaskedAssign(funcXiS, tmpM, tmpFuncXiS);
  }

  const Double_v varShat = varSprime / Math::Sqrt(funcXiS);

  for (int l = 0; l < kVecLenD; ++l) {
    double lFuncGS, lFuncPhiS, lVarShat;
    double lFuncXiS;
    lVarShat = Get(varShat, l);
    lFuncXiS = Get(funcXiS, l);
    GetLPMFunctions(lFuncGS, lFuncPhiS, lVarShat);
    Set(funcGS, l, lFuncGS);
    Set(funcPhiS, l, lFuncPhiS);

    if (lFuncXiS * lFuncPhiS > 1. || lVarShat > 0.57) {
      Set(funcXiS, l, 1. / lFuncPhiS);
    }
  }
  // MAKE SURE SUPPRESSION IS SMALLER THAN 1: due to Migdal's approximation on xi
}
}
