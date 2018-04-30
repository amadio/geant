
#include "Geant/VecRelativisticBremsModel.h"

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

VecRelativisticBremsModel::VecRelativisticBremsModel(const std::string &modelname) : RelativisticBremsModel(modelname)
{
}

void VecRelativisticBremsModel::Initialize()
{
  RelativisticBremsModel::Initialize();

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

void VecRelativisticBremsModel::SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td)
{
  const int N               = tracks.GetNtracks();
  double *gammaeEnergyArray = td->fPhysicsData->fPhysicsScratchpad.fEps;
  double *gammaCutArr       = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr;
  //  int *IZet                 = td->fPhysicsData->fPhysicsScratchpad.fIzet;
  double *Zet                 = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr2;
  double *densityCorrConstArr = td->fPhysicsData->fPhysicsScratchpad.fR0;
  double *lpmEnergyArray      = td->fPhysicsData->fPhysicsScratchpad.fR1;

  for (int i = 0; i < N; i += kVecLenD) {

    Double_v primEkin = tracks.GetKinEVec(i);
    Double_v gammaCut;
    IndexD_v mcIndxLocal; // Used by alias method
    Double_v densityCorConst;
    Double_v lpmEnergy;
    for (int l = 0; l < kVecLenD; ++l) {
      const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(tracks.GetMaterialCutCoupleIndex(i + l));
      Set(gammaCut, l, matCut->GetProductionCutsInEnergy()[0]);
      Set(densityCorConst, l,
          matCut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol() * gDensityFactor);
      Set(lpmEnergy, l, matCut->GetMaterial()->GetMaterialProperties()->GetRadiationLength() * gLPMFactor);
      if (GetUseSamplingTables()) {
        Set(mcIndxLocal, l, fGlobalMatGCutIndxToLocal[matCut->GetIndex()]);
      } else {
        // sample target element
        const Vector_t<Element *> &theElements = matCut->GetMaterial()->GetElementVector();
        double targetElemIndx                  = 0;
        if (theElements.size() > 1) {
          targetElemIndx = SampleTargetElementIndex(matCut, Get(primEkin, l), td->fRndm->uniform());
        }
        Zet[i + l] = theElements[targetElemIndx]->GetZ();
        //        IZet[i + l]      = (int) std::lrint(zet);//std::min((int)std::lrint(zet), fDCSMaxZet);
      }
    }
    Double_v totalEn    = primEkin + geant::units::kElectronMassC2;
    Double_v densityCor = densityCorConst * (totalEn * totalEn);
    assert((primEkin >= gammaCut).isFull()); // Cut filtering should be applied up the call chain.
    if (GetUseSamplingTables()) {
      Double_v r1          = td->fRndm->uniformV();
      Double_v r2          = td->fRndm->uniformV();
      Double_v r3          = td->fRndm->uniformV();
      Double_v primEkin    = tracks.GetKinEVec(i);
      Double_v gammaEnergy = SampleEnergyTransfer(gammaCut, densityCor, mcIndxLocal, fAliasData.fLogEmin.data(),
                                                  fAliasData.fILDelta.data(), primEkin, r1, r2, r3);
      vecCore::Store(gammaEnergy, gammaeEnergyArray + i);
    } else {
      vecCore::Store(gammaCut, gammaCutArr + i);
      vecCore::Store(densityCorConst, densityCorrConstArr + i);
      vecCore::Store(lpmEnergy, lpmEnergyArray + i);
    }
  }

  if (!GetUseSamplingTables()) {
    tracks.GetKinEArr()[N] = tracks.GetKinEArr()[N - 1];
    gammaCutArr[N]         = gammaCutArr[N - 1];
    //    IZet[N]                = IZet[N - 1];
    Zet[N]                 = Zet[N - 1];
    densityCorrConstArr[N] = densityCorrConstArr[N - 1];
    lpmEnergyArray[N]      = lpmEnergyArray[N - 1];

    SampleEnergyTransfer(tracks.GetKinEArr(), gammaCutArr, Zet, densityCorrConstArr, lpmEnergyArray, gammaeEnergyArray,
                         N, td);
  }

  for (int i = 0; i < N; i += kVecLenD) {
    Double_v ekin = tracks.GetKinEVec(i);
    Double_v gammaEnergy;
    vecCore::Load(gammaEnergy, gammaeEnergyArray + i);

    Double_v cosTheta = 1.0;
    Double_v sinTheta = 0.0;
    Double_v rnd0     = td->fRndm->uniformV();
    SamplePhotonDirection(ekin, sinTheta, cosTheta, rnd0);

    Double_v rnd1      = td->fRndm->uniformV();
    const Double_v phi = geant::units::kTwoPi * (rnd1);
    Double_v sinPhi, cosPhi;
    Math::SinCos(phi, sinPhi, cosPhi);
    // gamma direction in the scattering frame
    Double_v gamDirX = sinTheta * cosPhi;
    Double_v gamDirY = sinTheta * sinPhi;
    Double_v gamDirZ = cosTheta;
    // rotate gamma direction to the lab frame:
    RotateToLabFrame(gamDirX, gamDirY, gamDirZ, tracks.GetDirXVec(i), tracks.GetDirYVec(i), tracks.GetDirZVec(i));

    LightTrack_v &secondaries = td->fPhysicsData->GetSecondarySOA();
    for (int l = 0; l < kVecLenD; ++l) {
      int idx = secondaries.InsertTrack();
      secondaries.SetKinE(Get(gammaEnergy, l), idx);
      secondaries.SetDirX(Get(gamDirX, l), idx);
      secondaries.SetDirY(Get(gamDirY, l), idx);
      secondaries.SetDirZ(Get(gamDirZ, l), idx);
      secondaries.SetGVcode(fSecondaryInternalCode, idx);
      secondaries.SetMass(0.0, idx);
      secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx);
    }

    // compute the primary e-/e+ post interaction direction: from momentum vector conservation
    const Double_v elInitTotalMomentum = Math::Sqrt(ekin * (ekin + 2.0 * geant::units::kElectronMassC2));
    // final momentum of the e-/e+ in the lab frame
    Double_v elDirX = elInitTotalMomentum * tracks.GetDirXVec(i) - gammaEnergy * gamDirX;
    Double_v elDirY = elInitTotalMomentum * tracks.GetDirYVec(i) - gammaEnergy * gamDirY;
    Double_v elDirZ = elInitTotalMomentum * tracks.GetDirZVec(i) - gammaEnergy * gamDirZ;

    Double_v norm = 1.0 / Math::Sqrt(elDirX * elDirX + elDirY * elDirY + elDirZ * elDirZ);

    // update primary track direction
    for (int l = 0; l < kVecLenD; ++l) {
      tracks.SetDirX(Get(elDirX * norm, l), i + l);
      tracks.SetDirY(Get(elDirY * norm, l), i + l);
      tracks.SetDirZ(Get(elDirZ * norm, l), i + l);
      tracks.SetKinE(Get(ekin - gammaEnergy, l), i + l);
    }
  }
}

Double_v VecRelativisticBremsModel::SampleEnergyTransfer(Double_v gammaCut, Double_v densityCor, IndexD_v mcLocalIdx,
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
  if (!MaskEmpty(mask)) {
    vecCore::MaskedAssign(indxPrimEkin, mask, indxPrimEkin + 1);
  }

  Double_v egammaV;
  for (int l = 0; l < kVecLenD; ++l) {
    assert(indxPrimEkin[l] >= 0);
    assert(indxPrimEkin[l] <= fSamplingTables[mcLocalIdx[l]]->fAliasData.size() - 1);
    //    LinAliasCached& als = fAliasData.fTablesPerMatCut[mcLocalIdx[l]]->fAliasData[indxPrimEkin[l]];
    //    double xi = AliasTableAlternative::SampleLinear(als,fSTNumSamplingElecEnergies,r2[l],r3[l]);

    const LinAlias *als = fSamplingTables[Get(mcLocalIdx, l)]->fAliasData[Get(indxPrimEkin, l)];
    const double egamma =
        fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]), &(als->fAliasIndx[0]),
                                    fSTNumSamplingPhotEnergies, Get(r2, l), Get(r3, l));
    Set(egammaV, l, egamma);
  }

  const Double_v dum1 = gammaCut * gammaCut + densityCor;
  const Double_v dum2 = (primekin * primekin + densityCor) / dum1;

  return Math::Sqrt(dum1 * Math::Exp(egammaV * Math::Log(dum2)) - densityCor);
}

void VecRelativisticBremsModel::SampleEnergyTransfer(const double *eEkin, const double *gammaCut, const double *zetArr,
                                                     const double *densityCorConstArr, const double *lpmEnergyArr,
                                                     double *gammaEn, int N, const geant::TaskData *td)
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
    Double_v eekin            = vecCore::Gather<Double_v>(eEkin, idx);
    Double_v gcut             = vecCore::Gather<Double_v>(gammaCut, idx);
    Double_v densityCorrConst = vecCore::Gather<Double_v>(densityCorConstArr, idx);
    Double_v lpmEnergy        = vecCore::Gather<Double_v>(lpmEnergyArr, idx);
    Double_v zet              = vecCore::Gather<Double_v>(zetArr, idx);

    Double_v etot        = eekin + geant::units::kElectronMassC2;
    Double_v densityCorr = densityCorrConst * etot * etot;

    Double_v energyThLPM = Math::Sqrt(densityCorrConst) * lpmEnergy;
    MaskD_v isLPM        = MaskD_v(fIsUseLPM) && (etot > energyThLPM);

    Double_v minVal   = Math::Log(gcut * gcut + densityCorr);
    Double_v valRange = Math::Log(eekin * eekin + densityCorr) - minVal;
    Double_v funcMax;
    std::array<int, kVecLenD> izetV;
    for (int l = 0; l < kVecLenD; ++l) {
      const int izet = std::min(std::lrint(Get(zet, l)), gMaxZet - 1);
      Set(funcMax, l, gElementData[izet]->fZFactor1 + gElementData[izet]->fZFactor2);
      izetV[l] = izet;
    }

    Double_v egamma =
        Math::Sqrt(Math::Max(Math::Exp(minVal + td->fRndm->uniformV() * valRange) - densityCorr, (Double_v)0.));
    Double_v funcVal;
    if (!MaskEmpty(isLPM)) {
      vecCore::MaskedAssign(funcVal, isLPM, ComputeURelDXSecPerAtom(egamma, etot, lpmEnergy, densityCorr, izetV));
    }
    if (!MaskEmpty(!isLPM)) {
      vecCore::MaskedAssign(funcVal, !isLPM, ComputeDXSecPerAtom(egamma, etot, zet));
    }

    MaskD_v accepted = funcVal > td->fRndm->uniformV() * funcMax;

    if (!MaskEmpty(accepted)) {
      vecCore::Scatter(egamma, gammaEn, idx);
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

Double_v VecRelativisticBremsModel::ComputeURelDXSecPerAtom(Double_v egamma, Double_v etotal, Double_v lpmenergy,
                                                            Double_v densitycor, std::array<int, kVecLenD> izet)
{
  const Double_v y     = egamma / etotal;
  const Double_v onemy = 1. - y;
  const Double_v dum0  = 0.25 * y * y;
  Double_v funcGS, funcPhiS, funcXiS;
  ComputeLPMfunctions(funcXiS, funcGS, funcPhiS, lpmenergy, egamma, etotal, densitycor, izet);
  Double_v dcs;
  Double_v ZFactor1, ZFactor2;
  for (int l = 0; l < kVecLenD; ++l) {
    Set(ZFactor1, l, gElementData[izet[l]]->fZFactor1);
    Set(ZFactor2, l, gElementData[izet[l]]->fZFactor2);
  }
  dcs = funcXiS * (dum0 * funcGS + (onemy + 2.0 * dum0) * funcPhiS) * ZFactor1 + onemy * ZFactor2;
  return Math::Max(dcs, (Double_v)0.0);
}

void VecRelativisticBremsModel::ComputeLPMfunctions(Double_v &funcXiS, Double_v &funcGS, Double_v &funcPhiS,
                                                    const Double_v lpmenergy, const Double_v egamma,
                                                    const Double_v etot, const Double_v densitycor,
                                                    const std::array<int, kVecLenD> izet)
{
  static const Double_v sqrt2 = Math::Sqrt(2.);
  const Double_v redegamma    = egamma / etot;
  // const Double_v varSprime = std::sqrt(0.125*redegamma/(1.0-redegamma)*lpmenergy/etot);
  const Double_v varSprime = Math::Sqrt(0.125 * redegamma * lpmenergy / ((1.0 - redegamma) * etot));
  Double_v varS1, ILVarS1Cond, ILVarS1;
  for (int l = 0; l < kVecLenD; ++l) {
    Set(varS1, l, gElementData[izet[l]]->fVarS1);
    Set(ILVarS1Cond, l, gElementData[izet[l]]->fILVarS1Cond);
    Set(ILVarS1, l, gElementData[izet[l]]->fILVarS1);
  }
  const Double_v condition = sqrt2 * varS1;
  Double_v funcXiSprime    = 2.0;

  MaskD_v tmp1, tmp2;
  tmp1 = varSprime > 1.0;
  tmp2 = varSprime > condition;
  vecCore::MaskedAssign(funcXiSprime, tmp1, (Double_v)1.0);
  const Double_v funcHSprime = Math::Log(varSprime) * ILVarS1Cond;
  Double_v tmpFuncXiSprime =
      1.0 + funcHSprime - 0.08 * (1.0 - funcHSprime) * funcHSprime * (2.0 - funcHSprime) * ILVarS1Cond;
  vecCore::MaskedAssign(funcXiSprime, !tmp1 && tmp2, tmpFuncXiSprime);

  const Double_v varS = varSprime / Math::Sqrt(funcXiSprime);
  // - include dielectric suppression effect into s according to Migdal
  const Double_v varShat = varS * (1.0 + densitycor / (egamma * egamma));
  funcXiS                = 2.0;
  tmp1                   = varShat > 1.0;
  tmp2                   = varShat > varS1;
  vecCore::MaskedAssign(funcXiS, tmp1, (Double_v)1.0);
  vecCore::MaskedAssign(funcXiS, !tmp1 && tmp2, 1.0 + Math::Log(varShat) * ILVarS1);
  GetLPMFunctions(funcGS, funcPhiS, varShat);
  // ComputeLPMGsPhis(funcGS, funcPhiS, varShat);
  //
  // MAKE SURE SUPPRESSION IS SMALLER THAN 1: due to Migdal's approximation on xi
  tmp1 = funcXiS * funcPhiS > 1.0 || varShat > 0.57;
  vecCore::MaskedAssign(funcXiS, tmp1, 1.0 / funcPhiS);
}

void VecRelativisticBremsModel::GetLPMFunctions(Double_v &lpmGs, Double_v &lpmPhis, const Double_v s)
{
  MaskD_v tmp = s < gLPMFuncs.fSLimit;

  Double_v val  = s / gLPMFuncs.fSDelta;
  IndexD_v ilow = IndexD_v(val);
  val -= ilow;
  for (int l = 0; l < kVecLenD; ++l) {
    if (!Get(tmp, l)) Set(ilow, l, gLPMFuncs.fLPMFuncPhi.size() - 2); // above limit
  }
  Double_v lpmG_lowp1   = vecCore::Gather<Double_v>(gLPMFuncs.fLPMFuncG.data(), ilow + (IndexD_v)1);
  Double_v lpmG_low     = vecCore::Gather<Double_v>(gLPMFuncs.fLPMFuncG.data(), ilow);
  Double_v lpmPhi_lowp1 = vecCore::Gather<Double_v>(gLPMFuncs.fLPMFuncPhi.data(), ilow + (IndexD_v)1);
  Double_v lpmPhi_low   = vecCore::Gather<Double_v>(gLPMFuncs.fLPMFuncPhi.data(), ilow);
  lpmGs                 = (lpmG_lowp1 - lpmG_low) * val + lpmG_low;
  lpmPhis               = (lpmPhi_lowp1 - lpmPhi_low) * val + lpmPhi_low;

  if (!MaskEmpty(!tmp)) {
    Double_v ss = s * s;
    ss *= ss;
    vecCore::MaskedAssign(lpmPhis, !tmp, 1.0 - 0.01190476 / ss);
    vecCore::MaskedAssign(lpmGs, !tmp, 1.0 - 0.0230655 / ss);
  }
}

void VecRelativisticBremsModel::SamplePhotonDirection(Double_v elenergy, Double_v &sinTheta, Double_v &cosTheta,
                                                      Double_v rndm)
{
  const Double_v c = 4. - 8. * rndm;
  Double_v a       = c;
  Double_v signc   = 1.;
  MaskD_v tmp      = c < 0.0;
  vecCore::MaskedAssign(signc, tmp, (Double_v)-1.0);
  vecCore::MaskedAssign(a, tmp, -c);

  const Double_v delta = 0.5 * (Math::Sqrt(a * a + 4.) + a);
  //  delta += a;
  //  delta *= 0.5;

  const Double_v cofA = -signc * Math::Exp(Math::Log(delta) / 3.0);
  cosTheta            = cofA - 1. / cofA;

  const Double_v tau  = elenergy * geant::units::kInvElectronMassC2;
  const Double_v beta = Math::Sqrt(tau * (tau + 2.)) / (tau + 1.);

  cosTheta = (cosTheta + beta) / (1. + cosTheta * beta);
  // check cosTheta limit
  cosTheta = Math::Min((Double_v)1.0, cosTheta);
  // if (cosTheta>1.0) {
  //  cosTheta = 1.0;
  //}
  sinTheta = Math::Sqrt((1. - cosTheta) * (1. + cosTheta));
}

Double_v VecRelativisticBremsModel::ComputeDXSecPerAtom(Double_v egamma, Double_v etotal, Double_v zet)
{
  // constexpr Double_v factor =
  // 16.*geant::units::kFineStructConst*geant::units::kClassicElectronRadius*geant::units::kClassicElectronRadius/3.;
  Double_v dcs         = 0.;
  const Double_v y     = egamma / etotal;
  const Double_v onemy = 1. - y;
  Double_v ZFactor1;
  Double_v ZFactor2;
  Double_v Fz;
  Double_v LogZ;
  Double_v GammaFactor;
  Double_v EpsilonFactor;
  IndexD_v izetV;
  for (int l = 0; l < kVecLenD; ++l) {
    int izet      = std::lrint(Get(zet, l));
    ZFactor1      = gElementData[izet]->fZFactor1;
    ZFactor2      = gElementData[izet]->fZFactor2;
    Fz            = gElementData[izet]->fFz;
    LogZ          = gElementData[izet]->fLogZ;
    GammaFactor   = gElementData[izet]->fGammaFactor;
    EpsilonFactor = gElementData[izet]->fEpsilonFactor;
    Set(izetV, l, izet);
  }
  MaskD_v smallZ = izetV < 5;
  if (!MaskEmpty(smallZ)) {
    vecCore::MaskedAssign(dcs, smallZ, (onemy + 0.75 * y * y) * ZFactor1 + onemy * ZFactor2);
  }
  if (!MaskEmpty(!smallZ)) {
    // Tsai: screening from Thomas-Fermi model of atom; Tsai Eq.(3.82)
    // variables gamma and epsilon from Tsai Eq.(3.30) and Eq.(3.31)
    const Double_v invZ    = 1. / zet;
    const Double_v dum0    = y / (etotal - egamma);
    const Double_v gamma   = dum0 * GammaFactor;
    const Double_v epsilon = dum0 * EpsilonFactor;
    Double_v phi1, phi1m2, xsi1, xsi1m2;
    ComputeScreeningFunctions(phi1, phi1m2, xsi1, xsi1m2, gamma, epsilon);
    vecCore::MaskedAssign(dcs, !smallZ,
                          (onemy + 0.75 * y * y) * ((0.25 * phi1 - Fz) + (0.25 * xsi1 - 2. * LogZ / 3.) * invZ) +
                              0.125 * onemy * (phi1m2 + xsi1m2 * invZ));
    // dcs *= factor*zet*zet/egamma;
  }

  return Math::Max(dcs, (Double_v)0.0);
}

void VecRelativisticBremsModel::ComputeScreeningFunctions(Double_v &phi1, Double_v &phi1m2, Double_v &xsi1,
                                                          Double_v &xsi1m2, const Double_v gamma,
                                                          const Double_v epsilon)
{
  const Double_v gamma2 = gamma * gamma;
  phi1 =
      16.863 - 2.0 * Math::Log(1.0 + 0.311877 * gamma2) + 2.4 * Math::Exp(-0.9 * gamma) + 1.6 * Math::Exp(-1.5 * gamma);
  phi1m2                  = 2.0 / (3.0 + 19.5 * gamma + 18.0 * gamma2); // phi1-phi2
  const Double_v epsilon2 = epsilon * epsilon;
  xsi1                    = 24.34 - 2.0 * Math::Log(1.0 + 13.111641 * epsilon2) + 2.8 * Math::Exp(-8.0 * epsilon) +
         1.2 * Math::Exp(-29.2 * epsilon);
  xsi1m2 = 2.0 / (3.0 + 120.0 * epsilon + 1200.0 * epsilon2); // xsi1-xsi2
}

bool VecRelativisticBremsModel::IsModelUsable(const MaterialCuts *matCut, double ekin)
{
  const double gammaCut = matCut->GetProductionCutsInEnergy()[0];
  return ekin < GetHighEnergyUsageLimit() && ekin > GetLowEnergyUsageLimit() && ekin > gammaCut;
}

} // namespace geantphysics
