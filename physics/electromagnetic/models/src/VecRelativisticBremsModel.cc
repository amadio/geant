
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

  for (int i = 0; i < N; i += kPhysDVWidth) {

    PhysDV primEkin = tracks.GetKinEVec(i);
    PhysDV gammaCut;
    PhysDI mcIndxLocal; // Used by alias method
    PhysDV densityCorConst;
    PhysDV lpmEnergy;
    for (int l = 0; l < kPhysDVWidth; ++l) {
      const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(tracks.GetMaterialCutCoupleIndex(i + l));
      gammaCut[l]                = matCut->GetProductionCutsInEnergy()[0];
      densityCorConst[l] =
          matCut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol() * gDensityFactor;
      lpmEnergy[l] = matCut->GetMaterial()->GetMaterialProperties()->GetRadiationLength() * gLPMFactor;
      if (GetUseSamplingTables()) {
        mcIndxLocal[l] = fGlobalMatGCutIndxToLocal[matCut->GetIndex()];
      } else {
        // sample target element
        const Vector_t<Element *> &theElements = matCut->GetMaterial()->GetElementVector();
        double targetElemIndx                  = 0;
        if (theElements.size() > 1) {
          targetElemIndx = SampleTargetElementIndex(matCut, primEkin[l], td->fRndm->uniform());
        }
        Zet[i + l] = theElements[targetElemIndx]->GetZ();
        //        IZet[i + l]      = (int) std::lrint(zet);//std::min((int)std::lrint(zet), fDCSMaxZet);
      }
    }
    PhysDV totalEn    = primEkin + geant::units::kElectronMassC2;
    PhysDV densityCor = densityCorConst * (totalEn * totalEn);
    assert((primEkin >= gammaCut).isFull()); // Cut filtering should be applied up the call chain.
    if (GetUseSamplingTables()) {
      PhysDV r1          = td->fRndm->uniformV();
      PhysDV r2          = td->fRndm->uniformV();
      PhysDV r3          = td->fRndm->uniformV();
      PhysDV primEkin    = tracks.GetKinEVec(i);
      PhysDV gammaEnergy = SampleEnergyTransfer(gammaCut, densityCor, mcIndxLocal, fAliasData.fLogEmin.data(),
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

  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV ekin = tracks.GetKinEVec(i);
    PhysDV gammaEnergy;
    vecCore::Load(gammaEnergy, gammaeEnergyArray + i);

    PhysDV cosTheta = 1.0;
    PhysDV sinTheta = 0.0;
    PhysDV rnd0     = td->fRndm->uniformV();
    SamplePhotonDirection(ekin, sinTheta, cosTheta, rnd0);

    PhysDV rnd1      = td->fRndm->uniformV();
    const PhysDV phi = geant::units::kTwoPi * (rnd1);
    PhysDV sinPhi, cosPhi;
    SinCos(phi, &sinPhi, &cosPhi);
    // gamma direction in the scattering frame
    PhysDV gamDirX = sinTheta * cosPhi;
    PhysDV gamDirY = sinTheta * sinPhi;
    PhysDV gamDirZ = cosTheta;
    // rotate gamma direction to the lab frame:
    RotateToLabFrame(gamDirX, gamDirY, gamDirZ, tracks.GetDirXVec(i), tracks.GetDirYVec(i), tracks.GetDirZVec(i));

    LightTrack_v &secondaries = td->fPhysicsData->GetSecondarySOA();
    for (int l = 0; l < kPhysDVWidth; ++l) {
      int idx = secondaries.InsertTrack();
      secondaries.SetKinE(gammaEnergy[l], idx);
      secondaries.SetDirX(gamDirX[l], idx);
      secondaries.SetDirY(gamDirY[l], idx);
      secondaries.SetDirZ(gamDirZ[l], idx);
      secondaries.SetGVcode(fSecondaryInternalCode, idx);
      secondaries.SetMass(0.0, idx);
      secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx);
    }

    // compute the primary e-/e+ post interaction direction: from momentum vector conservation
    const PhysDV elInitTotalMomentum = Sqrt(ekin * (ekin + 2.0 * geant::units::kElectronMassC2));
    // final momentum of the e-/e+ in the lab frame
    PhysDV elDirX = elInitTotalMomentum * tracks.GetDirXVec(i) - gammaEnergy * gamDirX;
    PhysDV elDirY = elInitTotalMomentum * tracks.GetDirYVec(i) - gammaEnergy * gamDirY;
    PhysDV elDirZ = elInitTotalMomentum * tracks.GetDirZVec(i) - gammaEnergy * gamDirZ;

    PhysDV norm = 1.0 / Sqrt(elDirX * elDirX + elDirY * elDirY + elDirZ * elDirZ);

    // update primary track direction
    for (int l = 0; l < kPhysDVWidth; ++l) {
      tracks.SetDirX((elDirX * norm)[l], i + l);
      tracks.SetDirY((elDirY * norm)[l], i + l);
      tracks.SetDirZ((elDirZ * norm)[l], i + l);
      tracks.SetKinE((ekin - gammaEnergy)[l], i + l);
    }
  }
}

PhysDV VecRelativisticBremsModel::SampleEnergyTransfer(PhysDV gammaCut, PhysDV densityCor, PhysDI mcLocalIdx,
                                                       double *tableEmin, double *tableILDeta, PhysDV primekin,
                                                       PhysDV r1, PhysDV r2, PhysDV r3)
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

  PhysDV egammaV;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    assert(indxPrimEkin[l] >= 0);
    assert(indxPrimEkin[l] <= fSamplingTables[mcLocalIdx[l]]->fAliasData.size() - 1);
    //    LinAliasCached& als = fAliasData.fTablesPerMatCut[mcLocalIdx[l]]->fAliasData[indxPrimEkin[l]];
    //    double xi = AliasTableAlternative::SampleLinear(als,fSTNumSamplingElecEnergies,r2[l],r3[l]);

    const LinAlias *als = fSamplingTables[mcLocalIdx[l]]->fAliasData[indxPrimEkin[l]];
    const double egamma = fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]),
                                                      &(als->fAliasIndx[0]), fSTNumSamplingPhotEnergies, r2[l], r3[l]);
    egammaV[l] = egamma;
  }

  const PhysDV dum1 = gammaCut * gammaCut + densityCor;
  const PhysDV dum2 = (primekin * primekin + densityCor) / dum1;

  return Sqrt(dum1 * Exp(egammaV * Log(dum2)) - densityCor);
}

void VecRelativisticBremsModel::SampleEnergyTransfer(const double *eEkin, const double *gammaCut, const double *zetArr,
                                                     const double *densityCorConstArr, const double *lpmEnergyArr,
                                                     double *gammaEn, int N, const geant::TaskData *td)
{

  // assert(N>=kPhysDVWidth)
  int currN        = 0;
  PhysDM lanesDone = PhysDM::Zero();
  PhysDI idx;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    idx[l] = currN++;
  }

  while (currN < N || !lanesDone.isFull()) {
    PhysDV eekin            = vecCore::Gather<PhysDV>(eEkin, idx);
    PhysDV gcut             = vecCore::Gather<PhysDV>(gammaCut, idx);
    PhysDV densityCorrConst = vecCore::Gather<PhysDV>(densityCorConstArr, idx);
    PhysDV lpmEnergy        = vecCore::Gather<PhysDV>(lpmEnergyArr, idx);
    PhysDV zet              = vecCore::Gather<PhysDV>(zetArr, idx);

    PhysDV etot        = eekin + geant::units::kElectronMassC2;
    PhysDV densityCorr = densityCorrConst * etot * etot;

    PhysDV energyThLPM = Sqrt(densityCorrConst) * lpmEnergy;
    PhysDM isLPM       = PhysDM(fIsUseLPM) && (etot > energyThLPM);

    PhysDV minVal   = Log(gcut * gcut + densityCorr);
    PhysDV valRange = Log(eekin * eekin + densityCorr) - minVal;
    PhysDV funcMax;
    std::array<int, kPhysDVWidth> izetV;
    for (int l = 0; l < kPhysDVWidth; ++l) {
      const int izet = std::min(std::lrint(zet[l]), gMaxZet - 1);
      funcMax[l]     = gElementData[izet]->fZFactor1 + gElementData[izet]->fZFactor2;
      izetV[l]       = izet;
    }

    PhysDV egamma = Sqrt(Max(Exp(minVal + td->fRndm->uniformV() * valRange) - densityCorr, (PhysDV)0.));
    PhysDV funcVal;
    if (isLPM.isNotEmpty()) {
      vecCore::MaskedAssign(funcVal, isLPM, ComputeURelDXSecPerAtom(egamma, etot, lpmEnergy, densityCorr, izetV));
    }
    if ((!isLPM).isNotEmpty()) {
      vecCore::MaskedAssign(funcVal, !isLPM, ComputeDXSecPerAtom(egamma, etot, zet));
    }

    PhysDM accepted = funcVal > td->fRndm->uniformV() * funcMax;

    if (accepted.isNotEmpty()) {
      vecCore::Scatter(egamma, gammaEn, idx);
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

PhysDV VecRelativisticBremsModel::ComputeURelDXSecPerAtom(PhysDV egamma, PhysDV etotal, PhysDV lpmenergy,
                                                          PhysDV densitycor, std::array<int, kPhysDVWidth> izet)
{
  const PhysDV y     = egamma / etotal;
  const PhysDV onemy = 1. - y;
  const PhysDV dum0  = 0.25 * y * y;
  PhysDV funcGS, funcPhiS, funcXiS;
  ComputeLPMfunctions(funcXiS, funcGS, funcPhiS, lpmenergy, egamma, etotal, densitycor, izet);
  PhysDV dcs;
  PhysDV ZFactor1, ZFactor2;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    ZFactor1[l] = gElementData[izet[l]]->fZFactor1;
    ZFactor2[l] = gElementData[izet[l]]->fZFactor2;
  }
  dcs = funcXiS * (dum0 * funcGS + (onemy + 2.0 * dum0) * funcPhiS) * ZFactor1 + onemy * ZFactor2;
  return Max(dcs, (PhysDV)0.0);
}

void VecRelativisticBremsModel::ComputeLPMfunctions(PhysDV &funcXiS, PhysDV &funcGS, PhysDV &funcPhiS,
                                                    const PhysDV lpmenergy, const PhysDV egamma, const PhysDV etot,
                                                    const PhysDV densitycor, const std::array<int, kPhysDVWidth> izet)
{
  static const PhysDV sqrt2 = Sqrt(2.);
  const PhysDV redegamma    = egamma / etot;
  // const PhysDV varSprime = std::sqrt(0.125*redegamma/(1.0-redegamma)*lpmenergy/etot);
  const PhysDV varSprime = Sqrt(0.125 * redegamma * lpmenergy / ((1.0 - redegamma) * etot));
  PhysDV varS1, ILVarS1Cond, ILVarS1;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    varS1[l]       = gElementData[izet[l]]->fVarS1;
    ILVarS1Cond[l] = gElementData[izet[l]]->fILVarS1Cond;
    ILVarS1[l]     = gElementData[izet[l]]->fILVarS1;
  }
  const PhysDV condition = sqrt2 * varS1;
  PhysDV funcXiSprime    = 2.0;

  PhysDM tmp1, tmp2;
  tmp1 = varSprime > 1.0;
  tmp2 = varSprime > condition;
  vecCore::MaskedAssign(funcXiSprime, tmp1, (PhysDV)1.0);
  const PhysDV funcHSprime = Log(varSprime) * ILVarS1Cond;
  PhysDV tmpFuncXiSprime =
      1.0 + funcHSprime - 0.08 * (1.0 - funcHSprime) * funcHSprime * (2.0 - funcHSprime) * ILVarS1Cond;
  vecCore::MaskedAssign(funcXiSprime, !tmp1 && tmp2, tmpFuncXiSprime);

  const PhysDV varS = varSprime / Sqrt(funcXiSprime);
  // - include dielectric suppression effect into s according to Migdal
  const PhysDV varShat = varS * (1.0 + densitycor / (egamma * egamma));
  funcXiS              = 2.0;
  tmp1                 = varShat > 1.0;
  tmp2                 = varShat > varS1;
  vecCore::MaskedAssign(funcXiS, tmp1, (PhysDV)1.0);
  vecCore::MaskedAssign(funcXiS, !tmp1 && tmp2, 1.0 + Log(varShat) * ILVarS1);
  GetLPMFunctions(funcGS, funcPhiS, varShat);
  // ComputeLPMGsPhis(funcGS, funcPhiS, varShat);
  //
  // MAKE SURE SUPPRESSION IS SMALLER THAN 1: due to Migdal's approximation on xi
  tmp1 = funcXiS * funcPhiS > 1.0 || varShat > 0.57;
  vecCore::MaskedAssign(funcXiS, tmp1, 1.0 / funcPhiS);
}

void VecRelativisticBremsModel::GetLPMFunctions(PhysDV &lpmGs, PhysDV &lpmPhis, const PhysDV s)
{
  PhysDM tmp = s < gLPMFuncs.fSLimit;

  PhysDV val  = s / gLPMFuncs.fSDelta;
  PhysDI ilow = PhysDI(val);
  val -= ilow;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    if (!tmp[l]) ilow[l] = gLPMFuncs.fLPMFuncPhi.size() - 2; // above limit
  }
  PhysDV lpmG_lowp1   = vecCore::Gather<PhysDV>(gLPMFuncs.fLPMFuncG.data(), ilow + (PhysDI)1);
  PhysDV lpmG_low     = vecCore::Gather<PhysDV>(gLPMFuncs.fLPMFuncG.data(), ilow);
  PhysDV lpmPhi_lowp1 = vecCore::Gather<PhysDV>(gLPMFuncs.fLPMFuncPhi.data(), ilow + (PhysDI)1);
  PhysDV lpmPhi_low   = vecCore::Gather<PhysDV>(gLPMFuncs.fLPMFuncPhi.data(), ilow);
  lpmGs               = (lpmG_lowp1 - lpmG_low) * val + lpmG_low;
  lpmPhis             = (lpmPhi_lowp1 - lpmPhi_low) * val + lpmPhi_low;

  if ((!tmp).isNotEmpty()) {
    PhysDV ss = s * s;
    ss *= ss;
    vecCore::MaskedAssign(lpmPhis, !tmp, 1.0 - 0.01190476 / ss);
    vecCore::MaskedAssign(lpmGs, !tmp, 1.0 - 0.0230655 / ss);
  }
}

void VecRelativisticBremsModel::SamplePhotonDirection(PhysDV elenergy, PhysDV &sinTheta, PhysDV &cosTheta, PhysDV rndm)
{
  const PhysDV c = 4. - 8. * rndm;
  PhysDV a       = c;
  PhysDV signc   = 1.;
  PhysDM tmp     = c < 0.0;
  vecCore::MaskedAssign(signc, tmp, (PhysDV)-1.0);
  vecCore::MaskedAssign(a, tmp, -c);

  const PhysDV delta = 0.5 * (Sqrt(a * a + 4.) + a);
  //  delta += a;
  //  delta *= 0.5;

  const PhysDV cofA = -signc * Exp(Log(delta) / 3.0);
  cosTheta          = cofA - 1. / cofA;

  const PhysDV tau  = elenergy * geant::units::kInvElectronMassC2;
  const PhysDV beta = Sqrt(tau * (tau + 2.)) / (tau + 1.);

  cosTheta = (cosTheta + beta) / (1. + cosTheta * beta);
  // check cosTheta limit
  cosTheta = Min((PhysDV)1.0, cosTheta);
  // if (cosTheta>1.0) {
  //  cosTheta = 1.0;
  //}
  sinTheta = Sqrt((1. - cosTheta) * (1. + cosTheta));
}

PhysDV VecRelativisticBremsModel::ComputeDXSecPerAtom(PhysDV egamma, PhysDV etotal, PhysDV zet)
{
  // constexpr PhysDV factor =
  // 16.*geant::units::kFineStructConst*geant::units::kClassicElectronRadius*geant::units::kClassicElectronRadius/3.;
  PhysDV dcs         = 0.;
  const PhysDV y     = egamma / etotal;
  const PhysDV onemy = 1. - y;
  PhysDV ZFactor1;
  PhysDV ZFactor2;
  PhysDV Fz;
  PhysDV LogZ;
  PhysDV GammaFactor;
  PhysDV EpsilonFactor;
  PhysDI izetV;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    int izet      = std::lrint(zet[l]);
    ZFactor1      = gElementData[izet]->fZFactor1;
    ZFactor2      = gElementData[izet]->fZFactor2;
    Fz            = gElementData[izet]->fFz;
    LogZ          = gElementData[izet]->fLogZ;
    GammaFactor   = gElementData[izet]->fGammaFactor;
    EpsilonFactor = gElementData[izet]->fEpsilonFactor;
    izetV[l]      = izet;
  }
  PhysDM smallZ = izetV < 5;
  if (smallZ.isNotEmpty()) {
    vecCore::MaskedAssign(dcs, smallZ, (onemy + 0.75 * y * y) * ZFactor1 + onemy * ZFactor2);
  }
  if ((!smallZ).isNotEmpty()) {
    // Tsai: screening from Thomas-Fermi model of atom; Tsai Eq.(3.82)
    // variables gamma and epsilon from Tsai Eq.(3.30) and Eq.(3.31)
    const PhysDV invZ    = 1. / zet;
    const PhysDV dum0    = y / (etotal - egamma);
    const PhysDV gamma   = dum0 * GammaFactor;
    const PhysDV epsilon = dum0 * EpsilonFactor;
    PhysDV phi1, phi1m2, xsi1, xsi1m2;
    ComputeScreeningFunctions(phi1, phi1m2, xsi1, xsi1m2, gamma, epsilon);
    vecCore::MaskedAssign(dcs, !smallZ,
                          (onemy + 0.75 * y * y) * ((0.25 * phi1 - Fz) + (0.25 * xsi1 - 2. * LogZ / 3.) * invZ) +
                              0.125 * onemy * (phi1m2 + xsi1m2 * invZ));
    // dcs *= factor*zet*zet/egamma;
  }

  return Max(dcs, (PhysDV)0.0);
}

void VecRelativisticBremsModel::ComputeScreeningFunctions(PhysDV &phi1, PhysDV &phi1m2, PhysDV &xsi1, PhysDV &xsi1m2,
                                                          const PhysDV gamma, const PhysDV epsilon)
{
  const PhysDV gamma2 = gamma * gamma;
  phi1                = 16.863 - 2.0 * Log(1.0 + 0.311877 * gamma2) + 2.4 * Exp(-0.9 * gamma) + 1.6 * Exp(-1.5 * gamma);
  phi1m2              = 2.0 / (3.0 + 19.5 * gamma + 18.0 * gamma2); // phi1-phi2
  const PhysDV epsilon2 = epsilon * epsilon;
  xsi1   = 24.34 - 2.0 * Log(1.0 + 13.111641 * epsilon2) + 2.8 * Exp(-8.0 * epsilon) + 1.2 * Exp(-29.2 * epsilon);
  xsi1m2 = 2.0 / (3.0 + 120.0 * epsilon + 1200.0 * epsilon2); // xsi1-xsi2
}

bool VecRelativisticBremsModel::IsModelUsable(const MaterialCuts *matCut, double ekin)
{
  const double gammaCut = matCut->GetProductionCutsInEnergy()[0];
  return ekin < GetHighEnergyUsageLimit() && ekin > GetLowEnergyUsageLimit() && ekin > gammaCut;
}

} // namespace geantphysics
