
#include "Geant/VecSeltzerBergerBremsModel.h"

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

VecSeltzerBergerBremsModel::VecSeltzerBergerBremsModel(bool iselectron, const std::string &modelname)
    : SeltzerBergerBremsModel(iselectron, modelname)
{
}

void VecSeltzerBergerBremsModel::Initialize()
{
  SeltzerBergerBremsModel::Initialize();

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

void VecSeltzerBergerBremsModel::SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td)
{
  const int N               = tracks.GetNtracks();
  double *gammaeEnergyArray = td->fPhysicsData->fPhysicsScratchpad.fEps;
  double *gammaCutArr       = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr;
  int *IZet                 = td->fPhysicsData->fPhysicsScratchpad.fIzet;
  double *Zet               = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr2;
  double *densityCorrArr    = td->fPhysicsData->fPhysicsScratchpad.fR0;

  for (int i = 0; i < N; i += kVecLenD) {

    Double_v primEkin = tracks.GetKinEVec(i);
    Double_v gammaCut;
    IndexD_v mcIndxLocal; // Used by alias method
    Double_v densityCor;
    for (int l = 0; l < kVecLenD; ++l) {
      const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(tracks.GetMaterialCutCoupleIndex(i + l));
      Set(gammaCut, l, matCut->GetProductionCutsInEnergy()[0]);
      Set(densityCor, l, matCut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol());
      if (GetUseSamplingTables()) {
        Set(mcIndxLocal, l, fGlobalMatGCutIndxToLocal[matCut->GetIndex()]);
      } else {
        // sample target element
        const Vector_t<Element *> &theElements = matCut->GetMaterial()->GetElementVector();
        double targetElemIndx                  = 0;
        if (theElements.size() > 1) {
          targetElemIndx = SampleTargetElementIndex(matCut, Get(primEkin, l), td->fRndm->uniform());
        }
        const double zet = theElements[targetElemIndx]->GetZ();
        IZet[i + l]      = std::min((int)std::lrint(zet), fDCSMaxZet);
        Zet[i + l]       = zet;
      }
    }
    Double_v totalEn = primEkin + geant::units::kElectronMassC2;
    densityCor *= gMigdalConst * (totalEn * totalEn);
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
      vecCore::Store(densityCor, densityCorrArr + i);
    }
  }

  if (!GetUseSamplingTables()) {
    tracks.GetKinEArr()[N] = tracks.GetKinEArr()[N - 1];
    gammaCutArr[N]         = gammaCutArr[N - 1];
    IZet[N]                = IZet[N - 1];
    Zet[N]                 = Zet[N - 1];
    densityCorrArr[N]      = densityCorrArr[N - 1];

    SampleEnergyTransfer(tracks.GetKinEArr(), gammaCutArr, IZet, Zet, gammaeEnergyArray, densityCorrArr, N, td);
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
    Math::RotateToLabFrame(gamDirX, gamDirY, gamDirZ, tracks.GetDirXVec(i), tracks.GetDirYVec(i), tracks.GetDirZVec(i));

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

Double_v VecSeltzerBergerBremsModel::SampleEnergyTransfer(Double_v gammaCut, Double_v densityCor, IndexD_v mcLocalIdx,
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

void VecSeltzerBergerBremsModel::SampleEnergyTransfer(const double *eEkin, const double *gammaCut, const int *IZet,
                                                      const double *zet, double *gammaEn, const double *densityCorArr,
                                                      int N, const geant::TaskData *td)
{

  //  // assert(N>=kVecLenD)
  int currN = 0;
  MaskD_v lanesDone;
  IndexD_v idx;
  for (int l = 0; l < kVecLenD; ++l) {
    AssignMaskLane(lanesDone, l, false);
    Set(idx, l, currN++);
  }

  //
  while (currN < N || !MaskFull(lanesDone)) {
    Double_v eekin       = vecCore::Gather<Double_v>(eEkin, idx);
    Double_v gcut        = vecCore::Gather<Double_v>(gammaCut, idx);
    Double_v densityCorr = vecCore::Gather<Double_v>(densityCorArr, idx);

    const Double_v kappac = gcut / eekin;

    Double_v lekin  = Math::Log(eekin);
    Double_v eresid = (lekin - fLogLoadDCSMinElecEnergy) * fInvLogLoadDCSDeltaEnergy;
    IndexD_v ie     = (IndexD_v)eresid; // y1 index
    eresid -= (Double_v)ie;             // (y2-y1)*resid + y1
    assert((eekin < fLoadDCSElectronEnergyGrid[0]).isEmpty());
    assert((eekin >= fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies - 1]).isEmpty());

    Double_v vmax;
    for (int l = 0; l < kVecLenD; ++l) {
      Set(vmax, l, 1.02 * GetDXSECValue(IZet[Get(idx, l)], Get(ie, l), Get(eresid, l), Get(kappac, l)));
      //
      // majoranta corrected vmax for e-
      constexpr double epeaklimit = 300.0 * geant::units::MeV;
      constexpr double elowlimit  = 20.0 * geant::units::keV;
      if (fIsElectron && Get(kappac, l) < 0.97 && ((Get(eekin, l) > epeaklimit) || (Get(eekin, l) < elowlimit))) {
        vmax = std::max(vmax, std::min(fXsecLimits[IZet[Get(idx, l)] - 1],
                                       1.1 * GetDXSECValue(IZet[Get(idx, l)], Get(ie, l), Get(eresid, l), 0.97)));
      }
    }
    MaskD_v tmp = kappac < 0.05;
    vecCore::MaskedAssign(vmax, tmp, vmax * 1.2);

    const Double_v minXi = Math::Log(gcut * gcut + densityCorr);
    const Double_v delXi = Math::Log(eekin * eekin + densityCorr) - minXi;

    Double_v rnd0 = td->fRndm->uniformV();

    Double_v egamma      = Math::Sqrt(Math::Max(Math::Exp(minXi + rnd0 * delXi) - densityCorr, (Double_v)0.));
    const Double_v kappa = egamma / eekin;
    Double_v val         = 0.0;
    for (int l = 0; l < kVecLenD; ++l) {
      Set(val, l, GetDXSECValue(IZet[Get(idx, l)], Get(ie, l), Get(eresid, l), Get(kappa, l)));
    }

    if (!fIsElectron) {
      Double_v zetV = vecCore::Gather<Double_v>(zet, idx);
      val *= PositronCorrection1(eekin, kappa, gcut, zetV);
    }

    Double_v rnd1    = td->fRndm->uniformV();
    MaskD_v accepted = val > vmax * rnd1;

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

void VecSeltzerBergerBremsModel::SamplePhotonDirection(Double_v elenergy, Double_v &sinTheta, Double_v &cosTheta,
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

Double_v VecSeltzerBergerBremsModel::PositronCorrection1(Double_v ekinelectron, Double_v ephoton, Double_v gcutener,
                                                         Double_v z)
{
  const Double_v dum1 = geant::units::kTwoPi * geant::units::kFineStructConst;
  Double_v poscor     = 0.0;
  Double_v e1         = ekinelectron - gcutener; // here is the dif.
  Double_v ibeta1 = (e1 + geant::units::kElectronMassC2) / Math::Sqrt(e1 * (e1 + 2.0 * geant::units::kElectronMassC2));
  Double_v e2     = ekinelectron * (1.0 - ephoton);
  Double_v ibeta2 = (e2 + geant::units::kElectronMassC2) / Math::Sqrt(e2 * (e2 + 2.0 * geant::units::kElectronMassC2));
  Double_v ddum   = dum1 * z * (ibeta1 - ibeta2);
  MaskD_v tmp     = ddum < -12.0;
  vecCore::MaskedAssign(poscor, tmp, (Double_v)-12.0);
  vecCore::MaskedAssign(poscor, !tmp, Math::Exp(ddum));
  return poscor;
}

bool VecSeltzerBergerBremsModel::IsModelUsable(const MaterialCuts *matCut, double ekin)
{
  const double gammaCut = matCut->GetProductionCutsInEnergy()[0];
  return ekin < GetHighEnergyUsageLimit() && ekin > GetLowEnergyUsageLimit() && ekin > gammaCut;
}

} // namespace geantphysics
