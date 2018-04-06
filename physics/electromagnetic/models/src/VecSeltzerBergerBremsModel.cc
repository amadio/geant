
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

  for (int i = 0; i < N; i += kPhysDVWidth) {

    PhysDV primEkin = tracks.GetKinEVec(i);
    PhysDV gammaCut;
    PhysDI mcIndxLocal; // Used by alias method
    PhysDV densityCor;
    for (int l = 0; l < kPhysDVWidth; ++l) {
      const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(tracks.GetMaterialCutCoupleIndex(i + l));
      gammaCut[l]                = matCut->GetProductionCutsInEnergy()[0];
      densityCor[l]              = matCut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol();
      if (GetUseSamplingTables()) {
        mcIndxLocal[l] = fGlobalMatGCutIndxToLocal[matCut->GetIndex()];
      } else {
        // sample target element
        const Vector_t<Element *> &theElements = matCut->GetMaterial()->GetElementVector();
        double targetElemIndx                  = 0;
        if (theElements.size() > 1) {
          targetElemIndx = SampleTargetElementIndex(matCut, primEkin[l], td->fRndm->uniform());
        }
        const double zet = theElements[targetElemIndx]->GetZ();
        IZet[i + l]      = std::min((int)std::lrint(zet), fDCSMaxZet);
        Zet[i + l]       = zet;
      }
    }
    PhysDV totalEn = primEkin + geant::units::kElectronMassC2;
    densityCor *= gMigdalConst * (totalEn * totalEn);
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

PhysDV VecSeltzerBergerBremsModel::SampleEnergyTransfer(PhysDV gammaCut, PhysDV densityCor, PhysDI mcLocalIdx,
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

void VecSeltzerBergerBremsModel::SampleEnergyTransfer(const double *eEkin, const double *gammaCut, const int *IZet,
                                                      const double *zet, double *gammaEn, const double *densityCorArr,
                                                      int N, const geant::TaskData *td)
{

  //  // assert(N>=kPhysDVWidth)
  int currN        = 0;
  PhysDM lanesDone = PhysDM::Zero();
  PhysDI idx;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    idx[l] = currN++;
  }
  //
  while (currN < N || !lanesDone.isFull()) {
    PhysDV eekin       = vecCore::Gather<PhysDV>(eEkin, idx);
    PhysDV gcut        = vecCore::Gather<PhysDV>(gammaCut, idx);
    PhysDV densityCorr = vecCore::Gather<PhysDV>(densityCorArr, idx);

    const PhysDV kappac = gcut / eekin;

    PhysDV lekin  = Log(eekin);
    PhysDV eresid = (lekin - fLogLoadDCSMinElecEnergy) * fInvLogLoadDCSDeltaEnergy;
    PhysDI ie     = (PhysDI)eresid; // y1 index
    eresid -= (PhysDV)ie;           // (y2-y1)*resid + y1
    assert((eekin < fLoadDCSElectronEnergyGrid[0]).isEmpty());
    assert((eekin >= fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies - 1]).isEmpty());

    PhysDV vmax;
    for (int l = 0; l < kPhysDVWidth; ++l) {
      vmax[l] = 1.02 * GetDXSECValue(IZet[idx[l]], ie[l], eresid[l], kappac[l]);
      //
      // majoranta corrected vmax for e-
      constexpr double epeaklimit = 300.0 * geant::units::MeV;
      constexpr double elowlimit  = 20.0 * geant::units::keV;
      if (fIsElectron && kappac[l] < 0.97 && ((eekin[l] > epeaklimit) || (eekin[l] < elowlimit))) {
        vmax = std::max(
            vmax, std::min(fXsecLimits[IZet[idx[l]] - 1], 1.1 * GetDXSECValue(IZet[idx[l]], ie[l], eresid[l], 0.97)));
      }
    }
    PhysDM tmp = kappac < 0.05;
    vecCore::MaskedAssign(vmax, tmp, vmax * 1.2);

    const PhysDV minXi = Log(gcut * gcut + densityCorr);
    const PhysDV delXi = Log(eekin * eekin + densityCorr) - minXi;

    PhysDV rnd0 = td->fRndm->uniformV();

    PhysDV egamma      = Sqrt(Max(Exp(minXi + rnd0 * delXi) - densityCorr, (PhysDV)0.));
    const PhysDV kappa = egamma / eekin;
    PhysDV val         = 0.0;
    for (int l = 0; l < kPhysDVWidth; ++l) {
      val[l] = GetDXSECValue(IZet[idx[l]], ie[l], eresid[l], kappa[l]);
    }

    if (!fIsElectron) {
      PhysDV zetV = vecCore::Gather<PhysDV>(zet, idx);
      val *= PositronCorrection1(eekin, kappa, gcut, zetV);
    }

    PhysDV rnd1     = td->fRndm->uniformV();
    PhysDM accepted = val > vmax * rnd1;

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

void VecSeltzerBergerBremsModel::SamplePhotonDirection(PhysDV elenergy, PhysDV &sinTheta, PhysDV &cosTheta, PhysDV rndm)
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

PhysDV VecSeltzerBergerBremsModel::PositronCorrection1(PhysDV ekinelectron, PhysDV ephoton, PhysDV gcutener, PhysDV z)
{
  const PhysDV dum1 = geant::units::kTwoPi * geant::units::kFineStructConst;
  PhysDV poscor     = 0.0;
  PhysDV e1         = ekinelectron - gcutener; // here is the dif.
  PhysDV ibeta1     = (e1 + geant::units::kElectronMassC2) / Sqrt(e1 * (e1 + 2.0 * geant::units::kElectronMassC2));
  PhysDV e2         = ekinelectron * (1.0 - ephoton);
  PhysDV ibeta2     = (e2 + geant::units::kElectronMassC2) / Sqrt(e2 * (e2 + 2.0 * geant::units::kElectronMassC2));
  PhysDV ddum       = dum1 * z * (ibeta1 - ibeta2);
  PhysDM tmp        = ddum < -12.0;
  vecCore::MaskedAssign(poscor, tmp, (PhysDV)-12.0);
  vecCore::MaskedAssign(poscor, !tmp, Exp(ddum));
  return poscor;
}

bool VecSeltzerBergerBremsModel::IsModelUsable(const MaterialCuts *matCut, double ekin)
{
  const double gammaCut = matCut->GetProductionCutsInEnergy()[0];
  return ekin < GetHighEnergyUsageLimit() && ekin > GetLowEnergyUsageLimit() && ekin > gammaCut;
}

} // namespace geantphysics
