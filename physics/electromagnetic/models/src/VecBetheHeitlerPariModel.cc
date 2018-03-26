#include <Geant/SystemOfUnits.h>
#include <Geant/PhysicalConstants.h>
#include <Geant/TaskData.h>
#include <Geant/PhysicsData.h>
#include <Geant/MaterialCuts.h>
#include <TError.h>
#include "Geant/VecBetheHeitlerPairModel.h"

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
}

void VecBetheHeitlerPairModel::SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td)
{
  int N     = tracks.GetNtracks();
  int *IZet = td->fPhysicsData->fPhysicsScratchpad.fIzet;

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

  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV ekin;
    vecCore::Load(ekin, tracks.GetKinEVec(i));
    PhysDV eps0 = geant::units::kElectronMassC2 / ekin;

    PhysDV eps;
    PhysDM smallE = ekin < fGammaEneregyLimit;
    if (smallE.isNotEmpty()) {
      vecCore::MaskedAssign(eps, smallE, eps0 + (0.5 - eps0) * td->fRndm->uniformV());
    }

    if ((!smallE).isNotEmpty()) {
      if (GetUseSamplingTables()) {
        PhysDV r1      = td->fRndm->uniformV();
        PhysDV r2      = td->fRndm->uniformV();
        PhysDV r3      = td->fRndm->uniformV();
        PhysDV tmpEkin = ekin;
        vecCore::MaskedAssign(
            tmpEkin, smallE,
            (PhysDV)fMinimumPrimaryEnergy *
                (1.1)); // To be sure that energy that is too small for alias table will not slip in alias sampling
        PhysDV sampledEps = SampleTotalEnergyTransferAliasOneShot(tmpEkin, &IZet[i], r1, r2, r3);
        vecCore::MaskedAssign(eps, !smallE, sampledEps);
      } else {
        Error("BetheHeitlerVec", "RejSampling not implemented");
      }
    }

    PhysDV electronTotE, positronTotE;
    PhysDM tmpM = td->fRndm->uniformV() > 0.5;
    vecCore::MaskedAssign(electronTotE, tmpM, (1. - eps) * ekin);
    vecCore::MaskedAssign(positronTotE, !tmpM, (1. - eps) * ekin);
    vecCore::MaskedAssign(electronTotE, !tmpM, eps * ekin);
    vecCore::MaskedAssign(positronTotE, tmpM, eps * ekin);

    PhysDV r0 = td->fRndm->uniformV();
    PhysDV r1 = td->fRndm->uniformV();
    PhysDV r2 = td->fRndm->uniformV();

    PhysDV uvar  = -vecCore::math::Log(r0 * r1);
    PhysDM tmpM2 = 9. > 36. * r2;
    vecCore::MaskedAssign(uvar, tmpM2, uvar * 1.6);
    vecCore::MaskedAssign(uvar, !tmpM2, uvar * 0.53333);

    PhysDV r3                  = td->fRndm->uniformV();
    const PhysDV thetaElectron = uvar * geant::units::kElectronMassC2 / electronTotE;
    PhysDV sintEle, costEle;
    vecCore::math::SinCos(thetaElectron, &sintEle, &costEle);
    const PhysDV thetaPositron = uvar * geant::units::kElectronMassC2 / positronTotE;
    PhysDV sintPos, costPos;
    vecCore::math::SinCos(thetaPositron, &sintPos, &costPos);
    const PhysDV phi = geant::units::kTwoPi * td->fRndm->uniformV();
    PhysDV sinphi, cosphi;
    vecCore::math::SinCos(phi, &sinphi, &cosphi);

    PhysDV eleDirX = sintEle * cosphi;
    PhysDV eleDirY = sintEle * sinphi;
    PhysDV eleDirZ = costEle;

    PhysDV posDirX = sintPos * cosphi;
    PhysDV posDirY = sintPos * sinphi;
    PhysDV posDirZ = costPos;

    for (int l = 0; l < kPhysDVWidth; ++l) {
      tracks.SetEnergyDeposit(0.0, i + l);
      tracks.SetTrackStatus(LTrackStatus::kKill, i + l);
    }

    const PhysDV ekinElectron = vecCore::math::Max((electronTotE - geant::units::kElectronMassC2), (PhysDV)0.);
    const PhysDV ekinPositron = vecCore::math::Max((positronTotE - geant::units::kElectronMassC2), (PhysDV)0.);
    // 5. rotate direction back to the lab frame: current directions are relative to the photon dir as z-dir
    PhysDV gammaX, gammaY, gammaZ;
    vecCore::Load(gammaX, tracks.GetDirXV(i));
    vecCore::Load(gammaY, tracks.GetDirYV(i));
    vecCore::Load(gammaZ, tracks.GetDirZV(i));
    RotateToLabFrame(eleDirX, eleDirY, eleDirZ, gammaX, gammaY, gammaZ);
    RotateToLabFrame(posDirX, posDirY, posDirZ, gammaX, gammaY, gammaZ);

    for (int l = 0; l < kPhysDVWidth; ++l) {
      auto &secondarySoA = td->fPhysicsData->GetSecondarySOA();

      int idx = secondarySoA.InsertTrack();
      secondarySoA.SetDirX(vecCore::Get(eleDirX, l), idx);
      secondarySoA.SetDirY(vecCore::Get(eleDirY, l), idx);
      secondarySoA.SetDirZ(vecCore::Get(eleDirZ, l), idx);
      secondarySoA.SetKinE(vecCore::Get(ekinElectron, l), idx);
      secondarySoA.SetGVcode(fElectronInternalCode, idx);
      secondarySoA.SetMass(geant::units::kElectronMassC2, idx);
      secondarySoA.SetTrackIndex(tracks.GetTrackIndex(i + l), idx); // parent Track index
      // then set the e+
      idx = secondarySoA.InsertTrack();
      secondarySoA.SetDirX(vecCore::Get(posDirX, l), idx);
      secondarySoA.SetDirY(vecCore::Get(posDirY, l), idx);
      secondarySoA.SetDirZ(vecCore::Get(posDirZ, l), idx);
      secondarySoA.SetKinE(vecCore::Get(ekinPositron, l), idx);
      secondarySoA.SetGVcode(fPositronInternalCode, idx);
      secondarySoA.SetMass(geant::units::kElectronMassC2, idx);
      secondarySoA.SetTrackIndex(tracks.GetTrackIndex(i + l), idx); // parent Track index
    }
  }
}

void VecBetheHeitlerPairModel::SampleTotalEnergyTransferAliasVec(const double *egamma, const int *izet,
                                                                 const double *r1, const double *r2, const double *r3,
                                                                 int N, double *eps)
{

  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV egammaV;
    vecCore::Load(egammaV, &egamma[i]);
    const PhysDV legamma = vecCore::math::Log(egammaV);

    PhysDV val        = (legamma - fSTLogMinPhotonEnergy) * fSTILDeltaPhotonEnergy;
    PhysDI indxEgamma = (PhysDI)val; // lower electron energy bin index
    PhysDV pIndxHigh  = val - indxEgamma;
    PhysDV R1;
    vecCore::Load(R1, &r1[i]);
    PhysDM mask = R1 < pIndxHigh;
    if (!mask.isEmpty()) {
      vecCore::MaskedAssign(indxEgamma, mask, indxEgamma + 1);
    }

    PhysDV xiV;
    for (int l = 0; l < kPhysDVWidth; ++l) {
      int idx                  = (int)vecCore::Get(indxEgamma, l);
      RatinAliasDataTrans &als = fAliasTablesPerZ[izet[i + l]].fTablePerEn[idx];
      double xi = AliasTableAlternative::SampleRatin(als, fSTNumDiscreteEnergyTransferVals, r2[i + l], r3[i + l], 0);

      vecCore::Set(xiV, l, xi);
    }
    PhysDV deltaMax;
    PhysDV deltaFactor;
    for (int l = 0; l < kPhysDVWidth; ++l) {
      double deltaMaxScalar = egamma[i + l] > 50. * geant::units::MeV ? gElementData[izet[i + l]]->fDeltaMaxHighTsai
                                                                      : gElementData[izet[i + l]]->fDeltaMaxLowTsai;
      double deltaFactorScalar = gElementData[izet[i + l]]->fDeltaFactor;
      vecCore::Set(deltaMax, l, deltaMaxScalar);
      vecCore::Set(deltaFactor, l, deltaFactorScalar);
    }
    const PhysDV eps0   = geant::units::kElectronMassC2 / egammaV;
    const PhysDV epsp   = 0.5 - 0.5 * vecCore::math::Sqrt(1. - 4. * eps0 * deltaFactor / deltaMax);
    const PhysDV epsMin = vecCore::math::Max(eps0, epsp);
    const PhysDV epsV   = epsMin * vecCore::math::Exp(xiV * vecCore::math::Log(0.5 / epsMin));
    vecCore::Store(epsV, eps);
  }
}

PhysDV VecBetheHeitlerPairModel::SampleTotalEnergyTransferAliasOneShot(const PhysDV egamma, const int *izet,
                                                                       const PhysDV r1, const PhysDV r2,
                                                                       const PhysDV r3)
{
  const PhysDV legamma = vecCore::math::Log(egamma);

  PhysDV val        = (legamma - fSTLogMinPhotonEnergy) * fSTILDeltaPhotonEnergy;
  PhysDI indxEgamma = (PhysDI)val; // lower electron energy bin index
  PhysDV pIndxHigh  = val - indxEgamma;
  PhysDM mask       = r1 < pIndxHigh;
  if (!mask.isEmpty()) {
    vecCore::MaskedAssign(indxEgamma, mask, indxEgamma + 1);
  }

  PhysDV xiV;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    int idx                  = (int)vecCore::Get(indxEgamma, l);
    RatinAliasDataTrans &als = fAliasTablesPerZ[izet[l]].fTablePerEn[idx];
    double r2l               = vecCore::Get(r2, l);
    double r3l               = vecCore::Get(r3, l);
    double xi                = AliasTableAlternative::SampleRatin(als, fSTNumDiscreteEnergyTransferVals, r2l, r3l, 0);

    vecCore::Set(xiV, l, xi);
  }
  PhysDV deltaMax;
  PhysDV deltaFactor;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    double deltaMaxScalar = egamma[l] > 50. * geant::units::MeV ? gElementData[izet[l]]->fDeltaMaxHighTsai
                                                                : gElementData[izet[l]]->fDeltaMaxLowTsai;
    double deltaFactorScalar = gElementData[izet[l]]->fDeltaFactor;
    vecCore::Set(deltaMax, l, deltaMaxScalar);
    vecCore::Set(deltaFactor, l, deltaFactorScalar);
  }
  const PhysDV eps0   = geant::units::kElectronMassC2 / egamma;
  const PhysDV epsp   = 0.5 - 0.5 * vecCore::math::Sqrt(1. - 4. * eps0 * deltaFactor / deltaMax);
  const PhysDV epsMin = vecCore::math::Max(eps0, epsp);
  const PhysDV epsV   = epsMin * vecCore::math::Exp(xiV * vecCore::math::Log(0.5 / epsMin));
  return epsV;
}

void VecBetheHeitlerPairModel::SampleTotalEnergyTransferRejVec(const double *egamma, const int *izet,
                                                               geant::TaskData *td)
{
}
}
