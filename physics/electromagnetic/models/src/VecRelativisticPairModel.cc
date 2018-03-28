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
  if (GetUseSamplingTables() && fIsUseLPM) {
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
    for (int i = 0; i < N; i += kPhysDVWidth) {
      PhysDV ekin;
      vecCore::Load(ekin, tracks.GetKinEVec(i));
      PhysDV r1 = td->fRndm->uniformV();
      PhysDV r2 = td->fRndm->uniformV();
      PhysDV r3 = td->fRndm->uniformV();

      PhysDV sampledEps = SampleTotalEnergyTransferAliasOneShot(ekin, &IZet[i], r1, r2, r3);
      vecCore::Store(sampledEps, &td->fPhysicsData->fPhysicsScratchpad.fEps[i]);
    }
  } else {

    tracks.GetKinEVec()[N] = tracks.GetKinEVec()[N - 1];
    LPMEnergy[N]           = LPMEnergy[N - 1];
    IZet[N]                = IZet[N - 1];
    SampleTotalEnergyTransferRejVec(tracks.GetKinEVec(), LPMEnergy, IZet, td->fPhysicsData->fPhysicsScratchpad.fEps, N,
                                    td);
  }

  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV ekin;
    vecCore::Load(ekin, tracks.GetKinEVec(i));

    PhysDV eps;
    vecCore::Load(eps, &td->fPhysicsData->fPhysicsScratchpad.fEps[i]);

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

PhysDV VecRelativisticPairModel::SampleTotalEnergyTransferAliasOneShot(const PhysDV egamma, const int *matIDX,
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
    double r2l               = vecCore::Get(r2, l);
    double r3l               = vecCore::Get(r3, l);
    int idx                  = (int)vecCore::Get(indxEgamma, l);
    RatinAliasDataTrans &als = fAliasTablesPerMaterial[matIDX[l]].fTablePerEn[idx];
    double xi                = AliasTableAlternative::SampleRatin(als, fSTNumDiscreteEnergyTransferVals, r2l, r3l, 0);

    vecCore::Set(xiV, l, xi);
  }

  PhysDV deltaMax;
  PhysDV deltaFactor;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    int izet = fAliasTablesPerMaterial[matIDX[l]].fILowestZ;

    const double deltaMaxScalar = gElementData[izet]->fDeltaMaxTsai;
    double deltaFactorScalar    = gElementData[izet]->fDeltaFactor;
    vecCore::Set(deltaMax, l, deltaMaxScalar);
    vecCore::Set(deltaFactor, l, deltaFactorScalar);
  }
  const PhysDV eps0   = geant::units::kElectronMassC2 / egamma;
  const PhysDV epsp   = 0.5 - 0.5 * vecCore::math::Sqrt(1. - 4. * eps0 * deltaFactor / deltaMax);
  const PhysDV epsMin = vecCore::math::Max(eps0, epsp);
  const PhysDV epsV   = epsMin * vecCore::math::Exp(xiV * vecCore::math::Log(0.5 / epsMin));
  return epsV;
}

void VecRelativisticPairModel::SampleTotalEnergyTransferRejVec(const double *egamma, const double *lpmEnergy,
                                                               const int *izet, double *epsOut, int N,
                                                               geant::TaskData *td)
{
  // assert(N>=kPhysDVWidth)
  int currN        = 0;
  PhysDM lanesDone = PhysDM::Zero();
  PhysDI idx;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    vecCore::Set(idx, l, currN);
    ++currN;
  }

  while (currN < N || !lanesDone.isFull()) {
    PhysDV fz;
    PhysDV deltaMax;
    PhysDV deltaFac;
    PhysDV varS1Cond, ilVarS1Cond;

    for (int l = 0; l < kPhysDVWidth; ++l) {
      int lIdx = vecCore::Get(idx, l);
      int lZet = izet[lIdx];
      vecCore::Set(fz, l, gElementData[lZet]->fFz);
      if (fIsUseTsaisScreening) {
        vecCore::Set(deltaMax, l, gElementData[lZet]->fDeltaMaxTsai);
      } else {
        vecCore::Set(deltaMax, l, gElementData[lZet]->fDeltaMax);
      }
      vecCore::Set(deltaFac, l, gElementData[lZet]->fDeltaFactor);
      vecCore::Set(varS1Cond, l, gElementData[lZet]->fVarS1Cond);
      vecCore::Set(ilVarS1Cond, l, gElementData[lZet]->fILVarS1Cond);
    }

    PhysDV egammaV        = vecCore::Gather<PhysDV>(egamma, idx);
    const PhysDV eps0     = geant::units::kElectronMassC2 / egammaV;
    const PhysDV deltaMin = 4. * eps0 * deltaFac;
    const PhysDV epsp     = 0.5 - 0.5 * vecCore::math::Sqrt(1. - deltaMin / deltaMax);
    const PhysDV epsMin   = vecCore::math::Max(eps0, epsp);
    const PhysDV epsRange = 0.5 - epsMin;

    PhysDV F10, F20;
    ScreenFunction12(F10, F20, deltaMin, fIsUseTsaisScreening);
    F10 -= fz;
    F20 -= fz;

    const PhysDV NormF1   = vecCore::math::Max(F10 * epsRange * epsRange, (PhysDV)0.);
    const PhysDV NormF2   = vecCore::math::Max(1.5 * F20, (PhysDV)0.);
    const PhysDV NormCond = NormF1 / (NormF1 + NormF2);

    PhysDV rnd0 = td->fRndm->uniformV();
    PhysDV rnd1 = td->fRndm->uniformV();
    PhysDV rnd2 = td->fRndm->uniformV();

    PhysDV eps     = 0.0;
    PhysDV greject = 0.0;
    PhysDM cond1   = NormCond > rnd0;
    vecCore::MaskedAssign(eps, cond1, 0.5 - epsRange * vecCore::math::Pow(rnd1, (PhysDV)1. / 3.));
    vecCore::MaskedAssign(eps, !cond1, epsMin + epsRange * rnd1);
    const PhysDV delta = deltaFac * eps0 / (eps * (1. - eps));

    PhysDM lpmMask    = PhysDM(fIsUseLPM) && egammaV > fLPMEnergyLimit;
    PhysDV lpmEnergyV = vecCore::Gather<PhysDV>(lpmEnergy, idx);
    if (cond1.isNotEmpty()) {
      if (lpmMask.isNotEmpty()) {
        PhysDV lpmPhiS, lpmGS, lpmXiS, phi1, phi2;
        ComputeScreeningFunctions(phi1, phi2, delta, fIsUseTsaisScreening);
        ComputeLPMfunctions(lpmXiS, lpmGS, lpmPhiS, lpmEnergyV, eps, egammaV, varS1Cond, ilVarS1Cond);
        PhysDV tmpRej = lpmXiS * ((lpmGS + 2. * lpmPhiS) * phi1 - lpmGS * phi2 - lpmPhiS * fz) / F10;
        vecCore::MaskedAssign(greject, cond1 && lpmMask, tmpRej);
      }
      if ((!lpmMask).isNotEmpty()) {
        vecCore::MaskedAssign(greject, cond1 && !lpmMask, (ScreenFunction1(delta, fIsUseTsaisScreening) - fz) / F10);
      }
    }
    if ((!cond1).isNotEmpty()) {
      if (lpmMask.isNotEmpty()) {
        PhysDV lpmPhiS, lpmGS, lpmXiS, phi1, phi2;
        ComputeScreeningFunctions(phi1, phi2, delta, fIsUseTsaisScreening);
        ComputeLPMfunctions(lpmXiS, lpmGS, lpmPhiS, lpmEnergyV, eps, egammaV, varS1Cond, ilVarS1Cond);
        PhysDV tmpRej =
            lpmXiS * ((0.5 * lpmGS + lpmPhiS) * phi1 + 0.5 * lpmGS * phi2 - 0.5 * (lpmGS + lpmPhiS) * fz) / F20;
        vecCore::MaskedAssign(greject, !cond1 && lpmMask, tmpRej);
      }
      if ((!lpmMask).isNotEmpty()) {
        vecCore::MaskedAssign(greject, !cond1 && !lpmMask, (ScreenFunction2(delta, fIsUseTsaisScreening) - fz) / F20);
      }
    }

    PhysDM accepted = greject > rnd2;

    if (accepted.isNotEmpty()) {
      vecCore::Scatter(eps, epsOut, idx);
    }

    lanesDone = lanesDone || accepted;
    for (int l = 0; l < kPhysDVWidth; ++l) {
      auto laneDone = vecCore::Get(accepted, l);
      if (laneDone) {
        if (currN < N) {
          vecCore::Set(idx, l, currN++);
          vecCore::Set(lanesDone, l, false);
        } else {
          vecCore::Set(idx, l, N);
        }
      }
    }
  }
}

void VecRelativisticPairModel::ScreenFunction12(PhysDV &val1, PhysDV &val2, const PhysDV delta, const bool istsai)
{
  if (istsai) {
    const PhysDV gamma  = delta * 0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const PhysDV gamma2 = gamma * gamma;
    const PhysDV dum1   = 33.726 - 4. * vecCore::math::Log(1.0 + 0.311877 * gamma2) +
                        4.8 * vecCore::math::Exp(-0.9 * gamma) + 3.2 * vecCore::math::Exp(-1.5 * gamma);
    const PhysDV dum2 = 2. / (3. + 19.5 * gamma + 18. * gamma2);
    val1              = dum1 + dum2;
    val2              = dum1 - 0.5 * dum2;
  } else {
    PhysDM tmp = delta > 1.;
    vecCore::MaskedAssign(val1, tmp, 42.24 - 8.368 * vecCore::math::Log(delta + 0.952));
    vecCore::MaskedAssign(val2, tmp, val1);
    vecCore::MaskedAssign(val1, !tmp, 42.392 - delta * (7.796 - 1.961 * delta));
    vecCore::MaskedAssign(val2, !tmp, 41.405 - delta * (5.828 - 0.8945 * delta));
  }
}

// 3xPhi_1 - Phi_2: used in case of rejection (either istsai or not)
PhysDV VecRelativisticPairModel::ScreenFunction1(const PhysDV delta, const bool istsai)
{
  PhysDV val;
  if (istsai) {
    const PhysDV gamma  = delta * 0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const PhysDV gamma2 = gamma * gamma;
    val = 33.726 - 4. * vecCore::math::Log(1.0 + 0.311877 * gamma2) + 4.8 * vecCore::math::Exp(-0.9 * gamma) +
          3.2 * vecCore::math::Exp(-1.5 * gamma) + 2. / (3. + 19.5 * gamma + 18. * gamma2);
  } else {
    PhysDM tmp = delta > 1.;
    vecCore::MaskedAssign(val, tmp, 42.24 - 8.368 * vecCore::math::Log(delta + 0.952));
    vecCore::MaskedAssign(val, !tmp, 42.392 - delta * (7.796 - 1.961 * delta));
  }
  return val;
}

// 1.5*Phi_1 + 0.5*Phi_2: used in case of rejection (either istsai or not)
PhysDV VecRelativisticPairModel::ScreenFunction2(const PhysDV delta, const bool istsai)
{
  PhysDV val;
  if (istsai) {
    const PhysDV gamma  = delta * 0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const PhysDV gamma2 = gamma * gamma;
    val = 33.726 - 4. * vecCore::math::Log(1.0 + 0.311877 * gamma2) + 4.8 * vecCore::math::Exp(-0.9 * gamma) +
          3.2 * vecCore::math::Exp(-1.5 * gamma) - 1. / (3. + 19.5 * gamma + 18. * gamma2);
  } else {
    PhysDM tmp = delta > 1.;
    vecCore::MaskedAssign(val, tmp, 42.24 - 8.368 * vecCore::math::Log(delta + 0.952));
    vecCore::MaskedAssign(val, !tmp, 41.405 - delta * (5.828 - 0.8945 * delta));
  }
  return val;
}

void VecRelativisticPairModel::ComputeScreeningFunctions(PhysDV &phi1, PhysDV &phi2, const PhysDV delta,
                                                         const bool istsai)
{
  if (istsai) {
    const PhysDV gamma  = delta * 0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const PhysDV gamma2 = gamma * gamma;
    phi1 = 16.863 - 2.0 * vecCore::math::Log(1.0 + 0.311877 * gamma2) + 2.4 * vecCore::math::Exp(-0.9 * gamma) +
           1.6 * vecCore::math::Exp(-1.5 * gamma);
    phi2 = phi1 - 2.0 / (3.0 + 19.5 * gamma + 18.0 * gamma2); // phi1-phi2
  } else {
    PhysDM tmp = delta > 1.;
    vecCore::MaskedAssign(phi1, tmp, 21.12 - 4.184 * vecCore::math::Log(delta + 0.952));
    vecCore::MaskedAssign(phi2, tmp, phi1);

    PhysDV delta2 = delta * delta;
    vecCore::MaskedAssign(phi1, !tmp, 20.867 - 3.242 * delta + 0.625 * delta2);
    vecCore::MaskedAssign(phi2, !tmp, 20.209 - 1.93 * delta - 0.086 * delta2);
  }
}

void VecRelativisticPairModel::ComputeLPMfunctions(PhysDV &funcXiS, PhysDV &funcGS, PhysDV &funcPhiS, PhysDV lpmenergy,
                                                   PhysDV eps, PhysDV egamma, PhysDV varS1Cond, PhysDV ilVarS1Cond)
{

  //  1. y = E_+/E_{\gamma} with E_+ being the total energy transfered to one of the e-/e+ pair
  //     s'  = \sqrt{ \frac{1}{8} \frac{1}{y(1-y)}   \frac{E^{KL}_{LPM}}{E_{\gamma}}  }
  const PhysDV varSprime = vecCore::math::Sqrt(0.125 * lpmenergy / (eps * egamma * (1.0 - eps)));
  funcXiS                = 2.0;

  PhysDM tmpM = varSprime > 1.0;
  vecCore::MaskedAssign(funcXiS, tmpM, (PhysDV)1.0);
  tmpM = varSprime > varS1Cond;
  if (tmpM.isNotEmpty()) {
    const PhysDV funcHSprime = vecCore::math::Log(varSprime) * ilVarS1Cond;
    PhysDV tmpFuncXiS =
        1.0 + funcHSprime - 0.08 * (1.0 - funcHSprime) * funcHSprime * (2.0 - funcHSprime) * ilVarS1Cond;
    vecCore::MaskedAssign(funcXiS, tmpM, tmpFuncXiS);
  }

  const PhysDV varShat = varSprime / vecCore::math::Sqrt(funcXiS);

  for (int l = 0; l < kPhysDVWidth; ++l) {
    double lFuncGS, lFuncPhiS, lVarShat;
    double lFuncXiS;
    lVarShat = vecCore::Get(varShat, l);
    lFuncXiS = vecCore::Get(funcXiS, l);
    GetLPMFunctions(lFuncGS, lFuncPhiS, lVarShat);
    vecCore::Set(funcGS, l, lFuncGS);
    vecCore::Set(funcPhiS, l, lFuncPhiS);

    if (lFuncXiS * lFuncPhiS > 1. || lVarShat > 0.57) {
      vecCore::Set(funcXiS, l, 1. / lFuncPhiS);
    }
  }
  // MAKE SURE SUPPRESSION IS SMALLER THAN 1: due to Migdal's approximation on xi
}
}
