#include <Geant/PhysicalConstants.h>
#include <TString.h>
#include <Geant/TaskData.h>
#include <Geant/PhysicsData.h>
#include "Geant/VecKleinNishinaComptonModel.h"
#include "Geant/AliasTable.h"

namespace geantphysics {

void VecKleinNishinaComptonModel::Initialize()
{
  KleinNishinaComptonModel::Initialize();
  for (int i = 0; i < (int)fSamplingTables.size(); ++i) {
    LinAlias *alias = fSamplingTables[i];
    int size        = (int)alias->fAliasIndx.size();

    LinAliasCached aliasCached(size);
    for (int j = 0; j < size - 1; ++j) { // TODO: What if sampled index is really numdata - 1?
      aliasCached.fLinAliasData[j].fAliasIndx  = alias->fAliasIndx[j];
      aliasCached.fLinAliasData[j].fAliasW     = alias->fAliasW[j];
      aliasCached.fLinAliasData[j].fX          = alias->fXdata[j];
      aliasCached.fLinAliasData[j].fYdata      = alias->fYdata[j];
      aliasCached.fLinAliasData[j].fXdelta     = alias->fXdata[j + 1] - alias->fXdata[j];
      aliasCached.fLinAliasData[j].fYdataDelta = (alias->fYdata[j + 1] - alias->fYdata[j]) / alias->fYdata[j];
      aliasCached.fLinAliasData[j].fXdivYdelta =
          aliasCached.fLinAliasData[j].fXdelta / aliasCached.fLinAliasData[j].fYdataDelta;
    }
    fAliasTablePerGammaEnergy.push_back(aliasCached);
  }
}

void VecKleinNishinaComptonModel::SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td)
{
  int N                        = tracks.GetNtracks();
  PhysicsModelScratchpad &data = td->fPhysicsData->fPhysicsScratchpad;

  if (GetUseSamplingTables()) {
    for (int i = 0; i < N; i += kPhysDVWidth) {
      PhysDV ekin = tracks.GetKinEVec(i);
      PhysDV r1   = td->fRndm->uniformV();
      PhysDV r2   = td->fRndm->uniformV();
      PhysDV r3   = td->fRndm->uniformV();
      PhysDV eps  = SampleReducedPhotonEnergyVec(ekin, r1, r2, r3);
      vecCore::Store(eps, data.fEps + i);
    }
  } else {
    tracks.GetKinEArr()[N] = tracks.GetKinEArr()[N - 1];
    SampleReducedPhotonEnergyRej(tracks.GetKinEArr(), data.fDoubleArr, data.fDoubleArr2, data.fEps, N, td);
  }

  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV ekin = tracks.GetKinEVec(i);
    PhysDV eps;
    vecCore::Load(eps, &data.fEps[i]);

    PhysDV oneMinusCost;
    PhysDV sint2;
    if (GetUseSamplingTables()) {
      const PhysDV kappa = ekin / geant::units::kElectronMassC2;
      oneMinusCost       = (1. / eps - 1.) / kappa;
      sint2              = oneMinusCost * (2. - oneMinusCost);
    } else {
      vecCore::Load(oneMinusCost, &data.fDoubleArr[i]);
      vecCore::Load(sint2, &data.fDoubleArr2[i]);
    }

    sint2             = Max((PhysDV)0., sint2);
    const PhysDV cost = 1.0 - oneMinusCost;
    const PhysDV sint = Sqrt(sint2);

    PhysDV r3        = td->fRndm->uniformV();
    const PhysDV phi = geant::units::kTwoPi * (r3);
    // direction of the scattered gamma in the scattering frame
    PhysDV sinPhi, cosPhi;
    SinCos(phi, &sinPhi, &cosPhi);
    PhysDV dirX = sint * sinPhi;
    PhysDV dirY = sint * cosPhi;
    PhysDV dirZ = cost;

    PhysDV gammaX = tracks.GetDirXVec(i);
    PhysDV gammaY = tracks.GetDirYVec(i);
    PhysDV gammaZ = tracks.GetDirZVec(i);
    RotateToLabFrame(dirX, dirY, dirZ, gammaX, gammaY, gammaZ);

    PhysDV enDeposit = 0.0;

    PhysDV postGammaE = eps * ekin;
    PhysDM gammaAlive = postGammaE > GetLowestSecondaryEnergy();
    for (int l = 0; l < kPhysDVWidth; ++l) {
      bool alive = gammaAlive[l];
      if (alive) {
        tracks.SetKinE(postGammaE[l], i + l);
        tracks.SetDirX(dirX[l], i + l);
        tracks.SetDirY(dirY[l], i + l);
        tracks.SetDirZ(dirZ[l], i + l);
      } else {
        enDeposit[l] = postGammaE[l];
        tracks.SetKinE(0.0, i + l);
        tracks.SetTrackStatus(LTrackStatus::kKill, i + l);
      }
    }

    PhysDV elEnergy = ekin - postGammaE;
    PhysDM elAlive  = elEnergy > GetLowestSecondaryEnergy();

    PhysDV elDirX = 1.0, elDirY = 0.0, elDirZ = 0.0;
    if (!elAlive.isEmpty()) {
      elDirX            = ekin * gammaX - postGammaE * dirX;
      elDirY            = ekin * gammaY - postGammaE * dirY;
      elDirZ            = ekin * gammaZ - postGammaE * dirZ;
      const PhysDV norm = 1.0 / Sqrt(elDirX * elDirX + elDirY * elDirY + elDirZ * elDirZ);
      elDirX *= norm;
      elDirY *= norm;
      elDirZ *= norm;
    }

    for (int l = 0; l < kPhysDVWidth; ++l) {
      LightTrack_v &secondaries = td->fPhysicsData->GetSecondarySOA();

      bool alive = elAlive[l];
      if (alive) {
        int idx = secondaries.InsertTrack();
        secondaries.SetKinE(elEnergy[l], idx);
        secondaries.SetDirX(elDirX[l], idx);
        secondaries.SetDirY(elDirY[l], idx);
        secondaries.SetDirZ(elDirZ[l], idx);
        secondaries.SetGVcode(fSecondaryInternalCode, idx);
        secondaries.SetMass(geant::units::kElectronMassC2, idx);
        secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx);
      } else {
        enDeposit[l] = elEnergy[l];
      }
    }

    tracks.SetEnergyDepositVec(enDeposit, i);
  }
}

PhysDV VecKleinNishinaComptonModel::SampleReducedPhotonEnergyVec(PhysDV egamma, PhysDV r1, PhysDV r2, PhysDV r3)
{
  // determine electron energy lower grid point
  PhysDV legamma = Log(egamma);
  //
  PhysDV val        = (legamma - fSTLogMinPhotonEnergy) * fSTILDeltaPhotonEnergy;
  PhysDI indxEgamma = (PhysDI)val; // lower electron energy bin index
  PhysDV pIndxHigh  = val - indxEgamma;
  PhysDM mask       = r1 < pIndxHigh;
  if (!mask.isEmpty()) {
    vecCore::MaskedAssign(indxEgamma, mask, indxEgamma + 1);
  }

  PhysDV xiV;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    int idx = (int)indxEgamma[l];
    //    LinAliasCached &als = fAliasTablePerGammaEnergy[idx];
    //    double xi           = AliasTableAlternative::SampleLinear(als, fSTNumDiscreteEnergyTransferVals, r2[l],
    //    r3[l]);

    // Standard version:
    const LinAlias *als = fSamplingTables[idx];
    const double xi =
        fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]), &(als->fAliasIndx[0]),
                                    fSTNumDiscreteEnergyTransferVals, r2[l], r3[l]);
    vecCore::Set(xiV, l, xi);
    xiV[l] = xi;
  }
  // transform it back to eps = E_1/E_0
  // \epsion(\xi) = \exp[ \alpha(1-\xi) ] = \exp [\ln(1+2\kappa)(\xi-1)]
  PhysDV kappa = egamma * geant::units::kInvElectronMassC2;

  return Exp(Log(1. + 2. * kappa) * (xiV - 1.));
}

void VecKleinNishinaComptonModel::SampleReducedPhotonEnergyRej(const double *egamma, double *onemcostOut,
                                                               double *sint2Out, double *epsOut, int N,
                                                               const geant::TaskData *td)
{
  // assert(N>=kPhysDVWidth)
  int currN        = 0;
  PhysDM lanesDone = PhysDM::Zero();
  PhysDI idx;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    idx[l] = currN++;
  }

  while (currN < N || !lanesDone.isFull()) {

    PhysDV kappa = vecCore::Gather<PhysDV>(egamma, idx) / geant::units::kElectronMassC2;
    PhysDV eps0  = 1. / (1. + 2. * kappa);
    PhysDV eps02 = eps0 * eps0;
    PhysDV al1   = -Log(eps0);
    PhysDV cond  = al1 / (al1 + 0.5 * (1. - eps02));

    PhysDV rnd1 = td->fRndm->uniformV();
    PhysDV rnd2 = td->fRndm->uniformV();
    PhysDV rnd3 = td->fRndm->uniformV();
    // This part is left for future way of rng generating
    // for (int l = 0; l < kPhysDVWidth; ++l){
    //  if(!lanesDone[l]) {
    //    vecCore::Set(rnd1, l, td->fRndm->uniform());
    //    vecCore::Set(rnd2, l, td->fRndm->uniform());
    //    vecCore::Set(rnd3, l, td->fRndm->uniform());
    //  }
    //}

    PhysDV eps, eps2, gf;

    PhysDM cond1 = cond > rnd1;
    if (cond1.isNotEmpty()) {
      vecCore::MaskedAssign(eps, cond1, Exp(-al1 * rnd2));
      vecCore::MaskedAssign(eps2, cond1, eps * eps);
    }
    if ((!cond1).isNotEmpty()) {
      vecCore::MaskedAssign(eps2, !cond1, eps02 + (1.0 - eps02) * rnd2);
      vecCore::MaskedAssign(eps, !cond1, Sqrt(eps2));
    }

    PhysDV onemcost = (1. - eps) / (eps * kappa);
    PhysDV sint2    = onemcost * (2. - onemcost);
    gf              = 1. - eps * sint2 / (1. + eps2);

    PhysDM accepted = gf > rnd3;

    if (accepted.isNotEmpty()) {
      vecCore::Scatter(onemcost, onemcostOut, idx);
      vecCore::Scatter(sint2, sint2Out, idx);
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

bool VecKleinNishinaComptonModel::IsModelUsable(const MaterialCuts *, double ekin)
{
  return ekin < GetHighEnergyUsageLimit() && ekin > GetLowEnergyUsageLimit();
}
}
