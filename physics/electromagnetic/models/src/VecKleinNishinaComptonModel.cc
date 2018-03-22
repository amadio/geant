#include <Geant/PhysicalConstants.h>
#include <TString.h>
#include <Geant/TaskData.h>
#include <Geant/PhysicsData.h>
#include "Geant/VecKleinNishinaComptonModel.h"

namespace geantphysics {

void VecKleinNishinaComptonModel::Initialize()
{
  KleinNishinaComptonModel::Initialize();
  for (int i = 0; i < fSamplingTables.size(); ++i) {
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
          aliasCached.fLinAliasData[j].fX / aliasCached.fLinAliasData[j].fYdataDelta;
    }
    fCachedAliasTable.push_back(aliasCached);
  }
}

void VecKleinNishinaComptonModel::SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td)
{
  if (GetUseSamplingTables()) {
    SampleSecondariesVectorAlias(tracks, td);
  } else {
    SampleSecondariesVectorRej(tracks, td);
  }
}

void VecKleinNishinaComptonModel::SampleSecondariesVectorAlias(LightTrack_v &tracks, geant::TaskData *td)
{
  int N                  = tracks.GetNtracks();
  KleinNishinaData &data = td->fPhysicsData->fKleinNishinaData;

  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV r = td->fRndm->uniformV();
    vecCore::Store(r, &data.fR0[i]);
  }
  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV r = td->fRndm->uniformV();
    vecCore::Store(r, &data.fR1[i]);
  }
  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV r = td->fRndm->uniformV();
    vecCore::Store(r, &data.fR2[i]);
  }
  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV r = td->fRndm->uniformV();
    vecCore::Store(r, &data.fR3[i]);
  }

  SampleReducedPhotonEnergyVec(tracks.GetKinEVec(), data.fR0, data.fR1, data.fR2, data.fEps, N);

  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV ekin;
    vecCore::Load(ekin, &tracks.GetKinEVec()[i]);
    PhysDV eps;
    vecCore::Load(eps, &data.fEps[i]);

    const PhysDV kappa  = ekin / geant::units::kElectronMassC2;
    PhysDV oneMinusCost = (1. / eps - 1.) / kappa;
    PhysDV sint2        = oneMinusCost * (2. - oneMinusCost);

    sint2             = vecCore::math::Max((PhysDV)0., sint2);
    const PhysDV cost = 1.0 - oneMinusCost;
    const PhysDV sint = vecCore::math::Sqrt(sint2);

    PhysDV r3;
    vecCore::Load(r3, &data.fR3[i]);
    const PhysDV phi = geant::units::kTwoPi * (r3);
    // direction of the scattered gamma in the scattering frame
    PhysDV tempSin, tempCos;
    vecCore::math::SinCos(phi, &tempSin, &tempCos);
    PhysDV dirX = sint * tempCos;
    PhysDV dirY = sint * tempSin;
    PhysDV dirZ = cost;

    PhysDV gammaX, gammaY, gammaZ;
    vecCore::Load(gammaX, &tracks.GetDirXV()[i]);
    vecCore::Load(gammaY, &tracks.GetDirYV()[i]);
    vecCore::Load(gammaZ, &tracks.GetDirZV()[i]);
    RotateToLabFrame(dirX, dirY, dirZ, gammaX, gammaY, gammaZ);

    PhysDV enDeposit = 0.0;

    PhysDV postGammaE = eps * ekin;
    PhysDM gammaAlive = postGammaE > GetLowestSecondaryEnergy();
    for (int l = 0; l < kPhysDVWidth; ++l) {
      bool alive = vecCore::Get(gammaAlive, l);
      if (alive) {
        tracks.SetKinE(vecCore::Get(postGammaE, l), i + l);
        tracks.SetDirX(vecCore::Get(dirX, l), i + l);
        tracks.SetDirY(vecCore::Get(dirY, l), i + l);
        tracks.SetDirZ(vecCore::Get(dirZ, l), i + l);
      } else {
        vecCore::Set(enDeposit, l, vecCore::Get(postGammaE, l));
        vecCore::Set(enDeposit, l, vecCore::Get(postGammaE, l));
        tracks.SetKinE(0.0, i + l);
        tracks.SetTrackStatus(LTrackStatus::kKill, i + l);
      }
    }

    PhysDV elEnergy = ekin - postGammaE;
    PhysDM elAlive  = elEnergy > GetLowestSecondaryEnergy();

    PhysDV elDirX = 1.0, elDirY = 0.0, elDirZ = 0.0;
    if (!elAlive.isEmpty()) {
      elDirX = ekin * gammaX - postGammaE * dirX;
      elDirY = ekin * gammaY - postGammaE * dirY;
      elDirZ = ekin * gammaZ - postGammaE * dirZ;
      // normalisation factor
      const PhysDV norm = 1.0 / vecCore::math::Sqrt(elDirX * elDirX + elDirY * elDirY + elDirZ * elDirZ);
      elDirX *= norm;
      elDirY *= norm;
      elDirZ *= norm;
    }

    for (int l = 0; l < kPhysDVWidth; ++l) {
      LightTrack_v &secondaries = td->fPhysicsData->GetSecondarySOA();

      bool alive = vecCore::Get(elAlive, l);
      if (alive) {
        int idx = secondaries.InsertTrack();
        secondaries.SetKinE(vecCore::Get(elEnergy, l), idx);
        secondaries.SetDirX(vecCore::Get(elDirX, l), idx);
        secondaries.SetDirY(vecCore::Get(elDirY, l), idx);
        secondaries.SetDirZ(vecCore::Get(elDirZ, l), idx);
        secondaries.SetGVcode(fSecondaryInternalCode, idx);
        secondaries.SetMass(geant::units::kElectronMassC2, idx);
        secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx);
      } else {
        vecCore::Set(enDeposit, l, vecCore::Get(elEnergy, l));
      }
    }

    vecCore::Store(enDeposit, &tracks.GetEnergyDepositVec()[i]);
  }
}

void VecKleinNishinaComptonModel::SampleSecondariesVectorRej(LightTrack_v &tracks, geant::TaskData *td)
{
  int N                  = tracks.GetNtracks();
  KleinNishinaData &data = td->fPhysicsData->fKleinNishinaData;

  SampleReducedPhotonEnergyRej(tracks.GetKinEVec(), data.fOneMCos, data.fSin2t, data.fEps, N, td);

  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV ekin;
    vecCore::Load(ekin, &tracks.GetKinEVec()[i]);
    PhysDV eps;
    vecCore::Load(eps, &data.fEps[i]);

    PhysDV oneMinusCost;
    PhysDV sint2;
    vecCore::Load(oneMinusCost, &data.fOneMCos[i]);
    vecCore::Load(sint2, &data.fSin2t[i]);

    sint2             = vecCore::math::Max((PhysDV)0., sint2);
    const PhysDV cost = 1.0 - oneMinusCost;
    const PhysDV sint = vecCore::math::Sqrt(sint2);

    PhysDV r3        = td->fRndm->uniformV();
    const PhysDV phi = geant::units::kTwoPi * (r3);
    // direction of the scattered gamma in the scattering frame
    PhysDV tempSin, tempCos;
    vecCore::math::SinCos(phi, &tempSin, &tempCos);
    PhysDV dirX = sint * tempCos;
    PhysDV dirY = sint * tempSin;
    PhysDV dirZ = cost;

    PhysDV gammaX, gammaY, gammaZ;
    vecCore::Load(gammaX, &tracks.GetDirXV()[i]);
    vecCore::Load(gammaY, &tracks.GetDirYV()[i]);
    vecCore::Load(gammaZ, &tracks.GetDirZV()[i]);
    RotateToLabFrame(dirX, dirY, dirZ, gammaX, gammaY, gammaZ);

    PhysDV enDeposit = 0.0;

    PhysDV postGammaE = eps * ekin;
    PhysDM gammaAlive = postGammaE > GetLowestSecondaryEnergy();
    for (int l = 0; l < kPhysDVWidth; ++l) {
      bool alive = vecCore::Get(gammaAlive, l);
      if (alive) {
        tracks.SetKinE(vecCore::Get(postGammaE, l), i + l);
        tracks.SetDirX(vecCore::Get(dirX, l), i + l);
        tracks.SetDirY(vecCore::Get(dirY, l), i + l);
        tracks.SetDirZ(vecCore::Get(dirZ, l), i + l);
      } else {
        vecCore::Set(enDeposit, l, vecCore::Get(postGammaE, l));
        vecCore::Set(enDeposit, l, vecCore::Get(postGammaE, l));
        tracks.SetKinE(0.0, i + l);
        tracks.SetTrackStatus(LTrackStatus::kKill, i + l);
      }
    }

    PhysDV elEnergy = ekin - postGammaE;
    PhysDM elAlive  = elEnergy > GetLowestSecondaryEnergy();

    PhysDV elDirX = 1.0, elDirY = 0.0, elDirZ = 0.0;
    if (!elAlive.isEmpty()) {
      elDirX = ekin * gammaX - postGammaE * dirX;
      elDirY = ekin * gammaY - postGammaE * dirY;
      elDirZ = ekin * gammaZ - postGammaE * dirZ;
      // normalisation factor
      const PhysDV norm = 1.0 / vecCore::math::Sqrt(elDirX * elDirX + elDirY * elDirY + elDirZ * elDirZ);
      elDirX *= norm;
      elDirY *= norm;
      elDirZ *= norm;
    }

    for (int l = 0; l < kPhysDVWidth; ++l) {
      LightTrack_v &secondaries = td->fPhysicsData->GetSecondarySOA();

      bool alive = vecCore::Get(elAlive, l);
      if (alive) {
        int idx = secondaries.InsertTrack();
        secondaries.SetKinE(vecCore::Get(elEnergy, l), idx);
        secondaries.SetDirX(vecCore::Get(elDirX, l), idx);
        secondaries.SetDirY(vecCore::Get(elDirY, l), idx);
        secondaries.SetDirZ(vecCore::Get(elDirZ, l), idx);
        secondaries.SetGVcode(fSecondaryInternalCode, idx);
        secondaries.SetMass(geant::units::kElectronMassC2, idx);
        secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx);
      } else {
        vecCore::Set(enDeposit, l, vecCore::Get(elEnergy, l));
      }
    }

    vecCore::Store(enDeposit, &tracks.GetEnergyDepositVec()[i]);
  }
}

void VecKleinNishinaComptonModel::SampleReducedPhotonEnergyVec(const double *egamma, const double *r1, const double *r2,
                                                               const double *r3, double *out, int N)
{

  for (int i = 0; i < N; i += kPhysDVWidth) {
    // determine electron energy lower grid point
    PhysDV legamma;
    vecCore::Load(legamma, &egamma[i]);
    legamma = vecCore::math::Log(legamma);
    //
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
      int idx             = (int)vecCore::Get(indxEgamma, l);
      LinAliasCached &als = fCachedAliasTable[idx];
      double xi = AliasTableAlternative::SampleLinear(als, fSTNumDiscreteEnergyTransferVals, r2[i + l], r3[i + l]);

      // Standard version:
      // const LinAlias *als = fSamplingTables[idx];
      //      const double xi = fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]),
      //                                                    &(als->fAliasIndx[0]), fSTNumDiscreteEnergyTransferVals,
      //                                                    r2[i+l], r3[i+l]);
      vecCore::Set(xiV, l, xi);
    }
    // transform it back to eps = E_1/E_0
    // \epsion(\xi) = \exp[ \alpha(1-\xi) ] = \exp [\ln(1+2\kappa)(\xi-1)]
    PhysDV kappa;
    vecCore::Load(kappa, &egamma[i]);
    kappa = kappa * (geant::units::kInvElectronMassC2);

    PhysDV outV = vecCore::math::Exp(vecCore::math::Log(1. + 2. * kappa) * (xiV - 1.));
    vecCore::Store(outV, &out[i]);
  }
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
    vecCore::Set(idx, l, currN);
    ++currN;
  }

  while (currN < N || !lanesDone.isFull()) {

    PhysDV kappa = vecCore::Gather<PhysDV>(egamma, idx) / geant::units::kElectronMassC2;
    PhysDV eps0  = 1. / (1. + 2. * kappa);
    PhysDV eps02 = eps0 * eps0;
    PhysDV al1   = -vecCore::math::Log(eps0);
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
      vecCore::MaskedAssign(eps, cond1, vecCore::math::Exp(-al1 * rnd2));
      vecCore::MaskedAssign(eps2, cond1, eps * eps);
    }
    if ((!cond1).isNotEmpty()) {
      vecCore::MaskedAssign(eps2, !cond1, eps02 + (1.0 - eps02) * rnd2);
      vecCore::MaskedAssign(eps, !cond1, vecCore::math::Sqrt(eps2));
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
      auto laneDone = vecCore::Get(accepted, l);
      if (laneDone) {
        if (currN < N) {
          vecCore::Set(idx, l, currN++);
          vecCore::Set(lanesDone, l, false);
        } else {
          vecCore::Set(idx, l, N + 1);
        }
      }
    }
  }
}
}
