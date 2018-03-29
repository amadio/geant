#include <Geant/PhysicalConstants.h>
#include <TString.h>
#include <Geant/TaskData.h>
#include <Geant/PhysicsData.h>
#include "Geant/VecPositronTo2GammaModel.h"
#include "Geant/AliasTable.h"

namespace geantphysics {

void VecPositronTo2GammaModel::Initialize()
{
  PositronTo2GammaModel::Initialize();
  //  for (int i = 0; i < (int)fSamplingTables.size(); ++i) {
  //    LinAlias *alias = fSamplingTables[i];
  //    int size        = (int)alias->fAliasIndx.size();
  //
  //    LinAliasCached aliasCached(size);
  //    for (int j = 0; j < size - 1; ++j) {
  //      aliasCached.fLinAliasData[j].fAliasIndx  = alias->fAliasIndx[j];
  //      aliasCached.fLinAliasData[j].fAliasW     = alias->fAliasW[j];
  //      aliasCached.fLinAliasData[j].fX          = alias->fXdata[j];
  //      aliasCached.fLinAliasData[j].fYdata      = alias->fYdata[j];
  //      aliasCached.fLinAliasData[j].fXdelta     = alias->fXdata[j + 1] - alias->fXdata[j];
  //      aliasCached.fLinAliasData[j].fYdataDelta = (alias->fYdata[j + 1] - alias->fYdata[j]) / alias->fYdata[j];
  //      aliasCached.fLinAliasData[j].fXdivYdelta =
  //          aliasCached.fLinAliasData[j].fXdelta / aliasCached.fLinAliasData[j].fYdataDelta;
  //    }
  //    fCachedAliasTable.push_back(aliasCached);
  //  }
}

void VecPositronTo2GammaModel::SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td)
{
  int N                        = tracks.GetNtracks();
  PhysicsModelScratchpad &data = td->fPhysicsData->fPhysicsScratchpad;
  double *epsArr               = data.fEps;
  double *gamma                = data.fDoubleArr; // Used by rejection

  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV pekin  = tracks.GetKinEVec(i);
    PhysDV gammaV = pekin * geant::units::kInvElectronMassC2 + 1.0;
    if (GetUseSamplingTables()) {
      PhysDV rnd1 = td->fRndm->uniformV();
      PhysDV rnd2 = td->fRndm->uniformV();
      PhysDV rnd3 = td->fRndm->uniformV();
      PhysDV epsV = SampleEnergyTransferAlias(pekin, rnd1, rnd2, rnd3, gammaV);
      vecCore::Store(epsV, &epsArr[i]);
    } else {
      vecCore::Store(gammaV, &gamma[i]);
    }
  }
  if (!GetUseSamplingTables()) {
    gamma[N] = gamma[N - 1];
    SampleEnergyTransferRej(gamma, epsArr, N, td);
  }

  for (int i = 0; i < N; i += kPhysDVWidth) {
    PhysDV pekin = tracks.GetKinEVec(i);
    PhysDV tau   = pekin * geant::units::kInvElectronMassC2;
    PhysDV tau2  = tau + 2.0;
    PhysDV eps;
    vecCore::Load(eps, &epsArr[i]);
    // direction of the first gamma
    PhysDV ct         = (eps * tau2 - 1.) / (eps * Sqrt(tau * tau2));
    const PhysDV cost = Max(Min(ct, (PhysDV)1.), (PhysDV)-1.);
    const PhysDV sint = Sqrt((1. + cost) * (1. - cost));
    const PhysDV phi  = geant::units::kTwoPi * td->fRndm->uniformV();
    PhysDV sinPhi, cosPhi;
    SinCos(phi, &sinPhi, &cosPhi);
    PhysDV gamDirX = sint * cosPhi;
    PhysDV gamDirY = sint * sinPhi;
    PhysDV gamDirZ = cost;
    // rotate gamma direction to the lab frame:
    PhysDV posX = tracks.GetDirXVec(i);
    PhysDV posY = tracks.GetDirYVec(i);
    PhysDV posZ = tracks.GetDirZVec(i);
    RotateToLabFrame(gamDirX, gamDirY, gamDirZ, posX, posY, posZ);
    //

    LightTrack_v &secondaries = td->fPhysicsData->GetSecondarySOA();
    // kinematics of the first gamma

    const PhysDV tEnergy = pekin + 2 * geant::units::kElectronMassC2;
    const PhysDV gamEner = eps * tEnergy;
    for (int l = 0; l < kPhysDVWidth; ++l) {
      int idx = secondaries.InsertTrack();
      secondaries.SetKinE(gamEner[l], idx);
      secondaries.SetDirX(gamDirX[l], idx);
      secondaries.SetDirY(gamDirY[l], idx);
      secondaries.SetDirZ(gamDirZ[l], idx);
      secondaries.SetGVcode(fSecondaryInternalCode, idx);
      secondaries.SetMass(0.0, idx);
      secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx);
    }
    //
    // go for the second gamma properties
    const PhysDV posInitTotalMomentum = Sqrt(pekin * (pekin + 2.0 * geant::units::kElectronMassC2));
    // momentum of the second gamma in the lab frame (mom. cons.)
    gamDirX = posInitTotalMomentum * posX - gamEner * gamDirX;
    gamDirY = posInitTotalMomentum * posY - gamEner * gamDirY;
    gamDirZ = posInitTotalMomentum * posZ - gamEner * gamDirZ;
    // normalisation
    const PhysDV norm = 1.0 / Sqrt(gamDirX * gamDirX + gamDirY * gamDirY + gamDirZ * gamDirZ);
    // set up the second gamma track
    for (int l = 0; l < kPhysDVWidth; ++l) {
      int idx = secondaries.InsertTrack();
      secondaries.SetKinE((tEnergy - gamEner)[l], idx);
      secondaries.SetDirX((gamDirX * norm)[l], idx);
      secondaries.SetDirY((gamDirY * norm)[l], idx);
      secondaries.SetDirZ((gamDirZ * norm)[l], idx);
      secondaries.SetGVcode(fSecondaryInternalCode, idx);
      secondaries.SetMass(0.0, idx);
      secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx);
    }
    for (int l = 0; l < kPhysDVWidth; ++l) {
      tracks.SetKinE(0.0, i + l);
      tracks.SetTrackStatus(LTrackStatus::kKill, i + l);
    }
  }
}

PhysDV VecPositronTo2GammaModel::SampleEnergyTransferAlias(PhysDV pekin, PhysDV r1, PhysDV r2, PhysDV r3, PhysDV gamma)
{
  PhysDV lpekin = Log(pekin);
  //
  PhysDV val       = (lpekin - fSTLogMinPositronEnergy) * fSTILDeltaPositronEnergy;
  PhysDI indxPekin = (PhysDI)val; // lower electron energy bin index
  PhysDV pIndxHigh = val - indxPekin;
  PhysDM mask      = r1 < pIndxHigh;
  if (!mask.isEmpty()) {
    vecCore::MaskedAssign(indxPekin, mask, indxPekin + 1);
  }

  PhysDV xiV;
  for (int l = 0; l < kPhysDVWidth; ++l) {
    int idx = (int)indxPekin[l];

    // LinAliasCached &als = fCachedAliasTable[idx];
    //    double xi = AliasTableAlternative::SampleLinear(als, fSTNumDiscreteEnergyTransferVals, r2l, r3l);

    // Standard version:
    const LinAlias *als = fSamplingTables[idx];
    const double xi =
        fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]), &(als->fAliasIndx[0]),
                                    fSTNumDiscreteEnergyTransferVals, r2[l], r3[l]);
    xiV[l] = xi;
  }

  const PhysDV minEps = 0.5 * (1. - Sqrt((gamma - 1.) / (gamma + 1.)));
  const PhysDV maxEps = 0.5 * (1. + Sqrt((gamma - 1.) / (gamma + 1.)));
  const PhysDV eps    = Exp(xiV * Log(maxEps / minEps)) * minEps;

  return eps;
}

void VecPositronTo2GammaModel::SampleEnergyTransferRej(const double *gammaArr, double *epsOut, int N,
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
    PhysDV gamma = vecCore::Gather<PhysDV>(gammaArr, idx);

    const PhysDV minEps = 0.5 * (1. - Sqrt((gamma - 1.) / (gamma + 1.)));
    const PhysDV maxEps = 0.5 * (1. + Sqrt((gamma - 1.) / (gamma + 1.)));
    const PhysDV dum1   = Log(maxEps / minEps);
    const PhysDV dum2   = (gamma + 1.) * (gamma + 1.);

    PhysDV rnd0 = td->fRndm->uniformV();
    PhysDV rnd1 = td->fRndm->uniformV();

    PhysDV eps = minEps * Exp(dum1 * rnd0);

    PhysDV grej     = 1. - eps + (2. * gamma * eps - 1.) / (eps * dum2);
    PhysDM accepted = grej > rnd1;

    if (accepted.isNotEmpty()) {
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
}
