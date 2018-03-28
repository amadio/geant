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
    PhysDV pekin;
    vecCore::Load(pekin, tracks.GetKinEVec(i));
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
    PhysDV pekin;
    vecCore::Load(pekin, tracks.GetKinEVec(i));
    PhysDV tau  = pekin * geant::units::kInvElectronMassC2;
    PhysDV tau2 = tau + 2.0;
    PhysDV eps;
    vecCore::Load(eps, &epsArr[i]);
    // direction of the first gamma
    PhysDV ct         = (eps * tau2 - 1.) / (eps * vecCore::math::Sqrt(tau * tau2));
    const PhysDV cost = vecCore::math::Max(vecCore::math::Min(ct, (PhysDV)1.), (PhysDV)-1.);
    const PhysDV sint = vecCore::math::Sqrt((1. + cost) * (1. - cost));
    const PhysDV phi  = geant::units::kTwoPi * td->fRndm->uniformV();
    PhysDV sinPhi, cosPhi;
    vecCore::math::SinCos(phi, &sinPhi, &cosPhi);
    PhysDV gamDirX = sint * cosPhi;
    PhysDV gamDirY = sint * sinPhi;
    PhysDV gamDirZ = cost;
    // rotate gamma direction to the lab frame:
    PhysDV posX, posY, posZ;
    vecCore::Load(posX, tracks.GetDirXV(i));
    vecCore::Load(posY, tracks.GetDirYV(i));
    vecCore::Load(posZ, tracks.GetDirZV(i));
    RotateToLabFrame(gamDirX, gamDirY, gamDirZ, posX, posY, posZ);
    //

    LightTrack_v &secondaries = td->fPhysicsData->GetSecondarySOA();
    // kinematics of the first gamma

    const PhysDV tEnergy = pekin + 2 * geant::units::kElectronMassC2;
    const PhysDV gamEner = eps * tEnergy;
    for (int l = 0; l < kPhysDVWidth; ++l) {
      int idx = secondaries.InsertTrack();
      secondaries.SetKinE(vecCore::Get(gamEner, l), idx);
      secondaries.SetDirX(vecCore::Get(gamDirX, l), idx);
      secondaries.SetDirY(vecCore::Get(gamDirY, l), idx);
      secondaries.SetDirZ(vecCore::Get(gamDirZ, l), idx);
      secondaries.SetGVcode(fSecondaryInternalCode, idx);
      secondaries.SetMass(0.0, idx);
      secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx);
    }
    //
    // go for the second gamma properties
    const PhysDV posInitTotalMomentum = vecCore::math::Sqrt(pekin * (pekin + 2.0 * geant::units::kElectronMassC2));
    // momentum of the second gamma in the lab frame (mom. cons.)
    gamDirX = posInitTotalMomentum * posX - gamEner * gamDirX;
    gamDirY = posInitTotalMomentum * posY - gamEner * gamDirY;
    gamDirZ = posInitTotalMomentum * posZ - gamEner * gamDirZ;
    // normalisation
    const PhysDV norm = 1.0 / vecCore::math::Sqrt(gamDirX * gamDirX + gamDirY * gamDirY + gamDirZ * gamDirZ);
    // set up the second gamma track
    for (int l = 0; l < kPhysDVWidth; ++l) {
      int idx = secondaries.InsertTrack();
      secondaries.SetKinE(vecCore::Get(tEnergy - gamEner, l), idx);
      secondaries.SetDirX(vecCore::Get(gamDirX * norm, l), idx);
      secondaries.SetDirY(vecCore::Get(gamDirY * norm, l), idx);
      secondaries.SetDirZ(vecCore::Get(gamDirZ * norm, l), idx);
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
  PhysDV lpekin = vecCore::math::Log(pekin);
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
    int idx    = (int)vecCore::Get(indxPekin, l);
    double r2l = vecCore::Get(r2, l);
    double r3l = vecCore::Get(r3, l);

    // LinAliasCached &als = fCachedAliasTable[idx];
    //    double xi = AliasTableAlternative::SampleLinear(als, fSTNumDiscreteEnergyTransferVals, r2l, r3l);

    // Standard version:
    const LinAlias *als = fSamplingTables[idx];
    const double xi     = fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]),
                                                  &(als->fAliasIndx[0]), fSTNumDiscreteEnergyTransferVals, r2l, r3l);
    vecCore::Set(xiV, l, xi);
  }

  const PhysDV minEps = 0.5 * (1. - vecCore::math::Sqrt((gamma - 1.) / (gamma + 1.)));
  const PhysDV maxEps = 0.5 * (1. + vecCore::math::Sqrt((gamma - 1.) / (gamma + 1.)));
  const PhysDV eps    = vecCore::math::Exp(xiV * vecCore::math::Log(maxEps / minEps)) * minEps;

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
    vecCore::Set(idx, l, currN);
    ++currN;
  }

  while (currN < N || !lanesDone.isFull()) {
    PhysDV gamma = vecCore::Gather<PhysDV>(gammaArr, idx);

    const PhysDV minEps = 0.5 * (1. - vecCore::math::Sqrt((gamma - 1.) / (gamma + 1.)));
    const PhysDV maxEps = 0.5 * (1. + vecCore::math::Sqrt((gamma - 1.) / (gamma + 1.)));
    const PhysDV dum1   = vecCore::math::Log(maxEps / minEps);
    const PhysDV dum2   = (gamma + 1.) * (gamma + 1.);

    PhysDV rnd0 = td->fRndm->uniformV();
    PhysDV rnd1 = td->fRndm->uniformV();

    PhysDV eps = minEps * vecCore::math::Exp(dum1 * rnd0);

    PhysDV grej     = 1. - eps + (2. * gamma * eps - 1.) / (eps * dum2);
    PhysDM accepted = grej > rnd1;

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
          vecCore::Set(idx, l, N + 1);
        }
      }
    }
  }
}
}
