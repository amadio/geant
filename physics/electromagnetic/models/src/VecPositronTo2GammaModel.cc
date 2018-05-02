#include <Geant/PhysicalConstants.h>
#include <TString.h>
#include <Geant/TaskData.h>
#include <Geant/PhysicsData.h>
#include "Geant/VecPositronTo2GammaModel.h"
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

  for (int i = 0; i < N; i += kVecLenD) {
    Double_v pekin  = tracks.GetKinEVec(i);
    Double_v gammaV = pekin * geant::units::kInvElectronMassC2 + 1.0;
    if (GetUseSamplingTables()) {
      Double_v rnd1 = td->fRndm->uniformV();
      Double_v rnd2 = td->fRndm->uniformV();
      Double_v rnd3 = td->fRndm->uniformV();
      Double_v epsV = SampleEnergyTransferAlias(pekin, rnd1, rnd2, rnd3, gammaV);
      vecCore::Store(epsV, &epsArr[i]);
    } else {
      vecCore::Store(gammaV, &gamma[i]);
    }
  }
  if (!GetUseSamplingTables()) {
    gamma[N] = gamma[N - 1];
    SampleEnergyTransferRej(gamma, epsArr, N, td);
  }

  for (int i = 0; i < N; i += kVecLenD) {
    Double_v pekin = tracks.GetKinEVec(i);
    Double_v tau   = pekin * geant::units::kInvElectronMassC2;
    Double_v tau2  = tau + 2.0;
    Double_v eps;
    vecCore::Load(eps, &epsArr[i]);
    // direction of the first gamma
    Double_v ct         = (eps * tau2 - 1.) / (eps * Math::Sqrt(tau * tau2));
    const Double_v cost = Math::Max(Math::Min(ct, (Double_v)1.), (Double_v)-1.);
    const Double_v sint = Math::Sqrt((1. + cost) * (1. - cost));
    const Double_v phi  = geant::units::kTwoPi * td->fRndm->uniformV();
    Double_v sinPhi, cosPhi;
    Math::SinCos(phi, sinPhi, cosPhi);
    Double_v gamDirX = sint * cosPhi;
    Double_v gamDirY = sint * sinPhi;
    Double_v gamDirZ = cost;
    // rotate gamma direction to the lab frame:
    Double_v posX = tracks.GetDirXVec(i);
    Double_v posY = tracks.GetDirYVec(i);
    Double_v posZ = tracks.GetDirZVec(i);
    Math::RotateToLabFrame(gamDirX, gamDirY, gamDirZ, posX, posY, posZ);
    //

    LightTrack_v &secondaries = td->fPhysicsData->GetSecondarySOA();
    // kinematics of the first gamma

    const Double_v tEnergy = pekin + 2 * geant::units::kElectronMassC2;
    const Double_v gamEner = eps * tEnergy;
    for (int l = 0; l < kVecLenD; ++l) {
      int idx = secondaries.InsertTrack();
      secondaries.SetKinE(Get(gamEner, l), idx);
      secondaries.SetDirX(Get(gamDirX, l), idx);
      secondaries.SetDirY(Get(gamDirY, l), idx);
      secondaries.SetDirZ(Get(gamDirZ, l), idx);
      secondaries.SetGVcode(fSecondaryInternalCode, idx);
      secondaries.SetMass(0.0, idx);
      secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx);
    }
    //
    // go for the second gamma properties
    const Double_v posInitTotalMomentum = Math::Sqrt(pekin * (pekin + 2.0 * geant::units::kElectronMassC2));
    // momentum of the second gamma in the lab frame (mom. cons.)
    gamDirX = posInitTotalMomentum * posX - gamEner * gamDirX;
    gamDirY = posInitTotalMomentum * posY - gamEner * gamDirY;
    gamDirZ = posInitTotalMomentum * posZ - gamEner * gamDirZ;
    // normalisation
    const Double_v norm = 1.0 / Math::Sqrt(gamDirX * gamDirX + gamDirY * gamDirY + gamDirZ * gamDirZ);
    // set up the second gamma track
    for (int l = 0; l < kVecLenD; ++l) {
      int idx = secondaries.InsertTrack();
      secondaries.SetKinE(Get(tEnergy - gamEner, l), idx);
      secondaries.SetDirX(Get(gamDirX * norm, l), idx);
      secondaries.SetDirY(Get(gamDirY * norm, l), idx);
      secondaries.SetDirZ(Get(gamDirZ * norm, l), idx);
      secondaries.SetGVcode(fSecondaryInternalCode, idx);
      secondaries.SetMass(0.0, idx);
      secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx);
    }
    for (int l = 0; l < kVecLenD; ++l) {
      tracks.SetKinE(0.0, i + l);
      tracks.SetTrackStatus(LTrackStatus::kKill, i + l);
    }
  }
}

Double_v VecPositronTo2GammaModel::SampleEnergyTransferAlias(Double_v pekin, Double_v r1, Double_v r2, Double_v r3,
                                                             Double_v gamma)
{
  Double_v lpekin = Math::Log(pekin);
  //
  Double_v val       = (lpekin - fSTLogMinPositronEnergy) * fSTILDeltaPositronEnergy;
  IndexD_v indxPekin = (IndexD_v)val; // lower electron energy bin index
  Double_v pIndxHigh = val - indxPekin;
  MaskD_v mask       = r1 < pIndxHigh;
  if (!MaskEmpty(mask)) {
    vecCore::MaskedAssign(indxPekin, mask, indxPekin + 1);
  }

  Double_v xiV;
  for (int l = 0; l < kVecLenD; ++l) {
    int idx = (int)Get(indxPekin, l);

    // LinAliasCached &als = fCachedAliasTable[idx];
    //    double xi = AliasTableAlternative::SampleLinear(als, fSTNumDiscreteEnergyTransferVals, r2l, r3l);

    // Standard version:
    const LinAlias *als = fSamplingTables[idx];
    const double xi =
        fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]), &(als->fAliasIndx[0]),
                                    fSTNumDiscreteEnergyTransferVals, Get(r2, l), Get(r3, l));
    Set(xiV, l, xi);
  }

  const Double_v minEps = 0.5 * (1. - Math::Sqrt((gamma - 1.) / (gamma + 1.)));
  const Double_v maxEps = 0.5 * (1. + Math::Sqrt((gamma - 1.) / (gamma + 1.)));
  const Double_v eps    = Math::Exp(xiV * Math::Log(maxEps / minEps)) * minEps;

  return eps;
}

void VecPositronTo2GammaModel::SampleEnergyTransferRej(const double *gammaArr, double *epsOut, int N,
                                                       const geant::TaskData *td)
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
    Double_v gamma = vecCore::Gather<Double_v>(gammaArr, idx);

    const Double_v minEps = 0.5 * (1. - Math::Sqrt((gamma - 1.) / (gamma + 1.)));
    const Double_v maxEps = 0.5 * (1. + Math::Sqrt((gamma - 1.) / (gamma + 1.)));
    const Double_v dum1   = Math::Log(maxEps / minEps);
    const Double_v dum2   = (gamma + 1.) * (gamma + 1.);

    Double_v rnd0 = td->fRndm->uniformV();
    Double_v rnd1 = td->fRndm->uniformV();

    Double_v eps = minEps * Math::Exp(dum1 * rnd0);

    Double_v grej    = 1. - eps + (2. * gamma * eps - 1.) / (eps * dum2);
    MaskD_v accepted = grej > rnd1;

    if (!MaskEmpty(accepted)) {
      vecCore::Scatter(eps, epsOut, idx);
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
}
