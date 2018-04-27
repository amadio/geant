#include <Geant/PhysicalConstants.h>
#include <TString.h>
#include <Geant/TaskData.h>
#include <Geant/PhysicsData.h>
#include "Geant/VecKleinNishinaComptonModel.h"
#include "Geant/AliasTable.h"

namespace geantphysics {

using geant::Double_v;
using geant::IndexD_v;
using geant::kVecLenD;
using geant::MaskD_v;

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
    for (int i = 0; i < N; i += kVecLenD) {
      Double_v ekin = tracks.GetKinEVec(i);
      Double_v r1   = td->fRndm->uniformV();
      Double_v r2   = td->fRndm->uniformV();
      Double_v r3   = td->fRndm->uniformV();
      Double_v eps  = SampleReducedPhotonEnergyVec(ekin, r1, r2, r3);
      vecCore::Store(eps, data.fEps + i);
    }
  } else {
    tracks.GetKinEArr()[N] = tracks.GetKinEArr()[N - 1];
    SampleReducedPhotonEnergyRej(tracks.GetKinEArr(), data.fDoubleArr, data.fDoubleArr2, data.fEps, N, td);
  }

  for (int i = 0; i < N; i += kVecLenD) {
    Double_v ekin = tracks.GetKinEVec(i);
    Double_v eps;
    vecCore::Load(eps, &data.fEps[i]);

    Double_v oneMinusCost;
    Double_v sint2;
    if (GetUseSamplingTables()) {
      const Double_v kappa = ekin / geant::units::kElectronMassC2;
      oneMinusCost         = (1. / eps - 1.) / kappa;
      sint2                = oneMinusCost * (2. - oneMinusCost);
    } else {
      vecCore::Load(oneMinusCost, &data.fDoubleArr[i]);
      vecCore::Load(sint2, &data.fDoubleArr2[i]);
    }

    sint2               = Math::Max((Double_v)0., sint2);
    const Double_v cost = 1.0 - oneMinusCost;
    const Double_v sint = Math::Sqrt(sint2);

    Double_v r3        = td->fRndm->uniformV();
    const Double_v phi = geant::units::kTwoPi * (r3);
    // direction of the scattered gamma in the scattering frame
    Double_v sinPhi, cosPhi;
    Math::SinCos(phi, sinPhi, cosPhi);
    Double_v dirX = sint * sinPhi;
    Double_v dirY = sint * cosPhi;
    Double_v dirZ = cost;

    Double_v gammaX = tracks.GetDirXVec(i);
    Double_v gammaY = tracks.GetDirYVec(i);
    Double_v gammaZ = tracks.GetDirZVec(i);
    RotateToLabFrame(dirX, dirY, dirZ, gammaX, gammaY, gammaZ);

    Double_v enDeposit = 0.0;

    Double_v postGammaE = eps * ekin;
    MaskD_v gammaAlive  = postGammaE > GetLowestSecondaryEnergy();
    for (int l = 0; l < kVecLenD; ++l) {
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

    Double_v elEnergy = ekin - postGammaE;
    MaskD_v elAlive   = elEnergy > GetLowestSecondaryEnergy();

    Double_v elDirX = 1.0, elDirY = 0.0, elDirZ = 0.0;
    if (!elAlive.isEmpty()) {
      elDirX              = ekin * gammaX - postGammaE * dirX;
      elDirY              = ekin * gammaY - postGammaE * dirY;
      elDirZ              = ekin * gammaZ - postGammaE * dirZ;
      const Double_v norm = 1.0 / Math::Sqrt(elDirX * elDirX + elDirY * elDirY + elDirZ * elDirZ);
      elDirX *= norm;
      elDirY *= norm;
      elDirZ *= norm;
    }

    for (int l = 0; l < kVecLenD; ++l) {
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

Double_v VecKleinNishinaComptonModel::SampleReducedPhotonEnergyVec(Double_v egamma, Double_v r1, Double_v r2,
                                                                   Double_v r3)
{
  // determine electron energy lower grid point
  Double_v legamma = Math::Log(egamma);
  //
  Double_v val        = (legamma - fSTLogMinPhotonEnergy) * fSTILDeltaPhotonEnergy;
  IndexD_v indxEgamma = (IndexD_v)val; // lower electron energy bin index
  Double_v pIndxHigh  = val - indxEgamma;
  MaskD_v mask        = r1 < pIndxHigh;
  if (!mask.isEmpty()) {
    vecCore::MaskedAssign(indxEgamma, mask, indxEgamma + 1);
  }

  Double_v xiV;
  for (int l = 0; l < kVecLenD; ++l) {
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
  Double_v kappa = egamma * geant::units::kInvElectronMassC2;

  return Math::Exp(Math::Log(1. + 2. * kappa) * (xiV - 1.));
}

void VecKleinNishinaComptonModel::SampleReducedPhotonEnergyRej(const double *egamma, double *onemcostOut,
                                                               double *sint2Out, double *epsOut, int N,
                                                               const geant::TaskData *td)
{
  // assert(N>=kVecLenD)
  int currN         = 0;
  MaskD_v lanesDone = MaskD_v::Zero();
  IndexD_v idx;
  for (int l = 0; l < kVecLenD; ++l) {
    idx[l] = currN++;
  }

  while (currN < N || !lanesDone.isFull()) {

    Double_v kappa = vecCore::Gather<Double_v>(egamma, idx) / geant::units::kElectronMassC2;
    Double_v eps0  = 1. / (1. + 2. * kappa);
    Double_v eps02 = eps0 * eps0;
    Double_v al1   = -Math::Log(eps0);
    Double_v cond  = al1 / (al1 + 0.5 * (1. - eps02));

    Double_v rnd1 = td->fRndm->uniformV();
    Double_v rnd2 = td->fRndm->uniformV();
    Double_v rnd3 = td->fRndm->uniformV();
    // This part is left for future way of rng generating
    // for (int l = 0; l < kVecLenD; ++l){
    //  if(!lanesDone[l]) {
    //    vecCore::Set(rnd1, l, td->fRndm->uniform());
    //    vecCore::Set(rnd2, l, td->fRndm->uniform());
    //    vecCore::Set(rnd3, l, td->fRndm->uniform());
    //  }
    //}

    Double_v eps, eps2, gf;

    MaskD_v cond1 = cond > rnd1;
    if (cond1.isNotEmpty()) {
      vecCore::MaskedAssign(eps, cond1, Math::Exp(-al1 * rnd2));
      vecCore::MaskedAssign(eps2, cond1, eps * eps);
    }
    if ((!cond1).isNotEmpty()) {
      vecCore::MaskedAssign(eps2, !cond1, eps02 + (1.0 - eps02) * rnd2);
      vecCore::MaskedAssign(eps, !cond1, Math::Sqrt(eps2));
    }

    Double_v onemcost = (1. - eps) / (eps * kappa);
    Double_v sint2    = onemcost * (2. - onemcost);
    gf                = 1. - eps * sint2 / (1. + eps2);

    MaskD_v accepted = gf > rnd3;

    if (accepted.isNotEmpty()) {
      vecCore::Scatter(onemcost, onemcostOut, idx);
      vecCore::Scatter(sint2, sint2Out, idx);
      vecCore::Scatter(eps, epsOut, idx);
    }

    lanesDone = lanesDone || accepted;
    for (int l = 0; l < kVecLenD; ++l) {
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
