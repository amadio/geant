
#include "Geant/PositronTo2GammaModel.h"

// from amterial
#include "Geant/Types.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/Material.h"
#include "Geant/MaterialProperties.h"
#include "Geant/Element.h"

#include "Geant/MaterialCuts.h"

#include "Geant/AliasTable.h"

#include "Geant/PhysicsParameters.h"

#include "Geant/Gamma.h"

#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"

// from geantV
#include "Geant/TaskData.h"

#include <cmath>
#include "Geant/math_wrappers.h"

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

PositronTo2GammaModel::PositronTo2GammaModel(const std::string &modelname) : EMModel(modelname)
{
  fSecondaryInternalCode = -1.;

  fSTNumPositronEnergiesPerDecade  = 8;  // ST=>SamplingTables
  fSTNumDiscreteEnergyTransferVals = 60; // ST=>SamplingTables
  fSTNumPositronEnergies           = -1; // ST=>SamplingTables: will be set at init

  fSTLogMinPositronEnergy  = -1.; // ST=>SamplingTables: will be set at init
  fSTILDeltaPositronEnergy = -1.; // ST=>SamplingTables: will be set at init

  fAliasSampler = nullptr;
}

PositronTo2GammaModel::~PositronTo2GammaModel()
{
  if (GetUseSamplingTables()) {
    ClearSamplingTables();
  }
  if (fAliasSampler) {
    delete fAliasSampler;
  }
}

void PositronTo2GammaModel::Initialize()
{
  EMModel::Initialize();
  fSecondaryInternalCode = Gamma::Definition()->GetInternalCode();
  if (GetUseSamplingTables()) { // if sampling tables were requested
    InitSamplingTables();
  } // else => rejection
}

double PositronTo2GammaModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy,
                                                         const Particle * /*particle*/)
{
  double mxsec = 0.0;
  if (kinenergy < GetLowEnergyUsageLimit() || kinenergy > GetHighEnergyUsageLimit()) {
    return mxsec;
  }
  double elDensity = matcut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol();
  mxsec            = elDensity * ComputeXsectionPerElectron(kinenergy);
  return mxsec;
}

double PositronTo2GammaModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts * /*matcut*/,
                                                     double kinenergy, const Particle * /*particle*/)
{
  double xsec = 0.0;
  if (kinenergy < GetLowEnergyUsageLimit() || kinenergy > GetHighEnergyUsageLimit()) {
    return xsec;
  }
  xsec = elem->GetZ() * ComputeXsectionPerElectron(kinenergy);
  return xsec;
}

double PositronTo2GammaModel::ComputeXsectionPerElectron(double pekin)
{
  constexpr double factor =
      geant::units::kPi * geant::units::kClassicElectronRadius * geant::units::kClassicElectronRadius;
  //
  pekin        = std::max(1. * geant::units::eV, pekin);
  double tau   = pekin / geant::units::kElectronMassC2; // E_kin of the e+ in rest mass units
  double gamma = tau + 1.;
  double g2m1  = tau * (tau + 2.); // gamma^2-1
  double sqx   = std::sqrt(g2m1);  // sqrt(gamma^2-1)
  //
  return factor * ((gamma * (gamma + 4.) + 1.) * Math::Log(gamma + sqx) - (gamma + 3.) * sqx) / (g2m1 * (gamma + 1.));
}

//
// Samples the in-flight e+ e- annihilation (NOTE: the at rest annihilation is implemented into the process)
int PositronTo2GammaModel::SampleSecondaries(LightTrack &track, geant::TaskData *td)
{
  int numSecondaries = 0;
  // sample gamma energy
  const double pekin = track.GetKinE();

  const double tau   = pekin / geant::units::kElectronMassC2; // E_kin of the e+ in rest mass units
  const double tau2  = tau + 2.;
  const double gamma = tau + 1.;
  double eps         = 0.;
  if (GetUseSamplingTables()) {
    double *rndArray = td->fDblArray;
    td->fRndm->uniform_array(3, rndArray);
    eps = SampleEnergyTransfer(pekin, gamma, rndArray[0], rndArray[1], rndArray[2]);
  } else {
    eps = SampleEnergyTransfer(gamma, td);
  }
  //
  // direction of the first gamma
  double ct         = (eps * tau2 - 1.) / (eps * std::sqrt(tau * tau2));
  const double cost = std::max(std::min(ct, 1.), -1.);
  const double sint = std::sqrt((1. + cost) * (1. - cost));
  const double phi  = geant::units::kTwoPi * td->fRndm->uniform();
  double gamDirX    = sint * Math::Cos(phi);
  double gamDirY    = sint * Math::Sin(phi);
  double gamDirZ    = cost;
  // rotate gamma direction to the lab frame:
  Math::RotateToLabFrame(gamDirX, gamDirY, gamDirZ, track.GetDirX(), track.GetDirY(), track.GetDirZ());
  //
  // kinematics of the first gamma
  const double tEnergy = pekin + 2 * geant::units::kElectronMassC2;
  const double gamEner = eps * tEnergy;
  //
  // create the secondary partcile i.e. the gamma
  numSecondaries          = 2;
  LightTrack &gamma1Track = td->fPhysicsData->InsertSecondary();
  gamma1Track.SetDirX(gamDirX);
  gamma1Track.SetDirY(gamDirY);
  gamma1Track.SetDirZ(gamDirZ);
  gamma1Track.SetKinE(gamEner);
  gamma1Track.SetGVcode(fSecondaryInternalCode); // gamma GV code
  gamma1Track.SetMass(0.0);
  gamma1Track.SetTrackIndex(track.GetTrackIndex()); // parent Track index
  //
  // go for the second gamma properties
  const double posInitTotalMomentum = std::sqrt(pekin * (pekin + 2.0 * geant::units::kElectronMassC2));
  // momentum of the second gamma in the lab frame (mom. cons.)
  gamDirX = posInitTotalMomentum * track.GetDirX() - gamEner * gamDirX;
  gamDirY = posInitTotalMomentum * track.GetDirY() - gamEner * gamDirY;
  gamDirZ = posInitTotalMomentum * track.GetDirZ() - gamEner * gamDirZ;
  // normalisation
  const double norm = 1.0 / std::sqrt(gamDirX * gamDirX + gamDirY * gamDirY + gamDirZ * gamDirZ);
  // set up the second gamma track
  LightTrack &gamma2Track = td->fPhysicsData->InsertSecondary();
  gamma2Track.SetDirX(gamDirX * norm);
  gamma2Track.SetDirY(gamDirY * norm);
  gamma2Track.SetDirZ(gamDirZ * norm);
  gamma2Track.SetKinE(tEnergy - gamEner);
  gamma2Track.SetGVcode(fSecondaryInternalCode); // gamma GV code
  gamma2Track.SetMass(0.0);
  gamma2Track.SetTrackIndex(track.GetTrackIndex()); // parent Track index
  // kill the primary e+
  track.SetKinE(0.0);
  track.SetTrackStatus(LTrackStatus::kKill);
  //
  return numSecondaries;
}

double PositronTo2GammaModel::SampleEnergyTransfer(double pekin, double gamma, double r1, double r2, double r3)
{
  // determine electron energy lower grid point
  const double lpekin = Math::Log(pekin);
  //
  int indxPekin = fSTNumPositronEnergies - 1;
  if (pekin < GetHighEnergyUsageLimit()) {
    const double val       = (lpekin - fSTLogMinPositronEnergy) * fSTILDeltaPositronEnergy;
    indxPekin              = (int)val; // lower e+ energy bin index
    const double pIndxHigh = val - indxPekin;
    if (r1 < pIndxHigh) ++indxPekin;
  }
  // sample the transformed variable
  const LinAlias *als = fSamplingTables[indxPekin];
  // sample the transformed variable xi=ln(eps/eps_min) / ln(eps_max/eps_min) that is in [0,1] when eps in [eps_min,
  // eps_max]
  const double xi = fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]),
                                                &(als->fAliasIndx[0]), fSTNumDiscreteEnergyTransferVals, r2, r3);
  // transform it back xi to eps = eps_min epx[xi*ln(eps_max/eps_min)]
  const double minEps = 0.5 * (1. - std::sqrt((gamma - 1.) / (gamma + 1.)));
  const double maxEps = 0.5 * (1. + std::sqrt((gamma - 1.) / (gamma + 1.)));
  const double eps    = Math::Exp(xi * Math::Log(maxEps / minEps)) * minEps;
  return eps;
}

double PositronTo2GammaModel::SampleEnergyTransfer(double gamma, geant::TaskData *td)
{
  const double minEps = 0.5 * (1. - std::sqrt((gamma - 1.) / (gamma + 1.)));
  const double maxEps = 0.5 * (1. + std::sqrt((gamma - 1.) / (gamma + 1.)));
  const double dum1   = Math::Log(maxEps / minEps);
  const double dum2   = (gamma + 1.) * (gamma + 1.);
  double *rndArray    = td->fDblArray;
  double eps          = 0.;
  do {
    td->fRndm->uniform_array(2, rndArray);
    eps = minEps * Math::Exp(dum1 * rndArray[0]);
  } while (1. - eps + (2. * gamma * eps - 1.) / (eps * dum2) < rndArray[1]);
  return eps;
}

void PositronTo2GammaModel::InitSamplingTables()
{
  ClearSamplingTables();
  // set number of primary e+ energy grid points
  const double minEprim  = GetLowEnergyUsageLimit();
  const double maxEprim  = GetHighEnergyUsageLimit();
  fSTNumPositronEnergies = fSTNumPositronEnergiesPerDecade * std::lrint(Math::Log10(maxEprim / minEprim)) + 1;
  fSTNumPositronEnergies = std::max(fSTNumPositronEnergies, 3);
  // set up the initial gamma energy grid
  const double delta       = Math::Log(maxEprim / minEprim) / (fSTNumPositronEnergies - 1.0);
  fSTLogMinPositronEnergy  = Math::Log(minEprim);
  fSTILDeltaPositronEnergy = 1. / delta;
  std::vector<double> primEVect(fSTNumPositronEnergies);
  primEVect[0]                          = minEprim;
  primEVect[fSTNumPositronEnergies - 1] = maxEprim;
  for (int i = 1; i < fSTNumPositronEnergies - 1; ++i) {
    primEVect[i] = Math::Exp(fSTLogMinPositronEnergy + i * delta);
  }
  // 3. create an AliasTable object
  if (fAliasSampler) {
    delete fAliasSampler;
  }
  fAliasSampler = new AliasTable();
  // 4. set up the container that stores sampling tables for all the materials (no material neither Z dependence)
  fSamplingTables.resize(fSTNumPositronEnergies, nullptr);
  // 5. prepare sampling tables one-by-one
  for (int i = 0; i < fSTNumPositronEnergies; ++i) {
    const double pekin = primEVect[i];
    const double gamma = pekin / geant::units::kElectronMassC2 + 1.0;
    BuildOneLinAlias(i, gamma);
  }
  primEVect.clear();
}

void PositronTo2GammaModel::BuildOneLinAlias(int indx, double gamma)
{
  LinAlias *tb          = new LinAlias(fSTNumDiscreteEnergyTransferVals);
  fSamplingTables[indx] = tb;
  //
  const double minEps = 0.5 * (1. - std::sqrt((gamma - 1.) / (gamma + 1.)));
  const double maxEps = 0.5 * (1. + std::sqrt((gamma - 1.) / (gamma + 1.)));
  // note: the transformd variable (xi) is in [0,1] when eps = E_g1/(E_t(e+) + mc^2) i.e. fraction of total energy
  //       transfered to one of the gammas is \in [eps_min, eps_max] where
  //       eps_min = 0.5*(1.-std::sqrt((gamma-1.)/(gamma+1.)));
  //       eps_max = 0.5*(1.+std::sqrt((gamma-1.)/(gamma+1.)));
  // So xi = ln(eps/eps_min)/ln(eps_max/eps_min) ==> eps = eps_min exp[xi ln(eps_max/eps_min)]
  // so fill 3 initial values of xi:
  //  -  xi_0 = x_min = 0
  //  -  xi_1 = (x_max-x_min)/2 = 0.5
  //  -  xi_2 = x_max = 1
  // and the corresponding y(i.e.~PDF) values
  tb->fXdata[0]  = 0.0;
  tb->fXdata[1]  = 0.5;
  tb->fXdata[2]  = 1.0;
  tb->fYdata[0]  = ComputeTransfDXSec(tb->fXdata[0], gamma, minEps, maxEps);
  tb->fYdata[1]  = ComputeTransfDXSec(tb->fXdata[1], gamma, minEps, maxEps);
  tb->fYdata[2]  = ComputeTransfDXSec(tb->fXdata[2], gamma, minEps, maxEps);
  int curNumData = 3;
  // expand the data up to numdata points
  while (curNumData < fSTNumDiscreteEnergyTransferVals) {
    // find the lower index of the bin, where we have the biggest linear interp. error compared to spline
    double maxerr  = 0.0; // value of the current maximum error
    double thexval = 0.0;
    double theyval = 0.0;
    int maxerrindx = 0; // the lower index of the corresponding bin
    for (int i = 0; i < curNumData - 1; ++i) {
      const double xx  = 0.5 * (tb->fXdata[i] + tb->fXdata[i + 1]);     // mid x point
      const double yy  = 0.5 * (tb->fYdata[i] + tb->fYdata[i + 1]);     // lin. interpolated pdf value at the mid point
      const double val = ComputeTransfDXSec(xx, gamma, minEps, maxEps); // real pdf value at the mid point
      const double err = std::abs(1. - (yy / val));
      if (err > maxerr) {
        maxerr     = err;
        maxerrindx = i;
        thexval    = xx;
        theyval    = val;
      }
    }
    // extend x,y data by puting a new real value at the mid point of the highest error bin
    // first shift all values to the right
    for (int j = curNumData; j > maxerrindx + 1; --j) {
      tb->fXdata[j] = tb->fXdata[j - 1];
      tb->fYdata[j] = tb->fYdata[j - 1];
    }
    // fill x mid point
    tb->fXdata[maxerrindx + 1] = thexval;
    tb->fYdata[maxerrindx + 1] = theyval;
    // increase number of data
    ++curNumData;
  } // end while
  // prepare the alias data for this PDF(x,y)
  fAliasSampler->PreparLinearTable(&(tb->fXdata[0]), &(tb->fYdata[0]), &(tb->fAliasW[0]), &(tb->fAliasIndx[0]),
                                   fSTNumDiscreteEnergyTransferVals);
}

void PositronTo2GammaModel::ClearSamplingTables()
{
  size_t num = fSamplingTables.size();
  for (size_t itb = 0; itb < num; ++itb) {
    LinAlias *tb = fSamplingTables[itb];
    if (tb) {
      tb->fXdata.clear();
      tb->fYdata.clear();
      tb->fAliasW.clear();
      tb->fAliasIndx.clear();
    }
    delete tb;
  }
  fSamplingTables.clear();
}

// transformed pdf of xi(eps) = ln[eps/eps_min]/ln[eps_max/eps_min]
double PositronTo2GammaModel::ComputeTransfDXSec(double xi, double gamma, double mineps, double maxeps)
{
  double eps    = Math::Exp(xi * Math::Log(maxeps / mineps)) * mineps;
  double igp1sq = 1. / ((gamma + 1.) * (gamma + 1.));
  return 1. + 2 * gamma * igp1sq - igp1sq / eps - eps;
}

void PositronTo2GammaModel::SampleSecondaries(LightTrack_v &tracks, geant::TaskData *td)
{
  // Prepare temporary arrays for SIMD processing
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
      Double_v epsV = SampleEnergyTransferAlias(pekin, gammaV, rnd2, rnd3, rnd1);
      vecCore::Store(epsV, &epsArr[i]);
    } else {
      vecCore::Store(gammaV, &gamma[i]);
    }
  }
  if (!GetUseSamplingTables()) {
    // Always create fake particle at the end of the input arrays to vector rejection sampling method
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

geant::Double_v PositronTo2GammaModel::SampleEnergyTransferAlias(geant::Double_v pekin, geant::Double_v gamma,
                                                                 geant::Double_v r2, geant::Double_v r3,
                                                                 geant::Double_v r1)
{
  // Compute energy index
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

void PositronTo2GammaModel::SampleEnergyTransferRej(const double *gammaArr, double *epsOut, int N,
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

} // namespace geantphysics
