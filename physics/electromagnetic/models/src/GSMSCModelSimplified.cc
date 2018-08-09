
#include "Geant/GSMSCModelSimplified.h"

// from material
#include "Geant/Types.h"
#include "Geant/Material.h"
#include "Geant/MaterialProperties.h"
#include "Geant/MaterialCuts.h"
#include "Geant/Element.h"

#include "Geant/Region.h"
#include "Geant/ELossTableManager.h"

#include "Geant/GSMSCTableSimplified.h"
#include "Geant/GSPWACorrections.h"

#include "Geant/Particle.h"
#include "Geant/Electron.h"
#include "Geant/Positron.h"

// from geantV
#include "Geant/TaskData.h"
#include "Geant/Track.h"
#include "Geant/math_wrappers.h"

#include <cmath>
#include "Geant/PhysicsData.h"
#include "Geant/VectorTypes.h"

namespace geantphysics {

using geant::Double_v;
using geant::IndexD_v;
using geant::kVecLenD;
using geant::MaskD_v;
using vecCore::AssignMaskLane;
using vecCore::Get;
using vecCore::MaskEmpty;
using vecCore::MaskFull;
using vecCore::Set;

GSMSCModelSimplified::GSMSCModelSimplified(bool iselectron, const std::string &name)
    : MSCModel(name), fIsElectron(iselectron), fGSTable(iselectron), fPWACorrection(iselectron)
{
  fParticle = Electron::Definition();
  if (!fIsElectron) {
    fParticle = Positron::Definition();
  }
  fCharge = fParticle->GetPDGCharge();
}

void GSMSCModelSimplified::Initialize()
{
  // init
  fGSTable.Initialize(GetLowEnergyUsageLimit(), GetHighEnergyUsageLimit(), GetListActiveRegions());
  // create PWA corrections table if it was requested (and not disactivated because active Mott-correction)
  fPWACorrection.Initialise(GetListActiveRegions());
  // Register MSCData into the tracks data
  fMSCdata = geant::TrackDataMgr::GetInstance()->RegisterDataType<MSCdata>("MSCdata");
}

void GSMSCModelSimplified::StepLimit(geant::Track *gtrack, geant::TaskData *td)
{
  bool isOnBoundary          = gtrack->Boundary();
  const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
      (const_cast<vecgeom::LogicalVolume *>(gtrack->GetVolume())->GetMaterialCutsPtr()));
  double kineticEnergy = gtrack->Ekin();

  if (kineticEnergy < GetLowEnergyUsageLimit() || kineticEnergy > GetHighEnergyUsageLimit()) {
    return;
  }

  double range = ELossTableManager::Instance().GetRestrictedRange(matCut, fParticle, kineticEnergy, gtrack->LogEkin());
  double presafety = gtrack->GetSafety(); // pre-step point safety
  //
  // Compute elastic mfp, first transport mfp, screening parameter and G1
  double lambel;
  double lambtr1;
  double scra;
  double g1;
  double pMCtoQ1;
  double pMCtoG2PerG1;
  MSCdata &mscdata = fMSCdata.Data<MSCdata>(gtrack);
  ComputeParameters(matCut, kineticEnergy, lambel, lambtr1, scra, g1, pMCtoQ1, pMCtoG2PerG1);
  mscdata.fLambda0 = lambel;
  mscdata.fLambda1 = lambtr1;
  mscdata.fScrA    = scra;
  mscdata.fG1      = g1;
  mscdata.fRange   = range;
  //
  // Set initial values:
  //  : lengths are initialised to the current minimum physics step  which is the true, minimum
  //    step length from all other physics
  mscdata.fTheTrueStepLenght    = gtrack->GetPstep();
  mscdata.fTheTransportDistance = gtrack->GetPstep();
  mscdata.fTheZPathLenght       = gtrack->GetPstep();
  mscdata.SetDisplacement(0., 0., 0.);
  mscdata.SetNewDirectionMsc(0., 0., 1.);

  // Can everything be done in the step limit phase ?
  mscdata.fIsEverythingWasDone = false;
  // Single scattering needs to be sample ?
  mscdata.fIsSingleScattering = false;
  // Was zero deflection in multiple scattering sampling ?
  mscdata.fIsNoScatteringInMSC = false;
  // Do not care about displacement in MSC sampling
  // ( used only in the case of fgIsOptimizationOn = true)
  mscdata.fIsNoDisplace = false;

  // Zeff  = TotNbOfElectPerVolume/TotNbOfAtomsPerVolume
  double fZeff = matCut->GetMaterial()->GetMaterialProperties()->GetEffectiveZ();
  // matCut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()/
  // matCut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfAtomsPerVol();

  // distance will take into account max-fluct: we don't have fluctuation but use this for consistency
  double distance = range;
  distance *= (1.20 - fZeff * (1.62e-2 - 9.22e-5 * fZeff));

  // indicate that MSC needs to be done (always and always after transportation)
  mscdata.fIsMultipleSacettring = true;
  // far from boundary-> in optimized mode do not sample dispalcement. (safety is zero if we are on boundary)
  if (distance < presafety) {
    mscdata.fIsNoDisplace = true;
  } else {
    // Urban like
    double fr = GetRangeFactor();
    if (mscdata.fIsFirstStep || isOnBoundary || mscdata.fTheInitialRange > 1.e+20) { // NOTE:
      mscdata.fTheInitialRange = range;
      // We don't use this: we won't converge to the single scattering results with
      //                    decreasing range-factor.
      //              rangeinit = std::max(rangeinit, fLambda1);
      //              if(fLambda1 > lambdalimit) {
      //                fr *= (0.75+0.25*fLambda1/lambdalimit);
      //              }
    }
    // step limit
    double tlimit = std::max(fr * mscdata.fTheInitialRange, GetSafetyFactor() * presafety);
    // first step randomization
    if (mscdata.fIsFirstStep ||
        isOnBoundary) { // NOTE: we should probably randomize only if the step was limited by msc
      mscdata.fTheTrueStepLenght = std::min(mscdata.fTheTrueStepLenght, RandomizeTrueStepLength(td, tlimit));
    } else {
      mscdata.fTheTrueStepLenght = std::min(mscdata.fTheTrueStepLenght, tlimit);
    }
  }
  // msc step limit is done!
  //
  //
  // reset first step flag
  mscdata.fIsFirstStep = false;
  // performe single scattering, multiple scattering if this later can be done safely here
  // convert the true step length to geometrical one
  ConvertTrueToGeometricLength(gtrack, td);
}

bool GSMSCModelSimplified::SampleScattering(geant::Track *gtrack, geant::TaskData *td)
{
  SampleMSC(gtrack, td);

  MSCdata &mscdata = fMSCdata.Data<MSCdata>(gtrack);
  if (!mscdata.fIsNoScatteringInMSC) {
    Math::RotateToLabFrame(mscdata.fTheNewDirectionX, mscdata.fTheNewDirectionY, mscdata.fTheNewDirectionZ,
                           gtrack->Dx(), gtrack->Dy(), gtrack->Dz());
    if (!mscdata.fIsNoDisplace) {
      Math::RotateToLabFrame(mscdata.fTheDisplacementVectorX, mscdata.fTheDisplacementVectorY,
                             mscdata.fTheDisplacementVectorZ, gtrack->Dx(), gtrack->Dy(), gtrack->Dz());
    }
    return true; // new direction
  }
  return false;
}
void GSMSCModelSimplified::SampleScattering(geant::TrackVec_t &gtracks, std::vector<bool> &hasNewDir,
                                            geant::TaskData *td)
{
  SampleMSC(gtracks, td);
  for (size_t i = 0; i < gtracks.size(); ++i) {
    auto gtrack      = gtracks[i];
    MSCdata &mscdata = fMSCdata.Data<MSCdata>(gtrack);
    if (!mscdata.fIsNoScatteringInMSC) {
      Math::RotateToLabFrame(mscdata.fTheNewDirectionX, mscdata.fTheNewDirectionY, mscdata.fTheNewDirectionZ,
                             gtrack->Dx(), gtrack->Dy(), gtrack->Dz());
      if (!mscdata.fIsNoDisplace) {
        Math::RotateToLabFrame(mscdata.fTheDisplacementVectorX, mscdata.fTheDisplacementVectorY,
                               mscdata.fTheDisplacementVectorZ, gtrack->Dx(), gtrack->Dy(), gtrack->Dz());
      }
      hasNewDir[i] = true;
    } else {
      hasNewDir[i] = false;
    }
  }
}

void GSMSCModelSimplified::ConvertTrueToGeometricLength(geant::Track *gtrack, geant::TaskData * /*td*/)
{
  MSCdata &mscdata = fMSCdata.Data<MSCdata>(gtrack);
  mscdata.fPar1    = -1.;
  mscdata.fPar2    = 0.;
  mscdata.fPar3    = 0.;
  // if fIsEverythingWasDone = TRUE  => fTheZPathLenght is already set so return with the already known value
  // Otherwise:
  // this correction needed to run MSC with eIoni and eBrem inactivated
  // and makes no harm for a normal run
  mscdata.fTheTrueStepLenght = std::min(mscdata.fTheTrueStepLenght, mscdata.fRange);
  //  do the true -> geom transformation
  mscdata.fTheZPathLenght = mscdata.fTheTrueStepLenght;
  // z = t for very small true-path-length
  if (mscdata.fTheTrueStepLenght < fTLimitMinfix2) {
    return;
  }
  //
  double ekin = gtrack->Ekin();
  double tau  = mscdata.fTheTrueStepLenght / mscdata.fLambda1;
  if (tau <= fTauSmall) {
    mscdata.fTheZPathLenght = std::min(mscdata.fTheTrueStepLenght, mscdata.fLambda1);
  } else if (mscdata.fTheTrueStepLenght < mscdata.fRange * fDtrl) {
    if (tau < fTauLim) {
      mscdata.fTheZPathLenght = mscdata.fTheTrueStepLenght * (1. - 0.5 * tau);
    } else {
      mscdata.fTheZPathLenght = mscdata.fLambda1 * (1. - Math::Exp(-tau));
    }
  } else if (ekin < geant::units::kElectronMassC2 || mscdata.fTheTrueStepLenght == mscdata.fRange) {
    mscdata.fPar1 = 1. / mscdata.fRange;                     // alpha =1/range_init for Ekin<mass
    mscdata.fPar2 = 1. / (mscdata.fPar1 * mscdata.fLambda1); // 1/(alphaxlambda01)
    mscdata.fPar3 = 1. + mscdata.fPar2;
    if (mscdata.fTheTrueStepLenght < mscdata.fRange) {
      mscdata.fTheZPathLenght = 1. / (mscdata.fPar1 * mscdata.fPar3) *
                                (1. - Math::Pow(1. - mscdata.fPar1 * mscdata.fTheTrueStepLenght, mscdata.fPar3));
    } else {
      mscdata.fTheZPathLenght = 1. / (mscdata.fPar1 * mscdata.fPar3);
    }
  } else {
    //      int    matIndx              = gtrack->GetMaterial()->GetIndex();
    //      int    regIndx              =
    //      const_cast<vecgeom::LogicalVolume*>(gtrack->GetVolume())->GetRegion()->GetIndex();
    //      const  MaterialCuts *matCut = MaterialCuts::GetMaterialCut(regIndx,matIndx);
    const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
        (const_cast<vecgeom::LogicalVolume *>(gtrack->GetVolume())->GetMaterialCutsPtr()));
    double rfin             = std::max(mscdata.fRange - mscdata.fTheTrueStepLenght, 0.01 * mscdata.fRange);
    double T1               = ELossTableManager::Instance().GetEnergyForRestrictedRange(matCut, fParticle, rfin);
    double lambda1          = GetTransportMeanFreePathOnly(matCut, T1);
    mscdata.fPar1           = (mscdata.fLambda1 - lambda1) / (mscdata.fLambda1 * mscdata.fTheTrueStepLenght); // alpha
    mscdata.fPar2           = 1. / (mscdata.fPar1 * mscdata.fLambda1);
    mscdata.fPar3           = 1. + mscdata.fPar2;
    mscdata.fTheZPathLenght = 1. / (mscdata.fPar1 * mscdata.fPar3) *
                              (1. - Math::Pow(1. - mscdata.fPar1 * mscdata.fTheTrueStepLenght, mscdata.fPar3));
  }
  mscdata.fTheZPathLenght = std::min(mscdata.fTheZPathLenght, mscdata.fLambda1); // NOTE:
}

void GSMSCModelSimplified::ConvertGeometricToTrueLength(geant::Track *gtrack, geant::TaskData * /*td*/)
{
  // init
  MSCdata &mscdata             = fMSCdata.Data<MSCdata>(gtrack);
  mscdata.fIsEndedUpOnBoundary = false;
  // step was not defined by transportation: i.e. physics
  if (!gtrack->Boundary()) {
    //  if ( std::abs(gtrack->fStep-mscdata.fTheZPathLenght)<1.e-8) {
    return; // fTheTrueStepLenght is known because the particle went as far as we expected
  }
  // else ::
  // - set the flag that transportation was the winer so DoNothin in DOIT !!
  // - convert geom -> true by using the mean value
  mscdata.fIsEndedUpOnBoundary = true; // OR LAST STEP
  // get the geometrical step length
  mscdata.fTheZPathLenght = gtrack->GetStep();
  // was a short single scattering step
  // t = z for very small step
  if (gtrack->GetStep() < fTLimitMinfix2) {
    mscdata.fTheTrueStepLenght = gtrack->GetStep();
    // recalculation
  } else {
    double tlength = gtrack->GetStep();
    if (gtrack->GetStep() > mscdata.fLambda1 * fTauSmall) {
      if (mscdata.fPar1 < 0.) {
        tlength = -mscdata.fLambda1 * Math::Log(1. - gtrack->GetStep() / mscdata.fLambda1);
      } else {
        double dum = mscdata.fPar1 * mscdata.fPar3 * gtrack->GetStep();
        if (dum < 1.) {
          tlength = (1. - Math::Pow(1. - dum, 1. / mscdata.fPar3)) / mscdata.fPar1;
        } else {
          tlength = mscdata.fRange;
        }
      }
      if (tlength < gtrack->GetStep() || tlength > mscdata.fTheTrueStepLenght) {
        tlength = gtrack->GetStep();
      }
    }
    mscdata.fTheTrueStepLenght = tlength;
  }
}

void GSMSCModelSimplified::SampleMSC(geant::Track *gtrack, geant::TaskData *td)
{
  MSCdata &mscdata             = fMSCdata.Data<MSCdata>(gtrack);
  mscdata.fIsNoScatteringInMSC = false;
  //
  const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
      (const_cast<vecgeom::LogicalVolume *>(gtrack->GetVolume())->GetMaterialCutsPtr()));
  double kineticEnergy = gtrack->Ekin();
  double range         = mscdata.fRange;             // set in the step limit phase
  double trueStepL     = mscdata.fTheTrueStepLenght; // proposed by all other physics
  //
  // Energy loss correction (2 versions): accurate and approximate
  //
  // what will be the energy loss due to along step energy losses if the particle goes to the current physics step
  double eloss =
      kineticEnergy - ELossTableManager::Instance().GetEnergyForRestrictedRange(matCut, fParticle, range - trueStepL);

  //  if (fTheTrueStepLenght > currentRange*dtrl) {
  //    eloss = kineticEnergy-
  //      GetEnergy(particle,range-mscdata.fTheTrueStepLenght,currentCouple);
  //  } else {
  //    eloss = fTheTrueStepLenght*GetDEDX(particle,kineticEnergy,currentCouple);
  //  }

  double tau  = 0.; //    = kineticEnergy/electron_mass_c2; // where kinEnergy is the mean kinetic energy
  double tau2 = 0.; //   = tau*tau;
  double eps0 = 0.; //   = eloss/kineticEnergy0; // energy loss fraction to the begin step energy
  double epsm = 0.; //   = eloss/kineticEnergy;  // energy loss fraction to the mean step energy
  //
  // - init.
  double efEnergy = kineticEnergy;
  double efStep   = trueStepL;
  //
  double kineticEnergy0 = kineticEnergy;
  kineticEnergy -= 0.5 * eloss; // mean energy along the full step
  // other parameters for energy loss corrections
  tau  = kineticEnergy / geant::units::kElectronMassC2; // where kinEnergy is the mean kinetic energy
  tau2 = tau * tau;
  eps0 = eloss / kineticEnergy0; // energy loss fraction to the begin step energy
  epsm = eloss / kineticEnergy;  // energy loss fraction to the mean step energy

  efEnergy   = kineticEnergy * (1. - epsm * epsm * (6. + 10. * tau + 5. * tau2) / (24. * tau2 + 48. * tau + 72.));
  double dum = 0.166666 * (4. + tau * (6. + tau * (7. + tau * (4. + tau)))) * (epsm / ((tau + 1.) * (tau + 2.))) *
               (epsm / ((tau + 1.) * (tau + 2.)));
  efStep = trueStepL * (1. - dum);
  //
  // Get elastic mfp, first transport mfp, screening parameter, and G1 (were computed at pre-step step limit and stored)
  double lambel;
  double lambtr1;
  double scra;
  double g1;
  double pMCtoQ1;
  double pMCtoG2PerG1;
  ComputeParameters(matCut, efEnergy, lambel, lambtr1, scra, g1, pMCtoQ1, pMCtoG2PerG1);
  // s/lambda_el i.e. mean number of elastic scattering along the step
  double lambdan = 0.;
  if (lambel > 0.0) {
    lambdan = efStep / lambel;
  }
  if (lambdan <= 1.0e-12) {
    mscdata.fIsNoScatteringInMSC = true;
    return;
  }
  //
  double Qn1 = lambdan * g1; // 2.* lambdan *scrA*((1.+scrA)*log(1.+1./scrA)-1.);
  //
  // Sample scattering angles
  // new direction, relative to the orriginal one is in {uss,vss,wss}
  double cosTheta1 = 1.0, sinTheta1 = 0.0, cosTheta2 = 1.0, sinTheta2 = 0.0;
  double cosPhi1 = 1.0, sinPhi1 = 0.0, cosPhi2 = 1.0, sinPhi2 = 0.0;
  double uss = 0.0, vss = 0.0, wss = 1.0;
  double x_coord = 0.0, y_coord = 0.0, z_coord = 1.0;
  double u2 = 0.0, v2 = 0.0;
  // if we are above the upper grid limit with lambdaxG1=true-length/first-trans-mfp
  // => izotropic distribution: lambG1_max =7.99 but set it to 7
  double *rndArray = td->fDblArray;
  if (0.5 * Qn1 > 7.0) {
    // get 2 random numbers
    td->fRndm->uniform_array(2, rndArray);
    cosTheta1 = 1. - 2. * rndArray[0];
    sinTheta1 = Math::Sqrt((1. - cosTheta1) * (1. + cosTheta1));
    cosTheta2 = 1. - 2. * rndArray[1];
    sinTheta2 = Math::Sqrt((1. - cosTheta2) * (1. + cosTheta2));
  } else {
    // sample 2 scattering cost1, sint1, cost2 and sint2 for half path
    // backup GS angular dtr pointer (kinetic energy and delta index in case of Mott-correction)
    // if the first was an msc sampling (the same will be used if the second is also an msc step)
    GSMSCTableSimplified::GSMSCAngularDtr *gsDtr = nullptr;
    double transfPar                             = 0.;
    bool isMsc = fGSTable.Sampling(0.5 * lambdan, 0.5 * Qn1, scra, cosTheta1, sinTheta1, &gsDtr, transfPar, td, true);
    fGSTable.Sampling(0.5 * lambdan, 0.5 * Qn1, scra, cosTheta2, sinTheta2, &gsDtr, transfPar, td, !isMsc);
    if (cosTheta1 + cosTheta2 >= 2.) { // no scattering happened
      mscdata.fIsNoScatteringInMSC = true;
      return;
    }
  }
  // sample 2 azimuthal angles
  // get 2 random numbers
  td->fRndm->uniform_array(2, rndArray);
  double phi1 = geant::units::kTwoPi * rndArray[0];
  sinPhi1     = Math::Sin(phi1);
  cosPhi1     = Math::Cos(phi1);
  double phi2 = geant::units::kTwoPi * rndArray[1];
  sinPhi2     = Math::Sin(phi2);
  cosPhi2     = Math::Cos(phi2);
  // compute final direction realtive to z-dir
  u2         = sinTheta2 * cosPhi2;
  v2         = sinTheta2 * sinPhi2;
  double u2p = cosTheta1 * u2 + sinTheta1 * cosTheta2;
  uss        = u2p * cosPhi1 - v2 * sinPhi1;
  vss        = u2p * sinPhi1 + v2 * cosPhi1;
  wss        = cosTheta1 * cosTheta2 - sinTheta1 * u2;
  //
  // set the new direction proposed by msc: will be applied if the step doesn't end on boundary
  mscdata.SetNewDirectionMsc(uss, vss, wss);
  // set the fTheZPathLenght if we don't sample displacement and
  // we should do everything at the step-limit-phase before we return
  //
  // in optimized-mode if the current-safety > current-range we do not use dispalcement
  if (mscdata.fIsNoDisplace) {
    return;
  }
  //
  //
  // Compute final position
  Qn1 *= pMCtoQ1;
  // correction parameter
  double par = 1.;
  if (Qn1 < 0.7) {
    par = 1.;
  } else if (Qn1 < 7.0) {
    par = -0.031376 * Qn1 + 1.01356;
  } else {
    par = 0.79;
  }
  // Moments with energy loss correction
  // --first the uncorrected (for energy loss) values of gamma, eta, a1=a2=0.5*(1-eta), delta
  // gamma = G_2/G_1 based on G2 computed from A by using the Wentzel DCS form of G2
  double loga  = Math::Log(1.0 + 1.0 / scra);
  double gamma = 6.0 * scra * (1.0 + scra) * (loga * (1.0 + 2.0 * scra) - 2.0) / g1;
  gamma *= pMCtoG2PerG1;
  // sample eta from p(eta)=2*eta i.e. P(eta) = eta_square ;-> P(eta) = rand --> eta = sqrt(rand)
  double eta  = Math::Sqrt(td->fRndm->uniform());
  double eta1 = 0.5 * (1 - eta); // used  more than once
  // 0.5 +sqrt(6)/6 = 0.9082483;
  // 1/(4*sqrt(6))  = 0.1020621;
  // (4-sqrt(6)/(24*sqrt(6))) = 0.026374715
  // delta = 0.9082483-(0.1020621-0.0263747*gamma)*Qn1 without energy loss cor.
  double delta = 0.9082483 - (0.1020621 - 0.0263747 * gamma) * Qn1;
  //
  // compute alpha1 and alpha2 for energy loss correction
  double temp1 = 2.0 + tau;
  double temp  = (2.0 + tau * temp1) / ((tau + 1.0) * temp1);
  // take into account logarithmic dependence
  temp     = temp - (tau + 1.0) / ((tau + 2.0) * (loga * (1.0 + scra) - 1.0));
  temp     = temp * epsm;
  temp1    = 1.0 - temp;
  delta    = delta + 0.40824829 * (eps0 * (tau + 1.0) /
                                    ((tau + 2.0) * (loga * (1.0 + scra) - 1.0) * (loga * (1.0 + 2.0 * scra) - 2.0)) -
                                0.25 * temp * temp);
  double b = eta * delta;
  double c = eta * (1.0 - delta);
  //
  // calculate transport direction cosines:
  // ut,vt,wt is the final position divided by the true step length
  double w1v2 = cosTheta1 * v2;
  double ut   = b * sinTheta1 * cosPhi1 + c * (cosPhi1 * u2 - sinPhi1 * w1v2) + eta1 * uss * temp1;
  double vt   = b * sinTheta1 * sinPhi1 + c * (sinPhi1 * u2 + cosPhi1 * w1v2) + eta1 * vss * temp1;
  double wt   = eta1 * (1 + temp) + b * cosTheta1 + c * cosTheta2 + eta1 * wss * temp1;
  //
  // long step correction (only if not error-free stepping algorithm is used)
  ut *= par;
  vt *= par;
  wt *= par;
  //
  // final position relative to the pre-step point in the scattering frame
  // ut = x_f/s so needs to multiply by s
  x_coord = ut * trueStepL;
  y_coord = vt * trueStepL;
  z_coord = wt * trueStepL;
  //
  // else:: we sample in the post step point so
  //       the fTheZPathLenght was already set and was taken as transport along zet
  mscdata.SetDisplacement(x_coord, y_coord, z_coord - mscdata.fTheZPathLenght);
}

void GSMSCModelSimplified::SampleMSCp1(geant::TrackVec_t &gtracks, geant::TaskData *td)
{
  double *lambdaNArr      = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr;
  double *pMCtoQ1Arr      = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr2;
  double *pMCtoG2PerG1Arr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr3;
  double *scraArr         = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr4;
  double *qn1Arr          = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr5;
  double *g1Arr           = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr6;
  double *tauArr          = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr7;
  double *epsMArr         = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr8;
  double *eps0Arr         = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr9;

  /** not used (yet)
  double *sinTheta1Arr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr10;
  double *cosTheta1Arr = td->fPhysicsData->fPhysicsScratchpad.fR0;
  double *cosPhi1Arr   = td->fPhysicsData->fPhysicsScratchpad.fR1;
  double *sinPhi1Arr   = td->fPhysicsData->fPhysicsScratchpad.fR2;
  double *sinTheta2Arr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr11;
  double *cosTheta2Arr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr12;
  double *cosPhi2Arr   = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr13;
  double *sinPhi2Arr   = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr14;
  **/
  bool *maskArr = td->fPhysicsData->fPhysicsScratchpad.fBoolArr;

  for (size_t i = 0; i < gtracks.size(); ++i) {
    maskArr[i] = false;
  }

  for (size_t i = 0; i < gtracks.size(); i += kVecLenD) {

    Double_v kineticEnergy; // = gtrack->T();
    Double_v range;         // = mscdata.fRange;             // set in the step limit phase
    Double_v trueStepL;     // = mscdata.fTheTrueStepLenght; // proposed by all other physics
    Double_v eloss;
    std::array<int, kVecLenD> matIdx;
    for (int l = 0; l < kVecLenD; ++l) {
      auto gtrack                  = gtracks[i + l];
      MSCdata &mscdata             = fMSCdata.Data<MSCdata>(gtrack);
      mscdata.fIsNoScatteringInMSC = false;
      //
      const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
          (const_cast<vecgeom::LogicalVolume *>(gtrack->GetVolume())->GetMaterialCutsPtr()));
      matIdx[l] = matCut->GetMaterial()->GetIndex();
      Set(kineticEnergy, l, gtrack->Ekin());
      Set(range, l, mscdata.fRange);                 // set in the step limit phase
      Set(trueStepL, l, mscdata.fTheTrueStepLenght); // proposed by all other physics
      //
      // Energy loss correction (2 versions): accurate and approximate
      //
      // what will be the energy loss due to along step energy losses if the particle goes to the current physics step
      Set(eloss, l,
          Get(kineticEnergy, l) -
              ELossTableManager::Instance().GetEnergyForRestrictedRange(matCut, fParticle, Get(range - trueStepL, l)));
    }

    //  if (fTheTrueStepLenght > currentRange*dtrl) {
    //    eloss = kineticEnergy-
    //      GetEnergy(particle,range-mscdata.fTheTrueStepLenght,currentCouple);
    //  } else {
    //    eloss = fTheTrueStepLenght*GetDEDX(particle,kineticEnergy,currentCouple);
    //  }

    Double_v tau  = 0.; //    = kineticEnergy/electron_mass_c2; // where kinEnergy is the mean kinetic energy
    Double_v tau2 = 0.; //   = tau*tau;
    Double_v eps0 = 0.; //   = eloss/kineticEnergy0; // energy loss fraction to the begin step energy
    Double_v epsm = 0.; //   = eloss/kineticEnergy;  // energy loss fraction to the mean step energy
    //
    // - init.
    Double_v efEnergy = kineticEnergy;
    Double_v efStep   = trueStepL;
    //
    Double_v kineticEnergy0 = kineticEnergy;
    kineticEnergy -= 0.5 * eloss; // mean energy along the full step
    // other parameters for energy loss corrections
    tau  = kineticEnergy / geant::units::kElectronMassC2; // where kinEnergy is the mean kinetic energy
    tau2 = tau * tau;
    eps0 = eloss / kineticEnergy0; // energy loss fraction to the begin step energy
    epsm = eloss / kineticEnergy;  // energy loss fraction to the mean step energy

    efEnergy     = kineticEnergy * (1. - epsm * epsm * (6. + 10. * tau + 5. * tau2) / (24. * tau2 + 48. * tau + 72.));
    Double_v dum = 0.166666 * (4. + tau * (6. + tau * (7. + tau * (4. + tau)))) * (epsm / ((tau + 1.) * (tau + 2.))) *
                   (epsm / ((tau + 1.) * (tau + 2.)));
    efStep = trueStepL * (1. - dum);
    //
    // Get elastic mfp, first transport mfp, screening parameter, and G1 (were computed at pre-step step limit and
    // stored)
    Double_v lambel;
    Double_v lambtr1;
    Double_v scra;
    Double_v g1;
    Double_v pMCtoQ1;
    Double_v pMCtoG2PerG1;
    ComputeParameters(matIdx, efEnergy, lambel, lambtr1, scra, g1, pMCtoQ1, pMCtoG2PerG1);
    // s/lambda_el i.e. mean number of elastic scattering along the step
    Double_v lambdan = 0.;
    vecCore::MaskedAssign(lambdan, lambel > 0.0, efStep / lambel);
    MaskD_v smallL = lambdan <= 1.0e-12;
    if (!MaskEmpty(smallL)) {
      for (int l = 0; l < kVecLenD; ++l) {
        auto gtrack                  = gtracks[i + l];
        MSCdata &mscdata             = fMSCdata.Data<MSCdata>(gtrack);
        mscdata.fIsNoScatteringInMSC = true;
        maskArr[i + l]               = true;
      }
    }
    //
    Double_v Qn1 = lambdan * g1; // 2.* lambdan *scrA*((1.+scrA)*log(1.+1./scrA)-1.);

    vecCore::Store(lambdan, lambdaNArr + i);
    vecCore::Store(pMCtoQ1, pMCtoQ1Arr + i);
    vecCore::Store(pMCtoG2PerG1, pMCtoG2PerG1Arr + i);
    vecCore::Store(scra, scraArr + i);
    vecCore::Store(Qn1, qn1Arr + i);
    vecCore::Store(g1, g1Arr + i);
    vecCore::Store(tau, tauArr + i);
    vecCore::Store(epsm, epsMArr + i);
    vecCore::Store(eps0, eps0Arr + i);
  }
}

void GSMSCModelSimplified::SampleMSCp2(geant::TrackVec_t &gtracks, geant::TaskData *td)
{
  double *lambdaNArr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr;
  // double *pMCtoQ1Arr      = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr2;
  // double *pMCtoG2PerG1Arr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr3;
  double *scraArr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr4;
  double *qn1Arr  = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr5;
  // double *g1Arr           = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr6;
  // double *tauArr          = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr7;
  // double *epsMArr         = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr8;
  // double *eps0Arr         = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr9;

  double *sinTheta1Arr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr10;
  double *cosTheta1Arr = td->fPhysicsData->fPhysicsScratchpad.fR0;
  double *cosPhi1Arr   = td->fPhysicsData->fPhysicsScratchpad.fR1;
  double *sinPhi1Arr   = td->fPhysicsData->fPhysicsScratchpad.fR2;
  double *sinTheta2Arr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr11;
  double *cosTheta2Arr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr12;
  double *cosPhi2Arr   = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr13;
  double *sinPhi2Arr   = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr14;

  bool *maskArr = td->fPhysicsData->fPhysicsScratchpad.fBoolArr;

  fGSTable.SampleTheta12(lambdaNArr, qn1Arr, scraArr, cosTheta1Arr, cosTheta2Arr, gtracks.size(), td);
  for (size_t i = 0; i < gtracks.size(); i += kVecLenD) {
    Double_v cost1;
    vecCore::Load(cost1, cosTheta1Arr + i);
    Double_v cost2;
    vecCore::Load(cost2, cosTheta2Arr + i);
    Double_v sinTheta1 = Math::Sqrt((1. - cost1) * (1. + cost1));
    Double_v sinTheta2 = Math::Sqrt((1. - cost2) * (1. + cost2));

    for (int l = 0; l < kVecLenD; ++l) {
      if (Get(cost1, l) + Get(cost2, l) >= 2.) { // no scattering happened
        if (!maskArr[i + l]) {
          auto gtrack                  = gtracks[i + l];
          MSCdata &mscdata             = fMSCdata.Data<MSCdata>(gtrack);
          mscdata.fIsNoScatteringInMSC = true;
          maskArr[i + l]               = true;
        }
      }
    }

    vecCore::Store(sinTheta1, sinTheta1Arr + i);
    vecCore::Store(sinTheta2, sinTheta2Arr + i);
  }

  for (size_t i = 0; i < gtracks.size(); i += kVecLenD) {
    Double_v cosPhi1 = 1.0, sinPhi1 = 0.0, cosPhi2 = 1.0, sinPhi2 = 0.0;
    Double_v uss = 0.0, vss = 0.0, wss = 1.0;
    Double_v u2 = 0.0, v2 = 0.0;
    Double_v sinTheta2;
    vecCore::Load(sinTheta2, sinTheta2Arr + i);
    Double_v cosTheta2;
    vecCore::Load(cosTheta2, cosTheta2Arr + i);
    Double_v sinTheta1;
    vecCore::Load(sinTheta1, sinTheta1Arr + i);
    Double_v cosTheta1;
    vecCore::Load(cosTheta1, cosTheta1Arr + i);
    // sample 2 azimuthal angles
    // get 2 random numbers
    Double_v phi1 = geant::units::kTwoPi * td->fRndm->uniformV();
    Math::SinCos(phi1, sinPhi1, cosPhi1);
    Double_v phi2 = geant::units::kTwoPi * td->fRndm->uniformV();
    Math::SinCos(phi2, sinPhi2, cosPhi2);
    // compute final direction realtive to z-dir
    u2           = sinTheta2 * cosPhi2;
    v2           = sinTheta2 * sinPhi2;
    Double_v u2p = cosTheta1 * u2 + sinTheta1 * cosTheta2;
    uss          = u2p * cosPhi1 - v2 * sinPhi1;
    vss          = u2p * sinPhi1 + v2 * cosPhi1;
    wss          = cosTheta1 * cosTheta2 - sinTheta1 * u2;

    for (int l = 0; l < kVecLenD; ++l) {
      auto gtrack = gtracks[i + l];

      MSCdata &mscdata = fMSCdata.Data<MSCdata>(gtrack);
      //
      // set the new direction proposed by msc: will be applied if the step doesn't end on boundary
      if (!maskArr[i + l]) {
        mscdata.SetNewDirectionMsc(Get(uss, l), Get(vss, l), Get(wss, l));
      }
      // set the fTheZPathLenght if we don't sample displacement and
      // we should do everything at the step-limit-phase before we return
      //
      // in optimized-mode if the current-safety > current-range we do not use dispalcement
      if (mscdata.fIsNoDisplace) {
        maskArr[i + l] = true;
      }
    }

    vecCore::Store(sinPhi1, sinPhi1Arr + i);
    vecCore::Store(cosPhi1, cosPhi1Arr + i);
    vecCore::Store(sinPhi2, sinPhi2Arr + i);
    vecCore::Store(cosPhi2, cosPhi2Arr + i);
  }
}

void GSMSCModelSimplified::SampleMSCp3(geant::TrackVec_t &gtracks, geant::TaskData *td)
{
  // double *lambdaNArr      = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr;
  double *pMCtoQ1Arr      = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr2;
  double *pMCtoG2PerG1Arr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr3;
  double *scraArr         = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr4;
  double *qn1Arr          = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr5;
  double *g1Arr           = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr6;
  double *tauArr          = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr7;
  double *epsMArr         = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr8;
  double *eps0Arr         = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr9;

  double *sinTheta1Arr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr10;
  double *cosTheta1Arr = td->fPhysicsData->fPhysicsScratchpad.fR0;
  double *cosPhi1Arr   = td->fPhysicsData->fPhysicsScratchpad.fR1;
  double *sinPhi1Arr   = td->fPhysicsData->fPhysicsScratchpad.fR2;
  double *sinTheta2Arr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr11;
  double *cosTheta2Arr = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr12;
  double *cosPhi2Arr   = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr13;
  double *sinPhi2Arr   = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr14;

  bool *maskArr = td->fPhysicsData->fPhysicsScratchpad.fBoolArr;

  for (size_t i = 0; i < gtracks.size(); i += kVecLenD) {
    Double_v Qn1;
    vecCore::Load(Qn1, qn1Arr + i);
    Double_v pMCtoQ1;
    vecCore::Load(pMCtoQ1, pMCtoQ1Arr + i);
    Double_v pMCtoG2PerG1;
    vecCore::Load(pMCtoG2PerG1, pMCtoG2PerG1Arr + i);
    Double_v scra;
    vecCore::Load(scra, scraArr + i);
    Double_v g1;
    vecCore::Load(g1, g1Arr + i);
    Double_v tau;
    vecCore::Load(tau, tauArr + i);
    Double_v eps0;
    vecCore::Load(eps0, eps0Arr + i);
    Double_v epsm;
    vecCore::Load(epsm, epsMArr + i);
    // Compute final position
    Qn1 *= pMCtoQ1;
    // correction parameter
    Double_v par = 1.;
    vecCore::MaskedAssign(par, (Qn1 > 0.7) && (Qn1 < 7.0), -0.031376 * Qn1 + 1.01356);
    vecCore::MaskedAssign(par, (Qn1 >= 7.0), (Double_v)0.79);

    // Moments with energy loss correction
    // --first the uncorrected (for energy loss) values of gamma, eta, a1=a2=0.5*(1-eta), delta
    // gamma = G_2/G_1 based on G2 computed from A by using the Wentzel DCS form of G2
    Double_v loga  = Math::Log(1.0 + 1.0 / scra);
    Double_v gamma = 6.0 * scra * (1.0 + scra) * (loga * (1.0 + 2.0 * scra) - 2.0) / g1;
    gamma *= pMCtoG2PerG1;
    // sample eta from p(eta)=2*eta i.e. P(eta) = eta_square ;-> P(eta) = rand --> eta = sqrt(rand)
    Double_v eta  = Math::Sqrt(td->fRndm->uniformV());
    Double_v eta1 = 0.5 * (1 - eta); // used  more than once
    // 0.5 +sqrt(6)/6 = 0.9082483;
    // 1/(4*sqrt(6))  = 0.1020621;
    // (4-sqrt(6)/(24*sqrt(6))) = 0.026374715
    // delta = 0.9082483-(0.1020621-0.0263747*gamma)*Qn1 without energy loss cor.
    Double_v delta = 0.9082483 - (0.1020621 - 0.0263747 * gamma) * Qn1;
    //
    // compute alpha1 and alpha2 for energy loss correction
    Double_v temp1 = 2.0 + tau;
    Double_v temp  = (2.0 + tau * temp1) / ((tau + 1.0) * temp1);
    // take into account logarithmic dependence
    temp       = temp - (tau + 1.0) / ((tau + 2.0) * (loga * (1.0 + scra) - 1.0));
    temp       = temp * epsm;
    temp1      = 1.0 - temp;
    delta      = delta + 0.40824829 * (eps0 * (tau + 1.0) /
                                      ((tau + 2.0) * (loga * (1.0 + scra) - 1.0) * (loga * (1.0 + 2.0 * scra) - 2.0)) -
                                  0.25 * temp * temp);
    Double_v b = eta * delta;
    Double_v c = eta * (1.0 - delta);
    //
    // calculate transport direction cosines:
    // ut,vt,wt is the final position divided by the true step length

    Double_v cosTheta1;
    vecCore::Load(cosTheta1, cosTheta1Arr + i);
    Double_v sinTheta1;
    vecCore::Load(sinTheta1, sinTheta1Arr + i);
    Double_v cosTheta2;
    vecCore::Load(cosTheta2, cosTheta2Arr + i);
    Double_v sinTheta2;
    vecCore::Load(sinTheta2, sinTheta2Arr + i);
    Double_v cosPhi1;
    vecCore::Load(cosPhi1, cosPhi1Arr + i);
    Double_v sinPhi1;
    vecCore::Load(sinPhi1, sinPhi1Arr + i);
    Double_v cosPhi2;
    vecCore::Load(cosPhi2, cosPhi2Arr + i);
    Double_v sinPhi2;
    vecCore::Load(sinPhi2, sinPhi2Arr + i);

    Double_v u2  = sinTheta2 * cosPhi2;
    Double_v v2  = sinTheta2 * sinPhi2;
    Double_v u2p = cosTheta1 * u2 + sinTheta1 * cosTheta2;
    Double_v uss = u2p * cosPhi1 - v2 * sinPhi1;
    Double_v vss = u2p * sinPhi1 + v2 * cosPhi1;
    Double_v wss = cosTheta1 * cosTheta2 - sinTheta1 * u2;

    Double_v w1v2 = cosTheta1 * v2;
    Double_v ut   = b * sinTheta1 * cosPhi1 + c * (cosPhi1 * u2 - sinPhi1 * w1v2) + eta1 * uss * temp1;
    Double_v vt   = b * sinTheta1 * sinPhi1 + c * (sinPhi1 * u2 + cosPhi1 * w1v2) + eta1 * vss * temp1;
    Double_v wt   = eta1 * (1 + temp) + b * cosTheta1 + c * cosTheta2 + eta1 * wss * temp1;
    //
    // long step correction (only if not error-free stepping algorithm is used)
    ut *= par;
    vt *= par;
    wt *= par;

    //
    for (int l = 0; l < kVecLenD; ++l) {
      auto gtrack      = gtracks[i + l];
      MSCdata &mscdata = fMSCdata.Data<MSCdata>(gtrack);
      // final position relative to the pre-step point in the scattering frame
      // ut = x_f/s so needs to multiply by s
      double trueStepL = mscdata.fTheTrueStepLenght; // proposed by all other physics
      double x_coord   = Get(ut, l) * trueStepL;
      double y_coord   = Get(vt, l) * trueStepL;
      double z_coord   = Get(wt, l) * trueStepL;
      // else:: we sample in the post step point so
      //       the fTheZPathLenght was already set and was taken as transport along zet
      if (!maskArr[i + l]) {
        mscdata.SetDisplacement(x_coord, y_coord, z_coord - mscdata.fTheZPathLenght);
      }
    }
  }
}

void GSMSCModelSimplified::SampleMSC(geant::TrackVec_t &gtracks, geant::TaskData *td)
{
  assert(gtracks.size() != 0);
  SampleMSCp1(gtracks, td);
  SampleMSCp2(gtracks, td);
  SampleMSCp3(gtracks, td);
}

//
// computes the elastic and first transport mfp, the screening parameter and the first transport coeficient values
// GetTransportMeanFreePath
void GSMSCModelSimplified::ComputeParameters(const MaterialCuts *matcut, double ekin, double &lambel, double &lambtr1,
                                             double &scra, double &g1, double &pMCtoQ1, double &pMCtoG2PerG1)
{
  const Material *mat = matcut->GetMaterial();
  lambel              = 0.0; // elastic mean free path
  lambtr1             = 0.0; // first transport mean free path
  scra                = 0.0; // screening parameter
  g1                  = 0.0; // first transport coef.
  pMCtoQ1             = 1.0;
  pMCtoG2PerG1        = 1.0;
  //
  // use Moliere's screening (with Mott-corretion if it was requested)
  if (ekin < 10. * geant::units::eV) ekin = 10. * geant::units::eV;
  // total mometum square
  double pt2 = ekin * (ekin + 2.0 * geant::units::kElectronMassC2);
  // beta square
  double beta2 = pt2 / (pt2 + geant::units::kElectronMassC2 * geant::units::kElectronMassC2);
  // current material index
  int matindx = mat->GetIndex();
  // Moliere's b_c
  double bc = fGSTable.GetMoliereBc(matindx);
  // get the Mott-correcton factors if Mott-correcton was requested by the user
  double pMCtoScrA = 1.0;
  double scpCor    = 1.0;
  fPWACorrection.GetPWACorrectionFactors(Math::Log(ekin), beta2, matindx, pMCtoScrA, pMCtoQ1, pMCtoG2PerG1);
  // screening parameter:
  // - if Mott-corretioncorrection: the Screened-Rutherford times Mott-corretion DCS with this
  //   screening parameter gives back the (elsepa) PWA first transport cross section
  // - if PWA correction: he Screened-Rutherford DCS with this screening parameter
  //   gives back the (elsepa) PWA first transport cross section
  scra = fGSTable.GetMoliereXc2(matindx) / (4.0 * pt2 * bc) * pMCtoScrA;
  // elastic mean free path in Geant4 internal lenght units: the neglected (1+screening parameter) term is corrected
  // (if Mott-corretion: the corrected screening parameter is used for this (1+A) correction + Moliere b_c is also
  // corrected with the screening parameter correction)
  lambel = beta2 * (1. + scra) * pMCtoScrA / bc / scpCor;
  // first transport coefficient (if Mott-corretion: the corrected screening parameter is used (it will be fully
  // consistent with the one used during the pre-computation of the Mott-correted GS angular distributions))
  g1 = 2.0 * scra * ((1.0 + scra) * Math::Log(1.0 / scra + 1.0) - 1.0);
  // first transport mean free path
  lambtr1 = lambel / g1;
}

void GSMSCModelSimplified::ComputeParameters(std::array<int, geant::kVecLenD> matIdx, geant::Double_v ekin,
                                             geant::Double_v &lambel, geant::Double_v &lambtr1, geant::Double_v &scra,
                                             geant::Double_v &g1, geant::Double_v &pMCtoQ1,
                                             geant::Double_v &pMCtoG2PerG1)
{
  lambel       = 0.0; // elastic mean free path
  lambtr1      = 0.0; // first transport mean free path
  scra         = 0.0; // screening parameter
  g1           = 0.0; // first transport coef.
  pMCtoQ1      = 1.0;
  pMCtoG2PerG1 = 1.0;
  //
  ekin = Math::Max(ekin, (Double_v)10 * geant::units::eV);
  // total mometum square
  Double_v pt2 = ekin * (ekin + 2.0 * geant::units::kElectronMassC2);
  // beta square
  Double_v beta2 = pt2 / (pt2 + geant::units::kElectronMassC2 * geant::units::kElectronMassC2);
  // current material index
  // Moliere's b_c
  Double_v bc, xc2;
  for (int l = 0; l < kVecLenD; ++l) {
    Set(bc, l, fGSTable.GetMoliereBc(matIdx[l]));
    Set(xc2, l, fGSTable.GetMoliereXc2(matIdx[l]));
  }
  // get the Mott-correcton factors if Mott-correcton was requested by the user
  Double_v pMCtoScrA = 1.0;
  Double_v scpCor    = 1.0;
  Double_v lekin     = Math::Log(ekin);
  for (int l = 0; l < kVecLenD; ++l) {
    double pMCtoScrA_tmp, pMCtoQ1_tmp, pMCtoG2PerG1_tmp;
    fPWACorrection.GetPWACorrectionFactors(Get(lekin, l), Get(beta2, l), matIdx[l], pMCtoScrA_tmp, pMCtoQ1_tmp,
                                           pMCtoG2PerG1_tmp);
    Set(pMCtoScrA, l, pMCtoScrA_tmp);
    Set(pMCtoQ1, l, pMCtoQ1_tmp);
    Set(pMCtoG2PerG1, l, pMCtoG2PerG1_tmp);
  }
  // screening parameter:
  // - if Mott-corretioncorrection: the Screened-Rutherford times Mott-corretion DCS with this
  //   screening parameter gives back the (elsepa) PWA first transport cross section
  // - if PWA correction: he Screened-Rutherford DCS with this screening parameter
  //   gives back the (elsepa) PWA first transport cross section
  scra = xc2 / (4.0 * pt2 * bc) * pMCtoScrA;
  // elastic mean free path in Geant4 internal lenght units: the neglected (1+screening parameter) term is corrected
  // (if Mott-corretion: the corrected screening parameter is used for this (1+A) correction + Moliere b_c is also
  // corrected with the screening parameter correction)
  lambel = beta2 * (1. + scra) * pMCtoScrA / bc / scpCor;
  // first transport coefficient (if Mott-corretion: the corrected screening parameter is used (it will be fully
  // consistent with the one used during the pre-computation of the Mott-correted GS angular distributions))
  g1 = 2.0 * scra * ((1.0 + scra) * Math::Log(1.0 / scra + 1.0) - 1.0);
  // first transport mean free path
  lambtr1 = lambel / g1;
}

double GSMSCModelSimplified::GetTransportMeanFreePathOnly(const MaterialCuts *matcut, double ekin)
{
  const Material *mat = matcut->GetMaterial();
  //
  double lambda0 = 0.0; // elastc mean free path
  double lambda1 = 0.0; // first transport mean free path
  double scrA    = 0.0; // screening parametr
  double g1      = 0.0; // first transport mean free path
  //
  // use Moliere's screening (with Mott-corretion if it was requested)
  if (ekin < 10. * geant::units::eV) ekin = 10. * geant::units::eV;
  // total mometum square in Geant4 internal energy2 units which is MeV2
  double pt2   = ekin * (ekin + 2.0 * geant::units::kElectronMassC2);
  double beta2 = pt2 / (pt2 + geant::units::kElectronMassC2 * geant::units::kElectronMassC2);
  int matindx  = mat->GetIndex();
  double bc    = fGSTable.GetMoliereBc(matindx);
  // get the Mott-correcton factors if Mott-correcton was requested by the user
  double mctoScrA    = 1.0;
  double mctoQ1      = 1.0;
  double mctoG2PerG1 = 1.0;
  double scpCor      = 1.0;
  fPWACorrection.GetPWACorrectionFactors(Math::Log(ekin), beta2, matindx, mctoScrA, mctoQ1, mctoG2PerG1);
  // scpCor = fGSTable->ComputeScatteringPowerCorrection(matcut, ekin);
  scrA = fGSTable.GetMoliereXc2(matindx) / (4.0 * pt2 * bc) * mctoScrA;
  // total elastic mean free path in Geant4 internal lenght units
  lambda0 = beta2 * (1. + scrA) * mctoScrA / bc / scpCor;
  g1      = 2.0 * scrA * ((1.0 + scrA) * Math::Log(1.0 / scrA + 1.0) - 1.0);
  lambda1 = lambda0 / g1;
  return lambda1;
}

double GSMSCModelSimplified::RandomizeTrueStepLength(geant::TaskData *td, double tlimit)
{
  double tempTLimit = tlimit;
  do {
    tempTLimit = td->fRndm->Gauss(tlimit, 0.1 * tlimit);
  } while ((tempTLimit < 0.) || (tempTLimit > 2. * tlimit));
  return tempTLimit;
}

} // namespace geantphysics
