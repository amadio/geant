
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

namespace geantphysics {

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
  double kineticEnergy = gtrack->T();

  if (kineticEnergy < GetLowEnergyUsageLimit() || kineticEnergy > GetHighEnergyUsageLimit()) {
    return;
  }

  double range     = ELossTableManager::Instance().GetRestrictedRange(matCut, fParticle, kineticEnergy);
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
  MSCdata &mscdata = fMSCdata.Data<MSCdata>(gtrack);

  SampleMSC(gtrack, td);
  if (!mscdata.fIsNoScatteringInMSC) {
    RotateToLabFrame(mscdata.fTheNewDirectionX, mscdata.fTheNewDirectionY, mscdata.fTheNewDirectionZ, gtrack->Dx(),
                     gtrack->Dy(), gtrack->Dz());
    if (!mscdata.fIsNoDisplace) {
      RotateToLabFrame(mscdata.fTheDisplacementVectorX, mscdata.fTheDisplacementVectorY,
                       mscdata.fTheDisplacementVectorZ, gtrack->Dx(), gtrack->Dy(), gtrack->Dz());
    }
    return true; // new direction
  }
  return false;
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
  double ekin = gtrack->T();
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
  double kineticEnergy = gtrack->T();
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
  temp  = temp - (tau + 1.0) / ((tau + 2.0) * (loga * (1.0 + scra) - 1.0));
  temp  = temp * epsm;
  temp1 = 1.0 - temp;
  delta = delta +
          0.40824829 *
              (eps0 * (tau + 1.0) / ((tau + 2.0) * (loga * (1.0 + scra) - 1.0) * (loga * (1.0 + 2.0 * scra) - 2.0)) -
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
