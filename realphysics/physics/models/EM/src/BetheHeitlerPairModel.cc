
#include "BetheHeitlerPairModel.h"

#include "PhysicalConstants.h"
// for Vector_t
#include "Types.h"

#include "Material.h"
#include "MaterialProperties.h"
#include "MaterialCuts.h"
#include "Element.h"
#include "ElementProperties.h"

#include "LightTrack.h"
#include "PhysicsData.h"

#include "Particle.h"
#include "Gamma.h"
#include "Electron.h"
#include "Positron.h"

#include "AliasTable.h"
#include "GLIntegral.h"

// from geantV
#include "GeantTaskData.h"

namespace geantphysics {

std::vector<BetheHeitlerPairModel::ElementData*>  BetheHeitlerPairModel::gElementData(gMaxZet,nullptr);

BetheHeitlerPairModel::BetheHeitlerPairModel(const std::string &modelname) : EMModel(modelname) {
  fIsUseTsaisScreening              = false;

  fElectronInternalCode             = -1;      // will be set at init
  fPositronInternalCode             = -1;      // will be set at init

  fSTNumPhotonEnergiesPerDecade     = 12;    // ST=>SamplingTables
  fSTNumDiscreteEnergyTransferVals  = 54; // ST=>SamplingTables
  fSTNumPhotonEnergies              = -1;      // will be set at init.

  fSTLogMinPhotonEnergy             = -1.;     // will be set at init in case of sampling tables
  fSTILDeltaPhotonEnergy            = -1.;     // will be set at init in case of sampling tables

  fMinimumPrimaryEnergy             =  2.*geant::kElectronMassC2; // final value will be set at init.
  fGammaEneregyLimit                =  2.*geant::MeV; // use simplified sampling below this gamma energy

  fAliasSampler                     = nullptr;
}


BetheHeitlerPairModel::~BetheHeitlerPairModel() {
  // clear ElementData
  for (size_t i=0; i<gElementData.size(); ++i) {
    if (gElementData[i]) {
      delete gElementData[i];
    }
  }
  gElementData.clear();
  // clear sampling tables if any
  if (GetUseSamplingTables()) {
    ClearSamplingTables();
  }
  if (fAliasSampler) {
    delete fAliasSampler;
  }
}


void BetheHeitlerPairModel::Initialize() {
  EMModel::Initialize();  // will set the PhysicsParameters member
  fElectronInternalCode = Electron::Definition()->GetInternalCode();
  fPositronInternalCode = Positron::Definition()->GetInternalCode();
  fMinimumPrimaryEnergy = 2.*geant::kElectronMassC2; // will be used to build table in target element selector
  if (GetLowEnergyUsageLimit()>fMinimumPrimaryEnergy) {
    fMinimumPrimaryEnergy = GetLowEnergyUsageLimit();
  }
  InitialiseElementData();
  if (GetUseSamplingTables()) {
    InitSamplingTables();
  }
  // request to build target element selector for this model, for gamma particle and element selectors per material
  InitialiseElementSelectors(this,Gamma::Definition(),true);
}


double BetheHeitlerPairModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle*) {
  double xsec = 0.0;
  if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
    return xsec;
  }
  // compute the macroscopic cross section as the sum of the atomic cross sections weighted by the number of atoms in
  // in unit volume.
  const Material *mat =  matcut->GetMaterial();
  const double egamma = kinenergy;
  // we will need the element composition of this material
  const Vector_t<Element*> &theElements   = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  size_t numElems = theElements.size();
  for (size_t iel=0; iel<numElems; ++iel) {
    xsec += theAtomicNumDensityVector[iel]*ComputeAtomicCrossSection(theElements[iel]->GetZ(), egamma);
  }
  return std::max(xsec,0.0);
}


double BetheHeitlerPairModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts*, double kinenergy,
                                                     const Particle*) {
   double xsec  = 0.0;
   if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
     return xsec;
   }
   // compute the parametrized atomic cross section: depends only on target Z and gamma energy.
   xsec = ComputeAtomicCrossSection(elem->GetZ(), kinenergy);
   return std::max(xsec,0.0);
}


int BetheHeitlerPairModel::SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td) {
  int    numSecondaries      = 0;
  const double ekin          = track.GetKinE();
  const double eps0          = geant::kElectronMassC2/ekin;
  // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
  // check if kinetic energy is above fHighEnergyUsageLimit andd o nothing if yes
  if (ekin<GetLowEnergyUsageLimit() || ekin>GetHighEnergyUsageLimit() || eps0>0.5) {
    return numSecondaries;
  }
  // interaction is possible so sample target element: will be needed anyway for the direction sampling
  MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[track.GetMaterialCutCoupleIndex()];
  const Vector_t<Element*> &theElements = matCut->GetMaterial()->GetElementVector();
  double targetElemIndx = 0;
  if (theElements.size()>1) {
    targetElemIndx = SampleTargetElementIndex(matCut, ekin, td->fRndm->uniform());
  }
  const double zet = theElements[targetElemIndx]->GetZ();
  //
  // sample the reduced total energy transfered to one of the secondary particles i.e. either e- or e+
  double eps = 0.0;
  if (ekin<fGammaEneregyLimit) {   // sample eps from uniform on [eps_0,0.5]
    eps = eps0 + (0.5-eps0)*td->fRndm->uniform();
  } else {
    const int izet = std::min(std::lrint(zet),gMaxZet-1);
    // sample the reduced total energy transferd to one of the e-/e+ pair by using tables or rejection
    if (GetUseSamplingTables()) { // sample eps using the sampling tables
      double *rndArray  = td->fDblArray;
      td->fRndm->uniform_array(4, rndArray);
      eps = SampleTotalEnergyTransfer(ekin, izet, rndArray[0], rndArray[1], rndArray[2]);
    } else { // sample eps by rejection
      eps = SampleTotalEnergyTransfer(ekin, izet, td);
    }
  }
  //
  // create the secondary partcicles:
  // 1. the total elengy of e-/e+
  double electronTotE;
  double positronTotE;
  if (td->fRndm->uniform()>0.5) {
    electronTotE = (1.-eps)*ekin;
    positronTotE = eps*ekin;
  } else {
    electronTotE = eps*ekin;
    positronTotE = (1.-eps)*ekin;
  }
  //
  // 2. sample the direction: theta is sampled based on Laszlo's approximation to Tsai-s dcs
  // note: we should investigate the possibility to use the leading term of the dcs and sample cos(theta(-+))
  double *rndArray  = td->fDblArray;
  td->fRndm->uniform_array(4, rndArray);
  double uvar = -std::log(rndArray[0]*rndArray[1]);
  if (9.>36.*rndArray[2]) {
    uvar *= 1.6;
  } else {
    uvar *= 0.53333;
  }
  const double thetaElectron = uvar*geant::kElectronMassC2/electronTotE;
  const double sintEle       = std::sin(thetaElectron);
  const double thetaPositron = uvar*geant::kElectronMassC2/positronTotE;
  const double sintPos       = -std::sin(thetaPositron);
  const double phi           = geant::kTwoPi*rndArray[3];
  const double sinphi        = std::sin(phi);
  const double cosphi        = std::cos(phi);
  // e- direction
  double eleDirX = sintEle*cosphi;
  double eleDirY = sintEle*sinphi;
  double eleDirZ = std::cos(thetaElectron);
  // e+ direction
  double posDirX = sintPos*cosphi;
  double posDirY = sintPos*sinphi;
  double posDirZ = std::cos(thetaPositron);
  //
  // 3. kill the primary photon and create the secondaries
  track.SetKinE(0.0);
  track.SetTrackStatus(LTrackStatus::kKill);
  // 4. compute kinetic energy of e-/e+
  const double ekinElectron = std::max((electronTotE-geant::kElectronMassC2),0.);
  const double ekinPositron = std::max((positronTotE-geant::kElectronMassC2),0.);
  // 5. rotate direction back to the lab frame: current directions are relative to the photon dir as z-dir
  RotateToLabFrame(eleDirX, eleDirY, eleDirZ, track.GetDirX(), track.GetDirY(), track.GetDirZ());
  RotateToLabFrame(posDirX, posDirY, posDirZ, track.GetDirX(), track.GetDirY(), track.GetDirZ());
  //
  // 6. insert the secondary e-/e+ into the secondary list:
  numSecondaries = 2;
  // current capacity of secondary track container
  int curSizeOfSecList = td->fPhysicsData->GetSizeListOfSecondaries();
  // currently used secondary tracks in the container
  int curNumUsedSecs   = td->fPhysicsData->GetNumUsedSecondaries();
  if (curSizeOfSecList-curNumUsedSecs<numSecondaries) {
    td->fPhysicsData->SetSizeListOfSecondaries(2*curSizeOfSecList);
  }
  int secIndx     = curNumUsedSecs;
  curNumUsedSecs += numSecondaries;
  td->fPhysicsData->SetNumUsedSecondaries(curNumUsedSecs);
  std::vector<LightTrack>& sectracks = td->fPhysicsData->GetListOfSecondaries();
  // first set the e-
  sectracks[secIndx].SetDirX(eleDirX);
  sectracks[secIndx].SetDirY(eleDirY);
  sectracks[secIndx].SetDirZ(eleDirZ);
  sectracks[secIndx].SetKinE(ekinElectron);
  sectracks[secIndx].SetGVcode(fElectronInternalCode);
  sectracks[secIndx].SetMass(geant::kElectronMassC2);
  sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); // parent GeantTrack index
  // then set the e+
  ++secIndx;
  sectracks[secIndx].SetDirX(posDirX);
  sectracks[secIndx].SetDirY(posDirY);
  sectracks[secIndx].SetDirZ(posDirZ);
  sectracks[secIndx].SetKinE(ekinPositron);
  sectracks[secIndx].SetGVcode(fPositronInternalCode);
  sectracks[secIndx].SetMass(geant::kElectronMassC2);
  sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); // parent GeantTrack index

  return numSecondaries;
}



/**
 * @internal
 * The method computes atomic cross section of converion of gamma particle into e-/e+ pair based on the Geant4
 * \cite cirrone2010validation \cite agostinelli2003geant4 parametrization.
 *
 * According to the Geant4 documentation \cite g4physref, the numerical atomic cross sections \f$ \sigma(Z,E_{\gamma}) \f$
 * as a function of the target atomic numberg \f$ Z \f$ and photon energy \f$ E_{\gamma}\f$ given in \cite hubbell1980pair
 * where approximated by the following function
 * \f[
 *    \sigma(Z,E_{\gamma}) \approx \tilde{\sigma}(Z,E_{\gamma}) \equiv Z(Z+1) \left[
 *           + F_1(\kappa) + Z F_2(\kappa) + \frac{F_3(\kappa) }{Z}
      \right]
 * \f]
 * where \f$ \kappa = \ln[E_{\gamma}/(m_ec^2)] \f$ is the logarithm of the primary gamma energy in electron rest mass units
 * and \f$ F_i(\kappa) \equiv \sum_{j=0}^{5} a_{ij} \kappa^j\f$ with the parameter values
 * \f[
 *   \begin{array}{lcrr}
 *    a_{10} & = & +8.7842e+2  & \text{[barn]} \\
 *    a_{11} & = & -1.9625e+3  & \text{[barn]} \\
 *    a_{12} & = & +1.2949e+3  & \text{[barn]} \\
 *    a_{13} & = & -2.0028e+2  & \text{[barn]} \\
 *    a_{14} & = & +1.2575e+1  & \text{[barn]} \\
 *    a_{15} & = & -2.8333e-1  & \text{[barn]} \\
 *    &&&\\
 *    a_{20} & = & -1.0342e+1  & \text{[barn]} \\
 *    a_{21} & = & +1.7692e+1  & \text{[barn]} \\
 *    a_{22} & = & -8.2381e+0 & \text{[barn]} \\
 *    a_{23} & = & +1.3063e+0  & \text{[barn]} \\
 *    a_{24} & = & -9.0815e-2  & \text{[barn]} \\
 *    a_{25} & = & +2.3586e-3  & \text{[barn]} \\
 *    &&&\\
 *    a_{20} & = & -4.5263e+2  & \text{[barn]} \\
 *    a_{21} & = & +1.1161e+3  & \text{[barn]} \\
 *    a_{22} & = & -8.6749e+2 & \text{[barn]} \\
 *    a_{23} & = & +2.1773e+2  & \text{[barn]} \\
 *    a_{24} & = & -2.0467e+1  & \text{[barn]} \\
 *    a_{25} & = & +6.5372e-1  & \text{[barn]} \\
 *   \end{array}
 * \f]
 * that were determined through a fit of \f$ \tilde{\sigma}(Z,E_{\gamma}) \f$ to the numerical values
 * \f$ \sigma(Z,E_{\gamma}) \f$ given in \cite hubbell1980pair . Note, that these numerical cross sections include pair
 * production both in the nuclear (with screening, Coulomb and radiative corrections) and electron fields (with
 * approximate screening, radiative, exchange and retardation effects).
 *
 * Accodring to \cite cirrone2010validation, the accuracy of the above parametrization is estimated to be \f$ 5 % \f$
 * with a mean value of \f$ 2.2% \f$ for \f$ 1 \leq Z \leq 100 \f$ and
 * \f$ E_{\gamma} \in [E_{\gamma}^{\text{low}}\equiv 1.5 \text{[MeV]}, 100 \text{[GeV]}] \f$. The extrapolation:
 * \f[
 *    \sigma(E_{\gamma}) = \sigma(E_{\gamma}^{\text{low}})
 *                      \left[ \frac{E_{\gamma}-2m_ec^2}{E_{\gamma}^{\text{low}}-2m_ec^2} \right]^2
 * \f]
 *  is used if \f$ E_{\gamma} \leq E_{\gamma}^{\text{low}} \f$.
 *
 * @endinternal
 */
double BetheHeitlerPairModel::ComputeAtomicCrossSection(double z, double egamma) {
  double xsec = 0.0;
  if (z<0.9 || egamma<=2.0*geant::kElectronMassC2) {
    return xsec;
  }
  // limit for low energy extrapolation
  constexpr double egammaLimit = 1.5*geant::MeV;
  // parameter values
  // a
  constexpr double a0 =  8.7842e+2*geant::microbarn;
  constexpr double a1 = -1.9625e+3*geant::microbarn;
  constexpr double a2 =  1.2949e+3*geant::microbarn;
  constexpr double a3 = -2.0028e+2*geant::microbarn;
  constexpr double a4 =  1.2575e+1*geant::microbarn;
  constexpr double a5 = -2.8333e-1*geant::microbarn;
  // b
  constexpr double b0 = -1.0342e+1*geant::microbarn;
  constexpr double b1 =  1.7692e+1*geant::microbarn;
  constexpr double b2 = -8.2381   *geant::microbarn;
  constexpr double b3 =  1.3063   *geant::microbarn;
  constexpr double b4 = -9.0815e-2*geant::microbarn;
  constexpr double b5 =  2.3586e-3*geant::microbarn;
  // c
  constexpr double c0 = -4.5263e+2*geant::microbarn;
  constexpr double c1 =  1.1161e+3*geant::microbarn;
  constexpr double c2 = -8.6749e+2*geant::microbarn;
  constexpr double c3 =  2.1773e+2*geant::microbarn;
  constexpr double c4 = -2.0467e+1*geant::microbarn;
  constexpr double c5 =  6.5372e-1*geant::microbarn;
  //
  double egammaOrg = egamma;
  if (egamma<egammaLimit) {
    egamma = egammaLimit;
  }
  // log photon energy in electron rest mass units
  const double kappa  = std::log(egamma/geant::kElectronMassC2);
  const double kappa2 = kappa*kappa;
  const double kappa3 = kappa2*kappa;
  const double kappa4 = kappa2*kappa2;
  const double kappa5 = kappa4*kappa;
  //
  const double F1 = a0 + a1*kappa + a2*kappa2 + a3*kappa3 + a4*kappa4 + a5*kappa5;
  const double F2 = b0 + b1*kappa + b2*kappa2 + b3*kappa3 + b4*kappa4 + b5*kappa5;
  const double F3 = c0 + c1*kappa + c2*kappa2 + c3*kappa3 + c4*kappa4 + c5*kappa5;
  // compute cross section
  xsec = (z+1.)*(F1*z+F2*z*z+F3);
  // low energy correction
  if (egammaOrg<egammaLimit) {
    const double dum = (egammaOrg-2.*geant::kElectronMassC2)/(egammaLimit-2.*geant::kElectronMassC2);
    xsec *= dum*dum;
  }
  // protection against negative values
  return std::max(xsec, 0.);
}

/**
 * @internal
 *  One ElementData structure will be created for each target atom that the model needs to respond (i.e.
 *  for each elements that appears in material that belongs to a region in which the model is active).
 *  Pointers to these data structures are stored in the gElementData array indexed by atomic number (Z).
 *  Those indices, that corresponds to atomic numbers that the model do not need to respond, will remain
 *  nullptr-s.
 * @endinternal
 */
void BetheHeitlerPairModel::InitialiseElementData() {
  size_t numMatCuts = MaterialCuts::GetTheMaterialCutsTable().size();
  // get list of active region
  std::vector<bool> &isActiveInRegion = GetListActiveRegions();
  for (size_t imc=0; imc<numMatCuts; ++imc) {
    const MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[imc];
    // if this MaterialCuts belongs to a region where this model is active:
    if (isActiveInRegion[matCut->GetRegionIndex()]) {
      // get the list of elements
      const Vector_t<Element*> &theElemVect = matCut->GetMaterial()->GetElementVector();
      size_t numElems = theElemVect.size();
      for (size_t ie=0; ie<numElems; ++ie) {
        const Element *elem = theElemVect[ie];
        double zet = elem->GetZ();
        int   izet = std::min(std::lrint(zet),gMaxZet-1);
        if (!gElementData[izet]) {
          ElementData *elemData       = new ElementData();
          elemData->fDeltaFactor      = 136./elem->GetElementProperties()->GetZ13(); // 136/pow(z,1/3)
          elemData->fCoulombCor       = elem->GetElementProperties()->GetCoulombCorrection();
          elemData->fFzLow            = 8.*elem->GetElementProperties()->GetLogZ13(); //8.*std::log(z)/3.;
          elemData->fFzHigh           = elemData->fFzLow+8*elemData->fCoulombCor;
          elemData->fDeltaMaxLow      = std::exp((42.24-elemData->fFzLow)/8.368)-0.952;
          elemData->fDeltaMaxHigh     = std::exp((42.24-elemData->fFzHigh)/8.368)-0.952;
          elemData->fDeltaMaxLowTsai  = 1.36*std::sqrt( std::exp(0.5*16.863-0.25*elemData->fFzLow)-1. )/0.55846;
          elemData->fDeltaMaxHighTsai = 1.36*std::sqrt( std::exp(0.5*16.863-0.25*elemData->fFzHigh)-1. )/0.55846;
          gElementData[izet] = elemData;
        }
      }
    }
  }
}


// samples the transformed variable for the given target atom(zindx), gamma energy, and transforms it back
double BetheHeitlerPairModel::SampleTotalEnergyTransfer(const double egamma, const int izet, const double r1,
                                                         const double r2, const double r3) {
  // determine electron energy lower grid point
  const double legamma = std::log(egamma);
  //
  int indxEgamma = fSTNumPhotonEnergies-1;
  if (egamma<GetHighEnergyUsageLimit()) {
    const double val       = (legamma-fSTLogMinPhotonEnergy)*fSTILDeltaPhotonEnergy;
    indxEgamma             = (int)val;  // lower electron energy bin index
    const double pIndxHigh = val-indxEgamma;
    if (r1<pIndxHigh)
      ++indxEgamma;
  }
  // sample the transformed variable
  const RatinAliasData *als = fSamplingTables[izet]->fRatinAliasData[indxEgamma];
  // we could put a static assert here to make sure that als is not nullptr
  //
  // sample the transformed variable xi=... (such that linear aprx will be used in the first interval)
  const double  xi = fAliasSampler->SampleRatin(&(als->fXdata[0]), &(als->fCumulative[0]), &(als->fParaA[0]),
                                                &(als->fParaB[0]), &(als->fAliasW[0]), &(als->fAliasIndx[0]),
                                                fSTNumDiscreteEnergyTransferVals, r2, r3, 0);
  // transform back xi to eps = E_total_energy_transfer/E_{\gamma}
  //
  double deltaMax = gElementData[izet]->fDeltaMaxLowTsai;
  if (egamma>50.*geant::MeV) {
    deltaMax = gElementData[izet]->fDeltaMaxHighTsai;
  }
  const double eps0   = geant::kElectronMassC2/egamma;
  const double epsp   = 0.5-0.5*std::sqrt(1.-4.*eps0*gElementData[izet]->fDeltaFactor/deltaMax);
  const double epsMin = std::max(eps0,epsp);
  return epsMin*std::exp(xi*std::log(0.5/epsMin));
}


/**
 * @internal
 * The Geant4 rejection algorithm \cite g4physref is adopted. By introducing the decreasing functions of
 * \f$ \delta \equiv \delta(\epsilon) = 136 Z^{-1/3} \epsilon_0/(\epsilon(1-\epsilon)) \f$
 * \f[
 *  \begin{array}{l c l c l}
 *   F_1(\delta) & = & 3\Phi_1(\delta) - \Phi_2(\delta) - F(Z) & = &
 *     \begin{cases}
 *        42.24  - 8.368 \ln[\delta+0.952] -F(Z)     \quad & \text{if } \delta > 1 \\
 *        \\
 *        42.392 - \delta[7.796 - 1.961\delta] -F(Z) \quad & \text{otherwise}
 *     \end{cases} \\
 *   &&&&\\
 *   F_2(\delta) & = & \frac{3}{2} \Phi_1(\delta) + \frac{1}{2} \Phi_2(\delta) -F(Z) & = &
 *     \begin{cases}
 *        42.24  - 8.368 \ln[\delta+0.952] -F(Z)      \quad & \text{if } \delta > 1 \\
 *        \\
 *        41.405 - \delta[5.828 - 0.8945\delta] -F(Z) \quad & \text{otherwise}
 *     \end{cases} \\
 *  \end{array}
 * \f]
 * These functions reach their maximum values \f$ F_1^{\text{max}} = F_1(\delta_{\text{min}}),
 * F_2^{\text{max}} = F_2(\delta_{\text{min}}) \f$  at the minimum \f$ \delta \f$ value which is
 * \f$ \delta_{\text{min}} \equiv \delta(\epsilon_{\text{max}}=0.5) = 4x136 Z^{-1/3} \epsilon_0 \f$.
 * The differential cross section \f$ \frac{\mathrm{d}\sigma(Z,\epsilon)}{\mathrm{d}\epsilon} \f$ given at the
 * #ComputeDXSection() method can be written with \f$ F_1(\delta),F_2(\delta),F_1^{\text{max}}, F_2^{\text{max}} \f$ as
 *  \f[
 *    \frac{\mathrm{d}\sigma(Z,\epsilon)}{\mathrm{d}\epsilon} = \alpha r_0^2 Z[Z+\eta(Z)]\frac{2}{9}
 *       [ 0.5 -\epsilon_{\text{min}}]
 *       \left\{
 *          N_1 f_1(\epsilon) g_1(\epsilon) + N_2 f_2(\epsilon) g_2(\epsilon)
 *       \right\}
 *  \f]
 * where
 * \f[
 *  \begin{array} {l c l l c l l c l}
 *  N_1                    & \equiv & [0.5-\epsilon_{\text{min}} ] F_1^{\text{max}},\;   &
 *  f_1(\epsilon) & \equiv & \frac{3}{[0.5-\epsilon_{\text{min}}]^3} [0.5-\epsilon]^2,\; &
 *  g_1(\epsilon) & \equiv & F_1(\delta(\epsilon)) / F_1^{\text{max}} \\
 *  N_2                    & \equiv & 1.5 F_2^{\text{max}},\;   &
 *  f_2(\epsilon) & \equiv & \text{const} = [0.5-\epsilon_{\text{min}}]^{-1},\; &
 *  g_2(\epsilon) & \equiv & F_2(\delta(\epsilon)) / F_2^{\text{max}} \\
 *   \end{array}
 * \f]
 * where \f$ f_{1,2}(\epsilon)\f$ are properly normalized pdf on \f$ \epsilon \in [\epsilon_{\text{min}}, 0.5]\f$ and
 * and \f$ g_{1,2}(\epsilon) \in (0,1] \f$ are valid rejection functions.
 *
 * Having 3 uniformly distributed random numbers \f$ \{r_1,r_2,r_3 \} \f$ :
 *  - determine which decomposition is used:
 *      -# if \f$ r_1 < N_1/(N_1+N_2) \quad \text{use} \quad f_1(\epsilon)g_1(\epsilon)\f$
 *      -# and use \f$ f_2(\epsilon)g_2(\epsilon)\f$ otherwise
 *  - use \f$ r_2 \f$ to sample \f$ \epsilon \f$:
 *      -# from \f$ f_1(\epsilon) \to \epsilon = 0.5-[0.5-\epsilon_{text{min}}]r_1^{1/3}\f$ or
 *      -# from from \f$ f_2(\epsilon) \to \epsilon = \epsilon_{text{min}}+[0.5-\epsilon_{text{min}}]r_1 \f$
 *  - compute the appropriate rejection function \f$ g_1(\epsilon)\text{ or } g_2(\epsilon)\f$ and reject
 *    \f$ \epsilon \f$ if \f$ g_i(\epsilon) < r_3 \f$
 * @endinternal
 */
double BetheHeitlerPairModel::SampleTotalEnergyTransfer(const double egamma, const int izet, const Geant::GeantTaskData *td) {
    double fz  = gElementData[izet]->fFzLow;
    double deltaMax;
    if (fIsUseTsaisScreening) {
      deltaMax = gElementData[izet]->fDeltaMaxLowTsai;
      if (egamma>50.*geant::MeV) {
        fz       = gElementData[izet]->fFzHigh;
        deltaMax = gElementData[izet]->fDeltaMaxHighTsai;
      }
    } else {
      deltaMax = gElementData[izet]->fDeltaMaxLow;
      if (egamma>50.*geant::MeV) {
        fz       = gElementData[izet]->fFzHigh;
        deltaMax = gElementData[izet]->fDeltaMaxHigh;
      }
    }
    const double eps0     = geant::kElectronMassC2/egamma;
    const double deltaFac = gElementData[izet]->fDeltaFactor;
    const double deltaMin = 4.*eps0*deltaFac;
    const double eps1     = 0.5-0.5*std::sqrt(1.-deltaMin/deltaMax);
    const double epsMin   = std::max(eps0,eps1);
    const double epsRange = 0.5-epsMin;
    //
    double F10,F20;
    ScreenFunction12(F10,F20,deltaMin,fIsUseTsaisScreening);
    F10 -= fz;
    F20 -= fz;
    const double NormF1   = std::max(F10*epsRange*epsRange,0.);
    const double NormF2   = std::max(1.5*F20,0.);
    const double NormCond = NormF1/(NormF1+NormF2);
    //
    double *rndArray = td->fDblArray;
    double greject   = 0.0;
    double eps       = 0.0;
    do {
      td->fRndm->uniform_array(3, rndArray);
      if (NormCond>rndArray[0]) {
      	eps = 0.5 - epsRange*std::pow(rndArray[1],1./3.);
	      const double delta = deltaFac*eps0/(eps*(1.-eps));
	      greject = (ScreenFunction1(delta,fIsUseTsaisScreening) - fz)/F10;
      } else {
	      eps = epsMin + epsRange*rndArray[1];
        const double delta = deltaFac*eps0/(eps*(1.-eps));
	      greject = (ScreenFunction2(delta,fIsUseTsaisScreening) - fz)/F20;
      }
    } while (greject<rndArray[2]);
  return eps;
}


void BetheHeitlerPairModel::ComputeScreeningFunctions(double &phi1, double &phi2, const double delta, const bool istsai) {
  if (istsai) {
    const double gamma   = delta*0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const double gamma2  = gamma*gamma;
    phi1   = 16.863-2.0*std::log(1.0+0.311877*gamma2)+2.4*std::exp(-0.9*gamma)+1.6*std::exp(-1.5*gamma);
    phi2   = phi1-2.0/(3.0+19.5*gamma+18.0*gamma2);  // phi1-phi2
  } else {
    if (delta>1.) {
      phi1 = 21.12 - 4.184*std::log(delta+0.952);
      phi2 = phi1;
    } else {
      const double delta2 = delta*delta;
      phi1 = 20.867 - 3.242*delta + 0.625*delta2;
      phi2 = 20.209 - 1.93 *delta - 0.086*delta2;
    }
  }
}

// val1 =  3xPhi_1 - Phi_2: used in case of rejection (either istsai or not)
// val2 =  1.5*Phi_1 + 0.5*Phi_2: used in case of rejection (either istsai or not)
void BetheHeitlerPairModel::ScreenFunction12(double &val1, double &val2, const double delta, const bool istsai) {
  if (istsai) {
    const double gamma  = delta*0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const double gamma2 = gamma*gamma;
    const double dum1   = 33.726 - 4.*std::log(1.0+0.311877*gamma2) + 4.8*std::exp(-0.9*gamma) + 3.2*std::exp(-1.5*gamma);
    const double dum2   = 2./(3.+19.5*gamma+18.*gamma2);
    val1 = dum1 +     dum2;
    val2 = dum1 - 0.5*dum2;
  } else {
    if (delta>1.) {
      val1 = 42.24 - 8.368*std::log(delta+0.952);
      val2 = val1;
    } else {
      val1 = 42.392 - delta*(7.796 - 1.961*delta);
      val2 = 41.405 - delta*(5.828 - 0.8945*delta);
    }
  }
}


// 3xPhi_1 - Phi_2: used in case of rejection (either istsai or not)
double BetheHeitlerPairModel::ScreenFunction1(const double delta, const bool istsai) {
  double val;
  if (istsai) {
    const double gamma   = delta*0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const double gamma2  = gamma*gamma;
    val = 33.726 - 4.*std::log(1.0+0.311877*gamma2) + 4.8*std::exp(-0.9*gamma) + 3.2*std::exp(-1.5*gamma)
          + 2./(3.+19.5*gamma+18.*gamma2);
  } else {
    if (delta>1.) {
      val = 42.24 - 8.368*std::log(delta+0.952);
    } else {
      val = 42.392 - delta*(7.796 - 1.961*delta);
    }
  }
  return val;
}

// 1.5*Phi_1 + 0.5*Phi_2: used in case of rejection (either istsai or not)
double BetheHeitlerPairModel::ScreenFunction2(const double delta, const bool istsai) {
  double val;
  if (istsai) {
    const double gamma   = delta*0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    const double gamma2  = gamma*gamma;
    val = 33.726 - 4.*std::log(1.0+0.311877*gamma2) + 4.8*std::exp(-0.9*gamma) + 3.2*std::exp(-1.5*gamma)
          - 1./(3.+19.5*gamma+18.*gamma2);
  } else {
    if (delta>1.) {
      val = 42.24 - 8.368*std::log(delta+0.952);
    } else {
      val = 41.405 - delta*(5.828 - 0.8945*delta);
    }
  }
  return val;
}


void BetheHeitlerPairModel::ClearSamplingTables() {
  size_t num = fSamplingTables.size();
  for (size_t im=0; im<num; ++im) {
    RatinAliasDataPerElement *elData = fSamplingTables[im];
    if (elData) {
      size_t nTables = elData->fRatinAliasData.size();
      for (size_t it=0; it<nTables; ++it) {
        RatinAliasData *tb = elData->fRatinAliasData[it];
        tb->fXdata.clear();
        tb->fCumulative.clear();
        tb->fParaA.clear();
        tb->fParaB.clear();
        tb->fAliasW.clear();
        tb->fAliasIndx.clear();
        delete tb;
      }
      elData->fRatinAliasData.clear();
      delete elData;
    }
  }
  fSamplingTables.clear();
}


void BetheHeitlerPairModel::InitSamplingTables() {
  // 1. clear sampling tables if any
  ClearSamplingTables();
  // 2. generate primary gamma energy grid
  const double minEprim  = fMinimumPrimaryEnergy;
  const double maxEprim  = GetHighEnergyUsageLimit();
  fSTNumPhotonEnergies   = fSTNumPhotonEnergiesPerDecade*std::lrint(std::log10(maxEprim/minEprim))+1;
  fSTNumPhotonEnergies   = std::max(fSTNumPhotonEnergies,3);
  // set up the initial gamma energy grid
  const double delta     = std::log(maxEprim/minEprim)/(fSTNumPhotonEnergies-1.0);
  fSTLogMinPhotonEnergy  = std::log(minEprim);
  fSTILDeltaPhotonEnergy = 1./delta;
  std::vector<double> primEVect(fSTNumPhotonEnergies);
  primEVect[0]                      = minEprim;
  primEVect[fSTNumPhotonEnergies-1] = maxEprim;
  for (int i=1; i<fSTNumPhotonEnergies-1; ++i) {
    primEVect[i] = std::exp(fSTLogMinPhotonEnergy+i*delta);
  }
  // 3. create an AliasTable object
  if (fAliasSampler) {
    delete fAliasSampler;
  }
  fAliasSampler = new AliasTable();
  // 4. set up the container that stores sampling tables for all elements (init to nullptr-s)
  fSamplingTables.resize(gMaxZet-1,nullptr);
  // 5. build the sampling tables:
  // get all MaterialCuts; loop over them and if they belongs to a region where this model is active:
  //  - get the material that the MaterialCuts belongs to
  //  - loop over its elements and build tables for them if not have been built yet.
  int numMatCuts = MaterialCuts::GetTheMaterialCutsTable().size();
  // get list of active region
  std::vector<bool> &isActiveInRegion = GetListActiveRegions();
  for (int i=0; i<numMatCuts; ++i) {
    const MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[i];
    // if this MaterialCuts belongs to a region where this model is active:
    if (isActiveInRegion[matCut->GetRegionIndex()]) {
      // get the material and the element vector
      const Vector_t<Element*> &theElemVect = matCut->GetMaterial()->GetElementVector();
      size_t numElems = theElemVect.size();
      for (size_t iel=0; iel<numElems; ++iel) {
        const Element *elem = theElemVect[iel];
        const int izet      = std::min(std::lrint(elem->GetZ()),gMaxZet-1);
        // if RatinAliasDataPerElement has not been built for this element yet => build
        if (!fSamplingTables[izet]) {
          BuildSamplingTablesForElement(elem, primEVect);
        }
      }
    }
  }
  primEVect.clear();
}


void BetheHeitlerPairModel::BuildSamplingTablesForElement(const Element *elem, const std::vector<double> &primevect) {
  // prepare one RatinAliasDataPerElement structure
  RatinAliasDataPerElement *perElem = new RatinAliasDataPerElement(fSTNumPhotonEnergies);
  const int izet = std::min(std::lrint(elem->GetZ()),gMaxZet-1);
  fSamplingTables[izet] = perElem;
  // allocate an array that will be used temporary:
  double *thePdfdata = new double[fSTNumDiscreteEnergyTransferVals]();
  // loop over the photon energy grid and build one RatinAliasData structure at each photon energy for this element
  for (int i=0; i<fSTNumPhotonEnergies; ++i) {
    const double egamma = primevect[i];
    // build one table for the given: egamma, element
    BuildOneRatinAlias(egamma, izet, thePdfdata, i);
  }
  // delete the temporary array
  delete [] thePdfdata;
}


/**
 * @internal
 * This internal method is called from the #BuildSamplingTablesForElement() method at initialization for different gamma
 * energies and a given target element to prepare sampling table for run-time sampling of the transformed variable
 * \f$ \xi \f$ (see more at the #ComputeDXSection() method) at a given gamma energy and target atom pair. The pdf of the
 * transformed variable will be perpresented by #fNumSamplingEnergies discrete points. The locations of these discrete
 * sample points are determined by minimizing the approximation error comming from the discretization.
 * Tsai's screening function approximation will be used (always)!
 *
 * @endinternal
 */
void BetheHeitlerPairModel::BuildOneRatinAlias(const double egamma, const int izet, double *pdfarray, const int egammaindx) {
  // compute the theoretical minimum of the reduced total energy transfer to the e+ (or to the e-)
  const double eps0   = geant::kElectronMassC2/egamma;
  const double dFact  = gElementData[izet]->fDeltaFactor;
  const double dMin   = 4.*eps0*dFact;
  double FZ           = gElementData[izet]->fFzLow;
  double dMax         = gElementData[izet]->fDeltaMaxLowTsai;
  if (egamma>50.*geant::MeV) {
    FZ   = gElementData[izet]->fFzHigh;
    dMax = gElementData[izet]->fDeltaMaxHighTsai;
  }
  const double eps1   = 0.5-0.5*std::sqrt(1.-dMin/dMax);
  FZ *= 0.125; //this is how we use in our dxsec computation
  const double epsMin = std::max(eps1,eps0);
  // allocate one RatinAliasData structure for the current element,gamma-energy combination
  RatinAliasData *raData = new RatinAliasData(fSTNumDiscreteEnergyTransferVals);
  // fill in 3 values of the transformed variable
  int curNumData    = 3;
  raData->fXdata[0] = 0.0;
  raData->fXdata[1] = 0.5;
  raData->fXdata[2] = 1.0;
  int glnum         = 4;
  // fill in the discrete pdf array by inserting the remaining fSTNumDiscreteEnergyTransferVals-3 discrete points such
  // that the interpolation error is minimized: one new discrete point will be inserted at the end of each  iteration of
  // this while loop:
  while (curNumData<fSTNumDiscreteEnergyTransferVals) {
    // compute the pdf values at the current discrete sample points
    for (int i=0;i<curNumData;++i) {
      pdfarray[i] = ComputeDXSection(epsMin,eps0,dFact,FZ,raData->fXdata[i]);
    }
    double maxerr      = 0.0;
    double thexval     = 0.0;
    int    maxErrIndex = 0;
    // prepare data for the approximation of the pdf based on the current discretization
    const double norm  = fAliasSampler->PreparRatinForPDF(&(raData->fXdata[0]), pdfarray, &(raData->fCumulative[0]),
                                                          &(raData->fParaA[0]), &(raData->fParaB[0]), curNumData,
                                                          false, glnum);
    // find the interval with the highest approximation error to the real pdf based on the current discretization
    for (int i=0; i<curNumData-1; ++i) {
      const double xx  = 0.5*(raData->fXdata[i]+raData->fXdata[i+1]);  // mid xi-point of the current(i.e. i-th) interval
      double err = 0.0;
      // we shuld compute the integrated error (by transforming the integral from [xi[i],xi[i+1] to [0,1])
      // but we will use sample points to estimate the error: divide the interval into isub sub intervals and compute
      // some error measure at each edge points
      const int isub   = 5;
      const double dd  = (raData->fXdata[i+1]-raData->fXdata[i])/((double)isub);
      for (int j=1; j<isub; ++j) {
        const double xval = raData->fXdata[i]+j*dd;
        double valAprx = 0.0;
        // - get the approximated pdf value based on the current discretization
        // - unless if the current interval is the first because we use linear appoximation in the first
        //   interval <= the pdf(xi[0])=0
        if (i>0) {
          valAprx = fAliasSampler->GetRatinForPDF1(xval, &(raData->fXdata[0]), &(raData->fCumulative[0]),
                                                   &(raData->fParaA[0]), &(raData->fParaB[0]), i);
        } else { //linear aprx in the first interval because pdf[0] = 0
          valAprx = norm*pdfarray[i+1]/(raData->fXdata[i+1]-raData->fXdata[i])*(xval-raData->fXdata[i]);
        }
        // get the real value
        const double valReal = norm*ComputeDXSection(epsMin, eps0, dFact, FZ, xval);
        err += std::fabs((1.-valAprx/valReal)*(valAprx-valReal));
      }
      err *=(raData->fXdata[i+1]-raData->fXdata[i]);
      // if the current interval gives the highest approximation error so far then store some values :
      // i.e. current maximum error value, the mid-point of the current interval and the lower edge index
      if (err>maxerr) {
        maxerr      = err;
        thexval     = xx;
        maxErrIndex = i;
      }
    } // end of the for loop i.e. the approximation error was checked in each of the currently used intervals
    //
    // after checking the approximation error in each interval we want to insert a new point into the midle of that
    // interval that gave the maximum approximation error:
    // 1. shift to the right all points above the maxErrIndex
    for (int i=curNumData; i>maxErrIndex+1; --i) {
      raData->fXdata[i] = raData->fXdata[i-1];
    }
    // 2. then insert the new discrete xi point
    raData->fXdata[maxErrIndex+1] = thexval;
    // 3. increase the number of currently used discrete sample points
    ++curNumData;
  }
  //
  // if all the available discrete xi points was filled, compute the corresponding values of the pdf at these discrete
  // points
  for (int i=0; i<curNumData; ++i) {
    pdfarray[i] = ComputeDXSection(epsMin, eps0, dFact, FZ, raData->fXdata[i]);
  }
  //
  // and prepare the final ratin-alias sampling table structure:
  fAliasSampler->PreparRatinTable(&(raData->fXdata[0]), pdfarray, &(raData->fCumulative[0]), &(raData->fParaA[0]),
                                  &(raData->fParaB[0]), &(raData->fAliasW[0]), &(raData->fAliasIndx[0]),
                                  fSTNumDiscreteEnergyTransferVals, false, glnum);
  // fill in
  fSamplingTables[izet]->fRatinAliasData[egammaindx]=raData;
}


/**
 * @internal
 *  The final state sampling is based on the Bethe-Heitler \cite bethe1934stopping differential cross section(DCS)
 *  corrected for various effects (screening, Coulomb correction, conversion in the field of atomic electrons) that, for
 *  a given target atom with atomic number \f$ Z \f$ and photon energy \f$ E_{\gamma}\f$ to create an e-/e+ pair of
 *  which one has a total energy \f$ E \f$, can be writte as \cite motz1969pair \cite tsai1974pair
 *  \f[
 *    \frac{\mathrm{d}\sigma(Z,\epsilon)}{\mathrm{d}\epsilon} = \alpha r_0^2 Z[Z+\eta(Z)]
 *       \left\{
 *         \left[ (\epsilon^2+(1-\epsilon)^2) \right] \left[ \Phi_1(\delta(\epsilon)) - \frac{F(Z)}{2} \right]
 *         + \frac{2}{3}\epsilon(1-\epsilon) \left[ \Phi_2(\delta(\epsilon)) - \frac{F(Z)}{2} \right]
 *       \right\}
 *  \f]
 *  where \f$ \alpha \f$ is the fine structure constant, \f$ r_0 \f$ is the classical electron radius and
 *  \f$ \epsilon \equiv E/E_{\gamma} \f$ is the reduced total energy of one of the created e-/e+ pair.
 *
 *  The threshold energy for pair production (considering conversion only in the the field of the nucleus with infinite
 *  mass) is \f$ 2m_ec^2 \f$ that gives the kinematical limits of \f$ \epsilon \f$
 *  \f[
 *    \epsilon_0 \equiv \frac{m_ec^2}{E_{\gamma}} \leq \epsilon \leq 1-\epsilon
 *  \f]
 *  Since the DCS is symmetric under the exchange \f$ \epsilon \Leftrightarrow 1-\epsilon \f$, the kinematical range of
 * \f$ \epsilon \f$ can be restricted to \f$ \epsilon \in [\epsilon_0,0.5] \f$ (i.e. the DCS is symmetric to 0.5).
 *
 * \f$ \textbf{Screening correction:} \f$
 *
 * The Bethe-Heitler DCS derived in \cite bethe1934stopping describes e-/e+ pair production in the Coulomb field of a
 * point nucleus. The effect of screening this nuclear Coulomb field by the atomic electrons was also accounted in
 * \cite bethe1934stopping by introducing screening functions \f$ \Phi_1, \Phi_2 \f$ based on atomic form factors using
 * the Thomas-Fermi model of atom. An analytical approximation of these screening funtions were provided in
 * \cite butcher1960electron
 * \f[
 *  \begin{array}{l r c l r c l}
 *   \text{if } & \delta(\epsilon) & \leq 1 &
 *       \Phi_1(\delta(\epsilon)) & = & 20.867-3.242\delta(\epsilon)+0.625\delta^2(\epsilon)\\
 *   &&& \Phi_2(\delta(\epsilon)) & = & 20.209-1.930\delta(\epsilon)-0.086\delta^2(\epsilon)\\
 *   \text{if } & \delta(\epsilon) & > 1 &
 *       \Phi_1(\delta(\epsilon)) & = &  \Phi_1(\delta(\epsilon)) = 21.12-4.184\ln[\delta(\epsilon)+0.952]
 *  \end{array}
 * \f]
 * where
 * \f[
 *   \delta(\epsilon) = \frac{136Z^{-1/3}\epsilon_0}{\epsilon(1-\epsilon)}
 * \f]
 *
 * \f$ \textbf{Coulomb correction:} \f$
 * The DCS in \cite bethe1934stopping was derived under the first order(plane wave) Born approximation. Correction,
 * accounting the effects of the higher order terms in the Born series, were derived in \cite davies1954theory giving
 * the Coulomb correction function \f$ F(Z) \f$ as
 * \f[
 *  F(Z) =
 *  \begin{cases}
 *    \frac{8}{3} \ln(Z) & \quad \text{if } E_{\gamma} < 50 \text{[MeV]} \\
 *    \frac{8}{3} \ln(Z) + 8 f_c(Z) & \quad \text{if } E_{\gamma} \geq 50 \text{[MeV]}
 *  \end{cases}
 * \f]
 * with
 *  \f[
 *   f_c(\nu) = \nu^2 \sum_{n=1}^{\infty} \frac{1}{n(n^2+\nu^2)} = \nu^2 \left[  1/(1+\nu^2)  + 0.20206 - 0.0369\nu^2
 *            + 0.0083\nu^4 - 0.002\nu^6 \right]
 *  \f]
 * where \f$\nu=\alpha Z\f$. It should be noted, that this Coulomb correction can result in negative DCS at high
 * \f$ \delta(\epsilon) \f$ values when \f$ \Phi_{1,2} -F(Z)/2 \f$ becomes \f$ < 0\f$. Using the expression for
 * \f$ \Phi_{1,2}\f$ in case of \f$ \delta(\epsilon)>1 \f$ one can derive the maximum value \f$ \delta'(\epsilon) \f$ at
 * which the DCS is still non-negative by solving \f$ \Phi_{1,2} -F(Z)/2 = 0 \f$ that gives
 * \f[
 *   \delta'(\epsilon) = \exp \left[ \frac{42.24-F(Z)}{8.368} \right] - 0.952
 * \f]
 * which (according to the expression for \f$ \delta(\epsilon) \f$) corresponds to an
 * \f$ \epsilon' = 0.5-0.5\sqrt{1-4\alpha/\delta'}\f$ with \f$ \alpha = 136Z^{-1/3}\epsilon_0 \f$. Therefore, the
 * kinematically allowed \f$ \epsilon \in [\epsilon_0, 0.5] \f$ interval needs to be constrained into
 * \f$ \epsilon \in [\epsilon_{\text{min}}, 0.5] \f$ with \f$ \epsilon_{\text{min}} = \text{max}[\epsilon_0,\epsilon']\f$
 * in order to exclude the possible \f$ \epsilon \f$ interval with negative DCS values.
 *
 *\f$ \textbf{Pair production in the field of atomic electrons:} \f$
 *
 * Conversion into e-/e+ pair in the field of atomic electrons, that is proportional to the target atomic number
 * \f$Z\f$, is approximately taken into account by the \f$ \eta(Z) \f$ factor
 * \f[
 *    \eta(Z) = \frac{L_{\text{inel}}}{L_{\text{el}}-f_c(Z)}
 * \f]
 * by assuming the same dependence on \f$ \epsilon \f$ as the contribution of the conversion in filed of the nucleus.
 * \f$ L_{\text{el}}\f$ and \f$ L_{\text{inel}} \f$ are as given at RelativisticBremsModel::ComputeDXSecPerAtom() .
 * More detailed discussion and more accurate description can be found in \cite hubbell1980pair.
 *
 *
 * Note, that the model is based on a parametrization or the atomic cross sections given in \cite hubbell1980pair.
 * Therefore, only the \f$ \epsilon \f$ dependent part of the above DCS will be used by the model: for sampling the
 * reduced total energy transfered to one of the e-/e+ pair. The above DCS can be written in the form of
 * \f[
 *   \frac{\mathrm{d}\sigma(Z,\epsilon)}{\mathrm{d}\epsilon}
 *     = \kappa(Z) \left( \frac{\mathrm{d}\sigma(Z,\epsilon)}{\mathrm{d}\epsilon} \right)^*
 * \f]
 * with the \f$ \epsilon \f$ dependent
 * \f[
 *         \left( \frac{\mathrm{d}\sigma(Z,\epsilon)}{\mathrm{d}\epsilon} \right)^* =
 *           \left[ (\epsilon^2+(1-\epsilon)^2) \right] \left[ \Phi_1(\delta(\epsilon)) - \frac{F(Z)}{2} \right]
 *         + \frac{2}{3}\epsilon(1-\epsilon) \left[ \Phi_2(\delta(\epsilon)) - \frac{F(Z)}{2} \right]
 * \f]
 * and constant \f$ \kappa(Z) = \alpha r_0^2 Z[Z+\eta(Z)]\f$. This DCS can be transformed to
 * \f[
 *      \left( \frac{\mathrm{d}\sigma(Z,\epsilon(\xi))}{\mathrm{d}\xi} \right)^* = \epsilon(\xi) \left\{
 *      \left[ (\epsilon(\xi)^2+(1-\epsilon(\xi))^2) \right] \left[ \Phi_1(\delta(\epsilon(\xi)))-\frac{F(Z)}{2} \right]
 *    + \frac{2}{3}\epsilon(\xi)(1-\epsilon(\xi)) \left[ \Phi_2(\delta(\epsilon(\xi))) - \frac{F(Z)}{2} \right]
 *     \right\}
 * \f]
 * by using \f$ \xi(\epsilon) = \ln[\epsilon/\epsilon_{\text{min}}]/\ln[0.5/\epsilon_{\text{min}}] \f$ and dropping a
 * \f$ \xi \f$ independent \f$ \ln[0.5/\epsilon_{\text{min}}] \f$ factor. Note, that \f$ \xi \in [0,1] \f$ when
 * \f$ \epsilon \in [\epsilon_{\text{min}},0.5] \f$ and
 * \f$ \epsilon(\xi) = \epsilon_{\text{min}}\exp[\xi\ln(0.5/\epsilon_{\text{min}})] \f$.
 * Therefore, the transformed variable \f$ \xi \f$ lies in the same interval independently from the primary photon
 * energy \f$ E_{\gamma} \f$ and the low \f$ \epsilon \f$ part of the original DCS, that can be strongly nonlinear, is
 * enhanced thanks to the logarithmic dependence of \f$ \xi \f$ on \f$ \epsilon \f$.
 *
 * The value of \f$ \left( \mathrm{d}\sigma(Z,\epsilon(\xi))/\mathrm{d}\xi \right)^* \f$ is computed in this method.
 * When sampling tables are requested to be used, sampling tables are built to sample this transformed variable based on
 * this method. These tables are used at run-time to sample \f$ \xi \f$ and the sampled value is transformed back to the
 * reduced total energy transfered to one of the e-/e+ pair by using the realtion
 * \f$ \epsilon(\xi) = \epsilon_{\text{min}}\exp[\xi\ln(0.5/\epsilon_{\text{min}})] \f$.
 *
 * @endinternal
 */
double BetheHeitlerPairModel::ComputeDXSection(double epsmin, double eps0, double deltafactor, double fz, double xi) {
  const double lHalfPerEpsMin = std::log(0.5/epsmin);
  const double eps            = epsmin*std::exp(xi*lHalfPerEpsMin);
  const double meps           = 1.-eps;
  const double delta          = deltafactor*eps0/(eps*meps);
  double phi1, phi2;
  ComputeScreeningFunctions(phi1, phi2, delta, true);
  double dxsec = eps*( (eps*eps+meps*meps)*(0.25*phi1-fz) + 2.*eps*meps*(0.25*phi2-fz)/3.);
  return std::max(dxsec, 0.0);
}


}  // namespace geantphysics
