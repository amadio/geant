
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


namespace geantphysics {


BetheHeitlerPairModel::BetheHeitlerPairModel(const std::string &modelname) : EMModel(modelname) {
  fElectronInternalCode             = -1;      // will be set at init
  fPositronInternalCode             = -1;      // will be set at init
  fMinimumPrimaryEnergy             =  2.*geant::kElectronMassC2; // final value will be set at init.
  fGammaEneregyLimit                =  2.*geant::MeV; // use simplified sampling below this gamma energy

  fElementData                      = nullptr;

  fNumSamplingPrimEnergiesPerDecade = 12;    // should be set/get and must be done before init
  fNumSamplingEnergies              = 54;    // should be set/get and must be done before init
  fMinPrimEnergy                    = -1.;
  fMaxPrimEnergy                    = -1.;
  fNumSamplingPrimEnergies          = -1;      // will be set in InitSamplingTables if needed
  fPrimEnLMin                       =  0.;     // will be set in InitSamplingTables if needed
  fPrimEnILDelta                    =  0.;     // will be set in InitSamplingTables if needed
  fSamplingPrimEnergies             = nullptr; // will be set in InitSamplingTables if needed
  fLSamplingPrimEnergies            = nullptr; // will be set in InitSamplingTables if needed

  fRatinAliasDataForAllElements     = nullptr;
  fAliasSampler                     = nullptr;
}


BetheHeitlerPairModel::~BetheHeitlerPairModel() {
  if (fElementData) {
    for (int i=0; i<gMaxZet; ++i) {
      if (fElementData[i]) {
        delete fElementData[i];
      }
    }
    delete [] fElementData;
  }
  //
  if (fSamplingPrimEnergies) {
    delete [] fSamplingPrimEnergies;
    delete [] fLSamplingPrimEnergies;
  }
  //
  if (fRatinAliasDataForAllElements) {
    for (int i=0; i<gMaxZet; ++i) {
      RatinAliasDataPerElement *oneElem = fRatinAliasDataForAllElements[i];
      if (oneElem) {
        for (int j=0; j<fNumSamplingPrimEnergies; ++j) {
          RatinAliasData *oneRatinAlias = oneElem->fRatinAliasDataForOneElement[j];
          if (oneRatinAlias) {
            delete [] oneRatinAlias->fXdata;
            delete [] oneRatinAlias->fCumulative;
            delete [] oneRatinAlias->fParaA;
            delete [] oneRatinAlias->fParaB;
            delete [] oneRatinAlias->fAliasW;
            delete [] oneRatinAlias->fAliasIndx;
            delete oneRatinAlias;
          }
        }
        delete [] oneElem->fRatinAliasDataForOneElement;
        delete oneElem;
      }
    }
    delete [] fRatinAliasDataForAllElements;
  }
  //
  if (fAliasSampler) {
    delete fAliasSampler;
  }
}

/**
 * @internal
 *  Internal codes of secondary particles (i.e. e-/e+) will be set, the minimum value of the primary photon energy will
 *  be set (\f$ E_{\gamma}^{\text{min}} = \text{max}\{2m_ec^2,\text{min-energy-usage} \}\f$), the model is initialised
 *  by calling the InitialiseModel() method and a target element selector is required and initialised.
 * @endinternal
 */
void   BetheHeitlerPairModel::Initialize() {
  EMModel::Initialize();  // will set the PhysicsParameters member
  fElectronInternalCode = Electron::Definition()->GetInternalCode();
  fPositronInternalCode = Positron::Definition()->GetInternalCode();
  fMinimumPrimaryEnergy = 2.*geant::kElectronMassC2; // will be used to build table in target element selector
  if (GetLowEnergyUsageLimit()>fMinimumPrimaryEnergy) {
    fMinimumPrimaryEnergy = GetLowEnergyUsageLimit();
  }
  InitialiseModel();
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
  double       egamma = kinenergy;
  // we will need the element composition of this material
  const Vector_t<Element*> theElements    = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int   numElems = theElements.size();
  for (int iel=0; iel<numElems; ++iel) {
    xsec += theAtomicNumDensityVector[iel]*ComputeAtomicCrossSection(theElements[iel]->GetZ(), egamma);
  }
  return xsec;
}


double BetheHeitlerPairModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts*, double kinenergy,
                                                     const Particle*) {
   double xsec  = 0.0;
   if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
     return xsec;
   }
   // compute the parametrized atomic cross section: depends only on target Z and gamma energy.
   xsec = ComputeAtomicCrossSection(elem->GetZ(), kinenergy);
   return xsec;
}


int BetheHeitlerPairModel::SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td) {
  int    numSecondaries      = 0;
  double ekin                = track.GetKinE();
  double eps0                = geant::kElectronMassC2/ekin;
  // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
  // check if kinetic energy is above fHighEnergyUsageLimit andd o nothing if yes
  if (ekin<GetLowEnergyUsageLimit() || ekin>GetHighEnergyUsageLimit() || eps0>0.5) {
    return numSecondaries;
  }
  // interaction is possible so sample target element: will be needed anyway for the direction sampling
  MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[track.GetMaterialCutCoupleIndex()];
  const Vector_t<Element*> theElements = matCut->GetMaterial()->GetElementVector();
  double targetElemIndx = 0;
  if (theElements.size()>1) {
    targetElemIndx = SampleTargetElementIndex(matCut, ekin, td->fRndm->uniform());
  }
  double  zet  = theElements[targetElemIndx]->GetZ();
  //
  // sample the reduced total energy transfered to one of the secondary particles i.e. either e- or e+
  double eps = 0.0;
  if (ekin<fGammaEneregyLimit) {   // sample eps from uniform on [eps_0,0.5]
    eps = eps0 + (0.5-eps0)*td->fRndm->uniform();
  } else {
    int           izet = std::lrint(zet);
    double      epsMin = eps0;
    double          FZ = fElementData[izet]->fFzLow;      // used only in case of rejection
    double    deltaMax = fElementData[izet]->fDeltaMaxLow;
    double deltaFactor = fElementData[izet]->fDeltaFactor;
    if (ekin>50.*geant::MeV) {
      FZ       = fElementData[izet]->fFzHigh;             // used only in case of rejection
      deltaMax = fElementData[izet]->fDeltaMaxHigh;
    }
    double deltaMin = 4.*eps0*deltaFactor;
    double eps1     = 0.5-0.5*std::sqrt(1.-deltaMin/deltaMax);
    if (eps1>epsMin) {
      epsMin = eps1;
    }
    // sample the reduced total energy transferd to one of the e-/e+ pair by using tables or rejection
    if (GetUseSamplingTables()) {
      // sample eps using the sampling tables built at initialisation over a photon energy grid for each possible target
      double *rndArray  = td->fDblArray;
      td->fRndm->uniform_array(4, rndArray);
      eps = SampleTotalEnergyTransfer(ekin, epsMin, izet, rndArray[0], rndArray[1], rndArray[2]);
    } else {
      // sample eps by rejection
      eps = SampleTotalEnergyTransfer(epsMin, eps0, deltaMin, FZ, deltaFactor, td);
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
  double thetaElectron = uvar*geant::kElectronMassC2/electronTotE;
  double sintEle       = std::sin(thetaElectron);
  double thetaPositron = uvar*geant::kElectronMassC2/positronTotE;
  double sintPos       = -std::sin(thetaPositron);
  double phi           = geant::kTwoPi*rndArray[3];
  double sinphi        = std::sin(phi);
  double cosphi        = std::cos(phi);
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
  double ekinElectron = (electronTotE-geant::kElectronMassC2);
  // final check
  if (ekinElectron<0.0) {
    ekinElectron = 0.0;
  }
  double ekinPositron = (positronTotE-geant::kElectronMassC2);
  // final check
  if (ekinPositron<0.0) {
    ekinPositron = 0.0;
  }
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
 *  A collection of frequently used target atom dependent variables is set up bu calling the InitialiseElementData()
 *  method. If sampling tables were requested to be used, then the min/max of the primary photon energy grid (which
 *  sampling tables will be built over) are set and sampling tables are initialized by calling the InitSamplingTables()
 *  method.
 * @endinternal
 */
void BetheHeitlerPairModel::InitialiseModel() {
  InitialiseElementData();
  // everything from this line is for building sampling tables:
  fMinPrimEnergy   = fGammaEneregyLimit;   // fGammaEneregyLimit = 2.*geant::MeV
  if (fMinPrimEnergy<GetLowEnergyUsageLimit()) {
    fMinPrimEnergy = GetLowEnergyUsageLimit();
  }
  fMaxPrimEnergy = GetHighEnergyUsageLimit();
  if (GetUseSamplingTables()) {
//    std::cerr<< "  === BH pair model: building sampling tables"<< std::endl;
    InitSamplingTables();
  }
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
  double kappa  = std::log(egamma/geant::kElectronMassC2);
  double kappa2 = kappa*kappa;
  double kappa3 = kappa2*kappa;
  double kappa4 = kappa2*kappa2;
  double kappa5 = kappa4*kappa;
  //
  double F1 = a0 + a1*kappa + a2*kappa2 + a3*kappa3 + a4*kappa4 + a5*kappa5;
  double F2 = b0 + b1*kappa + b2*kappa2 + b3*kappa3 + b4*kappa4 + b5*kappa5;
  double F3 = c0 + c1*kappa + c2*kappa2 + c3*kappa3 + c4*kappa4 + c5*kappa5;
  // compute cross section
  xsec = (z+1.)*(F1*z+F2*z*z+F3);
  // low energy correction
  if (egammaOrg<egammaLimit) {
    double dum = (egammaOrg-2.*geant::kElectronMassC2)/(egammaLimit-2.*geant::kElectronMassC2);
    xsec *= dum*dum;
  }
  // protection against negative values
  xsec = std::max(xsec, 0.);
  return xsec;
}

/**
 * @internal
 *  One ElementData structure will be created for each target atom that the model needs to respond (i.e.
 *  for each elements that appears in material that belongs to a region in which the model is active).
 *  Pointers to these data structures are stored in the fElementData array indexed by atomic number (Z).
 *  Those indices, that corresponds to atomic numbers that the model do not need to respond, will remain
 *  nullptr-s.
 * @endinternal
 */
void BetheHeitlerPairModel::InitialiseElementData() {
  // clear element data if needed
  if (fElementData) {
    for (int i=0; i<gMaxZet; ++i) {
      if (fElementData[i]) {
        delete fElementData[i];
      }
    }
    delete [] fElementData;
  }
  // create the container and fill with default nullptr-s
  fElementData = new ElementData*[gMaxZet]();
  for (int i=0; i<gMaxZet; ++i) {
    fElementData[i] = nullptr;
  }
  //
  int numMatCuts = MaterialCuts::GetTheMaterialCutsTable().size();
  // get list of active region
  std::vector<bool> isActiveInRegion = GetListActiveRegions();
  for (int i=0; i<numMatCuts; ++i) {
    const MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[i];
    // if this MaterialCuts belongs to a region where this model is active:
    if (isActiveInRegion[matCut->GetRegionIndex()]) {
      // get the list of elements
      const Vector_t<Element*> theElements =  matCut->GetMaterial()->GetElementVector();
      int numElems = theElements.size();
      for (int j=0; j<numElems; ++j) {
        // if RatinAliasDataPerElement has not been built for this element yet => build
        double zet = theElements[j]->GetZ();
        int elementIndx = std::lrint(zet);
//        std::cerr<< " === Building ElementData for " << theElements[j]->GetName() << std::endl;
        if (!fElementData[elementIndx]) {
          Element         *elem   = theElements[j];
          ElementData *elemData   = new ElementData();
          elemData->fDeltaFactor  = 136./elem->GetElementProperties()->GetZ13(); // 136/pow(z,1/3)
          elemData->fCoulombCor   = elem->GetElementProperties()->GetCoulombCorrection();
          elemData->fFzLow        = 8.*elem->GetElementProperties()->GetLogZ13(); //8.*std::log(z)/3.;
          elemData->fFzHigh       = elemData->fFzLow+8*elemData->fCoulombCor;
          elemData->fDeltaMaxLow  = std::exp((42.24-elemData->fFzLow)/8.368)-0.952;
          elemData->fDeltaMaxHigh = std::exp((42.24-elemData->fFzHigh)/8.368)-0.952;
          fElementData[elementIndx] = elemData;
        }
      }
    }
  }
}


// samples the transformed variable for the given target atom(zindx), gamma energy, and transforms it back
double BetheHeitlerPairModel::SampleTotalEnergyTransfer(double primekin, double epsmin, int zindx, double r1, double r2,
                                                        double r3) {
  // determine primary energy lower grid point
  double lGammaEnergy  = std::log(primekin);
  int gammaEnergyIndx  = (int) ((lGammaEnergy-fPrimEnLMin)*fPrimEnILDelta);
  //
  if (gammaEnergyIndx>=fNumSamplingPrimEnergies-1)
    gammaEnergyIndx = fNumSamplingPrimEnergies-2;

  double pLowerGammaEner = (fLSamplingPrimEnergies[gammaEnergyIndx+1]-lGammaEnergy)*fPrimEnILDelta;
  if (r1>pLowerGammaEner) {
     ++gammaEnergyIndx;
  }
  // get the RatinAliasData object pointer for the given Z-index and gamma-energy-index(that is gammaEnergyIndx)
  // fRatinAliasDataForOneElement[zindx] should not be nullptr because it was checked in the caller
  RatinAliasData *raData = fRatinAliasDataForAllElements[zindx]->fRatinAliasDataForOneElement[gammaEnergyIndx];

  // we could put a static assert here to make sure that raData is not nullptr
  //
  // sample the transformed variable xi=... (such that linear aprx will be used in the first interval)
  double  xi = fAliasSampler->SampleRatin(raData->fXdata, raData->fCumulative, raData->fParaA, raData->fParaB,
                                          raData->fAliasW, raData->fAliasIndx, raData->fNumdata, r2, r3, 0);

  // transform back xi to eps = E_total_energy_transfer/E_{\gamma}
  double  lHalfPerEpsMin = std::log(0.5/epsmin);
  return epsmin*std::exp(xi*lHalfPerEpsMin);
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
double BetheHeitlerPairModel::SampleTotalEnergyTransfer(double epsmin, double eps0, double deltamin, double fz,
                                                        double deltafactor, Geant::GeantTaskData *td) {
    double eps      = 0.0;
    double epsRange = 0.5-epsmin;
    //
    double F10      = ScreenFunction1(deltamin) - fz;
    double F20      = ScreenFunction2(deltamin) - fz;
    double NormF1   = std::max(F10*epsRange*epsRange,0.);
    double NormF2   = std::max(1.5*F20,0.);

    double *rndArray = td->fDblArray;
    double greject   = 0.0;
    do {
      td->fRndm->uniform_array(3, rndArray);
      if (NormF1/(NormF1+NormF2)>rndArray[0]) {
      	eps            = 0.5 - epsRange*std::pow(rndArray[1],1./3.);
	      double delta   = deltafactor*eps0/(eps*(1.-eps));
	      greject = (ScreenFunction1(delta) - fz)/F10;
      } else {
	      eps            = epsmin + epsRange*rndArray[1];
        double delta   = deltafactor*eps0/(eps*(1.-eps));
	      greject = (ScreenFunction2(delta) - fz)/F20;
      }
    } while (greject<rndArray[2]);
  return eps;
}

// 3xPhi_1 - Phi_2
double BetheHeitlerPairModel::ScreenFunction1(double delta) {
  double val;
  if (delta>1.) {
    val = 42.24 - 8.368*std::log(delta+0.952);
  } else {
    val = 42.392 - delta*(7.796 - 1.961*delta);
  }
  return val;
}

// 1.5*Phi_1 + 0.5*Phi_2
double BetheHeitlerPairModel::ScreenFunction2(double delta) {
  double val;
  if (delta>1.) {
    val = 42.24 - 8.368*std::log(delta+0.952);
  } else {
    val = 41.405 - delta*(5.828 - 0.8945*delta);
  }
  return val;
}


void BetheHeitlerPairModel::InitSamplingTables() {
  // set number of primary gamma energy grid points
  // keep the prev. value of primary energy grid points.
  int oldNumGridPoints     = fNumSamplingPrimEnergies;
  fNumSamplingPrimEnergies = fNumSamplingPrimEnergiesPerDecade*std::lrint(std::log10(fMaxPrimEnergy/fMinPrimEnergy))+1;
  if (fNumSamplingPrimEnergies<2) {
    fNumSamplingPrimEnergies = 2;
  }
  // set up the initial gamma energy grid
  if (fSamplingPrimEnergies) {
    delete [] fSamplingPrimEnergies;
    delete [] fLSamplingPrimEnergies;
    fSamplingPrimEnergies  = nullptr;
    fLSamplingPrimEnergies = nullptr;
  }
  fSamplingPrimEnergies  = new double[fNumSamplingPrimEnergies];
  fLSamplingPrimEnergies = new double[fNumSamplingPrimEnergies];
  fPrimEnLMin    = std::log(fMinPrimEnergy);
  double delta   = std::log(fMaxPrimEnergy/fMinPrimEnergy)/(fNumSamplingPrimEnergies-1.0);
  fPrimEnILDelta = 1.0/delta;
  fSamplingPrimEnergies[0]  = fMinPrimEnergy;
  fLSamplingPrimEnergies[0] = fPrimEnLMin;
  fSamplingPrimEnergies[fNumSamplingPrimEnergies-1]  = fMaxPrimEnergy;
  fLSamplingPrimEnergies[fNumSamplingPrimEnergies-1] = std::log(fMaxPrimEnergy);
  for (int i=1; i<fNumSamplingPrimEnergies-1; ++i) {
    fLSamplingPrimEnergies[i] = fPrimEnLMin+i*delta;
    fSamplingPrimEnergies[i]  = std::exp(fPrimEnLMin+i*delta);
  }
  //
  // clear (if needed) and create the container that can store sampling tables for all(active) elements
  // 1. clear:
  if (fRatinAliasDataForAllElements) {
    for (int i=0; i<gMaxZet; ++i) {
      RatinAliasDataPerElement *oneElem = fRatinAliasDataForAllElements[i];
      if (oneElem) {
        for (int j=0; j<oldNumGridPoints; ++j) {
          RatinAliasData *oneRatinAlias = oneElem->fRatinAliasDataForOneElement[j];
          if (oneRatinAlias) {
            delete [] oneRatinAlias->fXdata;
            delete [] oneRatinAlias->fCumulative;
            delete [] oneRatinAlias->fParaA;
            delete [] oneRatinAlias->fParaB;
            delete [] oneRatinAlias->fAliasW;
            delete [] oneRatinAlias->fAliasIndx;
            delete oneRatinAlias;
          }
        }
        delete [] oneElem->fRatinAliasDataForOneElement;
        delete oneElem;
      }
    }
    delete [] fRatinAliasDataForAllElements;
  }
  //
  // 2. create an AliasTable object
  if (fAliasSampler) {
    delete fAliasSampler;
  }
  fAliasSampler = new AliasTable();
  //
  // 3. create and fill with nullptr-s
  fRatinAliasDataForAllElements = new RatinAliasDataPerElement*[gMaxZet];
  for (int i=0; i<gMaxZet; ++i) {
    fRatinAliasDataForAllElements[i] = nullptr;
  }
  //
  // get all MaterialCuts; loop over them and if they belongs to a region where this model is active:
  //  - get the list of elements of the material that the MaterialCuts belongs to
  //  - loop over the this list of elements and build RatinAliasDataPerElement for which it has not been done yet
  int numMatCuts = MaterialCuts::GetTheMaterialCutsTable().size();
  // get list of active region
  std::vector<bool> isActiveInRegion = GetListActiveRegions();
  for (int i=0; i<numMatCuts; ++i) {
    const MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[i];
    // if this MaterialCuts belongs to a region where this model is active:
    if (isActiveInRegion[matCut->GetRegionIndex()]) {
      // get the list of elements
      const Vector_t<Element*> theElements =  matCut->GetMaterial()->GetElementVector();
      int numElems = theElements.size();
      for (int j=0; j<numElems; ++j) {
        // if RatinAliasDataPerElement has not been built for this element yet => build
        double zet = theElements[j]->GetZ();
        int elementIndx = std::lrint(zet);
        if (!fRatinAliasDataForAllElements[elementIndx]) {
          BuildSamplingTablesForElement(theElements[j]);
        }
      }
    }
  }
}


void BetheHeitlerPairModel::BuildSamplingTablesForElement(const Element *elem) {
  // prepare one RatinAliasDataPerElement structure
  RatinAliasDataPerElement *perElem          = new RatinAliasDataPerElement();
  perElem->fRatinAliasDataForOneElement      = new RatinAliasData*[fNumSamplingPrimEnergies];
  int elementindx                            = std::lrint(elem->GetZ());
  fRatinAliasDataForAllElements[elementindx] = perElem;
  // allocate an array that will be used temporary:
  double *thePdfdata = new double[fNumSamplingEnergies]();
  // loop over the photon energy grid and build one RatinAliasData structure at each photon energy for this element
  for (int i=0; i<fNumSamplingPrimEnergies; ++i) {
    // note that egamma is for sure > 2mc^2,; moreover, for sure >= 2MeV
    double egamma = fSamplingPrimEnergies[i];
    // build one table for the given: egamma, element
//    std::cerr<<"  === Bulding table for:  elem = " << elem->GetName() << " egamma = " << egamma/geant::MeV << " [MeV]" << std::endl;
    BuildOneRatinAlias(egamma, elem, thePdfdata, i);
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
 *
 * Due to the different approximation of the screening function below and above \f$ \delta(\epsilon) =1 \f$, there is a
 * step (within 0.5 %) in the pdf (DCS) if the \f$ \delta \in [\delta(\epsilon_{max}=0.5),\delta(\epsilon_{min})] \f$
 * interval contains the value 1. The step is located at
 * \f$ \xi(\delta=1) = \ln[\epsilon(\delta=1)/\epsilon_{min}]/\ln[0.5/\epsilon_{min}] \f$
 * where \f$ \epsilon(\delta=1) = 0.5-0.5\sqrt{1.-4 136 Z^{-1/3} \epsilon_0}\f$ is the \f$ \epsilon \f$ value that
 * coresponds to \f$ \delta = 1 \f$ and \f$ \xi(\delta=1) \f$ is the corresponding transformed variable.A small interval
 * will be inserted around this \f$ \xi(\delta=1) \f$ value in the discretized pdf and this interval will be ignored in
 * the error computation. Thanks to this, the individual distributions will reproduce this step properly. However, the
 * position of this step depends on the gamma energy through \f$ \epsilon_0 \f$ and sampling tables are built only at
 * discrete gamma energies. A linear interpolation on log gamma energy scale is used at run-time to sample the energy
 * transfer related variable \f$ \xi \f$ from the pre-prepared tables. Therefore, the sampled distributions will contain
 * this step in an interval (between the position of the step at the upper and lower gamma energies) instead at a point.
 *
 * @endinternal
 */
void BetheHeitlerPairModel::BuildOneRatinAlias(double egamma, const Element *elem, double *pdfarray, int egammaindx) {
  // compute the theoretical minimum of the reduced total energy transfer to the e+ (or to the e-)
  double eps0        = geant::kElectronMassC2/egamma;
  // set the real minimum to the theoretical one: will be updated
  double epsMin      = eps0;
  int    izet        = std::lrint(elem->GetZ());
  double FZ          = fElementData[izet]->fFzLow;
  double deltaMax    = fElementData[izet]->fDeltaMaxLow;
  double deltaFactor = fElementData[izet]->fDeltaFactor;
  if (egamma>50.*geant::MeV) {
    FZ       = fElementData[izet]->fFzHigh;
    deltaMax = fElementData[izet]->fDeltaMaxHigh;
  }
  double deltaMin = 4.*eps0*deltaFactor;  // 4*136*Z^(-1/3)*eps0
  double eps1     = 0.5-0.5*std::sqrt(1.-deltaMin/deltaMax);
  FZ *= 0.125; //this is how we use in our dxsec computation
  if (eps1>epsMin) {
    epsMin = eps1;
  }
  // there is an artificial step in the DCS due to the 2 different approximation to the screening funtion below/above
  // delta = 1 => we need to locate this point and in the corresponding xi variable a small interval will be inserted
  // that will be excluded from the error computation.
  // compute delta(eps) at eps-limits
  double deltaAtEpsMin = deltaFactor*eps0/(epsMin*(1.-epsMin)); // 136*Z^(-1/3) eps0/(epsMin*(1-epsMin))
  double deltaAtEpsMax = 4.*deltaFactor*eps0;  // 4*136*Z^(-1/3) eps0
  // compute eps that is such that delta(eps) = 1 if delta cross the value 1 between [epsMin,0.5]
  double xiDeltaU  = -1.0;
  double xiDeltaD  = -1.0;
  double xiDelta1  = -1.0;
  if (deltaAtEpsMax<1. && deltaAtEpsMin>1.) {
    double epsDelta1 = 0.5-0.5*std::sqrt(1.-4*deltaFactor*eps0);
    // the corresponding xi(epsDelta1) transformed variable
    xiDelta1  = std::log(epsDelta1/epsMin)/std::log(0.5/epsMin);
    // make a small interval that covers this xi(epsDelta1) value
    xiDeltaU = xiDelta1+xiDelta1*(1.-xiDelta1)*1.e-3;
    xiDeltaD = xiDelta1-xiDelta1*(1.-xiDelta1)*1.e-3;
  }
//std::cerr<< " xiDeltaD = " << xiDeltaD << std::endl;
  //
  // allocate one RatinAliasData structure for the current element,gamma-energy combination
  RatinAliasData *raData = new RatinAliasData();
  raData->fNumdata       = fNumSamplingEnergies;
  raData->fXiDelta1      = xiDelta1;
  raData->fXdata         = new double[fNumSamplingEnergies]();
  raData->fCumulative    = new double[fNumSamplingEnergies]();
  raData->fParaA         = new double[fNumSamplingEnergies]();
  raData->fParaB         = new double[fNumSamplingEnergies]();
  raData->fAliasW        = new double[fNumSamplingEnergies]();
  raData->fAliasIndx     = new int   [fNumSamplingEnergies]();

  //
  // fill in raData->fNumdata discrete sample point from the transformed variable and the corresponding pdf:
  // 1.a. fill in xi[min], xi[midle], xi[max] if there is no step in the DCS on xi [0,1]
  // 1.b. fill in xi[min], xi[max], the lower and upper limit of the small intervals plus one value
  // 2. set the currently used data points to 3 or 5 depending on if 1.a. or 1.b. the case
  int curNumData    = 3;
  int skipIndx      = -1;
  if (xiDeltaD<0.) {
    raData->fXdata[0] = 0.0;
    raData->fXdata[1] = 0.5;
    raData->fXdata[2] = 1.0;
  } else {
    raData->fXdata[0] = 0.0;
    raData->fXdata[4] = 1.0;
    if (xiDeltaU<0.5) {
      raData->fXdata[1] = xiDeltaD;
      raData->fXdata[2] = xiDeltaU;
      raData->fXdata[3] = 0.5;
      skipIndx = 1;
    } else if (xiDeltaD>0.5) {
      raData->fXdata[1] = 0.5;
      raData->fXdata[2] = xiDeltaD;
      raData->fXdata[3] = xiDeltaU;
      skipIndx = 2;
    } else {
      raData->fXdata[1] = xiDeltaD*0.5;
      raData->fXdata[2] = xiDeltaD;
      raData->fXdata[3] = xiDeltaU;
      skipIndx = 2;
    }
    curNumData = 5;
  }
  // create a 16 point GL integral to compute the interpolation error at each interval
  int glnum         = 4;
/*
  GLIntegral  *gl   = new GLIntegral(glnum, 0.0, 1.0);
  std::vector<double> glx = gl->GetAbscissas();
  std::vector<double> glw = gl->GetWeights();
*/
  //
  // fill in the discrete pdf array by inserting raData->fNumdata discrete points such that the
  // interpolation error is minimized: one new discrete point will be inserted at the end of each
  // iteration of this while loop:
  while (curNumData<raData->fNumdata) {
    // compute the pdf values at the current discrete sample points
    for (int i=0;i<curNumData;++i) {
      pdfarray[i] = ComputeDXSection(epsMin,eps0,deltaFactor,FZ,raData->fXdata[i]);
    }
    double maxerr      = 0.0;
    double thexval     = 0.0;
    int    maxErrIndex = 0;
    // prepare data for the approximation of the pdf based on the current discretization
    double norm = fAliasSampler->PreparRatinForPDF(raData->fXdata, pdfarray, raData->fCumulative, raData->fParaA,
                                                   raData->fParaB, curNumData, true, glnum);
    //
    // find the interval with the highest approximation error to the real pdf based on the current discretization
    for (int i=0; i<curNumData-1; ++i) {
      // exclude the possible insterted small interval around the step in the DCS
      if (i==skipIndx) {
        continue;
      }
      double xx  = 0.5*(raData->fXdata[i]+raData->fXdata[i+1]);  // mid xi-point of the current(i.e. i-th) interval
      double err = 0.0;
      // we shuld compute the integrated error (by transforming the integral from [xi[i],xi[i+1] to [0,1])
      // but we will use sample points to estimate the error: divide the interval into isub sub intervals and compute
      // some error measure at each edge points
      int isub   = 10;
      double dd  = (raData->fXdata[i+1]-raData->fXdata[i])/((double)isub);
      for (int j=1; j<isub; ++j) {
        double    xval = raData->fXdata[i]+j*dd;
        double valAprx = 0.0;
        // - get the approximated pdf value based on the current discretization
        // - unless if the current interval is the first because we use linear appoximation in the first
        //   interval <= the pdf(xi[0])=0
        if (i>0) {
          valAprx = fAliasSampler->GetRatinForPDF1(xval, raData->fXdata, raData->fCumulative, raData->fParaA,
                                                   raData->fParaB, i);
        } else { //linear aprx in the first interval because pdf[0] = 0
          valAprx = norm*pdfarray[i+1]/(raData->fXdata[i+1]-raData->fXdata[i])*(xval-raData->fXdata[i]);
        }
        // get the real value
        double valReal = norm*ComputeDXSection(epsMin, eps0, deltaFactor, FZ, xval);
        err += std::fabs((1.-valAprx/valReal)*(valAprx-valReal));
      }
      err *=(raData->fXdata[i+1]-raData->fXdata[i]);
/*
      // compute the integrated error (by transforming the integral from [xi[i],xi[i+1] to [0,1])
      for (int j=0; j<glnum; ++j) {
        double xval    = glx[j]*(raData->fXdata[i+1]-raData->fXdata[i])+raData->fXdata[i];
        double valAprx = 0.0;
        // - get the approximated pdf value based on the current discretization
        // - unless if the current interval is the first because we use linear appoximation in the first
        //   interval <= the pdf(xi[0])=0
        if (i>0) {
          valAprx = fAliasSampler->GetRatinForPDF1(xval, raData->fXdata, raData->fCumulative, raData->fParaA,
                                                   raData->fParaB, i);
        } else { //linear aprx in the first interval because pdf[0] = 0
          valAprx = norm*pdfarray[i+1]/(raData->fXdata[i+1]-raData->fXdata[i])*(xval-raData->fXdata[i]);
        }
        // get the real value
        double valReal = norm*ComputeDXSection(epsMin, eps0, deltaFactor, FZ, xval);
        // compute the error and add-up to the integrated value
        err += glw[j]*std::fabs((valAprx-valReal)*(1.-valAprx/valReal));
      }
      // due to the integral transform from [xi[i],xi[i+1]] to [0,1]
      err *=(raData->fXdata[i+1]-raData->fXdata[i]);
*/
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
      // update the position of the interval that we excluded from the procedure
      if (i-1==skipIndx) {
        skipIndx = i;
      }
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
    pdfarray[i] = ComputeDXSection(epsMin, eps0, deltaFactor, FZ, raData->fXdata[i]);
  }
  //
  // and prepare the final ratin-alias sampling table structure:
  fAliasSampler->PreparRatinTable(raData->fXdata, pdfarray, raData->fCumulative, raData->fParaA, raData->fParaB,
                                  raData->fAliasW, raData->fAliasIndx, raData->fNumdata, true, glnum);
  // fill in
  int elementindx = std::lrint(elem->GetZ());
  fRatinAliasDataForAllElements[elementindx]->fRatinAliasDataForOneElement[egammaindx]=raData;
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
  double lHalfPerEpsMin = std::log(0.5/epsmin);
  double eps            = epsmin*std::exp(xi*lHalfPerEpsMin);
  double meps           = 1.-eps;
  double delta          = deltafactor*eps0/(eps*meps);
  double phi1           = 0.0;
  double phi2           = 0.0;
  if (delta>1.) {
    phi1 = 21.12 - 4.184*std::log(delta+0.952);
    phi2 = phi1;
  } else {
    double delta2 = delta*delta;
    phi1 = 20.867 - 3.242*delta + 0.625*delta2;
    phi2 = 20.209 - 1.93 *delta - 0.086*delta2;
  }
  double dxsec = (eps*eps+meps*meps)*(0.25*phi1-fz) + 2.*eps*meps*(0.25*phi2-fz)/3.;
  dxsec *= eps;
  if (dxsec<0.) {
    dxsec = 0.0;
  }
  return dxsec;
}

} // namespace geantphysics
