
#include "RelativisticPairModel.h"

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

// use these elastic and inelatic form factors for light elements instead of TFM
// under the complete screening approximation
// Tsai Table.B2.
const double RelativisticPairModel::gFelLowZet  [] = {0.0, 5.310, 4.790, 4.740, 4.710, 4.680, 4.620, 4.570};
const double RelativisticPairModel::gFinelLowZet[] = {0.0, 6.144, 5.621, 5.805, 5.924, 6.012, 5.891, 5.788};


RelativisticPairModel::RelativisticPairModel(const std::string &modelname) : EMModel(modelname) {
  fIsUseTsaisScreening              = false;   // for xsec computation and rejection sampling only
  fIsUseLPM                         = true;    // photon energies above fLPMEnergyLimit energy

  fElectronInternalCode             = -1;      // will be set at init
  fPositronInternalCode             = -1;      // will be set at init

  fLPMEnergyLimit                   = 100.0*geant::GeV;  // photon energy limit above which LPM is active

  fElementData                      = nullptr;

  fNumSamplingPrimEnergiesPerDecade = 12;    // should be set/get and must be done before init
  fNumSamplingEnergies              = 84;    // should be set/get and must be done before init
  fMinPrimEnergy                    = -1.;
  fMaxPrimEnergy                    = -1.;
  fNumSamplingPrimEnergies          = -1;      // will be set in InitSamplingTables if needed
  fPrimEnLMin                       =  0.;     // will be set in InitSamplingTables if needed
  fPrimEnILDelta                    =  0.;     // will be set in InitSamplingTables if needed
  fSamplingPrimEnergies             = nullptr; // will be set in InitSamplingTables if needed
  fLSamplingPrimEnergies            = nullptr; // will be set in InitSamplingTables if needed

  //fRatinAliasDataForAllMaterials;
  fAliasSampler                     = nullptr;
}


RelativisticPairModel::~RelativisticPairModel() {}

void   RelativisticPairModel::Initialize() {
  EMModel::Initialize();  // will set the PhysicsParameters member
  fElectronInternalCode = Electron::Definition()->GetInternalCode();
  fPositronInternalCode = Positron::Definition()->GetInternalCode();
  if (GetLowEnergyUsageLimit()<50.*geant::MeV) {
    std::cerr<< "  *** ERROR: RelativisticPairModel::Initialize()  \n"
             << "      The model should not be used below 50 [MeV]   "
             << std::endl;
    exit(-1);
  }
  InitialiseModel();
  // In case of rejection sampling request to build target element selector for this model (to sample target atom)
  // (for gamma particle and element selectors per material)
  if (!GetUseSamplingTables()) {
    InitialiseElementSelectors(this,Gamma::Definition(),true);
  }
}

/**
 * @internal
 *  A collection of frequently used target atom dependent variables is set up bu calling the InitialiseElementData()
 *  method. If sampling tables were requested to be used, then the min/max of the primary photon energy grid (which
 *  sampling tables will be built over) are set and sampling tables are initialized by calling the InitSamplingTables()
 *  method.
 * @endinternal
 */
void RelativisticPairModel::InitialiseModel() {
  InitialiseElementData();
  // everything from this line is for building sampling tables:
  fMinPrimEnergy = GetLowEnergyUsageLimit();
  fMaxPrimEnergy = GetHighEnergyUsageLimit();
  if (GetUseSamplingTables()) {
//    std::cerr<< "  === Rel. pair model: building sampling tables"<< std::endl;
    InitSamplingTables();
  }
}


double RelativisticPairModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *part) {
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
    xsec += theAtomicNumDensityVector[iel]*ComputeXSectionPerAtom(theElements[iel], matcut, egamma, part);
  }
  return xsec;
}


double RelativisticPairModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts* matcut, double kinenergy,
                                                     const Particle*) {
   double xsec  = 0.0;
   if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
     return xsec;
   }
   // compute the parametrized atomic cross section: depends only on target Z and gamma energy.
   xsec = ComputeAtomicCrossSection(elem, matcut->GetMaterial(), kinenergy);
   return xsec;
}


int    RelativisticPairModel::SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td) {
  int    numSecondaries      = 0;
  double ekin                = track.GetKinE();
  double eps0                = geant::kElectronMassC2/ekin;
  // check if kinetic energy is below fLowEnergyUsageLimit (its minimum is 50 [MeV]) and do nothing if yes;
  // check if kinetic energy is above fHighEnergyUsageLimit andd o nothing if yes
  if (ekin<GetLowEnergyUsageLimit() || ekin>GetHighEnergyUsageLimit()) {
    return numSecondaries;
  }
  // interaction is possible so sample the sample the reduced total energy transfered to one of the secondary
  // particles i.e. either e- or e+ using sampling tables (if it was requested before initialisation) or with
  // rejection
  const MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[track.GetMaterialCutCoupleIndex()];
  const Material        *mat = matCut->GetMaterial();
  double eps = 0.0;
  if (GetUseSamplingTables()) {  // Sampling eps from table
    double *rndArray  = td->fDblArray;
    td->fRndm->uniform_array(4, rndArray);
    eps = SampleTotalEnergyTransfer(ekin, mat->GetIndex(), rndArray[0], rndArray[1], rndArray[2]);
  } else {                       // Sampling eps with rejection
    constexpr double lpmConstant = 0.25*geant::kFineStructConst*geant::kElectronMassC2*geant::kElectronMassC2/(geant::kPi*geant::kHBarPlanckCLight);
    double lpmEnergy = mat->GetMaterialProperties()->GetRadiationLength()*lpmConstant;
    // sample target element
    const Vector_t<Element*> &theElements = mat->GetElementVector();
    double targetElemIndx = 0;
    if (theElements.size()>1) {
      targetElemIndx = SampleTargetElementIndex(matCut, ekin, td->fRndm->uniform());
    }
    double  zet     = theElements[targetElemIndx]->GetZ();
    int    izet     = std::lrint(zet);
    double deltaFac = fElementData[izet]->fDeltaFactor;
    double deltaMin = 4.*eps0*deltaFac;
    double deltaMax = fElementData[izet]->fDeltaMax;
    if (fIsUseTsaisScreening) {
      deltaMax = fElementData[izet]->fDeltaMaxTsai;
    }
    double epsp   = 0.5-0.5*std::sqrt(1.-deltaMin/deltaMax);
    double epsMin = std::max(eps0,epsp);
    //
    double fz     = fElementData[izet]->fFz;
    double z23    = theElements[targetElemIndx]->GetElementProperties()->GetZ23();
    eps = SampleTotalEnergyTransfer(epsMin, eps0, fz, z23, ekin, lpmEnergy, deltaFac, td);
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



double RelativisticPairModel::SampleTotalEnergyTransfer(double egamma, int matindx, double r1, double r2, double r3) {
  // determine primary energy lower grid point
  double lGammaEnergy  = std::log(egamma);
  int gammaEnergyIndx  = (int) ((lGammaEnergy-fPrimEnLMin)*fPrimEnILDelta);
  //
  if (gammaEnergyIndx>=fNumSamplingPrimEnergies-1)
    gammaEnergyIndx = fNumSamplingPrimEnergies-2;

  double pLowerGammaEner = (fLSamplingPrimEnergies[gammaEnergyIndx+1]-lGammaEnergy)*fPrimEnILDelta;
  if (r1>pLowerGammaEner) {
     ++gammaEnergyIndx;
  }
  // get the RatinAliasData object pointer for the given material index and gamma-energy-index(that is gammaEnergyIndx)
  // fRatinAliasDataForAllMaterials[matindx] should not be nullptr because it was checked in the caller
  RatinAliasData *raData = fRatinAliasDataForAllMaterials[matindx]->fRatinAliasDataForOneMaterial[gammaEnergyIndx];
  //
  int    izet     = fRatinAliasDataForAllMaterials[matindx]->fILowestZ;
  double eps0     = geant::kElectronMassC2/egamma;
  double deltaMax = fElementData[izet]->fDeltaMaxTsai;
  double epsp     = 0.5-0.5*std::sqrt(1.-4.*eps0*fElementData[izet]->fDeltaFactor/deltaMax);
  double epsMin   = std::max(eps0,epsp);
  // we could put a static assert here to make sure that raData is not nullptr
  //
  // sample the transformed variable xi=... (such that linear aprx will be used in the first interval)
  double  xi = fAliasSampler->SampleRatin(raData->fXdata, raData->fCumulative, raData->fParaA, raData->fParaB,
                                          raData->fAliasW, raData->fAliasIndx, raData->fNumdata, r2, r3, 0);
  // transform back xi to eps = E_total_energy_transfer/E_{\gamma}
  return epsMin*std::exp(xi*std::log(0.5/epsMin));
}



// Sample totel energy (fraction) transfered to one of the e-/e+ pair with REJECTION: i.e. first the target atom
// needs to be sampled and there is an option if Tsai's screening is used or not: if Tsai's screening is used then
// epsmin, and so on must be evaluated with Tsai's screening
double RelativisticPairModel::SampleTotalEnergyTransfer(double epsmin, double eps0, double fz, double z23, double egamma,
                                                        double lpmenergy, double deltafactor, Geant::GeantTaskData *td) {
    double eps       = 0.0;
    double epsRange  = 0.5-epsmin;
    double deltamin  = 4.*deltafactor*eps0;
    //
    double F10       = ScreenFunction1(deltamin,fIsUseTsaisScreening) - fz;
    double F20       = ScreenFunction2(deltamin,fIsUseTsaisScreening) - fz;
    double NormF1    = std::max(F10*epsRange*epsRange,0.);
    double NormF2    = std::max(1.5*F20,0.);
    //
    double *rndArray = td->fDblArray;
    double greject   = 0.0;
    do {
      td->fRndm->uniform_array(3, rndArray);
      if (NormF1/(NormF1+NormF2)>rndArray[0]) {
      	eps            = 0.5 - epsRange*std::pow(rndArray[1],1./3.);
	      double delta   = deltafactor*eps0/(eps*(1.-eps));
        if (fIsUseLPM && egamma>fLPMEnergyLimit) {
          double funcPhiS, funcGS, funcXiS, phi1, phi2;
          ComputeScreeningFunctions(phi1, phi2, delta, fIsUseTsaisScreening);
          ComputeLPMFunctions(funcPhiS, funcGS, funcXiS, z23, lpmenergy, eps, egamma);
          greject = funcXiS*((funcGS+2.*funcPhiS)*phi1 - funcGS*phi2 - funcPhiS*fz)/F10;
        } else {
	        greject = (ScreenFunction1(delta,fIsUseTsaisScreening) - fz)/F10;
        }
      } else {
	      eps            = epsmin + epsRange*rndArray[1];
        double delta   = deltafactor*eps0/(eps*(1.-eps));
        if (fIsUseLPM && egamma>fLPMEnergyLimit) {
          double funcPhiS, funcGS, funcXiS, phi1, phi2;
          ComputeScreeningFunctions(phi1, phi2, delta, fIsUseTsaisScreening);
          ComputeLPMFunctions(funcPhiS, funcGS, funcXiS, z23, lpmenergy, eps, egamma);
          greject = funcXiS*((0.5*funcGS+funcPhiS)*phi1 + 0.5*funcGS*phi2 - 0.5*(funcGS+funcPhiS)*fz)/F20;
        } else {
	        greject = (ScreenFunction2(delta,fIsUseTsaisScreening) - fz)/F20;
        }
      }
    } while (greject<rndArray[2]);
  return eps;
}


double RelativisticPairModel::ComputeAtomicCrossSection(const Element *elem, const Material *mat, double egamma) {
  constexpr double constFactor = geant::kFineStructConst*geant::kClassicElectronRadius*geant::kClassicElectronRadius;
  constexpr double lpmConstant = 0.25*geant::kFineStructConst*geant::kElectronMassC2*geant::kElectronMassC2/(geant::kPi*geant::kHBarPlanckCLight);
  double lpmEnergy   = mat->GetMaterialProperties()->GetRadiationLength()*lpmConstant;
  double eps0        = geant::kElectronMassC2/egamma;
  double zet         = elem->GetZ();
  int izet           = std::lrint(zet);

  double deltaMax = fElementData[izet]->fDeltaMax;
  if (fIsUseTsaisScreening) {
    deltaMax = fElementData[izet]->fDeltaMaxTsai;
  }
  double epsp   = 0.5-0.5*std::sqrt(1.-4.*fElementData[izet]->fDeltaFactor*eps0/deltaMax);
  double epsMin = std::max(eps0,epsp);//+1.e-14;

  double fz     = fElementData[izet]->fFz;
  double df     = fElementData[izet]->fDeltaFactor;
  double eta    = fElementData[izet]->fEtaValue;
  double z23    = elem->GetElementProperties()->GetZ23();

  int ngl = 64;
  GLIntegral *gl = new GLIntegral(ngl,0.0,1.0);
  const std::vector<double> glW = gl->GetWeights();
  const std::vector<double> glX = gl->GetAbscissas();
  double xsec      = 0.0;
  double factor    = zet*(zet+eta);
  for (int i=0; i<ngl; ++i) {
    double xi  = glX[i];
    if (fIsUseLPM && egamma>fLPMEnergyLimit) {
      xsec += glW[i]*ComputeLPMDXsectionPerAtom(epsMin,df,fz,lpmEnergy,z23,egamma,xi,fIsUseTsaisScreening); // istsai = false by def.
    } else {
      xsec += glW[i]*ComputeDXsectionPerAtom(epsMin,eps0,df,fz,xi,fIsUseTsaisScreening); // istsai = false by def.
    }
  }
  xsec = std::max(0.0,2.*constFactor*factor*xsec);
  return xsec;
}


// transformd no - LPM
double RelativisticPairModel::ComputeDXsectionPerAtom(double epsmin, double eps0, double df, double fz, double xi, bool istsai) {
  double lHalfPerEpsMin = std::log(0.5/epsmin);
  double eps            = epsmin*std::exp(xi*lHalfPerEpsMin);
  double meps           = 1.-eps;
  double delta          = df*eps0/(eps*meps);
  double phi1           = 0.0;
  double phi2           = 0.0;
  ComputeScreeningFunctions(phi1,phi2,delta,istsai);
  double dxsec = (eps*eps+meps*meps)*(phi1-0.5*fz) + 2.*eps*meps*(phi2-0.5*fz)/3.;
  dxsec *= (eps*lHalfPerEpsMin);
  dxsec = std::max(dxsec,0.0);
  return dxsec;
}

// transformd with - LPM
double RelativisticPairModel::ComputeLPMDXsectionPerAtom(double epsmin, double df, double fz, double lpmenergy, double z23, double egamma, double xi, bool istsai) {
  double lHalfPerEpsMin = std::log(0.5/epsmin);
  double eps            = epsmin*std::exp(xi*lHalfPerEpsMin);
  double meps           = 1.-eps;
  double delta          = df*geant::kElectronMassC2/(egamma*eps*meps);
  double phi1           = 0.0;
  double phi2           = 0.0;
  double lpmPhiS        = 0.0;
  double lpmGS          = 0.0;
  double lpmXiS         = 0.0;
  ComputeScreeningFunctions(phi1,phi2,delta,istsai);
  ComputeLPMFunctions(lpmPhiS, lpmGS, lpmXiS, z23, lpmenergy, eps, egamma);
  double dxsec = ((lpmGS + 2.*lpmPhiS)/3.)*(eps*eps+meps*meps)*(phi1-0.5*fz) + 2.*lpmGS*eps*meps*(phi2-0.5*fz)/3.;
  dxsec *= (eps*lHalfPerEpsMin*lpmXiS);
  dxsec = std::max(dxsec,0.0);
  return dxsec;
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
void RelativisticPairModel::InitialiseElementData() {
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
        double zet = theElements[j]->GetZ();
        int elementIndx = std::lrint(zet);
//        std::cerr<< " === Building ElementData for " << theElements[j]->GetName() << std::endl;
        if (!fElementData[elementIndx]) {
          Element         *elem   = theElements[j];
          ElementData *elemData   = new ElementData();
          elemData->fDeltaFactor  = 136./elem->GetElementProperties()->GetZ13(); // 136/pow(z,1/3)
          elemData->fCoulombCor   = elem->GetElementProperties()->GetCoulombCorrection();
          elemData->fFz           = 8.*elem->GetElementProperties()->GetLogZ13()+8*elemData->fCoulombCor;
          elemData->fDeltaMax     = std::exp((42.24-elemData->fFz)/8.368)-0.952;
          elemData->fDeltaMaxTsai = 1.36*std::sqrt( std::exp(0.5*20.863-2.-0.25*elemData->fFz)-1. )/0.55846;
          double Fel   = 0.;
          double Finel = 0.;
          if (elementIndx<5) {
            Fel   = gFelLowZet[elementIndx];
            Finel = gFinelLowZet[elementIndx];
          } else {
            Fel   = std::log(184.15)  -    elem->GetElementProperties()->GetLogZ()/3.;
            Finel = std::log(1194.)   - 2.*elem->GetElementProperties()->GetLogZ()/3.;
          }
          elemData->fEtaValue  = Finel/(Fel-elemData->fCoulombCor);
          fElementData[elementIndx] = elemData;
        }
      }
    }
  }
}


void RelativisticPairModel::InitSamplingTables() {
  // set number of primary gamma energy grid points
  // keep the prev. value of primary energy grid points.
  int oldNumGridPoints     = fNumSamplingPrimEnergies;
  fNumSamplingPrimEnergies = fNumSamplingPrimEnergiesPerDecade*std::lrint(std::log10(fMaxPrimEnergy/fMinPrimEnergy))+1;
  if (fNumSamplingPrimEnergies<3) {
    fNumSamplingPrimEnergies = 3;
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
  // clear (if needed) and create the container that can store sampling tables for all(active) materials
  // 1. clear:
  for (size_t i=0; i<fRatinAliasDataForAllMaterials.size(); ++i) {
    RatinAliasDataPerMaterial *oneMat = fRatinAliasDataForAllMaterials[i];
    if (oneMat) {
      for (int j=0; j<oldNumGridPoints; ++j) {
        RatinAliasData *oneRatinAlias = oneMat->fRatinAliasDataForOneMaterial[j];
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
      delete [] oneMat->fRatinAliasDataForOneMaterial;
      delete oneMat;
    }
  }
  fRatinAliasDataForAllMaterials.clear();
  //
  // 2. create an AliasTable object
  if (fAliasSampler) {
    delete fAliasSampler;
  }
  fAliasSampler = new AliasTable();
  //
  // 3. create and fill with nullptr-s
  fRatinAliasDataForAllMaterials.resize(Material::GetTheMaterialTable().size(),nullptr);
  //
  // get all MaterialCuts; loop over them and if they belongs to a region where this model is active:
  //  - get the material that the MaterialCuts belongs to
  //  - build RatinAliasDataPerMaterial for this material if it has not been done yet
  int numMatCuts = MaterialCuts::GetTheMaterialCutsTable().size();
  // get list of active region
  std::vector<bool> isActiveInRegion = GetListActiveRegions();
  for (int i=0; i<numMatCuts; ++i) {
    const MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[i];
    // if this MaterialCuts belongs to a region where this model is active:
    if (isActiveInRegion[matCut->GetRegionIndex()]) {
      // get the material
      const Material *theMaterial =  matCut->GetMaterial();
      // if RatinAliasDataPerMaterial has not been built for this material yet => build
      if (!fRatinAliasDataForAllMaterials[theMaterial->GetIndex()]) {
        BuildSamplingTablesForMaterial(theMaterial);
      }
    }
  }
}


void RelativisticPairModel::BuildSamplingTablesForMaterial(const Material *mat) {
  // prepare one RatinAliasDataPerMaterial structure
  RatinAliasDataPerMaterial *perMat          = new RatinAliasDataPerMaterial();
  perMat->fRatinAliasDataForOneMaterial      = new RatinAliasData*[fNumSamplingPrimEnergies];
  fRatinAliasDataForAllMaterials[mat->GetIndex()] = perMat;
  // set lowest Z
  const Vector_t<Element*> theElements = mat->GetElementVector();
  int    numElems = theElements.size();
  double lowZ    = 200.;
  for (int iel=0; iel<numElems; ++iel) {
    double zet = theElements[iel]->GetZ();
    if (zet<lowZ) {
      lowZ = zet;
    }
  }
  perMat->fILowestZ = std::lrint(lowZ);
  // allocate an array that will be used temporary:
  double *thePdfdata = new double[fNumSamplingEnergies]();
  // loop over the photon energy grid and build one RatinAliasData structure at each photon energy for this material
  for (int i=0; i<fNumSamplingPrimEnergies; ++i) {
    // note that egamma is for sure > 50 [MeV]
    double egamma = fSamplingPrimEnergies[i];
    // build one table for the given: egamma, material
//    std::cerr<<"  === Bulding table for:  material = " << mat->GetName() << " egamma = " << egamma/geant::MeV << " [MeV]" << std::endl;
    BuildOneRatinAlias(egamma, mat, thePdfdata, i, perMat->fILowestZ);
  }
  // delete the temporary array
  delete [] thePdfdata;
}


void RelativisticPairModel::BuildOneRatinAlias(double egamma, const Material *mat, double *pdfarray, int egammaindx, int izet) {
  // compute the theoretical minimum of the reduced total energy transfer to the e+ (or to the e-)
  double eps0    = geant::kElectronMassC2/egamma;
  // set the real minimum to the theoretical one: will be updated
  double epsMin  = eps0;
  // deterime the lowest eps': that will belong to the lowest Z of the material
  double deltaMax    = fElementData[izet]->fDeltaMaxTsai;  // Tsai's screening apr. will be used for the tables
  double deltaFactor = fElementData[izet]->fDeltaFactor;
  double deltaMin    = 4.*eps0*deltaFactor;  // 4*136*Z^(-1/3)*eps0
  double eps1        = 0.5-0.5*std::sqrt(1.-deltaMin/deltaMax);
  // take the higest value between eps0 and eps'=eps1: eps will be in [epsMin,0.5]
  if (eps1>epsMin) {
    epsMin = eps1;
  }
  //
  // allocate one RatinAliasData structure for the current element,gamma-energy combination
  RatinAliasData *raData = new RatinAliasData();
  raData->fNumdata       = fNumSamplingEnergies;
  raData->fXdata         = new double[fNumSamplingEnergies]();
  raData->fCumulative    = new double[fNumSamplingEnergies]();
  raData->fParaA         = new double[fNumSamplingEnergies]();
  raData->fParaB         = new double[fNumSamplingEnergies]();
  raData->fAliasW        = new double[fNumSamplingEnergies]();
  raData->fAliasIndx     = new int   [fNumSamplingEnergies]();
  // insert possible points: at xi(epsMin)=0, xi(epsMax=0.5)=1, at the midle xi=0.5
  // and at xi(eps'(Z_i)) points if these eps'(Z_i)> epsMin
  int curNumData    = 1;
  raData->fXdata[0] = 0.0;
  const Vector_t<Element*> theElements = mat->GetElementVector();
  int numElems = theElements.size();
  for (int iel=0; iel<numElems; ++iel) {
    izet        = std::lrint(theElements[iel]->GetZ());
    deltaMax    = fElementData[izet]->fDeltaMaxTsai;  // Tsai's screening apr. will be used for the tables
    deltaFactor = fElementData[izet]->fDeltaFactor;
    deltaMin    = 4.*eps0*deltaFactor;  // 4*136*Z^(-1/3)*eps0
    eps1        = 0.5-0.5*std::sqrt(1.-deltaMin/deltaMax);
    if (eps1>epsMin) {
      raData->fXdata[curNumData] = std::log(eps1/epsMin)/std::log(0.5/epsMin);
      ++curNumData;
    }
  }
  // reorder
  for (int i=1; i<curNumData; ++i) {
    for (int j=i; j<curNumData; ++j) {
      if (raData->fXdata[j]<raData->fXdata[i]) {
        double a = raData->fXdata[i];
        raData->fXdata[i] = raData->fXdata[j];
        raData->fXdata[j] = a;
      }
    }
  }
  raData->fXdata[curNumData] = 0.5;
  raData->fXdata[curNumData+1] = 1.0;
  curNumData += 2;
  //
  //
  // create a 4 point GL integral to compute the interpolation error at each interval
  int glnum         = 4;
  //
  // fill in the discrete pdf array by inserting raData->fNumdata discrete points such that the
  // interpolation error is minimized: one new discrete point will be inserted at the end of each
  // iteration of this while loop:
  while (curNumData<raData->fNumdata) {
    // compute the pdf values at the current discrete sample points
    for (int i=0;i<curNumData;++i) {
      pdfarray[i] = ComputeDXsection(mat, egamma, epsMin, raData->fXdata[i], true);
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
        double valReal = norm*ComputeDXsection(mat, egamma, epsMin, xval, true);
        double curerr  = std::fabs((1.-valAprx/valReal));//*(valAprx-valReal));
        if (raData->fNumdata-curNumData<5) {
          curerr  = std::fabs((1.-valAprx/valReal)*(valAprx-valReal));
        }
        err += curerr;
        //err += std::fabs((1.-valAprx/valReal));//*(valAprx-valReal));
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
    pdfarray[i] = ComputeDXsection(mat, egamma, epsMin, raData->fXdata[i], true);
  }
  //
  // and prepare the final ratin-alias sampling table structure:
  fAliasSampler->PreparRatinTable(raData->fXdata, pdfarray, raData->fCumulative, raData->fParaA, raData->fParaB,
                                  raData->fAliasW, raData->fAliasIndx, raData->fNumdata, true, glnum);
  // fill in
  int matindx = mat->GetIndex();
  fRatinAliasDataForAllMaterials[matindx]->fRatinAliasDataForOneMaterial[egammaindx]=raData;
}

double RelativisticPairModel::ComputeDXsection(const Material *mat, double egamma, double epsmin, double xi, bool istsai) {
  constexpr double lpmConstant = 0.25*geant::kFineStructConst*geant::kElectronMassC2*geant::kElectronMassC2/(geant::kPi*geant::kHBarPlanckCLight);
  double lpmEnergy   = mat->GetMaterialProperties()->GetRadiationLength()*lpmConstant;
  double eps0        = geant::kElectronMassC2/egamma;
//  double logFactor   = std::log(0.5/epsmin);
//  double eps         = epsmin*std::exp(xi*logFactor);
  double dxsec       = 0.0;
  //
  const Vector_t<Element*> theElements    = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int    numElems = theElements.size();
  for (int iel=0; iel<numElems; ++iel) {
    double zet  = theElements[iel]->GetZ();
    int izet    = std::lrint(zet);
    double deltaMax = fElementData[izet]->fDeltaMax;
    if (istsai) {
      deltaMax = fElementData[izet]->fDeltaMaxTsai;
    }
//    double epsp   = 0.5-0.5*std::sqrt(1.-4.*fElementData[izet]->fDeltaFactor*eps0/deltaMax);
    double xsec   = 0.0;
//    if (epsp>=epsmin) {
      double fz     = fElementData[izet]->fFz;
      double df     = fElementData[izet]->fDeltaFactor;
      double eta    = fElementData[izet]->fEtaValue;
      double z23    = theElements[iel]->GetElementProperties()->GetZ23();
      if (fIsUseLPM && egamma>fLPMEnergyLimit) {
        xsec = ComputeLPMDXsectionPerAtom(epsmin,df,fz,lpmEnergy, z23, egamma, xi, istsai);
      } else {
        xsec = ComputeDXsectionPerAtom(epsmin,eps0,df,fz,xi,istsai);
      }
      xsec *= zet*(zet+eta);
//    }
    dxsec += theAtomicNumDensityVector[iel]*xsec;
  }
  return dxsec;
}


void RelativisticPairModel::ComputeScreeningFunctions(double &phi1, double &phi2, double delta, bool istsai) {
  if (!istsai) {
    if (delta>1.) {
      phi1 = 21.12 - 4.184*std::log(delta+0.952);
      phi2 = phi1;
    } else {
      double delta2 = delta*delta;
      phi1 = 20.867 - 3.242*delta + 0.625*delta2;
      phi2 = 20.209 - 1.93 *delta - 0.086*delta2;
    }
  } else {
    double gamma   = delta*0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    double gamma2  = gamma*gamma;
    phi1   = 20.863 - 2.0*std::log(1.0+0.311877*gamma2) - 4.0*(1.0-0.6*std::exp(-0.9*gamma)-0.4*std::exp(-1.5*gamma));
    phi2   = phi1-2.0/(3.0+19.5*gamma+18.0*gamma2);  // phi1-phi2
  }
}

// 3xPhi_1 - Phi_2: used in case of rejection (either istsai or not)
double RelativisticPairModel::ScreenFunction1(double delta, bool istsai) {
  double val;
  if (istsai) {
    double gamma   = delta*0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    double gamma2  = gamma*gamma;
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
double RelativisticPairModel::ScreenFunction2(double delta, bool istsai) {
  double val;
  if (istsai) {
    double gamma   = delta*0.735294; // 0.735 = 1/1.36 <= gamma = delta/1.36
    double gamma2  = gamma*gamma;
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


void RelativisticPairModel::ComputeLPMFunctions(double &funcPhiS, double &funcGS, double &funcXiS, double z23, double lpmenergy, double eps, double egamma) {
    // And these are the constant parts so they should be pulled out
    const double fact        = 1./(184.1499*184.1499);

    // radiation length: is a material property so it should be set by the caller
    // LPM energy is a material dependent constant so it should be set and provided by the caller

    //
    // compute the LPM functions s', s, \xi(s), G(s), \phi(s)
    //
    // DO IT ACCORDING TO Geant4 model:
    // -compute s according to Stanev to avoid iterative solution: s1, s', \xi(s') then s
    //  1. s_1 \equiv [Z^{1/3}/\exp(20.863/4)]^2 = Z^{2/3}/184.1499^2
    //     std::pow(z,2/3) should be pulled out as the constant 184.1499^2 =
    double varS1     = z23*fact;
    //  2. y = E_+/E_{\gamma} with E_+ being the total energy transfered to one of the e-/e+ pair
    //     s'  = \sqrt{ \frac{1}{8} \frac{1}{y(1-y)}   \frac{E^{KL}_{LPM}}{E_{\gamma}}  }
    double varSprime = std::sqrt(0.125*lpmenergy/(eps*egamma*(1.0-eps)));
    //  3. \xi(s') =
    //      \begin{cases}
    //        2 & \quad \text{if}\; s' \leq \sqrt{2}s_1
    //        1+h-\frac{0.08(1-h)[1-(1-h)^2]}{\ln(\sqrt{2}s_1)} & \quad \text{if}\;  \sqrt{2}s_1 < s' < 1
    //        1 & \quad \text{if}\; s' \geq 1
    //      \end{cases}
    // where h(s') \equiv (\ln(s'))/\ln(\sqrt{2}s_1)
    funcXiS   = 2.0;
    // sqrt(2) constant should be pulled out
    double condition = std::sqrt(2.0)*varS1;
    if (varSprime>1.0) {
      funcXiS = 1.0;
    } else if (varSprime>condition) {
      double dum0        = 1.0/std::log(condition);
      double funcHSprime = std::log(varSprime)*dum0;
      funcXiS            = 1.0+funcHSprime-(0.08*(1.0-funcHSprime)*(2.0*funcHSprime-funcHSprime*funcHSprime))*dum0;
    }

    //  4. s=\frac{s'}{\sqrt{\xi(s')}}
    double varShat    = varSprime/std::sqrt(funcXiS);

    // - compute \phi(s) and G(s) suppression functions according to Stanev (Geant4 approximations are slightly different)
    //   by replacing s with \hat{s} = s*(1+k_p^2/k^2)
    //  \f[
    //  \phi(s) \approx
    //   \begin{cases}
    //    6s(1-\pi s) & \quad \text{if}\; s \leq 0.01 \\
    //    1-\exp \left\{  -6s[1+s(3-\pi)] + \frac{s^3}{0.623+0.796s+0.658s^2} \right\}  & \quad \text{if}\; 0.01 < s < 2 \\
    //    1 & \quad \text{if}\; s \geq 2
    //   \end{cases}
    //  \f]
    //
    //  \f[
    //  G(s) = 3 \psi(s) - 2 \phi(s)
    //  \f]
    //  where
    //  \f[
    //  \psi(s) \approx
    //   \begin{cases}
    //    4s  & \quad \text{if}\; s \leq 0.01 \\
    //    1-\exp \left\{  -4s  - \frac{8s^2}{1+3.936s+4.97s^2-0.05s^3+7.5s^4} \right\}  & \quad \text{if}\; 0.01 < s < 2 \\
    //    1 & \quad \text{if}\; s \geq 2
    //   \end{cases}
    //  \f]
    funcPhiS = 6.0*varShat*(1.0-geant::kPi*varShat);
    double funcPsiS = 4.0*varShat;
    if (varShat>=2.0) {
      funcPhiS = 1.0;
      funcPsiS = 1.0;
    } else if (varShat>0.01) {
      double varShat2 = varShat*varShat;
      double varShat3 = varShat*varShat2;
      double varShat4 = varShat2*varShat2;
      funcPhiS = 1.0-std::exp(-6.0*varShat*(1.0+varShat*(3.0-geant::kPi)) + varShat3/(0.623+0.796*varShat+0.658*varShat2));
      funcPsiS = 1.0-std::exp(-4.0*varShat - 8.0*varShat2/(1.0+3.936*varShat+4.97*varShat2-0.05*varShat3+7.5*varShat4));
    }
    funcGS   = 3.0*funcPsiS-2.0*funcPhiS;

/*
  //// begin according to the genat4 model
    // - compute \phi(s) and G(s) suppression functions according Geant4 model (slightly different than Stanev)
    double varShat2 = varShat*varShat;
    double varShat3 = varShat*varShat2;
    double varShat4 = varShat2*varShat2;
    funcPhiS = 1.0;
    funcGS   = 1.0;
    if (varShat<0.1) { // high suppression limit:
      // 24s^3\sum_{j=1}^{\infty} \frac{1}{(j+s)^2+s^2} = 12 x^2 i {DigammaFunc[(1-i)s+1] - DigammaFunc[(1+i)s+1]}
      // where DigammaFunc is the digamma function
      // using Taylor approximation of  6s(1-\pi s) + 12 x^2 i {DigammaFunc[(1-i)s+1] - DigammaFunc[(1+i)s+1]} around s=0
      // ==> 6s - 6\pi s^2 + 4\pi^2 s^3 - 48 Zeta(3) s^4
      // where Zeta is the Riemann zeta function Zeta[3] \approx 1.202056903159594 is Apery's constant
      funcPhiS = 6.0*varShat-18.84955592153876*varShat2+39.47841760435743*varShat3-57.69873135166053*varShat4;
      // 48s^3\sum_{j=0}^{\infty} \frac{1}{(j+s+1/2)^2+s^2} = 24 x^2 i {DigammaFunc[(1-i)s+0.5] - DigammaFunc[(1+i)s+0.5]}
      // using Taylor approximation of 12\pi s^2 - 24 x^2 i {DigammaFunc[(1-i)s+0.5] - DigammaFunc[(1+i)s+0.5]} aound s=0
      // ==> 12\pi s^2 -24\pi^2 s^3 - 48 [Digamma''(0.5)] s^4
      // where [Digamma''(0.5)] is the second derivative of the digamma function at 0.5 and 48*[Digamma''(0.5)] \approx  -807.782238923247359
      funcGS   = 37.69911184307752*varShat2-236.8705056261446*varShat3+807.78223892324736*varShat4;
    } else if (varShat<1.9516) { // intermediate suppression: use Stanev approximation for \phi(s)
      // 1-\exp \left\{  -6s[1+s(3-\pi)] + \frac{s^3}{0.623+0.796s+0.658s^2} \right\}
      funcPhiS = 1.0-std::exp(-6.0*varShat*(1.0+varShat*(3.0-geant::kPi)) + varShat3/(0.623+0.796*varShat+0.658*varShat2));
      if (varShat<0.415827397755) { // use Stanev approximation: for \psi(s) and compute G(s)
        // 1-\exp \left\{  -4s  - \frac{8s^2}{1+3.936s+4.97s^2-0.05s^3+7.5s^4} \right\}
        double funcPsiS = 1.0-std::exp(-4.0*varShat - 8.0*varShat2/(1.0+3.936*varShat+4.97*varShat2-0.05*varShat3+7.5*varShat4));
        // G(s) = 3 \psi(s) - 2 \phi(s)
        funcGS = 3.0*funcPsiS - 2.0*funcPhiS;
      } else { // use alternative parametrisation
        double dum0 = -0.16072300849123999+3.7550300067531581*varShat-1.7981383069010097*varShat2+0.67282686077812381*varShat3-0.1207722909879257*varShat4;
        funcGS = std::tanh(dum0);
      }
    } else { //low suppression limit:
      // 24s^3\sum_{j=1}^{\infty} \frac{1}{(j+s)^2+s^2} = 12 x^2 i {DigammaFunc[(1-i)s+1] - DigammaFunc[(1+i)s+1]}
      // usign Laurent series expansion of 6s(1-\pi s) + 12 x^2 i {DigammaFunc[(1-i)s+1] - DigammaFunc[(1+i)s+1]} around s=\infty
      // ==> 1-1/(84s^4)
      funcPhiS = 1.0-0.01190476/varShat4;
      // 48s^3\sum_{j=0}^{\infty} \frac{1}{(j+s+1/2)^2+s^2} = 24 x^2 i {DigammaFunc[(1-i)s+0.5] - DigammaFunc[(1+i)s+0.5]}
      // usign Laurent series expansion of 12\pi s^2 - 24 x^2 i {DigammaFunc[(1-i)s+0.5] - DigammaFunc[(1+i)s+0.5]} aound s=\infty
      // ==> 1-0.0230655/s^4
      funcGS   = 1.0-0.0230655/varShat4;
    }

*/

    //MAKE SURE SUPPRESSION IS SMALLER THAN 1: due to Migdal's approximation on xi
    if (funcXiS*funcPhiS>1. || varShat>0.57) {
      funcXiS=1./funcPhiS;
    }
}



}    // namespace geantphysics
