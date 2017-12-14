
#include "SeltzerBergerBremsModel.h"
// from material
#include "Types.h"


#include "PhysicalConstants.h"

#include "PhysicalConstants.h"
#include "Material.h"
#include "Element.h"
#include "MaterialProperties.h"

#include "MaterialCuts.h"

#include "Spline.h"
#include "GLIntegral.h"
#include "AliasTable.h"

#include "PhysicsParameters.h"

#include "Gamma.h"

#include "LightTrack.h"
#include "PhysicsData.h"

// from geantV
#include "GeantTaskData.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>

namespace geantphysics {

void SeltzerBergerBremsModel::Initialize() {
  EMModel::Initialize();
  Initialise();
  // if it needs element selector: particle is coded in the model in case of this model and due to the cut dependence
  // we we need element selectors per MaterialCuts
  //  InitialiseElementSelectors(this, nullptr, false);
}

double SeltzerBergerBremsModel::ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle * /*particle*/,
                                            bool istotal){
  const Material *mat =  matcut->GetMaterial();
  const double *cuts  =  matcut->GetProductionCutsInEnergy();
  double gammacut     =  cuts[0];
  if (istotal) {
    // for the total stopping power we just need a gamma production cut >=kinenergy
    gammacut = 1.01*kinenergy;
  }
  return ComputeDEDXPerVolume(mat, gammacut, kinenergy);
}


double SeltzerBergerBremsModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy,
                                                           const Particle * /*particle*/) {
  double xsec = 0.0;
  if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
    return xsec;
  }
  const Material *mat =  matcut->GetMaterial();
  const double *cuts  =  matcut->GetProductionCutsInEnergy();
  double gammacut     =  cuts[0];
  xsec = ComputeXSectionPerVolume(mat, gammacut, kinenergy);
  return xsec;
}


double SeltzerBergerBremsModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut, double kinenergy, const Particle*) {
  double xsec = 0.0;
  if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
    return xsec;
  }
  const Material *mat =  matcut->GetMaterial();
  const double *cuts  =  matcut->GetProductionCutsInEnergy();
  double gammacut     =  cuts[0];
  xsec = ComputeXSectionPerAtom(elem, mat, gammacut, kinenergy);
  return xsec;
}


int SeltzerBergerBremsModel::SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td) {
  int    numSecondaries      = 0;
  double ekin                = track.GetKinE();
  const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(track.GetMaterialCutCoupleIndex());
  const double *cuts         = matCut->GetProductionCutsInEnergy();
  double gammacut            = cuts[0];
  // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
  // check if kinetic energy is above fHighEnergyUsageLimit andd o nothing if yes
  // check if kinetic energy is below gamma production cut and do nothing if yes
  if (ekin<GetLowEnergyUsageLimit() || ekin>GetHighEnergyUsageLimit() || ekin<=gammacut) {
    return numSecondaries;
  }
  // sample gamma energy
  // here we need 3 random number + 2 later for photon direction theta phi sampling
  double *rndArray = td->fDblArray;
  td->fRndm->uniform_array(5, rndArray);
  double gammaEnergy  = SamplePhotonEnergy(matCut, ekin, rndArray[0], rndArray[1], rndArray[2]);
  // sample gamma scattering angle in the scattering frame i.e. which z-dir points to the orginal e-/e+ direction
  double cosTheta = 1.0;
  double sinTheta = 0.0;
  SamplePhotonDirection(ekin, sinTheta, cosTheta, rndArray[3]);
  double phi      = geant::kTwoPi*(rndArray[4]);
  // gamma direction in the scattering frame
  double gamDirX  = sinTheta*std::cos(phi);
  double gamDirY  = sinTheta*std::sin(phi);
  double gamDirZ  = cosTheta;
  // rotate gamma direction to the lab frame:
  RotateToLabFrame(gamDirX, gamDirY, gamDirZ, track.GetDirX(), track.GetDirY(), track.GetDirZ());
  // create the secondary partcile i.e. the gamma
  numSecondaries = 1;
  // NO it can be dropped if we make sure that these secondary vectors are at least size 2
//  PhysicsData *physData = td->fPhysicsData;
  // current capacity of secondary track container
  int curSizeOfSecList = td->fPhysicsData->GetSizeListOfSecondaries();
  // currently used secondary tracks in the container
  int curNumUsedSecs   = td->fPhysicsData->GetNumUsedSecondaries();
  if (curSizeOfSecList-curNumUsedSecs<numSecondaries) {
    td->fPhysicsData->SetSizeListOfSecondaries(2*curSizeOfSecList);
  }
  int secIndx = curNumUsedSecs;
  curNumUsedSecs +=numSecondaries;
  td->fPhysicsData->SetNumUsedSecondaries(curNumUsedSecs);
  // this is known since it is a secondary track
//  sectracks[secIndx].SetTrackStatus(LTrackStatus::kNew); // to kew
  std::vector<LightTrack>& sectracks = td->fPhysicsData->GetListOfSecondaries();
  // this is known since it is a secondary track
//  sectracks[secIndx].SetTrackStatus(LTrackStatus::kNew); // to kew
  sectracks[secIndx].SetDirX(gamDirX);
  sectracks[secIndx].SetDirY(gamDirY);
  sectracks[secIndx].SetDirZ(gamDirZ);
  sectracks[secIndx].SetKinE(gammaEnergy);
  sectracks[secIndx].SetGVcode(fSecondaryInternalCode);  // gamma GV code
  sectracks[secIndx].SetMass(0.0);
  sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); // parent GeantTrack index
// these are known from the parent GeantTrack
//  sectracks[secIndx].SetMaterialCutCoupleIndex(track.GetMaterialCutCoupleIndex());
//  sectracks[secIndx].SetNumOfInteractionLegthLeft(-1.); // i.e. need to sample in the step limit
//  sectracks[secIndx].SetInvTotalMFP(0.);
//  sectracks[secIndx].SetStepLength(0.);
//  sectracks[secIndx].SetEnergyDeposit(0.);
//  sectracks[secIndx].SetTime(??);
//  sectracks[secIndx].SetWeight(??);
//  sectracks[secIndx].SetProcessIndex(-1); // ???
//  sectracks[secIndx].SetTargetZ(-1);
//  sectracks[secIndx].SetTargetN(-1);
  //
  // compute the primary e-/e+ post interaction direction: from momentum vector conservation
//  double elInitTotalEnergy   = ekin+geant::kElectronMassC2;  // initial total energy of the e-/e+
  double elInitTotalMomentum = std::sqrt(ekin*(ekin+2.0*geant::kElectronMassC2));
  // final momentum of the e-/e+ in the lab frame
  double elDirX = elInitTotalMomentum*track.GetDirX() - gammaEnergy*gamDirX;
  double elDirY = elInitTotalMomentum*track.GetDirY() - gammaEnergy*gamDirY;
  double elDirZ = elInitTotalMomentum*track.GetDirZ() - gammaEnergy*gamDirZ;
  // normalisation
  double norm  = 1.0/std::sqrt(elDirX*elDirX + elDirY*elDirY + elDirZ*elDirZ);
  // update primary track direction
  track.SetDirX(elDirX*norm);
  track.SetDirY(elDirY*norm);
  track.SetDirZ(elDirZ*norm);
  // update primary track kinetic energy
  track.SetKinE(ekin-gammaEnergy);
  // return with number of secondaries i.e. 1 gamma
  return numSecondaries;
}


double SeltzerBergerBremsModel::MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle*) const {
  double mine = (matcut->GetProductionCutsInEnergy())[0]; // gamma production cut in the given material-cuts
  return mine;
}



SeltzerBergerBremsModel::SeltzerBergerBremsModel(bool iselectron, int datafileindx, const std::string &modelname)
: EMModel(modelname), fIsElectron(iselectron), fDataFileIndx(datafileindx) {
  fNGL                             = 64;
  fDCSMaxZet                       = 0;
  fLoadDCSNumElectronEnergies      = 0;
  fLoadDCSReducedPhotonEnergyGrid  = 0;
  fLoadDCSElectronEnergyGrid       = nullptr;
  fLoadDCSReducedPhotonEnergyGrid  = nullptr;
  fLoadDCSForElements              = nullptr;


  fNumSamplingElecEnergies = 71;// 171=>25 per decada//71; // between the min/max e-/e+ kinetic energies
  fNumSamplingPhotEnergies = 64; // at each energy grid points
  fMinElecEnergy           =  1.0*geant::keV; // minimum kinetic energy of the interacting e-/e+
  fMaxElecEnergy           = 10.0*geant::GeV; // maximum kinetic energy of the interacting e-/e+
  fElEnLMin                = 0.0;
  fElEnILDelta             = 1.0;
  fSamplingElecEnergies    = nullptr;
  fLSamplingElecEnergies   = nullptr;

  fNumDifferentMaterialGCuts = 0;
  fGlobalMatGCutIndxToLocal = nullptr;
  fAliasData                = nullptr; //alias data for each matrial-gammacut pairs
  fAliasSampler             = nullptr;

  fSecondaryInternalCode    = -1;

  fGL                       = nullptr;

  //Initialise();
}

SeltzerBergerBremsModel::~SeltzerBergerBremsModel() {
  if (fLoadDCSElectronEnergyGrid)
    delete [] fLoadDCSElectronEnergyGrid;

  if (fLoadDCSReducedPhotonEnergyGrid)
    delete [] fLoadDCSReducedPhotonEnergyGrid;

  if (fLoadDCSForElements) {
    for (int i=0; i<fDCSMaxZet; ++i)
      if (fLoadDCSForElements[i])
        delete [] fLoadDCSForElements[i];
    delete [] fLoadDCSForElements;
  }

  if (fSamplingElecEnergies)
    delete [] fSamplingElecEnergies;
  if (fLSamplingElecEnergies)
    delete [] fLSamplingElecEnergies;

  if (fGlobalMatGCutIndxToLocal)
    delete [] fGlobalMatGCutIndxToLocal;

  if (fAliasData) {
    for (int i=0; i<fNumDifferentMaterialGCuts; ++i) {
      for (int j=0; j<fNumSamplingElecEnergies; ++j) {
        int indx = i*fNumSamplingElecEnergies+j;
        if (fAliasData[indx]) {
          delete [] fAliasData[indx]->fXdata;
          delete [] fAliasData[indx]->fYdata;
          delete [] fAliasData[indx]->fAliasW;
          delete [] fAliasData[indx]->fAliasIndx;
          delete fAliasData[indx];
        }
      }
    }
    delete [] fAliasData;
  }
  if (fAliasSampler) {
    delete fAliasSampler;
  }
  if (fGL) {
    delete fGL;
  }

}



void SeltzerBergerBremsModel::Initialise() {
  LoadDCSData();
  if (!fGL) {
    fGL = new GLIntegral(fNGL,0.0,1.0);
  }
  InitSamplingTables();
  fAliasSampler          = new AliasTable();
  fSecondaryInternalCode = Gamma::Definition()->GetInternalCode();
//  std::cerr<<"  ----> SB model = " << GetName() << "  is initialized!"<< std::endl;
}

/**
  *  The sampling is based on the sampling tables prepared at initialisation. Statistical interpolation is used to
  *  select one of the incident particle kinetic energy grid points out of \f$ E_i \leq E_{kin} < E_{i+1}\f$ (linear
  *  interpolation in log kinetic energy) at a given primary particle kinetic energy \f$E_{kin}\f$. Then the transformed
  *  variable \f$\xi\in[0,1]\f$ is sampled from the sampling table (prepared at initialisation) that belongs to the
  *  selected incident particle kinetic energy grid point. The emitted photon energy \f$k\f$ then is obtained by
  *  applying the following transformation:
  *  \f[
  *     k = \sqrt{[k_c^2+k_p^2]\exp{ \left[ \xi\ln\frac{E_{kin}^2+k_p^2}{k_c^2+k_p^2}\right]}-k_p^2}
  *  \f]
  *  where \f$E_{kin}\f$ is the current incident particle (i.e. e-/e+) kinetic energy, \f$k_c\f$ is the current gamma
  *  particle kinetic energy production threshold and \f$ k_p = \hbar \omega_p \gamma= \hbar \omega_p E/(mc^2)\f$,
  *  (\f$E\f$ is the total energy of the incident particle)
  *  \f$\hbar \omega_p= \hbar c \sqrt{4\pi n_e r_e}\f$ (\f$n_e\f$ is the electron density of the current material and
  *  \f$r_e\f$ is the classical electron radius).
  */
double SeltzerBergerBremsModel::SamplePhotonEnergy(const MaterialCuts *matcut, double eekin, double r1, double r2, double r3) {
  constexpr double  mgdl  = 4.0*geant::kPi*geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght*geant::kRedElectronComptonWLenght;
  constexpr double  ddum0 = 2.0*geant::kElectronMassC2;
  constexpr double  ddum1 = geant::kElectronMassC2*geant::kElectronMassC2;

  double densityFactor = matcut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*mgdl;
  double ptot2         = eekin*(eekin+ddum0);
  double etot2         = ptot2+ddum1;
  double densityCor    = densityFactor*etot2;
  double gcut          = matcut->GetProductionCutsInEnergy()[0];

/*
// these are checked now in the caller
  if (fMaxElecEnergy<eekin){
    std::cerr<<" **** Electron energy = "<<eekin/geant::GeV<<" [GeV] > fMaxElecEnergy = "<<fMaxElecEnergy<<std::endl;
    exit(-1);
  }
  if (eekin<gcut) { // protection: no gamma production at that e- energy and x-section should be zero so it should never happen
    return 0.0;
  }
  // we could put a protection here to use the lowest electron/positron energy or highest if eekin is below/above
*/
  int mcindx       = matcut->GetIndex();
  int macindxlocal = fGlobalMatGCutIndxToLocal[mcindx];
  // the location of the first-electron-energy lin-alias data of this mat-cut
  int indxstart    = macindxlocal*fNumSamplingElecEnergies;
  // determine electron energy lower grid point
  double leenergy  = std::log(eekin);
  int eenergyindx  = (int) ((leenergy-fElEnLMin)*fElEnILDelta);

  if (eenergyindx>=fNumSamplingElecEnergies-1)
    eenergyindx = fNumSamplingElecEnergies-2;

  double ploweener = (fLSamplingElecEnergies[eenergyindx+1]-leenergy)*fElEnILDelta;

  indxstart +=eenergyindx;
  if (r1>ploweener|| !(fAliasData[indxstart]))
    ++indxstart;

  // sample the transformed variable
  double gammae = fAliasSampler->SampleLinear(fAliasData[indxstart]->fXdata, fAliasData[indxstart]->fYdata,
                                       fAliasData[indxstart]->fAliasW, fAliasData[indxstart]->fAliasIndx,
                                       fNumSamplingPhotEnergies,r2,r3);
  // transform it back to gamma-energy
  double dum0  = (gcut*gcut+densityCor);
  double dum1  = (eekin*eekin+densityCor)/dum0;
  double    u  = dum0*std::exp(gammae*std::log(dum1));

  return std::sqrt(u-densityCor);
}

// the simple DipBustgenerator
void SeltzerBergerBremsModel::SamplePhotonDirection(double elenergy, double &sinTheta, double &cosTheta, double rndm) {
  double c = 4. - 8.*rndm;
  double a = c;
  double signc = 1.;
  if (c<0.) {
    signc = -1.;
    a     = -c;
  }
  double delta  = std::sqrt(a*a+4.);
  delta += a;
  delta *= 0.5;

  double cofA = -signc*std::exp(std::log(delta)/3.0);
  cosTheta = cofA-1./cofA;

  double tau  = elenergy/geant::kElectronMassC2;
  double beta = std::sqrt(tau*(tau+2.))/(tau+1.);

  cosTheta = (cosTheta+beta)/(1.+cosTheta*beta);
  // check cosTheta limit
  if (cosTheta>1.0) {
    cosTheta = 1.0;
  }
  sinTheta = std::sqrt((1.-cosTheta)*(1.+cosTheta));
}


void SeltzerBergerBremsModel::InitSamplingTables() {
  // set up the common electron energy grid
  if (fSamplingElecEnergies) {
    delete [] fSamplingElecEnergies;
    delete [] fLSamplingElecEnergies;
    fSamplingElecEnergies  = nullptr;
    fLSamplingElecEnergies = nullptr;
  }
  fSamplingElecEnergies  = new double[fNumSamplingElecEnergies];
  fLSamplingElecEnergies = new double[fNumSamplingElecEnergies];

  fElEnLMin    = std::log(fMinElecEnergy);
  double delta = std::log(fMaxElecEnergy/fMinElecEnergy)/(fNumSamplingElecEnergies-1.0);
  fElEnILDelta = 1.0/delta;
  fSamplingElecEnergies[0]  = fMinElecEnergy;
  fLSamplingElecEnergies[0] = fElEnLMin;
  fSamplingElecEnergies[fNumSamplingElecEnergies-1]  = fMaxElecEnergy;
  fLSamplingElecEnergies[fNumSamplingElecEnergies-1] = std::log(fMaxElecEnergy);
  for (int i=1; i<fNumSamplingElecEnergies-1; ++i) {
    fLSamplingElecEnergies[i] = fElEnLMin+i*delta;
    fSamplingElecEnergies[i]  = std::exp(fElEnLMin+i*delta);
//    std::cerr<<" E("<<i<<") = "<<fSamplingElecEnergies[i]/geant::GeV<<std::endl;
  }  // fMinElecEnergy fMaxElecEnergy az fLoadDCSElectronEnergyGrid[0] es fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies-1]

  // - get number of different material-gammacut pairs
  // - allocate space and fill the global to local material-cut index map
  const std::vector<MaterialCuts*> theMaterialCutsTable = MaterialCuts::GetTheMaterialCutsTable();
  int numMaterialCuts = theMaterialCutsTable.size();
  if (fGlobalMatGCutIndxToLocal) {
    delete [] fGlobalMatGCutIndxToLocal;
    fGlobalMatGCutIndxToLocal = nullptr;
  }
  fGlobalMatGCutIndxToLocal = new int[numMaterialCuts];
  //std::cerr<<" === Number of global Material-Cuts = "<<numMaterialCuts<<std::endl;

  // count diffenet material-gammacut pairs and set to global to local mat-cut index map
  int oldnumDif = fNumDifferentMaterialGCuts;
  int oldnumSEE = fNumSamplingElecEnergies;
  fNumDifferentMaterialGCuts = 0;
  for (int i=0; i<numMaterialCuts; ++i) {
    // if the current MaterialCuts does not belong to the current active regions
    if (!IsActiveRegion(theMaterialCutsTable[i]->GetRegionIndex())) {
      fGlobalMatGCutIndxToLocal[i] = -1;
      continue;
    }
    bool isnew = true;
    int j = 0;
    for (; j<fNumDifferentMaterialGCuts; ++j) {
      if (theMaterialCutsTable[i]->GetMaterial()->GetIndex()==theMaterialCutsTable[j]->GetMaterial()->GetIndex() &&
          theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]==theMaterialCutsTable[j]->GetProductionCutsInEnergy()[0]) {
        isnew = false;
        break;
      }
    }
    if (isnew) {
      fGlobalMatGCutIndxToLocal[i] = fNumDifferentMaterialGCuts;
      ++fNumDifferentMaterialGCuts;
    } else {
      fGlobalMatGCutIndxToLocal[i] = fGlobalMatGCutIndxToLocal[j];
    }
  }
  //std::cerr<<" === Number of local Material-Cuts = "<<fNumDifferentMaterialGCuts<<std::endl;

  // allocate space for the matrial-gcut sampling tables and init these pointers to null
  if (fAliasData) {
    for (int i=0; i<oldnumDif; ++i) {
      for (int j=0; j<oldnumSEE; ++j) {
        int indx = i*oldnumSEE+j;
        if (fAliasData[indx]) {
          delete [] fAliasData[indx]->fXdata;
          delete [] fAliasData[indx]->fYdata;
          delete [] fAliasData[indx]->fAliasW;
          delete [] fAliasData[indx]->fAliasIndx;
          delete fAliasData[indx];
        }
      }
    }
    delete [] fAliasData;
  }


  int *isdone = new int[fNumDifferentMaterialGCuts]();
  int  idum   = fNumDifferentMaterialGCuts*fNumSamplingElecEnergies;
  fAliasData = new LinAlias*[idum];
  for (int i=0; i<idum; ++i)
    fAliasData[i] = nullptr;

  for (int i=0; i<numMaterialCuts; ++i) {
    //std::cerr<<"   See if Material-Gcut ==> " <<theMaterialCutsTable[i]->GetMaterial()->GetName()<<"  gCut = "<< theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]<<std::endl;
    int localindx = fGlobalMatGCutIndxToLocal[i];
    if (localindx<0) {
      continue;
    }
    int ialias    = localindx*fNumSamplingElecEnergies;
    if (!isdone[localindx]) { // init for this material-gammacut pair if it has not been done yet
//       std::cerr<<"   -> Will init for Material-Gcut ==> " <<theMaterialCutsTable[i]->GetMaterial()->GetName()<<"  gCut = "<< theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]<<std::endl;
       BuildOneLinAlias(ialias, theMaterialCutsTable[i]->GetMaterial(), theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]);
       isdone[localindx] = 1;
    }
  }
  delete [] isdone;
  // test
//  for (int i=0; i<numMaterialCuts; ++i)
//    std::cerr<<"     --> Global MatCut-indx = "<< i << " local indx = "<<fGlobalMatGCutIndxToLocal[i] <<std::endl;
}


/**
 *  The distribution of the energy of the emitted bremsstrahlung photons \f$ k \f$ is determined by the differential
 *  cross section
 *  \f[
 *    p(k) \propto \frac{\mathrm{d}\sigma}{\mathrm{d}k} \frac{1}{k\Gamma(k)}
 *  \f]
 *  where \f$ \Gamma(k) = (1+k_p^2/k^2) \f$ is the main factor related to dielectric suppression. The dielectric
 *  suppression removes the infrared divergence and strongly suppresses the low photon energy tail of the distribution
 *  which results in a non-ideal shape from sampling point of view. Furthermore, since sampling tables are built only
 *  at discrete \f$ E_{i}^{kin} \f$ incident particle kinetic energies and during the simulation we always have incident
 *  particle with \f$ E_{i}^{kin} < E^{kin} < E_{i+1}^{kin}\f$, we need to transform the prepared distributions to the
 *  common range i.e. [0,1] in order to guarante that the sampled photon energy is always within the proper kinematic
 *  limits i.e. \f$ k \in [k_c,E^{kin}] \f$ where \f$ k_c \f$ is the kinetic energy threshold for gamma particle
 *  production. So we apply the following variable transformations:
 *
 *  - since the Seltzer-Berger DCS are available in a "scalled" form i.e.
 *    \f[
 *         \frac{\mathrm{d}\sigma}{\mathrm{d}k} = \frac{Z^2}{\beta^2} \frac{1}{k}\chi(\kappa; E,Z)
 *    \f]
 *    we transform
 *    \f$ k \to \kappa(k) = k/E,\; \frac{\mathrm{d}k}{\mathrm{d}\kappa} = E \f$ and by using
 *    \f$p(k)\mathrm{d}k \propto p(\kappa)\mathrm{d}\kappa \f$ one can get
 *    \f$ p(\kappa) \propto p(k) \frac{\mathrm{d}k}{\mathrm{d}\kappa}
 *    \propto \chi(\kappa; E,Z) \frac{1}{\kappa} \frac{\kappa^2E^2}{\kappa^2E^2+k_p^2}\f$
 *  - then we can get rid of the dielectric suppression tail by transforming
 *    \f$ \kappa \to u(\kappa)=\ln(\kappa^2E^2+k_p^2),\; \kappa=\sqrt{\frac{e^u-k_p^2}{E^2}},\;
 *    \frac{\mathrm{d}\kappa}{\mathrm{d}u}=\frac{\kappa^2E^2+k_p^2}{2E^2\kappa} \f$ and one can use
 *    \f$ p(\kappa)\mathrm{d}\kappa \propto p(u)\mathrm{d}u \f$ to get
 *    \f$ p(u) \propto p(\kappa)\frac{\mathrm{d}\kappa}{\mathrm{d}u} \propto \chi(\kappa;E,Z) \f$
 *  - then we appaly the \f$ u \to \xi(u) = [u-\ln(k_c^2+k_p^2)]/\ln[(E^2+k_p^2)/(k_c^2+k_p^2)] \in [0,1]\f$
 *    transformation that by taking into account that \f$ u = \xi\ln[(E^2+k_p^2)/(k_c^2+k_p^2)]+\ln[k_c^2+k_p^2],\;
 *    \frac{\mathrm{d}\xi}{\mathrm{d}u} = \ln[(E^2+k_p^2)/(k_c^2+k_p^2)]\f$ and by using \f$ p(u)\mathrm{d}u \propto
 *    p(\xi)\mathrm{d}\xi\f$ one can get \f$ p(\xi) \propto p(u)\frac{\mathrm{d}u}{\mathrm{d}\xi} \propto
 *    \ln[(E^2+k_p^2)/(k_c^2+k_p^2)] \chi(\kappa;E,Z) \propto \chi(\kappa;E,Z)\f$ since
 *    \f$ \ln[(E^2+k_p^2)/(k_c^2+k_p^2)]\f$ is just a constant.
 *
 *  When this transformed p.d.f. are prepared the numerical "scalled" Seltzer-Berger DCS are interpolated similarly
 *  like in the case of SeltzerBergerBremsModel::ComputeXSectionPerAtom. The variable \f$ \phi=\ln[1-\kappa+10^{-12}]\f$
 *  can be obtained at a given \f$ \xi \f$ value by \f$ \kappa = \sqrt{\frac{e^u-k_p^2}{E^2}}
 *  =\frac{1}{E}\sqrt{[k_c^2+k_p^2]e^{\xi\ln[(E^2+k_p^2)/(k_c^2+k_p^2)]}-k_p^2} \f$ and during the sampling, the emitted
 *  gamma photon energy \f$k = \kappa E = \sqrt{[k_c^2+k_p^2]e^{\xi\ln[(E^2+k_p^2)/(k_c^2+k_p^2)]}-k_p^2} \f$  for a
 *  given sampled transformed variable value \f$\xi\f$ (where \f$E\f$ is the actual run-time primary kinetic energ)y so
 *  the sampled photon energy will always be within the actual kinematic limits i.e. \f$ k_c<k \leq E\f$) .
 */
void SeltzerBergerBremsModel::BuildOneLinAlias(int ialias, const Material *mat, double gcut){
  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
                        *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;
  // we will need the element composition of this material
  const Vector_t<Element*> theElements    = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
//  int   numElems = theElements.size();

  //
  // for the common particle kinetic energy grid points
  //
  for (int ieener=0; ieener<fNumSamplingElecEnergies; ++ieener) {
    double eener = fSamplingElecEnergies[ieener];
    if (eener>gcut) { // otherwise no gamma production at that e- energy so let the sampling table to be nullptr
//      std::cerr<<"         indx ="<<ialias<<std::endl;
      // find the e- energy index in the available e- energies grid such that we are just above
      int eenerindx = 0;
      for (; eenerindx<fLoadDCSNumElectronEnergies; ++eenerindx)
        if (eener<fLoadDCSElectronEnergyGrid[eenerindx])
          break;
      if (eenerindx==fLoadDCSNumElectronEnergies)
        --eenerindx;
      //eenerindx is the upper index now

      double ptot2   = eener*(eener+2.0*geant::kElectronMassC2);
      double etot2   = ptot2+geant::kElectronMassC2*geant::kElectronMassC2;
//      double ibeta2  = etot2/ptot2;

      double densityCor = densityFactor*etot2;

      double *theDCS = new double[fLoadDCSNumReducedPhotonEnergies];
      double *logReducedPhotonEnergyGrid = new double[fLoadDCSNumReducedPhotonEnergies];
       // ln(x)-ln(x1)
       double dum0 = std::log(eener/fLoadDCSElectronEnergyGrid[eenerindx-1]);
       // ln(x2)-ln(x1)
       double dum1 = std::log(fLoadDCSElectronEnergyGrid[eenerindx]/fLoadDCSElectronEnergyGrid[eenerindx-1]);

       for (int irpener=0; irpener<fLoadDCSNumReducedPhotonEnergies; ++irpener) {
        logReducedPhotonEnergyGrid[irpener] = std::log(1.0-fLoadDCSReducedPhotonEnergyGrid[irpener]+1.0e-12);
         int indxdcsh = eenerindx*fLoadDCSNumReducedPhotonEnergies + irpener;
         int indxdcsl = (eenerindx-1)*fLoadDCSNumReducedPhotonEnergies + irpener;
         theDCS[irpener] = 0.0;
         for (unsigned long ielem=0; ielem<theElements.size(); ++ielem) {
           double zet  = theElements[ielem]->GetZ();
           int    izet = lrint(zet);
           // ln(y2) -ln(y1)
           //           double dum2 = std::log(fLoadDCSForElements[izet-1][indxdcsh]/fLoadDCSForElements[izet-1][indxdcsl]);
           //           double dcs  = dum2/dum1*dum0+std::log(fLoadDCSForElements[izet-1][indxdcsl]); //this is ln(dcs)
           //           dcs = std::exp(dcs);
           double dum2 = fLoadDCSForElements[izet-1][indxdcsh]-fLoadDCSForElements[izet-1][indxdcsl];
           double dcs  = dum2/dum1*dum0+fLoadDCSForElements[izet-1][indxdcsl];
           dcs *= theAtomicNumDensityVector[ielem]*zet*zet;
           // correction for positrons
           if (!fIsElectron) {
             dcs *= PositronCorrection1(eener, fLoadDCSReducedPhotonEnergyGrid[irpener], gcut, zet);
//             dcs *= PositronCorrection(eener, ibeta2, fLoadDCSReducedPhotonEnergyGrid[irpener], zet);
           }
           theDCS[irpener] += dcs;//dcs/(fLoadDCSReducedPhotonEnergyGrid[irpener]+1e-12);
         }
//         if (theDCS[irpener]<=0.0) // e+ dcs is zero at kappa=1
//           theDCS[irpener] = 1.0e-13;
//         theDCS[irpener] = std::log(theDCS[irpener]);
       }

       // set up a spline on the log(1-kappa)::log(dcs)
       Spline     *sp = new Spline();
       sp->SetUpSpline(logReducedPhotonEnergyGrid, theDCS, fLoadDCSNumReducedPhotonEnergies);

       //
       // fill in the initial x,y data
       //
       // create the alias data struct
       fAliasData[ialias] = new LinAlias();
       fAliasData[ialias]->fXdata     = new double[fNumSamplingPhotEnergies]();
       fAliasData[ialias]->fYdata     = new double[fNumSamplingPhotEnergies]();
       fAliasData[ialias]->fAliasW    = new double[fNumSamplingPhotEnergies]();
       fAliasData[ialias]->fAliasIndx = new int[fNumSamplingPhotEnergies]();
       // find reduced foton energy index such that we are just above
       double kappac = gcut/eener;
       int kappaindx = 0;
       for (; kappaindx<fLoadDCSNumReducedPhotonEnergies; ++kappaindx)
         if (kappac<fLoadDCSReducedPhotonEnergyGrid[kappaindx])
           break;
       if (std::abs(1.0-kappac/fLoadDCSReducedPhotonEnergyGrid[kappaindx])<1.e-12) // handle some possible rounding problems: if kappac is very close
         ++kappaindx;

       if (kappaindx>=fLoadDCSNumReducedPhotonEnergies)
         kappaindx = fLoadDCSNumReducedPhotonEnergies-1;

       // fill the first initial value; convert kappa scale to 0-1
       int numdata = 1;
       fAliasData[ialias]->fXdata[0] = 0.0;
       double kappaconv = std::log(1.0-kappac+1.0e-12);
//       fAliasData[ialias]->fYdata[0] = std::exp(sp->GetValueAt(kappaconv));
       fAliasData[ialias]->fYdata[0] = sp->GetValueAt(kappaconv);


       // fill possible values if any
       for (int k=kappaindx; k<fLoadDCSNumReducedPhotonEnergies-1; ++k) {
         double thekappa = fLoadDCSReducedPhotonEnergyGrid[k];
         kappaconv = std::log( (thekappa*thekappa*eener*eener+densityCor)/(gcut*gcut+densityCor) )/
                     std::log( (eener*eener+densityCor)/(gcut*gcut+densityCor)); // coverted to 0-1
         fAliasData[ialias]->fXdata[numdata] = kappaconv;  // xi
//         fAliasData[ialias]->fYdata[numdata] = std::exp(theDCS[k]);
         fAliasData[ialias]->fYdata[numdata] = theDCS[k];

         ++numdata;
       }
       // and the last point
       fAliasData[ialias]->fXdata[numdata] = 1.0;
//       fAliasData[ialias]->fYdata[numdata] = std::exp(theDCS[fLoadDCSNumReducedPhotonEnergies-1]);
       fAliasData[ialias]->fYdata[numdata] = theDCS[fLoadDCSNumReducedPhotonEnergies-1];
       ++numdata;

       // expand the data up to maximum
       while(numdata<fNumSamplingPhotEnergies) {
         // find the lower index of the bin, where we have the biggest linear interp. error compared to spline
         double maxerr     = 0.0; // value of the current maximum error
         double thexval    = 0.0;
         double theyval    = 0.0;
         int    maxerrindx = 0;   // the lower index of the corresponding bin

         for (int i=0; i<numdata-1; ++i) {
           double xx = 0.5*(fAliasData[ialias]->fXdata[i]+fAliasData[ialias]->fXdata[i+1]); // mid point
           double yy = 0.5*(fAliasData[ialias]->fYdata[i]+fAliasData[ialias]->fYdata[i+1]); // lin func val at the mid point

           dum0  = (gcut*gcut+densityCor);
           dum1  = (eener*eener+densityCor)/dum0;
           double    u  = dum0*std::exp(xx*std::log(dum1));
           double thekappa = std::sqrt( u-densityCor)/eener; // kappa

           double conv  = std::log(1.0-thekappa+1.0e-12);
//           double spval = std::exp(sp->GetValueAt(conv)); // spline intp. val. at mid point
           double spval = sp->GetValueAt(conv); // spline intp. val. at mid point
           double err   = std::fabs(yy-spval);// works better than fabs(1-yy/spval) might try constrained spline?
           if (err>maxerr) {
             maxerr     = err;
             maxerrindx = i;
             thexval    = xx;
             theyval    = spval;
           }
         }
         // extend x,y data by puting a spline interp.ted value at the mid point of the highest error bin
         // first shift all values to the right
         for (int j=numdata; j>maxerrindx+1; --j) {
           fAliasData[ialias]->fXdata[j] = fAliasData[ialias]->fXdata[j-1];
           fAliasData[ialias]->fYdata[j] = fAliasData[ialias]->fYdata[j-1];
         }
         // fill x mid point
         fAliasData[ialias]->fXdata[maxerrindx+1] = thexval;
         fAliasData[ialias]->fYdata[maxerrindx+1] = theyval;
         // increase number of data
         ++numdata;
       }

       // set up a linear alias smapler on this data
       AliasTable *alst = new AliasTable();
       alst->PreparLinearTable(fAliasData[ialias]->fXdata, fAliasData[ialias]->fYdata,
                               fAliasData[ialias]->fAliasW, fAliasData[ialias]->fAliasIndx,
                               fNumSamplingPhotEnergies);

      delete alst;
      delete sp;
      delete [] theDCS;
      delete [] logReducedPhotonEnergyGrid;
    }
    ++ialias;
  }
}

/**
 * The restricted atomic cross section for bremsstrahlung photon emission for target element with atomic number
 * \f$Z\f$, gamma photon production energy threshold \f$k_c\f$ is computed for e-/e+ kinetic energy \f$E\f$ according to
 * \f[
 *   \sigma(E;k_c,Z) = \int_{k_c}^{E} \frac{1}{\Gamma}\frac{\mathrm{d}\sigma}{\mathrm{d}k}\mathrm{d}k
 * \f]
 * if \f$E>k_c\f$ and immediate return with \f$0\f$ otherwise. (The \f$1/\Gamma\f$ factor is the main dielectric
 * suppression factor and \f$\Gamma = (1+k_p^2/k^2)\f$ where
 * \f$ k_p = \hbar \omega_p \gamma= \hbar \omega_p E_t/(mc^2)\f$, \f$\hbar \omega_p= \hbar c \sqrt{4\pi n_e r_e}\f$
 * (\f$n_e\f$ is the electron density) see more details at RelativisticBremsModel::ComputeURelDXSecPerAtom ).
 * The Seltzer-Berger numerical DCS are available in the from of "scalled" DCS as
 * \f[
 *      \chi(\kappa;E,Z) =  \frac{\beta^2}{Z^2}k \frac{\mathrm{d}\sigma}{\mathrm{d}k}
 * \f]
 * where \f$k\f$ is the emitted photon energy, \f$\kappa=k/E\f$ is the reduced photon energy. The above integral can be
 * written now with the "scalled" DCS as
 * \f[
 *   \sigma(E;k_c,Z) = \int_{k_c}^{E} \frac{1}{\Gamma}\frac{\mathrm{d}\sigma}{\mathrm{d}k}\mathrm{d}k =
 *                     \frac{Z^2}{\beta^2}\int_{\kappa_c}^{1} \frac{1}{\kappa\Gamma}\chi(\kappa;E,Z)\mathrm{d}\kappa
 * \f]
 * where \f$\kappa_c=k_c/E\f$.
 * Since the "scalled" DCS are available at fixed electron kinetic energies \f$\{E_i\}_i^N\f$, linear interpolation in
 * log electron energies is applied to obtain \f$\chi(\kappa;E,Z)\f$ from \f$\chi(\kappa;E_i,Z)\f$ and
 * \f$\chi(\kappa;E_{i+1},Z)\f$ such that \f$E_i \leq E < E_{i+1}\f$ over the available fixed set of reduced photon
 * energies \f$\{\kappa_j\}_j^M\f$ where \f$ \kappa_j \in [0,1]\; \forall\; j\f$. During the interpolation, the reduced
 * photon energy grid i.e. \f$\{\kappa_j\}_j^M\f$ is transformed to \f$ \phi = \ln[1-\kappa+10^{-12}] \f$ for getting a
 * more accurate interpolation later when the integral is computed (with 64-points Gauss-Legendre quadrature using cubic
 * spline interpolation of DCS values over the \f$\phi\f$ grid).
 *
 * The integral is computed by 64-points Gauss-Legendre quadrature after applying the following transformations
 * - first the reduced photon energy is transformed to \f$\kappa \to u=\ln(\kappa) \f$
 * - then we apply the following transformation \f$ u \to \xi = [u-\ln(\kappa_c)]/\ln(1/\kappa_c) \in [0,1] \f$
 *
 * The transformed integral
 * \f[
 *   \sigma(E;k_c,Z) = \frac{Z^2}{\beta^2}\int_{\kappa_c}^{1} \frac{1}{\kappa\Gamma}\chi(\kappa;E,Z)\mathrm{d}\kappa
 *                   = \frac{Z^2}{\beta^2}\int_{\ln(\kappa_c)}^{0} \frac{1}{\Gamma}\chi(\kappa;E,Z)\mathrm{d}u
 *                   = \frac{Z^2}{\beta^2}\ln\frac{1}{\kappa_c}\int_{0}^{1} \frac{1}{\Gamma}\chi(\phi;E,Z)\mathrm{d}\xi
 * \f]
 * where \f$\Gamma\f$ must be evaluated at \f$k=E\kappa_c e^{\xi\ln(1/\kappa_c)}\f$ and \f$ \chi(\phi;E,Z) \f$ must be
 * evaluated (interpolated) at \f$\phi =\ln[1-\kappa_c e^{\xi\ln(1/\kappa_c)}+10^{-12}]\f$ at a given value of \f$\xi\f$.
 */
double SeltzerBergerBremsModel::ComputeXSectionPerAtom(const Element *elem, const Material *mat, double gammaprodcutenergy,
                                                       double electronekin){
  double xsec = 0.0;
  if (electronekin<=gammaprodcutenergy || electronekin<fLoadDCSElectronEnergyGrid[0]
      || electronekin>fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies-1])
    return xsec;

  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
                        *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;

  double zet     = elem->GetZ();
  int    izet    = std::lrint(zet);
  double ptot2   = electronekin*(electronekin+2.0*geant::kElectronMassC2);
  double etot2   = ptot2+geant::kElectronMassC2*geant::kElectronMassC2;
  double ibeta2  = etot2/ptot2;
  double densityCor = densityFactor*etot2;//electronekin*electronekin;  // this is k_p^2
  double kappacr    = gammaprodcutenergy/electronekin;
  double logikappacr = std::log(1./kappacr);

  // find the electron energy index that
  int ieener = 0;
  for (ieener=0; ieener<fLoadDCSNumElectronEnergies; ++ieener)
    if (electronekin<fLoadDCSElectronEnergyGrid[ieener])
      break;
  if (ieener==fLoadDCSNumElectronEnergies)
    --ieener;

  double *theDCS = new double[fLoadDCSNumReducedPhotonEnergies];
  double *logReducedPhotonEnergyGrid = new double[fLoadDCSNumReducedPhotonEnergies];
  // ln(x)-ln(x1)
  double dum0 = std::log(electronekin/fLoadDCSElectronEnergyGrid[ieener-1]);
  // ln(x2)-ln(x1)
  double dum1 = std::log(fLoadDCSElectronEnergyGrid[ieener]/fLoadDCSElectronEnergyGrid[ieener-1]);
//  std::cerr<<" -->"<<electronekin<<" "<<ieener<<std::endl;

  for (int irpener=0; irpener<fLoadDCSNumReducedPhotonEnergies; ++irpener) {
    logReducedPhotonEnergyGrid[irpener] = std::log(1.0-fLoadDCSReducedPhotonEnergyGrid[irpener]+1.0e-12);
    int indxdcsh = ieener*fLoadDCSNumReducedPhotonEnergies + irpener;
    int indxdcsl = (ieener-1)*fLoadDCSNumReducedPhotonEnergies + irpener;
    double dcsl = fLoadDCSForElements[izet-1][indxdcsl];
    double dcsh = fLoadDCSForElements[izet-1][indxdcsh];
    // ln(y2) -ln(y1)
    //    double dum2 = std::log(dcsh/dcsl);
    //    double dcs  = dum2/dum1*dum0+std::log(dcsl);
    //    dcs  = std::exp(dcs);
    double dcs  = (dcsh-dcsl)/dum1*dum0+dcsl;
    // correction for positrons
//    if (!fIsElectron)
//      dcs *= PositronCorrection(electronekin, ibeta2, fLoadDCSReducedPhotonEnergyGrid[irpener], zet);
    theDCS[irpener] = dcs;
  }

  // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
  const std::vector<double> &glW = fGL->GetWeights();
  const std::vector<double> &glX = fGL->GetAbscissas();
  // we will need a natural cubic spile for the integral
  Spline     *sp = new Spline();
  sp->SetUpSpline(logReducedPhotonEnergyGrid, theDCS, fLoadDCSNumReducedPhotonEnergies);

  double integral = 0.0;
  for (int i=0; i<fNGL; ++i) {
    double dumx = 1.0-std::exp(glX[i]*logikappacr)*kappacr; // ez a felso x -nek az exp(x)-e
    double x = std::log(dumx+1.0e-12);
    double egamma = (1.0-dumx)*electronekin;
    double poscor = 1.0;
    if (!fIsElectron)
      poscor *= PositronCorrection(electronekin, ibeta2, 1.0-dumx, zet);

    integral+= glW[i]*poscor*sp->GetValueAt(x)/(1.+densityCor/(egamma*egamma));
  }

  delete [] theDCS;
  delete [] logReducedPhotonEnergyGrid;
  delete sp;

  return logikappacr*zet*zet*ibeta2*integral;
}

/**
 *   The restricted macroscopic cross section for bremsstrahlung photon emission for the given target material,
 *   gamma photon production energy threshold \f$k_c\f$ is computed for e-/e+ kinetic energy \f$E\f$ according to
 *  \f[
 *      \Sigma(E;k_c,\mathrm{material}) = \sum_i n_i \sigma_i(E;k_c,Z_i)
 *  \f]
 *  if \f$ E>k_c\f$ otherwise immediate return with \f$0\f$. The summation goes over the elements the matrial is
 *  composed from. \f$\sigma_i(E;k_c,Z_i)\f$ is the restricted atomic cross
 *  secion for the \f$i\f$-th element of the material with atomic number of \f$Z_i \f$ (computed similarly like
 *  SeltzerBergerBremsModel::ComputeXSectionPerAtom()) and \f$n_i\f$ is the number of atoms per unit volume of
 *  \f$i \f$-th element of the material that is \f$ n_i = \mathcal{N}\rho w_i/A_i \f$ where \f$\mathcal{N}\f$ is the
 *  Avogadro number, \f$\rho\f$ is the material density, \f$w_i\f$ is the proportion by mass (or mass fraction) of the
 *  \f$i \f$-th element and \f$A_i\f$ is the molar mass of the \f$i \f$-th element. The corresponding mean free path
 *  is \f$\lambda = 1/\Sigma \f$.
 */
double SeltzerBergerBremsModel::ComputeXSectionPerVolume(const Material *mat, double gammaprodcutenergy, double electronekin){
  double xsec = 0.0;
  if (electronekin<=gammaprodcutenergy || electronekin<fLoadDCSElectronEnergyGrid[0]
      || electronekin>fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies-1])
    return xsec;

  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
                        *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;

//  double zet     = elem->GetZ();
//  int    izet    = std::lrint(zet);
  double ptot2   = electronekin*(electronekin+2.0*geant::kElectronMassC2);
  double etot2   = ptot2+geant::kElectronMassC2*geant::kElectronMassC2;
  double ibeta2  = etot2/ptot2;
  double densityCor = densityFactor*etot2;//electronekin*electronekin;  // this is k_p^2
  double kappacr    = gammaprodcutenergy/electronekin;
  double logikappacr = std::log(1./kappacr);

  // find the electron energy index that
  int ieener = 0;
  for (ieener=0; ieener<fLoadDCSNumElectronEnergies; ++ieener)
    if (electronekin<fLoadDCSElectronEnergyGrid[ieener])
      break;
  if (ieener==fLoadDCSNumElectronEnergies)
    --ieener;


  // we will need the element composition of this material
  const Vector_t<Element*> theElements    = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int   numElems = theElements.size();

  double *theDCS = new double[fLoadDCSNumReducedPhotonEnergies*numElems];
  double *logReducedPhotonEnergyGrid = new double[fLoadDCSNumReducedPhotonEnergies];
  // ln(x)-ln(x1)
  double dum0 = std::log(electronekin/fLoadDCSElectronEnergyGrid[ieener-1]);
  // ln(x2)-ln(x1)
  double dum1 = std::log(fLoadDCSElectronEnergyGrid[ieener]/fLoadDCSElectronEnergyGrid[ieener-1]);

  for (int irpener=0; irpener<fLoadDCSNumReducedPhotonEnergies; ++irpener) {
    logReducedPhotonEnergyGrid[irpener] = std::log(1.0-fLoadDCSReducedPhotonEnergyGrid[irpener]+1.0e-12);
    int indxdcsh = ieener*fLoadDCSNumReducedPhotonEnergies + irpener;
    int indxdcsl = (ieener-1)*fLoadDCSNumReducedPhotonEnergies + irpener;
    for (unsigned long ielem=0; ielem<theElements.size(); ++ielem) {
      double zet  = theElements[ielem]->GetZ();
      int    izet = lrint(zet);
      // ln(y2) -ln(y1)
      //      double dum2 = std::log(fLoadDCSForElements[izet-1][indxdcsh]/fLoadDCSForElements[izet-1][indxdcsl]);
      //      double dcs  = dum2/dum1*dum0+std::log(fLoadDCSForElements[izet-1][indxdcsl]); //this is ln(dcs)
      //      dcs  = std::exp(dcs);
      double dum2 = fLoadDCSForElements[izet-1][indxdcsh]-fLoadDCSForElements[izet-1][indxdcsl];
      double dcs  = dum2/dum1*dum0+fLoadDCSForElements[izet-1][indxdcsl];
      dcs *= theAtomicNumDensityVector[ielem]*zet*zet;

      // correction for positrons
//      if (!fIsElectron) {
//        dcs *= PositronCorrection(electronekin, ibeta2, fLoadDCSReducedPhotonEnergyGrid[irpener], zet);
//      }
      theDCS[ielem*fLoadDCSNumReducedPhotonEnergies+irpener] = dcs;
    }
//    theDCS[irpener] = std::log(theDCS[irpener]);
  }

  // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
  const std::vector<double> &glW = fGL->GetWeights();
  const std::vector<double> &glX = fGL->GetAbscissas();
  // we will need a natural cubic spile for the integral
//  Spline     *sp = new Spline();
//  sp->SetUpSpline(logReducedPhotonEnergyGrid, theDCS, fLoadDCSNumReducedPhotonEnergies);
  // we will need as many natural cubic spile as elements
  Spline  **sp = new Spline*[numElems];
  for (int i=0; i<numElems; ++i) {
    sp[i] = new Spline();
    sp[i]->SetUpSpline(logReducedPhotonEnergyGrid, &theDCS[i*fLoadDCSNumReducedPhotonEnergies], fLoadDCSNumReducedPhotonEnergies);
  }

  double integral = 0.0;
  for (int i=0; i<fNGL; ++i) {
    double dumx = 1.0-std::exp(glX[i]*logikappacr)*kappacr; // ez a felso x -nek az exp(x)-e
    double x = std::log(dumx+1.0e-12);
    double egamma = (1.0-dumx)*electronekin;
    double sum = 0.0;
    for (int ielem=0; ielem<numElems; ++ielem) {
//      double val = std::exp(sp[ielem]->GetValueAt(x));
      double val = sp[ielem]->GetValueAt(x);
      if (!fIsElectron) {
        double zet  = theElements[ielem]->GetZ();
        val *= PositronCorrection(electronekin, ibeta2, (1.0-dumx), zet);
      }
      sum += val;
    }
    integral+= glW[i]*sum/(1.+densityCor/(egamma*egamma));
  }

  delete [] theDCS;
  delete [] logReducedPhotonEnergyGrid;
//  delete sp;
  for (int i=0; i<numElems; ++i)
    delete sp[i];
  delete [] sp;
  return logikappacr*ibeta2*integral;
}

/**
 *  The stopping power, i.e. average energy loss per unit path length, from bremsstrahlung photon emission is computed
 *  for the given e-/e+ kinetic energy \f$E\f$, the given material and gamma production threshold energy \f$k_c\f$
 *  \f[
 *      S(E;k_c,\mathrm{material})=\int_{0}^{\eta} k \sum_i n_i \frac{\mathrm{d}\sigma_i}{\mathrm{d}k} \mathrm{d}k
 *  \f]
 *  where
 *  \f[
 *     \eta =
 *      \begin{cases}
 *            E    & \quad \mathrm{if}\; E<k_c \\
 *            k_c  & \quad \mathrm{if}\; E \geq k_c
 *      \end{cases}
 *  \f]
 *  the summation goes over the elements the matrial is composed from. \f$ \mathrm{d}\sigma_i \mathrm{d}k\f$ is the
 *  differential cross section for bremsstrahlung photon emission for for the \f$i\f$-th element of the material with
 *  atomic number of \f$Z_i\f$ and \f$n_i\f$ is the number of atoms per unit volume of \f$i\f$-th element of the
 *  material that is \f$n_i=\mathcal{N}\rho w_i/A_i\f$ where \f$\mathcal{N}\f$ is the Avogadro number, \f$\rho\f$ is
 *  the material density, \f$w_i\f$ is the proportion by mass (or mass fraction) of the \f$i\f$-th element and \f$A_i\f$
 *  is the molar mass of the \f$i\f$-th element.
 *
 *  The Seltzer-Berger atomic DCS are used and interpolated similarly like in case of atomic cross section computation
 *  (SeltzerBergerBremsModel::ComputeXSectionPerAtom).  The Seltzer-Berger numerical atomic DCS are available in
 *  the from of "scalled" DCS as
 *  \f[
 *      \chi(\kappa;E,Z) =  \frac{\beta^2}{Z^2}k \frac{\mathrm{d}\sigma}{\mathrm{d}k}
 *  \f]
 *  where \f$k\f$ is the emitted photon energy, \f$\kappa=k/E \in [0,1]\f$ is the reduced photon energy. The above
 *  integral can be written now with the "scalled" DCS as
 *  \f[
 *      S(E;k_c,\mathrm{material})= \int_{0}^{\eta} k \sum_i n_i \frac{\mathrm{d}\sigma_i}{\mathrm{d}k} \mathrm{d}k
 *    = \frac{E}{\beta^2}\int_{0}^{\eta/E}\frac{1}{\Gamma}\sum_i n_i Z_i^2 \chi(\kappa;E,Z_i) \mathrm{d}\kappa
 *  \f]
 *  (The \f$1/\Gamma\f$ factor is the main dielectric
 *  suppression factor and \f$\Gamma = (1+k_p^2/k^2)\f$ where
 *  \f$ k_p = \hbar \omega_p \gamma= \hbar \omega_p E_t/(mc^2)\f$, \f$\hbar \omega_p= \hbar c \sqrt{4\pi n_e r_e}\f$
 *  (\f$n_e\f$ is the electron density) see more details at RelativisticBremsModel::ComputeURelDXSecPerAtom ).
 *
 *  The integral is computed by 64-points Gauss-Legendre quadrature after the following transformation
 *  - the reduced photon energy is transformed \f$ \kappa \to \xi = \kappa/(\eta/E) \in [0,1] \f$.
 *
 *  The integral then becomes
 *  \f[
 *   S(E;k_c,\mathrm{material})
 *    = \frac{E}{\beta^2}\int_{0}^{\eta/E}\frac{1}{\Gamma}\sum_i n_i Z_i^2 \chi(\kappa;E,Z_i) \mathrm{d}\kappa
 *    = \frac{\eta}{\beta^2}\int_{0}^{1}\frac{1}{\Gamma}\sum_i n_i Z_i^2 \chi(\kappa;E,Z_i) \mathrm{d}\xi
 *  \f]
 *  where \f$\Gamma\f$ must be evaluated at \f$k=\xi\eta\f$ and \f$\chi\f$ at \f$\kappa=\xi\eta/E\f$ at a given value of
 *  \f$\xi\f$.
 */
double SeltzerBergerBremsModel::ComputeDEDXPerVolume(const Material *mat, double gammaprodcutenergy, double electronekin){
  double dedx = 0.0;

  if (electronekin<fLoadDCSElectronEnergyGrid[0] || electronekin>fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies-1])
    return dedx;

  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
                        *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;

  double ptot2   = electronekin*(electronekin+2.0*geant::kElectronMassC2);
  double etot2   = ptot2+geant::kElectronMassC2*geant::kElectronMassC2;
  double ibeta2  = etot2/ptot2;
  double densityCor = densityFactor*etot2;//*electronekin*electronekin;  // this is k_p^2

  // find the index of the electron energy that is just below the requested one
  int ieener = 0;
  for (ieener=0; ieener<fLoadDCSNumElectronEnergies; ++ieener)
    if (electronekin<fLoadDCSElectronEnergyGrid[ieener])
      break;
  // handle the case when electronekin=fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies-1]
  if (ieener==fLoadDCSNumElectronEnergies) {
    --ieener;
  }

  // we will need the element composition of this material
  const Vector_t<Element*> theElements    = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int   numElems = theElements.size();

  double *theDCS = new double[fLoadDCSNumReducedPhotonEnergies*numElems];
  // ln(x)-ln(x1)
  double dum0 = std::log(electronekin/fLoadDCSElectronEnergyGrid[ieener-1]);
  // ln(x2)-ln(x1)
  double dum1 = std::log(fLoadDCSElectronEnergyGrid[ieener]/fLoadDCSElectronEnergyGrid[ieener-1]);
  for (int irpener=0; irpener<fLoadDCSNumReducedPhotonEnergies; ++irpener) {
    int indxdcsh = ieener*fLoadDCSNumReducedPhotonEnergies + irpener;
    int indxdcsl = (ieener-1)*fLoadDCSNumReducedPhotonEnergies + irpener;
    for (unsigned long ielem=0; ielem<theElements.size(); ++ielem) {
      double zet  = theElements[ielem]->GetZ();
      int    izet = lrint(zet);
      // ln(y2) -ln(y1)
      //      double dum2 = std::log(fLoadDCSForElements[izet-1][indxdcsh]/fLoadDCSForElements[izet-1][indxdcsl]);
      //      double dcs  = dum2/dum1*dum0+std::log(fLoadDCSForElements[izet-1][indxdcsl]); //this is ln(dcs)
      //      dcs  = std::exp(dcs);
      double dum2 = fLoadDCSForElements[izet-1][indxdcsh]-fLoadDCSForElements[izet-1][indxdcsl];
      double dcs  = dum2/dum1*dum0+fLoadDCSForElements[izet-1][indxdcsl];

      dcs *= theAtomicNumDensityVector[ielem]*zet*zet;
//      dcs += std::log(theAtomicNumDensityVector[ielem]*zet*zet);
      theDCS[ielem*fLoadDCSNumReducedPhotonEnergies+irpener] = dcs;
    }
  }

  // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
  const std::vector<double> &glW = fGL->GetWeights();
  const std::vector<double> &glX = fGL->GetAbscissas();
  // we will need a natural cubic spile for the integral
//  Spline     *sp = new Spline();
//  sp->SetUpSpline(fLoadDCSReducedPhotonEnergyGrid, theDCS, fLoadDCSNumReducedPhotonEnergies);
  // we will need as many natural cubic spile as elements
  Spline  **sp = new Spline*[numElems];
  for (int i=0; i<numElems; ++i) {
    sp[i] = new Spline();
    sp[i]->SetUpSpline(fLoadDCSReducedPhotonEnergyGrid, &theDCS[i*fLoadDCSNumReducedPhotonEnergies], fLoadDCSNumReducedPhotonEnergies);
  }

  // integrate
  double integral    = 0.0;
  double upperlimit  = gammaprodcutenergy;
  if (upperlimit>electronekin)
    upperlimit = electronekin;
  double kappacr = upperlimit/electronekin;
  for (int i=0; i<fNGL; ++i) {
    double t = glX[i]*kappacr;          // kappa
    double egamma = glX[i]*upperlimit;
    double sum = 0.0;
    for (int ielem=0; ielem<numElems; ++ielem) {
//      double val = std::exp(sp[ielem]->GetValueAt(t));
      double val = sp[ielem]->GetValueAt(t);
      if (!fIsElectron) {
        double zet  = theElements[ielem]->GetZ();
        val *= PositronCorrection(electronekin, ibeta2, egamma/electronekin, zet);
      }
      sum += val;
    }

    integral += glW[i]*sum/(1.+densityCor/(egamma*egamma));
    // x 1/(1+k_p^2/k^2) i.e. density effect correction
  }

  delete [] theDCS;
  //  delete sp;
  for (int i=0; i<numElems; ++i)
    delete sp[i];
  delete [] sp;

  return upperlimit*ibeta2*integral;
}


// correction for positrons : DCS must be multiplied by this for positrons
// ephoton is the reduced photon energy
double SeltzerBergerBremsModel::PositronCorrection(double ekinelectron, double ibeta2electron,
                                                   double ephoton, double z) {
    using geant::kElectronMassC2;
    constexpr double dum1 = geant::kTwoPi*geant::kFineStructConst;

    double poscor = 0.0;
    double ibeta1   = std::sqrt(ibeta2electron);
    double e2       = ekinelectron * (1.0 - ephoton);
    if (e2 > 0.0) {
      double ibeta2 = (e2 + kElectronMassC2)/std::sqrt(e2*(e2+2.0*kElectronMassC2));
      double dum0   = dum1*z*(ibeta1-ibeta2);
      if (dum0<-12.0) {
        poscor = 0.0;
      } else {
        poscor = std::exp(dum0);
      }
    } else {
      poscor = 0.0;
    }
 return poscor;
}



//ephoton is the reduced photon energy
double SeltzerBergerBremsModel::PositronCorrection1(double ekinelectron, double ephoton, double gcutener, double z) {
    using geant::kElectronMassC2;
    constexpr double dum1 = geant::kTwoPi*geant::kFineStructConst;

    double poscor = 0.0;
    double e1     = ekinelectron-gcutener;  // here is the dif.
    double ibeta1 = (e1+kElectronMassC2)/std::sqrt(e1*(e1+2.0*kElectronMassC2));
    double e2     = ekinelectron * (1.0 - ephoton);
    double ibeta2 = (e2+kElectronMassC2)/std::sqrt(e2*(e2+2.0*kElectronMassC2));
    double ddum   = dum1*z*(ibeta1-ibeta2);
    if (ddum<-12.0) {
      poscor = 0.0;
    } else {
      poscor = std::exp(ddum);
    }
 return poscor;
}

void SeltzerBergerBremsModel::LoadDCSData() {
   using geant::MeV;
   using geant::millibarn;

   // get the path to the main physics data directory
   char *path = std::getenv("GEANT_PHYSICS_DATA");
   if (!path) {
     std::cerr<<"******   ERROR in SeltzerBergerBremsModel::LoadDCSData() \n"
              <<"         GEANT_PHYSICS_DATA is not defined! Set the GEANT_PHYSICS_DATA\n"
              <<"         environmental variable to the location of Geant data directory!\n"
              <<std::endl;
     exit(1);
   }

   char baseFilename[512];
   switch (fDataFileIndx) {
     case 0: sprintf(baseFilename,"%s/brems/NIST_BREM1/nist_brems_",path);
             break;
     case 1: sprintf(baseFilename,"%s/brems/NIST_BREM/nist_brems_",path);
             break;
     case 2: sprintf(baseFilename,"%s/brems/NRC_BREM/nist_brems_",path);
             break;
     default: sprintf(baseFilename,"%s/brems/NIST_BREM1/nist_brems_",path);
   }

   FILE *f = nullptr;
   char filename[512];
   sprintf(filename,"%sgrid",baseFilename);
   f = fopen(filename,"r");
   if (!f) {
     std::cerr<<"******   ERROR in SeltzerBergerBremsModel::LoadDCSData() \n"
              <<"         "<< filename << " could not be found!\n"
              <<std::endl;
     exit(1);
   }
   // before we take the fDCSMaxZet make sure that the fDCSForElements if free
   if (fLoadDCSForElements) {
     for (int i=0; i<fDCSMaxZet; ++i)
       if (fLoadDCSForElements[i]) {
         delete [] fLoadDCSForElements[i];
         fLoadDCSForElements[i] = nullptr;
       }
     delete [] fLoadDCSForElements;
   }

   int nmatches = fscanf(f,"%d%d%d",&fDCSMaxZet,&fLoadDCSNumElectronEnergies,&fLoadDCSNumReducedPhotonEnergies);
   (void)nmatches;
   // allocate space for the elemental DCS data, for the electron energy and reduced photon energy grids and load them
   fLoadDCSForElements = new double*[fDCSMaxZet];
   for (int i=0; i<fDCSMaxZet; ++i)
     fLoadDCSForElements[i] = nullptr;
   if (fLoadDCSElectronEnergyGrid) {
     delete [] fLoadDCSElectronEnergyGrid;
     fLoadDCSElectronEnergyGrid = nullptr;
   }
   fLoadDCSElectronEnergyGrid = new double[fLoadDCSNumElectronEnergies];
   if (fLoadDCSReducedPhotonEnergyGrid) {
     delete [] fLoadDCSReducedPhotonEnergyGrid;
     fLoadDCSReducedPhotonEnergyGrid = nullptr;
   }
   fLoadDCSReducedPhotonEnergyGrid = new double[fLoadDCSNumReducedPhotonEnergies];
   for (int i=0;i<fLoadDCSNumElectronEnergies;++i) {
     nmatches = fscanf(f,"%lf",&(fLoadDCSElectronEnergyGrid[i]));
     fLoadDCSElectronEnergyGrid[i] *= MeV;   // change to internal energy units
   }
   for (int i=0;i<fLoadDCSNumReducedPhotonEnergies;++i)
     nmatches = fscanf(f,"%lf",&(fLoadDCSReducedPhotonEnergyGrid[i]));
   fclose(f);
   //
   // for (int i=0;i<fDCSNumElectronEnergies;++i)
   //   printf("%.6e\n",fDCSElectronEnergyGrid[i]);
   // for (int i=0;i<fDCSNumReducedPhotonEnergies;++i)
   //   printf("%.6e\n",fDCSReducedPhotonEnergyGrid[i]);

   // now go for each element that we have in the global element table and load data for them
   int numDCSdataPerElement = fLoadDCSNumElectronEnergies*fLoadDCSNumReducedPhotonEnergies;
   const Vector_t<Element*> theElements = Element::GetTheElementTable();
   //std::cout<<theElements;
   int numElements = theElements.size();
   for (int i=0; i<numElements; ++i) {
     int zet = std::lrint(theElements[i]->GetZ());
     sprintf(filename,"%s%d",baseFilename,zet);
     f = fopen(filename,"r");
     if (!f) {
       std::cerr<<"******   ERROR in SeltzerBergerBremsModel::LoadDCSData() \n"
                <<"         "<< filename << " could not be found!\n"
                <<std::endl;
       exit(1);
     }
     // allocate space for this elemental DCS
      fLoadDCSForElements[zet-1] = new double[numDCSdataPerElement];
      for (int j=0; j<numDCSdataPerElement; ++j) {
        nmatches = fscanf(f,"%lf",&(fLoadDCSForElements[zet-1][j]));
        fLoadDCSForElements[zet-1][j] *= millibarn;   // change to internal units
      }
     fclose(f);
   }
}

}   // namespace geantphysics
