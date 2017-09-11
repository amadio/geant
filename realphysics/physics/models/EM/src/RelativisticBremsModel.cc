
#include "RelativisticBremsModel.h"
// from amterial
#include "Types.h"


#include "PhysicalConstants.h"
#include "Material.h"
#include "Element.h"
#include "MaterialProperties.h"
#include "MaterialCuts.h"

#include "GLIntegral.h"
#include "AliasTable.h"

#include "PhysicsParameters.h"

#include "Gamma.h"

#include "LightTrack.h"
#include "PhysicsData.h"

// from geantV
#include "GeantTaskData.h"

#include <cmath>

namespace geantphysics {


void RelativisticBremsModel::Initialize() {
  EMModel::Initialize();
  Initialise();
  // if it needs element selector: particle is coded in the model in case of this model and due to the cut dependence
  // we we need element selectors per MaterialCuts
  // InitialiseElementSelectors(this, nullptr, false);
}

double RelativisticBremsModel::ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle*, bool istotal){
  const Material *mat =  matcut->GetMaterial();
  const double *cuts  =  matcut->GetProductionCutsInEnergy();
  double gammacut     =  cuts[0];
  if (istotal) {
    // for the total stopping power we just need a gamma production cut >=kinenergy
    gammacut = 1.01*kinenergy;
  }
  return ComputeDEDXPerVolume(mat, gammacut, kinenergy);
}

double RelativisticBremsModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle*) {
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


double RelativisticBremsModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut, double kinenergy, const Particle*) {
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


int RelativisticBremsModel::SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td) {
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
//  double elInitTotalEnergy   = ekin+geant::kElectronMassC2;
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


double RelativisticBremsModel::MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle*) const {
  double mine = (matcut->GetProductionCutsInEnergy())[0]; // gamma production cut in the given material-cuts
  return mine;
}

// use these elastic and inelatic form factors for light elements instead of TFM
// under the complete screening approximation
// Tsai Table.B2.
const double RelativisticBremsModel::gFelLowZet  [] = {0.0, 5.310, 4.790, 4.740, 4.710, 4.680, 4.620, 4.570};
const double RelativisticBremsModel::gFinelLowZet[] = {0.0, 6.144, 5.621, 5.805, 5.924, 6.012, 5.891, 5.788};

RelativisticBremsModel::RelativisticBremsModel(bool isuselpm, const std::string &modelname)
 : EMModel(modelname),fIsUseLPM(isuselpm) {
   // set default values:
   //
   // the following variables(+the fIsUseLPM) can be changed through public setters:
   //   ! All changes must be done before initialisation
   fIsUseLinearSamplingTable  =     true;     // use linear approximation of the p.d.f. within bins
   fNumSamplingElecEnergies   =       41;     // e-/e+ kinetic energy bin numbers
   fNumSamplingPhotEnergies   =      100;     // number of emitted photon energy related transformed variables
   fMinElecEnergy             =      1.0*geant::GeV; // minimum e-/e+ kinetic energy
   fMaxElecEnergy             =    100.0*geant::TeV; // maximum e-/e+ kinetic energy

   // the following variables will be set and used internally by the model
   fElEnLMin                  =      0.0;     // will be set at initialisation
   fElEnILDelta               =      1.0;     // will be set at initialisation
   fSamplingElecEnergies      =  nullptr;     // will be set at initialisation
   fLSamplingElecEnergies     =  nullptr;     // will be set at initialisation

   fNumDifferentMaterialGCuts =        0;      // will be set at initialisation
   fGlobalMatGCutIndxToLocal  =  nullptr;      // will be set at initialisation
   fAliasData                 =  nullptr;      // will be set at initialisation
   fRatinAliasData            =  nullptr;      // will be set at initialisation
   fAliasSampler              =  nullptr;      // will be set at initialisation

   fAliasSampler              =  new AliasTable(); // create an alias sampler util used for run time sampling

   fSecondaryInternalCode     = -1;

}

void RelativisticBremsModel::Initialise() {
   if (fIsUseLinearSamplingTable) {
     InitSamplingTables();  // build the linear approximation(pdf) based alias sampling tables
   } else {
     InitSamplingTables1(); // build the rational interpolation(cdf) based alias sampling tables
   }
   fSecondaryInternalCode = Gamma::Definition()->GetInternalCode();
}


RelativisticBremsModel::~RelativisticBremsModel() {
   if (fSamplingElecEnergies) {
     delete [] fSamplingElecEnergies;
     fSamplingElecEnergies = nullptr;
   }
   if (fLSamplingElecEnergies) {
     delete [] fLSamplingElecEnergies;
     fLSamplingElecEnergies = nullptr;
   }
   if (fGlobalMatGCutIndxToLocal) {
     delete [] fGlobalMatGCutIndxToLocal;
     fGlobalMatGCutIndxToLocal = nullptr;
   }

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

   if (fRatinAliasData) {
     for (int i=0; i<fNumDifferentMaterialGCuts; ++i) {
       for (int j=0; j<fNumSamplingElecEnergies; ++j) {
         int indx = i*fNumSamplingElecEnergies+j;
         if (fRatinAliasData[indx]) {
           delete [] fRatinAliasData[indx]->fXdata;
           delete [] fRatinAliasData[indx]->fAliasW;
           delete [] fRatinAliasData[indx]->fA;
           delete [] fRatinAliasData[indx]->fB;
           delete [] fRatinAliasData[indx]->fAliasIndx;
           delete fRatinAliasData[indx];
         }
       }
     }
     delete [] fRatinAliasData;
   }

   if (fAliasSampler) {
     delete fAliasSampler;
   }
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
double RelativisticBremsModel::SamplePhotonEnergy(const MaterialCuts *matcut, double eekin, double r1, double r2, double r3) {
   constexpr double  mgdl  = 4.0*geant::kPi*geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght*geant::kRedElectronComptonWLenght;
   // gamma production cut in energy
   double gcut = matcut->GetProductionCutsInEnergy()[0];
/*
// this is checked now in the caller
   // protection: make sure that post interaction particle kinetic energy is within the range of the model
   // should be nicer by using std::max and std::min
   if (fMaxElecEnergy<eekin){
     std::cerr<< "\n **** WARNING: RelativisticBremsModel::SamplePhotonEnergy()\n"
              << "          Particle energy = " << eekin/geant::GeV << " [GeV] > fMaxElecEnergy = " << fMaxElecEnergy << " [GeV]"
              << std::endl;
     eekin = fMaxElecEnergy;
   } else if (fMinElecEnergy>eekin) {
     std::cerr<< "\n **** WARNING: RelativisticBremsModel::SamplePhotonEnergy()\n"
              << "          Particle energy = " << eekin/geant::GeV << " [GeV] < fMinElecEnergy = " << fMinElecEnergy << " [GeV]"
              << std::endl;
     eekin = fMinElecEnergy;
   }

   // protection: no gamma production when the post-interaction particle kinetic energy is below the gamma production cut
   // (at that primary energy the x-section should be zero so it should never happen)
   if (eekin<gcut) {
     return 0.0;
   }
*/
   double densityFactor = matcut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*mgdl;
   double etot          = eekin+geant::kElectronMassC2;
   double densityCor    = densityFactor*etot*etot; // this is k_p^2

   int    mcindx        = matcut->GetIndex();
   int    macindxlocal  = fGlobalMatGCutIndxToLocal[mcindx];
   // the location of the first-electron-energy sampling table data structure of this material-cut pair
   int    indxstart     = macindxlocal*fNumSamplingElecEnergies;
   // determine electron energy lower grid point
   double leenergy      = std::log(eekin);
   int    eenergyindx   = (int) ((leenergy-fElEnLMin)*fElEnILDelta);

   if (eenergyindx>=fNumSamplingElecEnergies-1)
     eenergyindx = fNumSamplingElecEnergies-2;

   // prob of taking the lower energy edge of the bin
   double ploweener = (fLSamplingElecEnergies[eenergyindx+1]-leenergy)*fElEnILDelta;

   // the location of the lower kinetic energy sampling table data structure
   indxstart += eenergyindx;
   bool   isAtLimit = false;
   if (fIsUseLinearSamplingTable && !(fAliasData[indxstart]))
     isAtLimit = true;
   if (!fIsUseLinearSamplingTable && !(fRatinAliasData[indxstart]))
     isAtLimit = true;

   if (r1>ploweener || isAtLimit)
     ++indxstart;

   // sample the transformed variable
   double gammae = 0.0;
   if (fIsUseLinearSamplingTable) {
     gammae = fAliasSampler->SampleLinear(fAliasData[indxstart]->fXdata, fAliasData[indxstart]->fYdata,
                                          fAliasData[indxstart]->fAliasW, fAliasData[indxstart]->fAliasIndx,
                                          fNumSamplingPhotEnergies,r2,r3);
   } else {
     gammae = fAliasSampler->SampleRatin(fRatinAliasData[indxstart]->fXdata, fRatinAliasData[indxstart]->fC,
                                         fRatinAliasData[indxstart]->fA,fRatinAliasData[indxstart]->fB,
                                         fRatinAliasData[indxstart]->fAliasW, fRatinAliasData[indxstart]->fAliasIndx,
                                         fRatinAliasData[indxstart]->fNumdata,r2,r3);
   }
   // transform back to gamma energy
   double dum0  = (gcut*gcut+densityCor);
   double dum1  = (eekin*eekin+densityCor)/dum0;
   double    u  = dum0*std::exp(gammae*std::log(dum1));

  return std::sqrt(u-densityCor);
}


// the simple DipBustgenerator
void RelativisticBremsModel::SamplePhotonDirection(double elenergy, double &sinTheta, double &cosTheta, double rndm) {
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


void RelativisticBremsModel::InitSamplingTables() {
   // set up the common electron/positron energy grid
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
   }

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
   fAliasData  = new LinAlias*[idum];
   for (int i=0; i<idum; ++i)
     fAliasData[i] = nullptr;

   for (int i=0; i<numMaterialCuts; ++i) {
     //std::cerr<<"   See if Material-Gcut ==> " <<theMaterialCutsTable[i]->GetMaterial()->GetName()<<"  gCut = "<< theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]<<std::endl;
     int localindx = fGlobalMatGCutIndxToLocal[i];
     if (localindx<0) {
       continue;
     }
     int ialias    = localindx*fNumSamplingElecEnergies;
     if (!isdone[localindx]) { // init for this material-gammacut pair is it has not been done yet i.e. the last elec-energy bin sampling data is still null.
//       std::cerr<<"   -> Will init for Material-Gcut ==> " <<theMaterialCutsTable[i]->GetMaterial()->GetName()<<"  gCut = "<< theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]<<std::endl;
       BuildOneLinAlias(ialias, theMaterialCutsTable[i]->GetMaterial(), theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]);
       isdone[localindx] = 1;
     }
   }
   delete [] isdone;

   // print
//   for (int i=0; i<numMaterialCuts; ++i)
//     std::cerr<<"     --> Global MatCut-indx = "<< i << " local indx = "<<fGlobalMatGCutIndxToLocal[i] <<std::endl;
}

/**
 *  The distribution of the energy of the emitted bremsstrahlung photons \f$ k \f$ is determined by the differential
 *  cross section (see the documentation of RelativisticBremsModel::ComputeDXSecPerAtom and
 *  RelativisticBremsModel::ComputeURelDXSecPerAtom).
 *  So
 *  \f[
 *    p(k) \propto \frac{\mathrm{d}\sigma}{\mathrm{d}k}
 *         \propto \left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \frac{1}{k\Gamma(k)}
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
 *  - first \f$ k \to u(k)=\ln(k^2+k_p^2)\f$ that removes the strong low-energy \f$ k \f$ dependence introduced by the
 *    dielectric suppression factor:
 *
 *    \f$ k=\sqrt{\exp(u)-k_p^2} \f$ and \f$ \frac{\mathrm{d}k}{\mathrm{d}u}= \frac{\exp(u)}{2\sqrt{\exp(u)-k_p^2}}
 *    \left( = \frac{k^2+k_p^2}{2k} \right)\f$ so from the condition \f$ p(k)\mathrm{d}k \propto p(u)\mathrm{d}u\f$
 *    one can write \f$ p(u) \propto p(k)\frac{\mathrm{d}k}{\mathrm{d}u} = p(k)\frac{k^2+k_p^2}{2k} \propto
 *    \left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \frac{1}{k}\frac{k^2}{k^2+k_p^2}\frac{k^2+k_p^2}{2k} \propto
 *    \left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \f$
 *
 *  - then we apply the transformation \f$ u \to \xi=\frac{u-\ln(k_c^2+k_p^2)}{\ln[E^2+k_p^2]-\ln(k_c^2+k_p^2)}\;
 *    \in [0,1]\f$ where \f$ E \f$ is the kinetic energy of the incident particle:
 *
 *    \f$ u = \ln\frac{E^2+k_p^2}{k_c^2+k_p^2}\xi + \ln(k_c^2+k_p^2) \f$ and \f$ \frac{\mathrm{d}u}{\mathrm{d}\xi} =
 *    \ln\frac{E^2+k_p^2}{k_c^2+k_p^2} \f$ so from the condition \f$ p(u)\mathrm{d}u \propto p(\xi)\mathrm{d}\xi\f$
 *    one can write \f$ p(\xi) \propto p(u)\frac{\mathrm{d}u}{\mathrm{d}\xi} =  p(u) \ln\frac{E^2+k_p^2}{k_c^2+k_p^2}
 *    \propto \left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \ln\frac{E^2+k_p^2}{k_c^2+k_p^2} \propto
 *    \left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \f$ (the last factor was dropped since it is just a constant).
 *
 *  - therefore, insted of building tables to sample directly the emitted photon energy, we build tables to sample
 *    \f$ \xi \in [0,1]\f$ at some fixed, discrete values of incident particle kinetic energies
 *    \f$\left\{ E_{i}^{kin} \right\}_{i=1}^{N}\f$ (only at the kinematically allowed combinations of
 *    \f$(E_{i}^{kin},k_c)\f$) and at run-time for a given incident particle kinetic enegy \f$ E^{kin} \f$ the
 *    corresponding photon energy \f$ k \f$, distributed according to \f$ p(k) \f$, can be obtained by sampling
 *    \f$ \xi \f$ from the appropriate table and transforming back as
 *    \f[
 *     k=\sqrt{\exp(u)-k_p^2}=\sqrt{[k_c^2+k_p^2]\exp\left[ \xi \ln\frac{E^2+k_p^2}{k_c^2+k_p^2}
 *     \right] - k_p^2}
 *    \f]
 *    where \f$ E \f$ is the actual kinetic energy of the incident particle i.e.
 *    \f$ E_{i}^{kin} < E \equiv E^{kin} < E_{i+1}^{kin} \f$ and it is guaranted that \f$ k \f$ properly lies in the
 *    actual kinematically allowed energy range i.e. \f$ k\in [k_c,E^{kin}]\f$.
 */
void RelativisticBremsModel::BuildOneLinAlias(int ialias, const Material *mat, double gcut) {
  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  // material dependent constant
  double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
                        *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;
  // target i.e. Z dependent constant
  constexpr double lpmConstant = 0.25*geant::kFineStructConst*geant::kElectronMassC2*geant::kElectronMassC2/(geant::kPi*geant::kHBarPlanckCLight);
  // material dependent constant
  double lpmEnergy   = mat->GetMaterialProperties()->GetRadiationLength()*lpmConstant;
  // material dependent constant: the e-/e+ energy at which the corresponding k_LPM=k_p
  double energyThLPM = std::sqrt(densityFactor)*lpmEnergy;

  // we will need the element composition of this material
  const Vector_t<Element*> theElements    = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int   numElems = theElements.size();

  for (int ieener=0; ieener<fNumSamplingElecEnergies; ++ieener) {
    double eener = fSamplingElecEnergies[ieener];
    if (eener>gcut) { // otherwise no gamma production at that e- energy so let the sampling table to be nullptr
      // particle is always e-/e+
      double totalenergy = eener+geant::kElectronMassC2;
      double densityCor  = densityFactor*totalenergy*totalenergy; // this is k_p^2

      // create the alias data struct
      fAliasData[ialias] = new LinAlias();
      fAliasData[ialias]->fXdata     = new double[fNumSamplingPhotEnergies]();
      fAliasData[ialias]->fYdata     = new double[fNumSamplingPhotEnergies]();
      fAliasData[ialias]->fAliasW    = new double[fNumSamplingPhotEnergies]();
      fAliasData[ialias]->fAliasIndx = new int[fNumSamplingPhotEnergies]();

      // fill 3 values at 0,0.8,1 of the transformed variable
      int numdata = 3;
      fAliasData[ialias]->fXdata[0] = 0.0;
      fAliasData[ialias]->fXdata[1] = 0.8;
      fAliasData[ialias]->fXdata[2] = 1.0;

      for (int i=0; i<numdata; ++i) {
        double egamma = gcut;
        if (i==0) {
          egamma = gcut;
        } else if (i==2) {
          egamma = eener;
        } else {
          egamma = std::sqrt( (gcut*gcut+densityCor)*std::exp(fAliasData[ialias]->fXdata[i]*std::log((eener*eener+densityCor)/(gcut*gcut+densityCor) )) -densityCor);
        }
        double dcross = 0.0;
        for (int ielem=0; ielem<numElems; ++ielem) {
          double zet = theElements[ielem]->GetZ();
          double val = 0.0;
          if (totalenergy>energyThLPM)
            val = ComputeURelDXSecPerAtom(egamma, totalenergy, zet, lpmEnergy, densityCor);
          else
            val = ComputeDXSecPerAtom(egamma, totalenergy, zet);
          dcross += theAtomicNumDensityVector[ielem]*zet*zet*val;
        }
        fAliasData[ialias]->fYdata[i] = dcross;
      }

//      for (int i=0; i<numdata; ++i) {
//        std::cerr<< "  mat = "<<mat->GetName() <<" xi = "<<fAliasData[ialias]->fXdata[i] << " pdf = "<< fAliasData[ialias]->fYdata[i] << std::endl;
//      }

      // expand the data up to maximum
      while(numdata<fNumSamplingPhotEnergies) {
        // find the lower index of the bin, where we have the biggest linear interp. error compared to computed value
        double maxerr     = 0.0; // value of the current maximum error
        double thexval    = 0.0;
        double theyval    = 0.0;
        int    maxerrindx = 0;   // the lower index of the corresponding bin
//std::cerr<<" numdata = "<<numdata<<std::endl;
        for (int i=0; i<numdata-1; ++i) {
//          std::cerr<< "   "<<i<<"-th = "<<fAliasData[ialias]->fXdata[i+1]<<std::endl;
          double xx = 0.5*(fAliasData[ialias]->fXdata[i]+fAliasData[ialias]->fXdata[i+1]); // mid point
          double yy = 0.5*(fAliasData[ialias]->fYdata[i]+fAliasData[ialias]->fYdata[i+1]); // lin func val at the mid point

          double dum0  = (gcut*gcut+densityCor);
          double dum1  = (eener*eener+densityCor)/dum0;
          double    u  = dum0*std::exp(xx*std::log(dum1));
          double egamma = std::sqrt( u-densityCor);

          double dcross = 0.0;
          for (int ielem=0; ielem<numElems; ++ielem) {
            double zet = theElements[ielem]->GetZ();
            double val = 0.0;
            if (totalenergy>energyThLPM)
              val = ComputeURelDXSecPerAtom(egamma, totalenergy, zet, lpmEnergy, densityCor);
            else
              val = ComputeDXSecPerAtom(egamma, totalenergy, zet);
            dcross += theAtomicNumDensityVector[ielem]*zet*zet*val;
          }

//          double err   = std::fabs(yy-dcross);
//          double err   = std::fabs(1.0-yy/dcross);
           double err   = std::fabs((yy-dcross)*(1.0-yy/dcross));

          if (err>maxerr) {
            maxerr     = err;
            maxerrindx = i;
            thexval    = xx;
            theyval    = dcross;
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

      // set up a linear alias sampler on this data
      AliasTable *alst = new AliasTable();
      alst->PreparLinearTable(fAliasData[ialias]->fXdata, fAliasData[ialias]->fYdata,
                              fAliasData[ialias]->fAliasW, fAliasData[ialias]->fAliasIndx,
                              fNumSamplingPhotEnergies);
      delete alst;
    }
    ++ialias;
  }
}


/**
 *   The restricted atomic cross section for bremsstrahlung photon emeission for target element with atomic number
 *   \f$Z\f$, gamma photon production energy threshold \f$k_c\f$ is computed for e-/e+ kinetic energy \f$E\f$ according to
 *   \f[
 *     \sigma(E;k_c,Z) = \int_{k_c}^{E} \frac{\mathrm{d}\sigma}{\mathrm{d}k}\mathrm{d}k =
 *        \frac{16 \alpha r_e^2 Z^2}{3} \int_{k_c}^{E} \frac{1}{k\Gamma}
 *        \left( \frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \mathrm{d}k
 *   \f]
 *   if \f$E>k_c\f$ and immediate return with \f$0\f$ otherwise.
 *   At e-/e+ total energies \f$E_t\f$such that the corresponding \f$ k_{LPM} = E_t^2/E_{LPM} < k_p\f$ dielectric
 *   suppression overwhelms LPM suppression and only LPM suppression is observable we turn off LPM suppression i.e.
 *   \f$ \left( \frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \f$ is computed by
 *   RelativisticBremsModel::ComputeDXSecPerAtom(). Otherwise, when both LPM and dielectric suppression is active it is
 *   computed by RelativisticBremsModel::ComputeURelDXSecPerAtom(). The \f$1/\Gamma\f$ factor comes from the fact that
 *   dielectric suppression is always considered.
 *
 *   The integral is computed by 64-points Gauss-Legendre quadrature after applying the following transformations:
 *   - first the emitted photon energy is transformed \f$k\to u=\ln(k/E)\f$
 *   - then the following transformation is applied \f$u\to \xi = (u-\ln(k_c/E)/(\ln(E/k_c))) \in [0,1] \f$
 *
 *   The transformed integral
 *   \f[
 *     \int_{k_c}^{E} \frac{1}{k\Gamma} \left( \frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \mathrm{d}k =
 *     \int_{\ln(k_c/E)}^{0} \frac{1}{\exp(u)E\Gamma} \left( \frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^*
 *     \exp(u)E\mathrm{d}u
 *   = \ln\frac{E}{k_c} \int_{0}^{1} \frac{1}{\Gamma} \left( \frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \mathrm{d}\xi
 *   \f]
 *   where \f$\Gamma \f$ and \f$ (\mathrm{d}\sigma/\mathrm{d}k)^*\f$ must be evaluated at
 *   \f$ k = \exp(u)E= k_c\exp(\xi\ln(E/k_c))\f$ for a given \f$ \xi \f$.
 */
double RelativisticBremsModel::ComputeXSectionPerAtom(const Element *elem, const Material *mat, double gammaprodcutenergy, double particleekin) {
  double xsec = 0.0;
  if (particleekin<=gammaprodcutenergy) {
//    std::cerr<<" particleekin = "<<particleekin/geant::MeV<< " gammaprodcutenergy = "<<gammaprodcutenergy/geant::MeV << std::endl;
    return xsec;
  }
  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  // material dependent constant
  double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
                        *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;
  // particle mass should be set at initialisation
  double totalenergy = particleekin+geant::kElectronMassC2; // it is always e-/e+
  //
  double densityCor  = densityFactor*totalenergy*totalenergy; // this is k_p^2
  // target i.e. Z dependent constant
  constexpr double lpmConstant = 0.25*geant::kFineStructConst*geant::kElectronMassC2*geant::kElectronMassC2/(geant::kPi*geant::kHBarPlanckCLight);
  //
  constexpr double factor      = 16.0*geant::kFineStructConst*geant::kClassicElectronRadius*geant::kClassicElectronRadius/3.0;
  // material dependent constant
  double lpmEnergy   = mat->GetMaterialProperties()->GetRadiationLength()*lpmConstant;
  // material dependent constant: the e-/e+ energy at which the corresponding k_LPM=k_p
  double energyThLPM = std::sqrt(densityFactor)*lpmEnergy;

  double lKappaPrimePerCr = std::log(particleekin/gammaprodcutenergy);

  // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
  int ngl = 64;
  GLIntegral *gl = new GLIntegral(ngl,0.0,1.0);
  const std::vector<double> glW = gl->GetWeights();
  const std::vector<double> glX = gl->GetAbscissas();
  // do the integration on reduced photon energy transformd to log(kappa) that transformed to integral between 0-1
  double zet      = elem->GetZ();
  double zet2     = zet*zet;
  double integral = 0.0;
  for (int i=0;i<ngl;++i) {
    double egamma = std::exp(glX[i]*lKappaPrimePerCr)*gammaprodcutenergy;
    double val    = 0.0;
    if (totalenergy>energyThLPM)
      val = ComputeURelDXSecPerAtom(egamma, totalenergy, zet, lpmEnergy, densityCor);
    else
      val = ComputeDXSecPerAtom(egamma, totalenergy, zet);
    integral += glW[i]*val/(1.+densityCor/(egamma*egamma));
  }

  delete gl;

  return lKappaPrimePerCr*factor*zet2*integral;
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
 *  RelativisticBremsModel::ComputeXSectionPerAtom()) and \f$n_i\f$ is the number of atoms per unit volume of
 *  \f$i \f$-th element of the material that is \f$ n_i = \mathcal{N}\rho w_i/A_i \f$ where \f$\mathcal{N}\f$ is the
 *  Avogadro number, \f$\rho\f$ is the material density, \f$w_i\f$ is the proportion by mass (or mass fraction) of the
 *  \f$i \f$-th element and \f$A_i\f$ is the molar mass of the \f$i \f$-th element. The corresponding mean free path
 *  is \f$\lambda = 1/\Sigma \f$.
 */
double RelativisticBremsModel::ComputeXSectionPerVolume(const Material *mat, double gammaprodcutenergy, double particleekin) {
  double xsec = 0.0;
  if (particleekin<=gammaprodcutenergy) {
//    std::cerr<<" particleekin = "<<particleekin/geant::MeV<< " gammaprodcutenergy = "<<gammaprodcutenergy/geant::MeV << std::endl;
    return xsec;
  }
  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  // material dependent constant
  double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
                        *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;
  // particle mass should be set at initialisation
  double totalenergy = particleekin+geant::kElectronMassC2; // it's always e-/e+

  double densityCor  = densityFactor*totalenergy*totalenergy; // this is k_p^2
  // target i.e. Z dependent constant
  constexpr double lpmConstant = 0.25*geant::kFineStructConst*geant::kElectronMassC2*geant::kElectronMassC2/(geant::kPi*geant::kHBarPlanckCLight);
  // material dependent constant
  double lpmEnergy   = mat->GetMaterialProperties()->GetRadiationLength()*lpmConstant;
  // material dependent constant: the e-/e+ energy at which the corresponding k_LPM=k_p
  double energyThLPM = std::sqrt(densityFactor)*lpmEnergy;
  //
  double factor      = 16.0*geant::kFineStructConst*geant::kClassicElectronRadius*geant::kClassicElectronRadius/3.0;

  double lKappaPrimePerCr = std::log(particleekin/gammaprodcutenergy);

  // we will need the element composition of this material
  const Vector_t<Element*> theElements    = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int   numElems = theElements.size();

  // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
  int ngl = 64;
  GLIntegral *gl = new GLIntegral(ngl,0.0,1.0);
  const std::vector<double> glW = gl->GetWeights();
  const std::vector<double> glX = gl->GetAbscissas();
  // do the integration on reduced photon energy transformd to log(kappa) that transformed to integral between 0-1
  double integral = 0.0;
  for (int i=0;i<ngl;++i) {
    double egamma = std::exp(glX[i]*lKappaPrimePerCr)*gammaprodcutenergy;
    double sum    = 0.0;
    for (int ielem=0; ielem<numElems; ++ielem) {
      double zet  = theElements[ielem]->GetZ();
      double val  = 0.0;
      if (totalenergy>energyThLPM)
        val = ComputeURelDXSecPerAtom(egamma, totalenergy, zet, lpmEnergy, densityCor);
      else
        val = ComputeDXSecPerAtom(egamma, totalenergy, zet);
      sum  += theAtomicNumDensityVector[ielem]*zet*zet*val;
    }
    integral += glW[i]*sum/(1.+densityCor/(egamma*egamma));
  }

  delete gl;

  return lKappaPrimePerCr*factor*integral;
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
 *  Using the form of the applied differential cross section (see RelativisticBremsModel::ComputeDXSecPerAtom and
 *  RelativisticBremsModel::ComputeURelDXSecPerAtom one can write
 *  \f[
 *      S(E;k_c,\mathrm{material})=\frac{16\alpha r_e^2}{3}\int_{0}^{\eta} \frac{1}{\Gamma} \sum_i n_i
 *        Z_i^2\left( \frac{\mathrm{d}\sigma_i}{\mathrm{d}k} \right)^* \mathrm{d}k
 *  \f]
 *  The integral is computed by 64-points Gauss-Legendre quadrature after applying the following transformations:
 *  - first the emitted photon energy is transformed \f$ k \to u=\ln[1-\eta/E_t] \f$ where \f$E_t\f$ is the total energy
 *    of the e-/e+.
 *  - then the following transformation is applied \f$ u \to \xi = u/\ln[1-\eta/E_t] \in [0,1]\f$
 *
 *  The transformed integral
 *  \f[
 *  \int_{0}^{\eta} \frac{1}{\Gamma} \sum_i n_i Z_i^2\left( \frac{\mathrm{d}\sigma_i}{\mathrm{d}k} \right)^* \mathrm{d}k
 *  = -E_t \int_{0}^{\ln[1-\eta/E_t]} \frac{e^u}{\Gamma}
 *     \sum_i n_i Z_i^2\left( \frac{\mathrm{d}\sigma_i}{\mathrm{d}k} \right)^* \mathrm{d}u
 *  = -E_t \ln[1-\eta/E_t] \int_{0}^{1} \frac{e^{\xi\ln[1-\eta/E_t]}}{\Gamma}
 *     \sum_i n_i Z_i^2\left( \frac{\mathrm{d}\sigma_i}{\mathrm{d}k} \right)^* \mathrm{d}\xi
 *  \f]
 *  where \f$\Gamma\f$ and \f$ \left(\mathrm{d}\sigma_i / \mathrm{d}k \right)^* \f$ must be evaluated at
 *  \f$ k=E_t[1-e^{\xi\ln[1-\eta/E_t]}] \f$ for a given \f$\xi\f$.
 */
double RelativisticBremsModel::ComputeDEDXPerVolume(const Material *mat, double gammaprodcutenergy, double particleekin) {
  double dedx = 0.0;
  // TODO
  // HERE WE need to check if the particleekin is below the minimum energy that the model is set and we return zero in this case
  // if (GetLowEnergyLimit())


  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  // material dependent constant
  double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
                        *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;
  // particle mass should be set at initialisation
  double totalenergy = particleekin+geant::kElectronMassC2;  // it's always e-/e+

  double densityCor  = densityFactor*totalenergy*totalenergy; // this is k_p^2
  // target i.e. Z dependent constant
  constexpr double lpmConstant = 0.25*geant::kFineStructConst*geant::kElectronMassC2*geant::kElectronMassC2/(geant::kPi*geant::kHBarPlanckCLight);
  // material dependent constant
  double lpmEnergy   = mat->GetMaterialProperties()->GetRadiationLength()*lpmConstant;
  // material dependent constant: the e-/e+ energy at which the corresponding k_LPM=k_p
  double energyThLPM = std::sqrt(densityFactor)*lpmEnergy;
  // the constant factor
  double factor      = 16.0*geant::kFineStructConst*geant::kClassicElectronRadius*geant::kClassicElectronRadius/3.0;

  double upperlimit  = gammaprodcutenergy;
  if (upperlimit>particleekin)
    upperlimit = particleekin;
  double kappacr      = upperlimit/totalenergy;
  double log1mKappaCr = std::log(1.0-kappacr);

  // we will need the element composition of this material
  // const std::vector<Element*> theElements = mat->GetElementVector();
  const Vector_t<Element*> theElements    = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int   numElems = theElements.size();

  // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
  int ngl = 64;
  GLIntegral *gl = new GLIntegral(ngl,0.0,1.0);
  const std::vector<double> glW = gl->GetWeights();
  const std::vector<double> glX = gl->GetAbscissas();
  // do the integration on reduced photon energy transformd to log(kappa) that transformed to integral between 0-1
  double integral = 0.0;
  for (int i=0;i<ngl;++i) {
    double dumx   = 1.0-std::exp(glX[i]*log1mKappaCr);
    double egamma = totalenergy*dumx;
    double sum    = 0.0;
    for (int ielem=0; ielem<numElems; ++ielem) {
      double zet  = theElements[ielem]->GetZ();
      //int    izet = std::lrint(zet);
      double val  = 0.0;
      if (totalenergy>energyThLPM)
        val = ComputeURelDXSecPerAtom(egamma, totalenergy, zet, lpmEnergy, densityCor);
      else
        val = ComputeDXSecPerAtom(egamma, totalenergy, zet);
      sum  += theAtomicNumDensityVector[ielem]*zet*zet*val;
    }
    integral += glW[i]*sum*(egamma-totalenergy)/(1.+densityCor/(egamma*egamma));
  }

  delete gl;
  dedx = log1mKappaCr*factor*integral;
  return dedx;
}

/*
// Transformation: k \to u=ln[1-k/E_kin+1e-6] then u to xi = u-ln[1+1e-6]/ln((1+1e-6-upperlimit/e_kin)(1+1e-6))
//
double RelativisticBremsModel::ComputeDEDXPerVolume1(Material *mat, double gammaprodcutenergy, double particleekin) {
  double dedx = 0.0;
  // TODO
  // HERE WE need to check if the particleekin is below the minimum energy that the model is set and we return zero in this case

  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  // material dependent constant
  double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
                        *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;
  // particle mass should be set at initialisation
  double totalenergy = particleekin+geant::kElectronMassC2;  // it's always e-/e+

  double densityCor  = densityFactor*totalenergy*totalenergy; // this is k_p^2
  // target i.e. Z dependent constant
  constexpr double lpmConstant = 0.25*geant::kFineStructConst*geant::kElectronMassC2*geant::kElectronMassC2/(geant::kPi*geant::kHBarPlanckCLight);
  // material dependent constant
  double lpmEnergy   = mat->GetMaterialProperties()->GetRadiationLength()*lpmConstant;
  // material dependent constant: the e-/e+ energy at which the corresponding k_LPM=k_p
  double energyThLPM = std::sqrt(densityFactor)*lpmEnergy;
  //
  double factor      = 16.0*geant::kFineStructConst*geant::kClassicElectronRadius*geant::kClassicElectronRadius/3.0;

  double upperlimit  = gammaprodcutenergy;
  if (upperlimit>particleekin)
    upperlimit = particleekin;
  double kappacr      = upperlimit/particleekin;
  const  double phi   = 1e-6;
  const  double onePlusPhi = 1.0+phi;
  double log1mKappaCr = std::log(1.0-kappacr/onePlusPhi);

  // we will need the element composition of this material
  const std::vector<Element*> theElements = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int   numElems = theElements.size();

  // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
  int ngl = 64;
  GLIntegral *gl = new GLIntegral(ngl,0.0,1.0);
  const std::vector<double> glW = gl->GetWeights();
  const std::vector<double> glX = gl->GetAbscissas();
  // do the integration on reduced photon energy transformd to log(kappa) that transformed to integral between 0-1
  double integral = 0.0;
  for (int i=0;i<ngl;++i) {
    double dumx     = onePlusPhi*std::exp(glX[i]*log1mKappaCr);
    double egamma   = (onePlusPhi-dumx)*particleekin;
//    double egamma = totalenergy*dumx;
    double sum    = 0.0;
    for (int ielem=0; ielem<numElems; ++ielem) {
      double zet  = theElements[ielem]->GetZ();
      int    izet = std::lrint(zet);
      double val  = 0.0;
      if (totalenergy>energyThLPM)
        val = ComputeURelDXSecPerAtom(egamma, totalenergy, zet, lpmEnergy, densityCor);
      else
        val = ComputeDXSecPerAtom(egamma, totalenergy, zet);
      sum  += theAtomicNumDensityVector[ielem]*zet*zet*val;
    }
    integral -= glW[i]*sum*dumx/(1.+densityCor/(egamma*egamma));
  }

  delete gl;

  return log1mKappaCr*particleekin*factor*integral;
}
*/

// compute differential cross section per atom: Tsai, complete screening with LPM plus dielectric suppressions
// NOTE: to get the differential cross section the return value must be multiplied by (16 alpha r^2 Z^2)/(3*k)
//       where alpha is the fine structure constant, r is the classical electron radius, Z is the target atomic number and k is
//       the emitted photon energy.
/**
 *  When LPM efect is active, the atomic cross section, differential in emitted photon energy, is based on the
 *  <em>complete screening</em> form of Tsai's differential cross section \cite tsai1974pair [Eq. (3.83)] (see more
 *  details at RelativisticBremsModel::ComputeDXSecPerAtom):
 *  \f[
 *   \frac{\mathrm{d}\sigma}{\mathrm{d}k} = \frac{4\alpha r_e^2}{k}
 *   \left\{
 *   \left( \frac{4}{3} -\frac{4}{3}y+y^2 \right)
 *     \left[  Z^2\left(L_{el} -f \right) + ZL_{inel} \right] +
 *     \frac{1}{9}(1-y) \left[ Z^2+Z \right]
 *   \right\}
 *  \f]
 * where \f$\alpha\f$ is the fine structure constant, \f$r_e\f$ is the classical electron radius,
 * \f$k\f$ is the emitted photon energy, \f$y=k/E_t\f$ is the emitted photon energy in pre-interaction \f$e^-/e^+\f$ total
 *  energy (\f$ E_t\f$) units, \f$Z\f$ is the target atomic number, \f$f\f$ is the Coulomb correction \cite davies1954theory
 * [Eqs.(36-38)]
 *  \f[
 *   f(\nu) = \nu^2 \sum_{n=1}^{\infty} \frac{1}{n(n^2+\nu^2)} = \nu^2 \left[  1/(1+\nu^2)  + 0.20206 - 0.0369\nu^2
 *            + 0.0083\nu^4 - 0.002\nu^6 \right]
 *  \f]
 * where \f$\nu=\alpha Z\f$.  \f$L_{el}, L_{inel}\f$ are elastic and inelastic atomic form factor related variables
 * for which we use values from table for \f$Z<5\f$ and
 * \f$ L_{el} = \ln \left[ 184.1499 Z^{-1/3} \right] \f$ and  \f$ L_{inel} = \ln \left[ 1193.923 Z^{-2/3} \right] \f$
 * if \f$Z \geq 5\f$ (see more details on this at RelativisticBremsModel::ComputeDXSecPerAtom).
 * By pulling out \f$1/3\f$, taking into account \f$ 4-4y+3y^2=y^2+2[1+(1-y)^2]\f$, considering only elastic form factor
 * (interaction of electrons with the target nucleus and neglecting the interactions with the atomic electrons),
 * neglecting both the Coulomb correction and the last term proportional to \f$(1-y)\f$ the above differential cross
 * section transforms to the Bethe-Heitler one \cite bethe1934stopping
 * \f[
 *    \frac{\mathrm{d}\sigma_{BH}}{\mathrm{d}k} = \frac{4\alpha r_e^2}{3k}
 *     \left\{y^2 +2[1+(1-y)^2] \right\} Z^2 L_{el}
 * \f]
 * Migdal, based on quantum mechanical calculations, derived the following expression for the differential cross section
 * including LPM suppression \cite migdal1956bremsstrahlung
 * \f[
 *   \frac{\mathrm{d}\sigma_{M}^{LPM}}{\mathrm{d}k} = \frac{4\alpha r_e^2 \xi(s)}{3k}
 *    \left\{
 *     y^2 G(s) +2[1+(1-y)^2]\phi(s) Z^2 L_{el}
 *    \right\}
 * \f]
 * where \cite migdal1956bremsstrahlung[Eq.(46)]
 * \f[
 *    G(s) \equiv 24\left[  \frac{\pi}{2} - \int_{0}^{\infty} \exp(-st) \frac{\sin(st)}{\sinh(t/2)} \mathrm{d}t \right]
 *    = 12 \pi s^2 - 48 s^3\sum_{j=0}^{\infty} \frac{1}{(j+s+1/2)^2+s^2}
 * \f]
 * is the electron spin flip and \cite migdal1956bremsstrahlung[Eq.(47)]
 * \f[
 *   \phi(s) \equiv 12s^2\int_{0}^{\infty} \exp(-st) \coth(t/2)\sin(st) \mathrm{d}t - 6\pi s^2
 *    = 6s(1-\pi s) + 24s^3\sum_{j=1}^{\infty} \frac{1}{(j+s)^2+s^2}
 * \f]
 * is the electron no spin flip suppression function. Approximate expressions (within $0.15\%$) for the slowly
 * convergent series were derived in \cite stanev1982development [Eqs.(13-15)] (note that we use slightly different
 * approxiamtions both at the low and high \f$s\f$ cases)
 * \f[
 * \phi(s) \approx
 *  \begin{cases}
 *   6s-6\pi s^2 +4\pi^2 s^3 -48s^4\zeta(3) & \quad \text{if}\; s \leq 0.1 \\
 *   1-\exp \left\{  -6s[1+s(3-\pi)] + \frac{s^3}{0.623+0.796s+0.658s^2} \right\}  & \quad \text{if}\; 0.1 < s < 2 \\
 *   1-\frac{1}{84s^4} & \quad \text{if}\; s \geq 2
 *  \end{cases}
 * \f]
 * where instead of \cite stanev1982development we used the series expansion around \f$s\f$ zero and
 * infinity respectively
 * \f[
 *  24s^3\sum_{j=1}^{\infty} \frac{1}{(j+s)^2+s^2} = 12 s^2 i \left\{\psi^{(0)}[(1-i)s+1] - \psi^{(0)}[(1+i)s+1]\right\}  \\
 *  \begin{cases}
 *   \approx 4\pi^2 s^3-48\zeta(3) s^4 & \quad \mathrm{around}\; s=0 \\
 *   \approx 6s-6\pi s^2+1-1/(84s^4)   & \quad \mathrm{around}\; s=\infty
 *  \end{cases}
 * \f]
 * where \f$\psi^{(0)}(x)\f$ is the digamma function and \f$ \zeta(x) \f$ is the Riemann zeta function that
 * \f$ \zeta(3) \approx 1.202056903159594 \f$ is Apery's constant. And
 * \f[
 * G(s) \approx
 *  \begin{cases}
 *   12\pi s^2-24\pi^2 s^3 - 48 \psi^{2}(0.5) s^4  & \quad \text{if}\; s \leq 0.1 \\
 *   3\left\{1-\exp \left[-4s  - \frac{8s^2}{1+3.936s+4.97s^2-0.05s^3+7.5s^4} \right] \right\} -2\phi(s)  & \quad \text{if}\; 0.1 < s < 2 \\
 *   1 - 0.0230655/s^4  & \quad \text{if}\; s \geq 2
 *  \end{cases}
 * \f]
 * where instead of \cite stanev1982development we used the series expansion around \f$s\f$ zero and
 * infinity respectively
 * \f[
 *  48s^3\sum_{j=0}^{\infty} \frac{1}{(j+s+1/2)^2+s^2} = 24 s^2 i \left\{ \psi^{0}((1-i)s+0.5) - \psi^{0}((1+i)s+0.5) \right\} \approx \\
 *  \begin{cases}
 *   \approx 24\pi^2 s^3 + 48 \psi^{2}(0.5) s^4 & \quad \mathrm{around}\; s=0 \\
 *   \approx 12\pi s^2 -1 + 0.0230655/s^4   & \quad \mathrm{around}\; s=\infty
 *  \end{cases}
 * \f]
 * where \f$\psi^{2}(x)\f$ is the second derivative of the digamma function at 0.5 and \f$ 48 \psi^{2}(0.5) \approx
 * -807.782238923247359\f$.
 *
 * The variable \f$s\f$, (can be obtained from the LPM condition i.e. the relation of average photon emission angle to the
 * average multiple scattering deflection along the formation length) given by Migdal
 * \cite migdal1956bremsstrahlung [Eq.60.]
 * \f[
 *   s  = \frac{1}{2} \sqrt{ \frac{y}{1-y}   \frac{(mc^2)^2 \alpha X_0}{8\hbar c \pi} \frac{1}{E\xi(s)}  }
 * \f]
 * One can introduce the material dependent LPM energy
 * \f[
 *    E_{LPM}=\frac{ X_0\alpha m^2c^3}{ 4\pi \hbar} \left(=  \frac{X_0 ( mc^2)^4}{\hbar c E_s^2} \right)
 * \f]
 * and \f$s\f$ can be written as
 * \f[
 *    s  = \sqrt{\frac{1}{8} \frac{y}{1-y}   \frac{E_{LPM}}{E\xi(s)}  }
 * \f]
 * Migdal gave the following approximate expressions for \f$\xi(s)\f$ \cite migdal1956bremsstrahlung
 * \f[
 *  \xi(s) =
 *  \begin{cases}
 *  2 & \quad \text{if}\; s \leq s_1 \\
 *  1+\frac{\ln(s)}{\ln(s_1)} & \quad \text{if}\; s_1 < s < 1 \\
 *  1 & \quad \text{if}\; s \geq 1
 *  \end{cases}
 * \f]
 * where \f$ s_1 \equiv [Z^{1/3}/\exp(20.863/4)]^2 = Z^{2/3}/184.1499^2 \f$. Since \f$s\f$ depends on \f$\xi\f$ and
 * \f$\xi\f$ on \f$s\f$, an iterative procedure would be required to compute them. However, an approximate procedure
 * was suggested in \cite stanev1982development by introducing
 * \f[
 *  s'  = \sqrt{ \frac{1}{8} \frac{y}{1-y}   \frac{E_{LPM}}{E}  }
 * \f]
 * and \f$ h(s') \equiv (\ln(s'))/\ln(\sqrt{2}s_1) \f$
 * \f[
 * \xi(s') =
 * \begin{cases}
 *   2 & \quad \text{if}\; s' \leq \sqrt{2}s_1 \\
 *   1+h-\frac{0.08(1-h)[1-(1-h)^2]}{\ln(\sqrt{2}s_1)} & \quad \text{if}\;  \sqrt{2}s_1 < s' < 1 \\
 *   1 & \quad \text{if}\; s' \geq 1 \\
 * \end{cases}
 * \f]
 * and after determining \f$s'\f$ and \f$ \xi(s') \f$,  \f$ s\f$ can be obtained as
 * \f[
 *  s=\frac{s'}{\sqrt{\xi(s')}}
 * \f]
 * then \f$ \xi(s)\f$ can be computed by using Migdals's equations and the approximate expressions for
 * \f$G(s),\phi(s)\f$ can be used to compute them.
 *
 * According to Migdal, when both LPM and dielectric suppressions are active, one can include both effects by
 * replacing \f$\phi(s) \to \phi(s\Gamma)/\Gamma\f$ with \f$\Gamma = (1+k_p^2/k^2)\f$
 * where \f$ k_p = \hbar \omega_p \gamma= \hbar \omega_p E/(mc^2)\f$, \f$\hbar \omega_p= \hbar c \sqrt{4\pi n_e r_e}\f$
 * (\f$n_e\f$ is the electron density) and \f$k\f$ is the photon energy, and also replace \f$s\f$ with \f$s\Gamma\f$ in
 * \f$\xi(s) \to \xi(s\Gamma)\f$ and \f$G(s) \to G(s\Gamma)\f$. And the final form of the differential cross section
 * of bremsstrahlung photon emission
 * \f[
 * \frac{\mathrm{d}\sigma}{\mathrm{d}k} =  \frac{16 Z^2 \alpha r_e^2}{3k\Gamma}
 * \left\{ \xi(s\Gamma)
 *   \left[ \frac{y^2}{4}G(s\Gamma) +(1-y+\frac{y^2}{2}) \phi(s\Gamma) \right]
 *   \left[ (L_{el} -f)  +\frac{L_{inel}}{Z} \right]
 *  + \frac{1}{12}(1-y)\left[  1+\frac{1}{Z}  \right]
 * \right\}
 * \f]
 * What is computed below is \f$\left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^*\f$ part of the above expressions
 * such as
 *  \f[
 *    \frac{\mathrm{d}\sigma}{\mathrm{d}k} = \frac{16 \alpha r_e^2 Z^2}{3k\Gamma}
 *                                           \left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^*
 *  \f]
 *
 */
double RelativisticBremsModel::ComputeURelDXSecPerAtom(double gammaenergy, double totalenergy, double z, double lpmenergy, double densitycor) {
  //return zero if gamma energy is < 0
  if (gammaenergy<0.0) {
    return 0.0;
  }
  double y = gammaenergy/totalenergy; // reduced photon energy
  double firstTerm   = 0.0;
  double secondTerm  = 0.0;

  // Coulomb correction from Davis,Bethe,Maximom PRL 1954 Eqs.(36-38)
  // should be pulled out because it is constant for a given Z
  double mu  = z*geant::kFineStructConst;
  double mu2 = mu*mu;
  double mu4 = mu2*mu2;
  double mu6 = mu2*mu4;
  double coulombCorrection = mu2*(1.0/(1.0+mu2) + 0.20206 - 0.0369*mu2 + 0.0083*mu4 - 0.002*mu6);

  int    izet        = std::lrint(z);

  // elastic and inelastic screening function values in the complete screening case i.e. \gamma=\epsilon=0; used only for Z>=5
  // L_{el}   = [\ln (\exp(20.863/4)Z^{-1/3}) ] = \ln[184.1499 Z^{-1/3}]
  // L_{inel} = [\ln (\exp(28.340/4)Z^{-2/3}) ] = \ln[1193.923 Z^{-2/3}]
  // And these are the constant parts so they should be pulled out
  const double factorFel   = std::log(184.1499);
  const double factorFinel = std::log(1193.923);

  // radiation length: is a material property so it should be set by the caller
  // LPM energy is a material dependent constant so it should be set and provided by the caller

  //
  // compute the LPM functions s', s, \xi(s), G(s), \phi(s)
  //
  // DO IT ACCORDING TO Geant4 model:
  // -compute s according to Stanev to avoid iterative solution: s1, s', \xi(s') then s
  //  1. s_1 \equiv [Z^{1/3}/\exp(20.863/4)]^2 = Z^{2/3}/184.1499^2
  //     std::pow(z,2/3) should be pulled out as the constant 184.1499^2 =
  double varS1     = std::pow(z,2./3.)/(184.1499*184.1499);
  //  2. in Geant4 s' is computed according to Klein (!not like E_LPM or k_LPM bevause they computed like Anthony!)
  //     s'  = \sqrt{ \frac{1}{8} \frac{y}{1-y}   \frac{E^{KL}_{LPM}}{E}  }
  //     In Geant4 E^{AN}_{LPM} is used instead of E^{KL}_{LPM}.
  double varSprime = std::sqrt(0.125*y/(1.0-y)*lpmenergy/totalenergy);
  //  3. \xi(s') =
  //      \begin{cases}
  //        2 & \quad \text{if}\; s' \leq \sqrt{2}s_1
  //        1+h-\frac{0.08(1-h)[1-(1-h)^2]}{\ln(\sqrt{2}s_1)} & \quad \text{if}\;  \sqrt{2}s_1 < s' < 1
  //        1 & \quad \text{if}\; s' \geq 1
  //      \end{cases}
  // where h(s') \equiv (\ln(s'))/\ln(\sqrt{2}s_1)
  double funcXiSprime = 2.0;
  // sqrt(2) constant should be pulled out
  double condition    = std::sqrt(2.0)*varS1;
  if (varSprime>1.0) {
    funcXiSprime = 1.0;
  } else if (varSprime>condition) {
    double dum0        = 1.0/std::log(condition);
    double funcHSprime = std::log(varSprime)*dum0;
    funcXiSprime       = 1.0+funcHSprime-(0.08*(1.0-funcHSprime)*(2.0*funcHSprime-funcHSprime*funcHSprime))*dum0;
  }

  //  4. s=\frac{s'}{\sqrt{\xi(s')}}
  double varS    = varSprime/std::sqrt(funcXiSprime);
  // - include dielectric suppression effect into s according to Migdal i.e.
  //   change s to \hat{s}=s\Gamma  where Gamma = (1+k_p^2/k^2) where k_p = \hbar \omega_p \gamma= \hbar \omega_p E/(mc^2)
  //   where \hbar \omega_p= \hbar c \sqrt{4\pi n_e r_e}\f$ (\f$n_e\f$ is the electron density) and \f$k\f$ is the photon energy
  //   NOTE: that k_p^2 is our usual densityCor variable (where we compute the plasma energy so it is not taken from e.g. NIST data)
  //         densityFactor that is k_p^2/(E_t^2) is material dependent parameter!!! so should be provided by the caller
  double varShat = varS*(1.0+densitycor/(gammaenergy*gammaenergy));
  // - since we already know s compute \xi(s) according to Migdal
  //   by replacing s with \hat{s}
  // \xi(s) =
  //   \begin{cases}
  //    2 & \quad \text{if}\; s \leq s_1
  //    1+\frac{\ln(s)}{\ln(s_1)} & \quad \text{if}\; s_1 < s < 1
  //    1 & \quad \text{if}\; s \geq 1
  //    \end{cases}
  double funcXiS = 2.0;
  if (varShat>1.0) {
    funcXiS = 1.0;
  } else if (varShat>varS1) {
    funcXiS = 1.0+std::log(varShat)/std::log(varS1);
  }
/*
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
  double funcPhiS = 6.0*varShat*(1.0-geant::kPi*varShat);
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
  double funcGS   = 3.0*funcPsiS-2.0*funcPhiS;
*/

//// begin according to the genat4 model
  // - compute \phi(s) and G(s) suppression functions according Geant4 model (slightly different than Stanev)
  double varShat2 = varShat*varShat;
  double varShat3 = varShat*varShat2;
  double varShat4 = varShat2*varShat2;
  double funcPhiS = 1.0;
  double funcGS   = 1.0;
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


  //MAKE SURE SUPPRESSION IS SMALLER THAN 1: due to Migdal's approximation on xi
  if (funcXiS*funcPhiS>1. || varShat>0.57) {
    funcXiS=1./funcPhiS;
  }

  // -compute the dxsection: Tsai, complete screening, LMP according to Migdal
  //  NOTE: - what is computed is not the expression below! one needs to multiply this by (16 alpha r^2 Z^2)/(3k\Gamma) to get it
  //        - s\Gamma is varShat, \xi(s\Gamma) is funcXiS, G(s\Gamma) is funcGS, \phi(s\Gamma) is funcPhiS and all are
  //          computed above
  // \begin{split}
  //   \frac{\mathrm{d}\sigma}{\mathrm{d}k} = & \frac{16 Z^2 \alpha r_e^2}{3k\Gamma}
  //   \left\{ \xi(s\Gamma)
  //   \left[ \frac{y^2}{4}G(s\Gamma) +(1-y+\frac{y^2}{4}) \phi(s\Gamma) \right]
  //   \left[ (L_{el} -f)  +\frac{L_{inel}}{Z} \right]
  //  \right. \\ & \left.
  //   + \frac{1}{12}(1-y)\left[  1+\frac{1}{Z}  \right]
  //   \right\}
  //  \end{split}

  // don't forget: complete screening is used (fine because LPM is active at high energies and at high energy screening
  //               functions becoms independent on gamma and epsilon)
  double dum0 = 0.25*y*y;
  firstTerm   = funcXiS*(dum0*funcGS+(1.0-y+2.0*dum0)*funcPhiS);
  secondTerm  = (1.0-y)/12.0*(1.0+1.0/z);
  if (izet<5) {  // use complete screening (Tsai Eq.(3.83)) and L_el and L_inel from Tsai Table B2 for Z<5
    firstTerm  *= (gFelLowZet[izet]-coulombCorrection+gFinelLowZet[izet]/z);
  } else {       // use complete screening form of screening functions for Z>=5  case
    double dum1 = std::log(z)/3.0;
    firstTerm  *= (factorFel-dum1-coulombCorrection + (factorFinel-2.0*dum1)/z);
  }

  return (firstTerm+secondTerm);
}


/**
 *
 * The atomic cross section, differential in emitted photon energy, for \f$Z\geq5\f$ is computed according to Tsai
 * \cite tsai1974pair  [Eq. 3.82]
 *  \f[
 *    \begin{split}
 *    \frac{\mathrm{d}\sigma}{\mathrm{d}k} = \frac{\alpha r_e^2}{k} &
 *    \left\{
 *    \left(\frac{4}{3} -\frac{4}{3}y+y^2\right)
 *      \left[  Z^2\left(\varphi_1(\gamma) -\frac{4}{3}\ln(Z) -4f \right) + Z\left(\psi_1(\epsilon) -\frac{8}{3}\ln(Z) \right)  \right]
 *    \right. \\ & \left. \quad
 *      +\frac{2}{3}(1-y) \left[ Z^2(\varphi_1(\gamma) - \varphi_2(\gamma)) + Z(\psi_1(\epsilon) - \psi_2(\epsilon)) \right]
 *    \right\}
 *    \end{split}
 *  \f]
 * that we use in a slightly different form by pulling out \f$ \frac{16Z^2}{3} \f$
 *  \f[
 *    \frac{\mathrm{d}\sigma}{\mathrm{d}k} = \frac{16 \alpha r_e^2 Z^2}{3k} \left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^*
 *  \f]
 * and we compute
 *  \f[
 *    \begin{split}
 *    \left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* = &
 *    \left(1 - y + \frac{3}{4}y^2\right)
 *      \left[  \left(\frac{1}{4}\varphi_1(\gamma) -\frac{1}{3}\ln(Z) -f \right) + \frac{1}{Z}\left(\frac{1}{4}\psi_1(\epsilon) -\frac{2}{3}\ln(Z) \right)  \right]
 *    \\ & \quad
 *     + \frac{1}{8}(1-y) \left[ (\varphi_1(\gamma) - \varphi_2(\gamma)) + \frac{1}{Z}(\psi_1(\epsilon) - \psi_2(\epsilon)) \right]
 *    \end{split}
 *  \f]
 * in this method instead of \f$ \frac{\mathrm{d}\sigma}{\mathrm{d}k} \f$.
 * In the above equations, \f$\alpha\f$ is the fine structure constant, \f$r_e\f$ is the classical electron radius,
 * \f$k\f$ is the emitted photon energy, \f$y=k/E_t\f$ is the emitted photon energy in pre-interaction \f$e^-/e^+\f$ total
 *  energy (\f$ E_t\f$) units, \f$Z\f$ is the target atomic number, \f$f\f$ is the Coulomb correction \cite davies1954theory
 * [Eqs.(36-38)] (accounts that the interaction with electrons takes place in the field of the nucleus)
 *  \f[
 *   f(\nu) = \nu^2 \sum_{n=1}^{\infty} \frac{1}{n(n^2+\nu^2)} = \nu^2 \left[  1/(1+\nu^2)  + 0.20206 - 0.0369\nu^2
 *            + 0.0083\nu^4 - 0.002\nu^6 \right]
 *  \f]
 * where \f$\nu=\alpha Z\f$.
 * The elastic \f$\varphi_1, \varphi_2\f$ and inelastic \f$\psi_1, \psi_2\f$ screening functions are given in
 * \cite tsai1974pair [Eq.(3.14)-(3.17)]. Since these functions involve integrals that can be performed only numerically,
 * approximate expressions are given in \cite tsai1974pair [Eq. (3.38)-(3.41)]: the Thomas-Fermi model(TFM) of the atom
 * was adopted; the variables [Eq. (3.30) and (3.31)]
 *  \f[
 *   \begin{split}
 *   \gamma &=100 m_0c^2\frac{k}{E} \frac{Z^{-1/3}}{E'}
 *   \\ \epsilon &=100 m_0c^2\frac{k}{E} \frac{Z^{-2/3}}{E'}
 *   \end{split}
 *  \f]
 * with \f$E'=E_t-k\f$ were introduced since  \f$\varphi_1(\gamma), \varphi_2(\gamma), \psi_1(\epsilon), \psi_2(\epsilon)\f$
 * are universal i.e. \f$Z\f$ independent functions when the TFM is used; then  [Eqs.(3.14)-(3.17)] were computed
 * numerically by using the TFM model of atom; based on the numerical results, approximate analytical expressions were
 * obtained for \f$\varphi_1(\gamma), \varphi_2(\gamma)\f$ and \f$\psi_1(\epsilon), \psi_2(\epsilon)\f$ that reproduce
 * the numerical results with \f$0.5\f$\% \cite tsai1974pair [Eq. (3.38)-(3.41)]
 *  \f[
 *   \begin{split}
 *   \varphi_1(\gamma) = 20.863 -2\ln\left[ 1+(0.55846\gamma)^2  \right] - 4[ 1-0.6\exp(-0.9\gamma) - 0.4\exp(-1.5\gamma)]
 *   \end{split}
 *  \f]
 *  \f[
 *   \begin{split}
 *   \varphi_2(\gamma) = \varphi_1(\gamma) - \frac{2}{3(1+6.5\gamma+6\gamma^2)}
 *   \end{split}
 *  \f]
 *  \f[
 *   \begin{split}
 *   \psi_1(\epsilon) = 28.34 - 2\ln\left[ 1+(3.621\epsilon)^2  \right] - 4[ 1-0.7\exp(-8\epsilon) - 0.3\exp(-29.2\epsilon)]
 *   \end{split}
 *  \f]
 *  \f[
 *   \begin{split}
 *   \psi_2(\epsilon) = \psi_1(\epsilon) - \frac{2}{3(1+40\epsilon+400\epsilon^2)}
 *   \end{split}
 *  \f]
 * These expressions can be used for \f$Z\geq 5\f$ due to the limitation of the TFM model at low \f$Z\f$.
 * For \f$Z<5\f$ the <em>complete screening</em> i.e. \f$\gamma=\epsilon\approx 0\f$ form of the differential cross
 * section can be used \cite tsai1974pair (this will be valid at high energies \f$ > 10 \f$ [GeV] and when the detailed
 * shape of the bremsstrahlung spectrum at the high energy tip is not very important).
 * Under this conditions the differential cross section can be written \cite tsai1974pair [Eq. (3.83)]
 *  \f[
 *   \frac{\mathrm{d}\sigma}{\mathrm{d}k} = \frac{4\alpha r_e^2}{k}
 *   \left\{
 *   \left( \frac{4}{3} -\frac{4}{3}y+y^2 \right)
 *     \left[  Z^2\left(L_{el} -f \right) + ZL_{inel} \right] +
 *     \frac{1}{9}(1-y) \left[ Z^2+Z \right]
 *   \right\}
 *  \f]
 * that results in
 *  \f[
 *   \left( \frac{\mathrm{d}\sigma}{\mathrm{d}k} \right)^* =
 *     \left( 1 - y+\frac{3}{4}y^2 \right)
 *     \left[ \left(L_{el} -f \right) + \frac{1}{Z}L_{inel} \right] +
 *     \frac{1}{12}(1-y) \left[ 1+ \frac{1}{Z} \right]
 *  \f]
 * where
 *  \f[
 *   \begin{matrix}
 *    4L_{el}    & = & \left[\varphi_1(\gamma=0) -\frac{4}{3}\ln(Z)  \right]    & = & 4\left[ \ln \left(\exp(\varphi_1(\gamma=0)/4)Z^{-1/3}  \right)  \right] \\
 *    4L_{inel}  & = & \left[\psi_1(\epsilon=0)  -\frac{8}{3}\ln(Z)  \right]    & = & 4\left[  \ln \left( \exp(\psi_1(\epsilon=0)/4)Z^{-2/3} \right) \right]  \\
 *    \varphi_1(\gamma) - \varphi_2(\gamma)  & = & 2/3 && \\
 *    \psi_1(\epsilon) - \psi_2(\epsilon)    & =  & 2/3   &&
 *  \end{matrix}
 *  \f]
 * that in general gives \f$ 4L_{el} = 4\left[ \ln \left(\exp(20.863/4)Z^{-1/3}  \right)  \right] \f$ and
 * \f$ 4L_{inel} = 4\left[  \ln \left( \exp(\psi_1(\epsilon=0)/4)Z^{-2/3} \right) \right] \f$ when the above forms
 * of \f$ \varphi_1(\gamma),  \psi_1(\epsilon)\f$ (derived as approximations under the TFM) are used
 * i.e. \f$ L_{el} = \ln \left[ 184.1499 Z^{-1/3} \right] \f$ and  \f$ L_{inel} = \ln \left[ 1193.923 Z^{-2/3} \right] \f$.
 * Since the TFM model of atom is supposed to be good only for high \f$ Z \f$, Tsai \cite tsai1974pair [TABLE B.2.] gave the
 * best estimated values of \f$ L_{el} \f$ and \f$ L_{inel} \f$ for \f$ Z<5\f$
 *   <center>
 *   <table>
 *   <caption id="TsaiLowZ">Best estimated values for low Z elements </caption>
 *    <tr><th> Z                <th> 1     <th> 2     <th>  3     <th> 4     <th> 5     <th> 6     <th> 7
 *    <tr><td> L<sub>el</sub>   <td> 5.31  <td> 4.79  <td>  4.74  <td> 4.71  <td> 4.68  <td> 4.62  <td> 4.57
 *    <tr><td> L<sub>inel</sub> <td> 6.144 <td> 5.621 <td>  5.805 <td> 5.924 <td> 6.012 <td> 5.891 <td> 5.788
 *   </table>
 *   </center>
 * These values are used in the complete screening approximation to compute
 * \f$ \left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \f$ for \f$ Z<5 \f$.
 */
double RelativisticBremsModel::ComputeDXSecPerAtom(double gammaenergy, double totalenergy, double z) {
  //return zero if gamma energy is < 0
  if (gammaenergy<0.0) {
    return 0.0;
  }

  double y = gammaenergy/totalenergy; // reduced photon energy
  double firstTerm   = 0.0;
  double secondTerm  = 0.0;

  // Coulomb correction from Davis,Bethe,Maximom PRL 1954 Eqs.(36-38)
  // should be pulled out because it is constant for a given Z
  double mu  = z*geant::kFineStructConst;
  double mu2 = mu*mu;
  double mu4 = mu2*mu2;
  double mu6 = mu2*mu4;
  double coulombCorrection = mu2*(1.0/(1.0+mu2) + 0.20206 - 0.0369*mu2 + 0.0083*mu4 - 0.002*mu6);

  int izet = std::lrint(z);

  // 4/3Z^2 is pulled out in both cases so the multiplicative factor is (16 alpha r^2 Z^2)/(3*k)
  if (izet<5) {  // use complete screening (Tsai Eq.(3.83)) and L_el and L_inel from Tsai Table B2 for Z<5
    firstTerm  = (1.0-y+3./4.*y*y)*(gFelLowZet[izet]-coulombCorrection+gFinelLowZet[izet]/z);
    secondTerm = (1.0-y)/12.0*(1.0+1.0/z);
  } else { // Tsai: screening from Thomas-Fermi model of atom; Tsai Eq.(3.82)
    // variables gamma and epsilon from Tsai Eq.(3.30) and Eq.(3.31)
    // std::pow(...) should be pulled out because it is constant for a given Z
    double dum0    = 100.0*geant::kElectronMassC2*y/(totalenergy-gammaenergy);
    double gamma   = dum0/std::pow(z,1./3.);
    double epsilon = dum0/std::pow(z,2./3.);
    // compute TFM elastic and inelastic screening functions based on Tsai approximate expressions; Tsai Eqs.(3.38-3.41)
    double phi1   = 20.863 - 2.0*std::log(1.0+(0.55846*gamma)*(0.55846*gamma)) - 4.0*(1.0-0.6*std::exp(-0.9*gamma)-0.4*std::exp(-1.5*gamma));
    double phi1m2 = 2.0/(3.0+19.5*gamma+18.0*gamma*gamma);  // phi1-phi2
    double psi1   = 28.34  - 2.0*std::log(1.0+(3.621*epsilon)*(3.621*epsilon)) - 4.0*(1.0-0.7*std::exp(-8.0*epsilon)-0.3*std::exp(-29.2*epsilon));
    double psi1m2 = 2.0/(3.0+120.0*epsilon+1200.0*epsilon*epsilon); // psi1-psi2
    // std::log(z) should be pulled out because it is constant for a given Z
    firstTerm  = (1.0-y+3./4.*y*y)*( (0.25*phi1-1./3.*std::log(z)-coulombCorrection) + (0.25*psi1-2./3.*std::log(z))/z );
    secondTerm = 0.125*(1.0-y)*(phi1m2+psi1m2/z);
  }

  // Note: to get the differential cross section this must be multiplied by (16 alpha r^2 Z^2)/(3*k)
  // However, we will integrate by the transformed variable i.e. the (k * ds/dk) will be needed:
  //\sigma(Z,E|W_{c}) =  \int_{\ln(\kappa_c)}^{0}  k \frac{\mathrm{d}\sigma}{\mathrm{d}k}  \mathrm{d}u
  return firstTerm+secondTerm;
}









// Ratin Alias
void RelativisticBremsModel::InitSamplingTables1() {
  // set up the common electron energy grid
  if (fSamplingElecEnergies) {
    delete [] fSamplingElecEnergies;
    delete [] fLSamplingElecEnergies;
    fSamplingElecEnergies  = nullptr;
    fLSamplingElecEnergies = nullptr;
  }
  fSamplingElecEnergies = new double[fNumSamplingElecEnergies];
  fLSamplingElecEnergies = new double[fNumSamplingElecEnergies];
  double lMin  = std::log(fMinElecEnergy);
  fElEnLMin    = lMin;
  double lMax  = std::log(fMaxElecEnergy);
  double delta = (lMax-lMin)/(fNumSamplingElecEnergies-1.0);
  fElEnILDelta = 1.0/delta;
  double dumle = 0.0;
  for(int i = 0; i<fNumSamplingElecEnergies; ++i){
    double ddum = lMin+dumle;
    fLSamplingElecEnergies[i] = ddum;
    fSamplingElecEnergies[i]  = std::exp(lMin+dumle);
    std::cerr<<" E("<<i<<") = "<<fSamplingElecEnergies[i]/geant::MeV<<std::endl;
    dumle+=delta;
  }

  // - get number of different material-gammacut pairs
  // - allocate space and fill the global to local material-cut index map
  const std::vector<MaterialCuts*> theMaterialCutsTable = MaterialCuts::GetTheMaterialCutsTable();
  int numMaterialCuts = theMaterialCutsTable.size();
  if (fGlobalMatGCutIndxToLocal) {
    delete [] fGlobalMatGCutIndxToLocal;
    fGlobalMatGCutIndxToLocal = nullptr;
  }
  fGlobalMatGCutIndxToLocal = new int[numMaterialCuts];
  std::cerr<<" === Number of global Material-Cuts = "<<numMaterialCuts<<std::endl;
  // count diffenet material-gammacut pairs and set to global to local mat-cut index map
  int oldnumDif = fNumDifferentMaterialGCuts;
  int oldnumSEE = fNumSamplingElecEnergies;
  fNumDifferentMaterialGCuts = 0;
  for (int i=0; i<numMaterialCuts; ++i) {
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
      fGlobalMatGCutIndxToLocal[i] = j;
    }
  }
  std::cerr<<" === Number of local Material-Cuts = "<<fNumDifferentMaterialGCuts<<std::endl;
  // allocate space for the matrial-gcut sampling tables and init these pointers to null
  if (fRatinAliasData) {
    for (int i=0; i<oldnumDif; ++i) {
      for (int j=0; j<oldnumSEE; ++j) {
        int indx = i*oldnumSEE+j;
        if (fRatinAliasData[indx]) {
          delete [] fRatinAliasData[indx]->fXdata;
          delete [] fRatinAliasData[indx]->fAliasW;
          delete [] fRatinAliasData[indx]->fC;
          delete [] fRatinAliasData[indx]->fA;
          delete [] fRatinAliasData[indx]->fB;
          delete [] fRatinAliasData[indx]->fAliasIndx;
          delete fRatinAliasData[indx];
        }
      }
    }
    delete [] fRatinAliasData;
  }

  int *isdone = new int[fNumDifferentMaterialGCuts]();
  int  idum   = fNumDifferentMaterialGCuts*fNumSamplingElecEnergies;
  fRatinAliasData = new RatinAlias*[idum];
  for (int i=0; i<idum; ++i)
    fRatinAliasData[i] = nullptr;

  for (int i=0; i<numMaterialCuts; ++i) {
    std::cerr<<"   See if Material-Gcut ==> " <<theMaterialCutsTable[i]->GetMaterial()->GetName()<<"  gCut = "<< theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]<<std::endl;
    int localindx = fGlobalMatGCutIndxToLocal[i];
    int ialias    = localindx*fNumSamplingElecEnergies;
    if (!isdone[localindx]) { // init for this material-gammacut pair is it has not been done yet i.e. the last elec-energy bin sampling data is still null.
       std::cerr<<"   -> Will init for Material-Gcut ==> " <<theMaterialCutsTable[i]->GetMaterial()->GetName()<<"  gCut = "<< theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]<<std::endl;
       BuildOneRatinAlias1(ialias, theMaterialCutsTable[i]->GetMaterial(), theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0]);
       isdone[localindx] = 1;
    }
  }
  delete [] isdone;
  // test
  for (int i=0; i<numMaterialCuts; ++i)
    std::cerr<<"     --> Global MatCut-indx = "<< i << " local indx = "<<fGlobalMatGCutIndxToLocal[i] <<std::endl;

//
//  for (int i=0; i<fNumSamplingPhotEnergies; ++i) {
//    for (int j=0; j<fNumSamplingElecEnergies; ++j)
//      std::cout<< (datx[j])[i] <<"  "<<(daty[j])[i] <<" ";
//    std::cout<<std::endl;
//  }

}


void RelativisticBremsModel::BuildOneRatinAlias1(int ialias, const Material *mat, double gcut) {
  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  // material dependent constant
  double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol()*4.0*geant::kPi
                        *geant::kClassicElectronRadius*geant::kRedElectronComptonWLenght *geant::kRedElectronComptonWLenght;
  // target i.e. Z dependent constant
  constexpr double lpmConstant = 0.25*geant::kFineStructConst*geant::kElectronMassC2*geant::kElectronMassC2/(geant::kPi*geant::kHBarPlanckCLight);
  // material dependent constant
  double lpmEnergy   = mat->GetMaterialProperties()->GetRadiationLength()*lpmConstant;
  // material dependent constant: the e-/e+ energy at which the corresponding k_LPM=k_p
  double energyThLPM = std::sqrt(densityFactor)*lpmEnergy;

  // we will need the element composition of this material
  const Vector_t<Element*> theElements    = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int   numElems = theElements.size();

  // pdf is needed only locally for the preparation
  double *thePDF      = new double[fNumSamplingPhotEnergies]();
  // create a local alias data structure with maximum number of points
  RatinAlias *localRA = new RatinAlias();
  localRA->fXdata     = new double[fNumSamplingPhotEnergies]();
  localRA->fAliasW    = new double[fNumSamplingPhotEnergies]();
  localRA->fC         = new double[fNumSamplingPhotEnergies]();
  localRA->fA         = new double[fNumSamplingPhotEnergies]();
  localRA->fB         = new double[fNumSamplingPhotEnergies]();
  localRA->fAliasIndx = new    int[fNumSamplingPhotEnergies]();
  // a local AliasTable helper
  AliasTable *alst = new AliasTable();

  for (int ieener=0; ieener<fNumSamplingElecEnergies; ++ieener) {
    double eener = fSamplingElecEnergies[ieener];
    if (eener>gcut) { // otherwise no gamma production at that e- energy so let the sampling table to be nullptr
      // particle is always e-/e+
      double totalenergy = eener+geant::kElectronMassC2;
      double densityCor  = densityFactor*totalenergy*totalenergy; // this is k_p^2
      double dum0  = gcut*gcut+densityCor;
      double dum1  = (eener*eener+densityCor)/dum0;
      double ldum1 = std::log(dum1);

      // fill 3 initial values
      int numdata = 3;
      localRA->fXdata[0] = 0;
      localRA->fXdata[1] = 0.5;
      localRA->fXdata[2] = 1.0;

      //
      // expand the data up to maximum such that the error between the Rita interpolated and real PDF is minimised
      //
      double isDone = false;
      while(!isDone && (numdata<fNumSamplingPhotEnergies)) {
        // first get the distribution at the current x-points
        for (int i=0; i<numdata; ++i) {
          double egamma = gcut;
          if (i==0) {
            egamma = gcut;
          } else if (i==numdata-1) {
            egamma = eener;
          } else {
            egamma = std::sqrt(dum0*std::exp(localRA->fXdata[i]*ldum1)-densityCor);
          }
          double dcross = 0.0;
          for (int ielem=0; ielem<numElems; ++ielem) {
            double zet = theElements[ielem]->GetZ();
            double val = 0.0;
            if (totalenergy>energyThLPM)
              val = ComputeURelDXSecPerAtom(egamma, totalenergy, zet, lpmEnergy, densityCor);
            else
              val = ComputeDXSecPerAtom(egamma, totalenergy, zet);
            dcross += theAtomicNumDensityVector[ielem]*zet*zet*val;
          }
          thePDF[i] = dcross;
//          std::cerr << localRA->fXdata[i] << " "<<dcross<<std::endl;
        }
        // - set up a Rita table for the current pdf data (thePDF will be normaised so we need to recompute each time)
        double normFactor = alst->PreparRatinForPDF(localRA->fXdata, thePDF, localRA->fC, localRA->fA, localRA->fB, numdata);

        // - find the lower index of the bin, where we have the biggest interp. error compared to computed value
        //   i.e. compute the integral of |p(x)-phat(x)| at each interval
        double maxerr     = 0.0; // value of the current maximum error
        int    maxerrindx = 0;   // the lower index of the corresponding bin
        for (int i=0; i<numdata-1; ++i) {
          double xx = 0.5*(localRA->fXdata[i]+localRA->fXdata[i+1]);// current midpoint
          if (xx-localRA->fXdata[i]<1.e-6)
            continue;
          double  u = dum0*std::exp(xx*ldum1);
          double egamma = std::sqrt(u-densityCor);
          double dcross = 0.0;
          for (int ielem=0; ielem<numElems; ++ielem) {
            double zet = theElements[ielem]->GetZ();
            double val = 0.0;
            if (totalenergy>energyThLPM)
              val = ComputeURelDXSecPerAtom(egamma, totalenergy, zet, lpmEnergy, densityCor);
            else
              val = ComputeDXSecPerAtom(egamma, totalenergy, zet);
            dcross += theAtomicNumDensityVector[ielem]*zet*zet*val;
          }
          double curPx    = dcross*normFactor;
          double curAprPx = alst->GetRatinForPDF1(xx, localRA->fXdata, localRA->fC, localRA->fA, localRA->fB, i);
          double curError = std::fabs(1.0-curPx/curAprPx);

          if (curError>maxerr) {
            maxerr     = curError;
            maxerrindx = i;
          }
        }

        // check stopping
        if (maxerr<1.e-3)
          isDone = true;

//if (eener>1.0)
 //if(isDone)
 //std::cerr<<"Eener = "<<eener<<" numdata = "<<numdata << " maxerror = "<<maxerr  << "  xmin = "<< localRA->fXdata[maxerrindx]<<std::endl;
        // expand x,y data by puting an additional value at the mid point of the highest error bin
        // -first shift all values to the right
        for (int j=numdata; j>maxerrindx+1; --j) {
          localRA->fXdata[j] = localRA->fXdata[j-1];
        }
        // fill x mid point
        localRA->fXdata[maxerrindx+1] = 0.5*(localRA->fXdata[maxerrindx]+localRA->fXdata[maxerrindx+1]);
        // increase number of data
        ++numdata;
      }

      // compute the PDF at the final x-points
      for (int i=0; i<numdata; ++i) {
        double egamma = gcut;
        if (i==0) {
          egamma = gcut;
        } else if (i==numdata-1) {
          egamma = eener;
        } else {
          egamma = std::sqrt(dum0*std::exp(localRA->fXdata[i]*ldum1)-densityCor);
        }
        double dcross = 0.0;
        for (int ielem=0; ielem<numElems; ++ielem) {
          double zet = theElements[ielem]->GetZ();
          double val = 0.0;
          if (totalenergy>energyThLPM)
            val = ComputeURelDXSecPerAtom(egamma, totalenergy, zet, lpmEnergy, densityCor);
          else
            val = ComputeDXSecPerAtom(egamma, totalenergy, zet);
          dcross += theAtomicNumDensityVector[ielem]*zet*zet*val;
        }
        thePDF[i] = dcross;
      }

      // create the final alias struct
      fRatinAliasData[ialias]             = new RatinAlias();
      fRatinAliasData[ialias]->fNumdata   = numdata;
      fRatinAliasData[ialias]->fXdata     = new double[numdata]();
      fRatinAliasData[ialias]->fAliasW    = new double[numdata]();
      fRatinAliasData[ialias]->fC         = new double[numdata]();
      fRatinAliasData[ialias]->fA         = new double[numdata]();
      fRatinAliasData[ialias]->fB         = new double[numdata]();
      fRatinAliasData[ialias]->fAliasIndx = new    int[numdata]();
      for (int i=0;i<numdata;++i)
        fRatinAliasData[ialias]->fXdata[i] = localRA->fXdata[i];

      // set up the rational function aproxiamtion based alias sampler
      alst->PreparRatinTable(fRatinAliasData[ialias]->fXdata,thePDF, fRatinAliasData[ialias]->fC,
                            fRatinAliasData[ialias]->fA,fRatinAliasData[ialias]->fB,
                            fRatinAliasData[ialias]->fAliasW, fRatinAliasData[ialias]->fAliasIndx,
                            fRatinAliasData[ialias]->fNumdata);
      std::cerr<<"  -- Material = "<<mat->GetName()<<  "  Eel = " <<eener << " [GeV] " <<" #data = "<<numdata <<std::endl;
/*
if (eener>9.8e+4 || 1) {
//double nnorm =    alst->PreparRatinForPDF(fRatinAliasData[ialias]->fXdata, thePDF, fRatinAliasData[ialias]->fC,
//                           fRatinAliasData[ialias]->fA, fRatinAliasData[ialias]->fB, numdata);

  for (int i=0; i<1000001; ++i) {
    double vval = alst->GetRatinForPDF(i*1.e-6, fRatinAliasData[ialias]->fXdata, fRatinAliasData[ialias]->fC, fRatinAliasData[ialias]->fA,
                                           fRatinAliasData[ialias]->fB, numdata);
// for (int i=0; i<numdata; ++i) {
   double gg = std::sqrt( (gcut*gcut+densityCor)*std::exp(i*1.e-6*std::log((eener*eener+densityCor)/(gcut*gcut+densityCor) )) -densityCor);
//double gg = std::sqrt( (gcut*gcut+densityCor)*std::exp(fRatinAliasData[ialias]->fXdata[i]*std::log((eener*eener+densityCor)/(gcut*gcut+densityCor) )) -densityCor);


//   std::cout<<i*1.e-7<<" "<<vval /(1.0+densityCor/(gg*gg))<<std::endl;

//std::cout<<i*1.e-6<<" "<<vval <<std::endl;
std::cout<<std::log10(gg/eener)<<" "<<vval <<std::endl;
//std::cout<<std::log10(gg/eener)<<" "<<thePDF[i] <<std::endl;

 }

//for (int i=0; i<numdata; ++i)
//   std::cout<<i<<" "<<fRatinAliasData[ialias]->fXdata[i]<< " "<<thePDF[i]/nnorm<<std::endl;

exit(1);
}
*/

    }
    ++ialias;
  }
  delete [] thePDF;
  delete [] localRA->fXdata;
  delete [] localRA->fAliasW;
  delete [] localRA->fC;
  delete [] localRA->fA;
  delete [] localRA->fB;
  delete [] localRA->fAliasIndx;
  delete localRA;
  delete alst;
}



} //namespace geantphysics
