
#include "MollerBhabhaIonizationModel.h"

#include "PhysicalConstants.h"
#include "Material.h"
#include "Element.h"
#include "MaterialProperties.h"

#include "MaterialCuts.h"
#include "AliasTable.h"

#include "Electron.h"

#include "LightTrack.h"
#include "PhysicsData.h"

// from geantV
#include "GeantTaskData.h"

#include <iostream>
#include <cmath>

namespace geantphysics {

MollerBhabhaIonizationModel::MollerBhabhaIonizationModel(bool iselectron, const std::string &modelname)
: EMModel(modelname), fIsElectron(iselectron) {
  fSecondaryInternalCode         =  -1;

  fSTNumPrimaryEnergyPerDecade   =   8;   // ST=> sampling tables
  fSTNumSamplingElecEnergies     = 101;   // ST=> sampling tables

  fAliasSampler                  = nullptr;
}


MollerBhabhaIonizationModel::~MollerBhabhaIonizationModel() {
  if (GetUseSamplingTables()) {
    ClearSamplingTables();
  }
  if (fAliasSampler) {
    delete fAliasSampler;
  }
}


void MollerBhabhaIonizationModel::Initialize() {
  EMModel::Initialize();
  if (GetUseSamplingTables()) {
    InitSamplingTables();
  }
  fSecondaryInternalCode = Electron::Definition()->GetInternalCode();
}


double MollerBhabhaIonizationModel::ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle*, bool istotal){
  const Material *mat =  matcut->GetMaterial();
  double elCut        =  matcut->GetProductionCutsInEnergy()[1];
  if (istotal) {
    // normaly 2xEkin would be fine to get the total stopping power due to indistinguishability of the 2 final state
    // electrons but due to the low energy threshold in the Moller model kinenergy will make sure to be higher than a
    // threshold energy so we just set it to 1000xEkin
    elCut = 1000.0*kinenergy;
  }
  return ComputeDEDXPerVolume(mat, elCut, kinenergy);
}


double MollerBhabhaIonizationModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle*) {
  double xsec = 0.0;
  if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
    return xsec;
  }
  const Material *mat =  matcut->GetMaterial();
  const double elCut  =  matcut->GetProductionCutsInEnergy()[1];
  xsec = ComputeXSectionPerVolume(mat, elCut, kinenergy);
  return xsec;
}


double MollerBhabhaIonizationModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut,
                                                           double kinenergy, const Particle*) {
  double xsec = 0.0;
  if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
    return xsec;
  }
  const Material *mat =  matcut->GetMaterial();
  const double elCut  =  matcut->GetProductionCutsInEnergy()[1];
  xsec = ComputeXSectionPerAtom(elem, mat, elCut, kinenergy);
  return xsec;
}


int MollerBhabhaIonizationModel::SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td) {
  int    numSecondaries      = 0;
  const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(track.GetMaterialCutCoupleIndex());
  const double electronCut   = matCut->GetProductionCutsInEnergy()[1];
  const double ekin          = track.GetKinE();
  //
  double maxETransfer = ekin;
  if (fIsElectron) {
    maxETransfer *= 0.5;
  }
  // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
  // check if kinetic energy is above fHighEnergyUsageLimit andd o nothing if yes
  // check if max energy transfer is below e- production cut and do nothing if yes
  if (ekin<GetLowEnergyUsageLimit() || ekin>GetHighEnergyUsageLimit() || maxETransfer<=electronCut) {
    return numSecondaries;
  }
  //
  // sample energy transfer and compute direction
  const double elInitTotalEnergy   = ekin+geant::kElectronMassC2;  // initial total energy of the e-/e+
  const double elInitTotalMomentum = std::sqrt(ekin*(elInitTotalEnergy+geant::kElectronMassC2));
  double deltaKinEnergy = 0.;
  if (GetUseSamplingTables()) {
    double *rndArray = td->fDblArray;
    td->fRndm->uniform_array(3, rndArray);
    deltaKinEnergy = SampleEnergyTransfer(matCut, ekin, rndArray[0], rndArray[1], rndArray[2]);
  } else {
    deltaKinEnergy = SampleEnergyTransfer(matCut, ekin, td);
  }
  const double deltaTotalMomentum  = std::sqrt(deltaKinEnergy*(deltaKinEnergy+2.0*geant::kElectronMassC2));
  const double cost                = deltaKinEnergy*(elInitTotalEnergy+geant::kElectronMassC2)
                                    /(deltaTotalMomentum*elInitTotalMomentum);
  // check cosTheta limit
  const double cosTheta = std::min(cost,1.0);
  const double sinTheta = std::sqrt((1.0-cosTheta)*(1.0+cosTheta));
  const double phi       = geant::kTwoPi*td->fRndm->uniform();
  // direction of the delta e- in the scattering frame
  double deltaDirX = sinTheta*std::cos(phi);
  double deltaDirY = sinTheta*std::sin(phi);
  double deltaDirZ = cosTheta;
  // rotate back to lab frame
  RotateToLabFrame(deltaDirX, deltaDirY, deltaDirZ, track.GetDirX(), track.GetDirY(), track.GetDirZ());
  // create the secondary partcile i.e. the delta e-
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
  std::vector<LightTrack>& sectracks = td->fPhysicsData->GetListOfSecondaries();
  sectracks[secIndx].SetDirX(deltaDirX);
  sectracks[secIndx].SetDirY(deltaDirY);
  sectracks[secIndx].SetDirZ(deltaDirZ);
  sectracks[secIndx].SetKinE(deltaKinEnergy);
  sectracks[secIndx].SetGVcode(fSecondaryInternalCode);
  sectracks[secIndx].SetMass(geant::kElectronMassC2);
  sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); // parent GeantTrack index
  //
  // compute the primary e-/e+ post interaction kinetic energy and direction: from momentum vector conservation
  // final momentum of the primary e-/e+ in the lab frame
  double elDirX = elInitTotalMomentum*track.GetDirX() - deltaTotalMomentum*deltaDirX;
  double elDirY = elInitTotalMomentum*track.GetDirY() - deltaTotalMomentum*deltaDirY;
  double elDirZ = elInitTotalMomentum*track.GetDirZ() - deltaTotalMomentum*deltaDirZ;
  // normalisation
  const double norm  = 1.0/std::sqrt(elDirX*elDirX + elDirY*elDirY + elDirZ*elDirZ);
  // update primary track direction
  track.SetDirX(elDirX*norm);
  track.SetDirY(elDirY*norm);
  track.SetDirZ(elDirZ*norm);
  // update primary track kinetic energy
  track.SetKinE(ekin-deltaKinEnergy);
  //
  // return with number of secondaries i.e. 1 e-
  return numSecondaries;
}


double MollerBhabhaIonizationModel::MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle*) const {
  double mine = (matcut->GetProductionCutsInEnergy())[1]; // e- production cut in the material-cut
  // in case of e- e- interaction the minimum primary energy is 2x the e- production cut
  if (fIsElectron) {
    mine += mine;
  }
  return mine;
}





/**
 * Ionization part of the (restricted) dE/dx computed based on the formula given by Berger and Seltzer
 * \cite berger1964tables \cite crawford1970electron.
 * \f[
 *   \frac{\mathrm{d}E}{\mathrm{d}x} = \frac{2\pi r_{e}^{2} m_0c^2 n_{el}}{\beta^2}
 *   \left[ \ln \frac{2(\tau+2)}{(Im_0c^2)^2} + G^{\pm}(\tau,\tau_{up}) - \delta
 *   \right]
 * \f]
 * where
 * \f[
 *  \begin{array}{lcl}
 *  r_e       &\to&  \textrm{classical electron radius} \\
 *  m_0c^2    &\to&  \textrm{electron rest mass energy} \\
 *  n_{el}    &\to&  \textrm{electron density of the material} \\
 *  I         &\to&  \textrm{mean excitation energy of the material}(^*) \\
 *  \gamma    &\to&  \textrm{Lorentz factor} =E_t/(m_0c^2) \textrm{ i.e. total energy of the incident particle (e-/e+)
 *                    in rest mass energy unit} \\
 *  \beta     &\to&  \textrm{ratio of the speed of the incident particle to } c;\quad = P_t/E_t \textrm{ with } P_t=pc
 *                   \textrm{ i.e. total momentum of the incident particle (e-/e+) in rest mass energy unit} \\
 *  \tau      &\to&  \textrm{kinetic energy of the incident particle(e-/e+) in rest mass energy unit} \\
 *  \tau_{up} &\to&  \textrm{restricted maximum energy transfer in incident particle rest mass energy unit i.e. }
 *                    \min[\tau_{max},\tau_{c}] \textrm{ where} \\
 *            &&       \tau_{max} \to \textrm{maximum possible energy transfer in incident particle (e-/e+) rest mass energy unit i.e.}\\
 *            &&       \quad \tau_{max} = \begin{cases}
 *                                         \tau           & \quad \textrm{incident particle is e+} \\
 *                                         \frac{\tau}{2} & \quad \textrm{incident particle is e-}
 *                                        \end{cases}\\
 *            &&       \tau_{c} \to \textrm{e- production cut (kinetic) energy in rest mass energy unit}\\
 *  \delta    &\to&  \textrm{density effect correction}(^*) \\
 *  G^{\pm}(\tau,\tau_{up}) &\to& \textrm{represents the term that the integration for the Moller } [e^-+e^-\to e^-+e^-]
 *                                \textrm{ and Bhabha }[e^++e^-\to e^++e^-] \textrm{ differential cross sections yield
 *                                different results i.e.}\\
 *            &&      \textrm{for e-:}\\
 *            &&      \quad G^{-}(\tau,\tau_{up}) =-1-\beta^2+\ln[(\tau-\tau_{up})\tau_{up}]+\frac{\tau}{\tau-\tau_{up}}
 *                                                 +\frac{1}{\gamma^2}
 *                                     \left[
 *                                      \frac{\tau_{up}^2}{2} + (2\tau+1)\ln\left(1-\frac{\tau_{up}}{\tau}\right)
 *                                     \right] \\
 *            &&      \textrm{for e+:}\\
 *            &&      \quad G^{+}(\tau,\tau_{up}) = \ln(\tau\tau_{up})-\frac{\beta^2}{\tau}
 *                          \left[
 *                             \tau+2\tau_{up}-y\frac{3\tau_{up}^2}{2}-y^2\left( \tau_{up}-\frac{\tau_{up}^3}{3} \right)
 *                             -y^3\left( \frac{\tau_{up}^2}{2}-\tau\frac{\tau_{up}^3}{3}+\frac{\tau_{up}^4}{4} \right)
 *                          \right]\\
 *           &&       \quad \textrm{with } y=1/(\gamma+1)
 *  \end{array}
 * \f]
 *
 * (\f$^*\f$see more about \f$I\f$ and \f$\delta\f$ at MaterialProperties)
 */
double MollerBhabhaIonizationModel::ComputeDEDXPerVolume(const Material *mat, const double pcutenergy, const double primekin) {
  constexpr double factor     = geant::kTwoPi*geant::kClassicElectronRadius*geant::kClassicElectronRadius*geant::kElectronMassC2;
  const double twolog10inv    = 1.0/(2.0*std::log(10.0));
  // get the material properties
  MaterialProperties *matProp = mat->GetMaterialProperties();
  // get the electron denisty of the material
  const double elDensity      = matProp->GetTotalNumOfElectronsPerVol();
  // effective atomic number for the kinetic energy threshold computation
  const double effZ           = matProp->GetTotalNumOfElectronsPerVol()/matProp->GetTotalNumOfAtomsPerVol();
  // compute the kinetic energy threshold
  const double kineTh         = 0.25*std::sqrt(effZ)*geant::keV;
  // get the mean excitation energy
  double meanExcEnergy  = matProp->GetMeanExcitationEnergy();
  // set kinetic energy
  double kineticEnergy = primekin;
  if (kineticEnergy<kineTh) {
    kineticEnergy = kineTh;
  }
  // set other parameters
  const double tau       = kineticEnergy/geant::kElectronMassC2; // i.e. E_kin in (mc^2) units
  const double gamma     = tau + 1.0;       // = E_t/(mc^2)  i.e. E_t in (mc^2) units
  const double gamma2    = gamma*gamma;     // \gamma^2 = [E_t/(mc^2)]^2
  const double betagama2 = tau*(tau+2.0);   // = (\beta * \gamma)^2 = (P_t/(mc^2))^2 i.e. [P_t in (mc^2) units]^2
  const double beta2     = betagama2/gamma2;// \beta2 i.e. [P_t/E_t]^2
  meanExcEnergy   /= geant::kElectronMassC2;  // mean excitation energy in (mc^2) units
  meanExcEnergy   *= meanExcEnergy;
  // set maximum kinetic energy (in mc^2 units) that can be transformed to a free electron
  double taumax    = tau;  // for e+ : E_kin is the maximum kinetic energy transfer
  if (fIsElectron) {       // for e- : E_kin/2 is the maximum kinetic energy transfer
    taumax *= 0.5;
  }
  // set upper limit of tau: upper limit of the integral corresponding to the continuous part i.e.
  // min(e-/e+ production cut energy, maximum kinetic energy transfer) in mc^2 units
  double tauUpLim  = pcutenergy/geant::kElectronMassC2;
  if (tauUpLim>taumax) {
    tauUpLim = taumax;
  }
  //
  // compute dE/dx
  double dedx = 0.0;
  // first compute G^{-/+} that is different for e- (Moller) and e+ (Bhabha)
  if (fIsElectron) {
    dedx = std::log((tau-tauUpLim)*tauUpLim) + tau/(tau-tauUpLim)
           + (0.5*tauUpLim*tauUpLim + (2.0*tau + 1.0)*std::log(1.0-tauUpLim/tau))/gamma2 - 1.0 - beta2;
  } else {
    const double tauUpLim2 = tauUpLim*tauUpLim*0.5;   // \tau_up^2/2
    const double tauUpLim3 = tauUpLim2*tauUpLim/1.5;  // \tau_up^3/3
    const double tauUpLim4 = tauUpLim3*tauUpLim*0.75; // \tau_up^4/4
    const double y         = 1.0/(1.0 + gamma);
    dedx = std::log(tau*tauUpLim) - beta2*(tau + 2.0*tauUpLim - y*(3.0*tauUpLim2 + y*(tauUpLim - tauUpLim3
           + y*(tauUpLim2 - tau*tauUpLim3 + tauUpLim4))))/tau;
  }
  // add the common term
  dedx        += std::log(2.0*(tau + 2.0)/meanExcEnergy);
  // get the density effect correction term
  double dumx  = std::log(betagama2)*twolog10inv;
  dedx        -= matProp->GetDensityEffectFunctionValue(dumx);
  // apply the multiplicative factor to get the final dE/dx
  dedx        *= factor*elDensity/beta2;
  // very low energy extrapolation
  if (primekin<kineTh) {
    dumx = primekin/kineTh;
    if (dumx>0.25) {
      dedx /= std::sqrt(dumx);
    } else {
      dedx *= 1.4*std::sqrt(dumx)/(dumx+0.1);
    }
  }
  return std::max(dedx,0.0);
}


double MollerBhabhaIonizationModel::ComputeXSectionPerAtom(const Element *elem, const Material* /*mat*/,
                                                           const double pcutenergy, const double primekin) {
  const double xsec = elem->GetZ()*ComputeXSectionPerElectron(pcutenergy, primekin);
  return std::max(xsec,0.0);
}


double MollerBhabhaIonizationModel::ComputeXSectionPerVolume(const Material *mat, const double pcutenergy,
                                                             const double primekin) {
  const double elDensity = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol();
  const double xsec      = elDensity*ComputeXSectionPerElectron(pcutenergy, primekin);
  return std::max(xsec,0.0);
}


/**
  *  The sampling is based on the sampling tables prepared at initialisation. Statistical interpolation is used to
  *  select one of the primary particle kinetic energy grid points out of \f$ E_i \leq E_{kin} < E_{i+1}\f$ (linear
  *  interpolation in log kinetic energy) at a given primary particle kinetic energy \f$E_{kin}\f$. Then the transformed
  *  variable \f$\xi\in[0,1]\f$ is sampled from the sampling table (prepared at initialisation) that belongs to the
  *  selected primary particle kinetic energy grid point. The kinetic energy transfered to the electron \f$T\f$ then is
  *  obtained by applying the following transformation:
  *  \f[
  *     T =
  *     \begin{cases}
  *      T_{cut}^{e-}e^{\xi\ln(0.5E_{kin}/T_{cut}^{e-})} & \textrm{in case of Moller scattering }[e^-+e^-\to e^-+e^-]\\
  *      T_{cut}^{e-}e^{\xi\ln(E_{kin}/T_{cut}^{e-})}    & \textrm{in case of Bhabha scattering }[e^++e^-\to e^++e^-]
  *     \end{cases}
  *  \f]
  *  where \f$E_{kin}\f$ is the current primary particle (i.e. e-/e+) kinetic energy and \f$T_{cut}^{e-}\f$ is the
  *  current electron kinetic energy production threshold.
  */
double MollerBhabhaIonizationModel::SampleEnergyTransfer(const MaterialCuts *matcut, const double primekin,
                                                         const double r1, const double r2, const double r3) {
  const double elProdCut = matcut->GetProductionCutsInEnergy()[1]; // e- production cut
  const int  mcIndxLocal = fGlobalMatECutIndxToLocal[matcut->GetIndex()];
  // determine electron energy lower grid point
  const double lPrimEkin = std::log(primekin);
  //
  int indxPrimEkin = fSamplingTables[mcIndxLocal]->fNData-1;
  if (primekin<GetHighEnergyUsageLimit()) {
    const double val       = (lPrimEkin-fSamplingTables[mcIndxLocal]->fLogEmin)*fSamplingTables[mcIndxLocal]->fILDelta;
    indxPrimEkin           = (int)val;  // lower primary electron/prositron energy bin index
    const double pIndxHigh = val-indxPrimEkin;
    if (r1<pIndxHigh)
      ++indxPrimEkin;
  }
  // sample the transformed variable
  const LinAlias *als = fSamplingTables[mcIndxLocal]->fAliasData[indxPrimEkin];
  const double xi     = fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]),
                                                    &(als->fAliasIndx[0]), fSTNumSamplingElecEnergies, r2, r3);
  double dum1 = std::log(primekin/elProdCut);
  if (fIsElectron) {
    dum1 -= 0.693147180559945; // dum1 = dum1 + log(0.5)
  }
  // return with the sampled kinetic energy transfered to the electron
  return std::exp(xi*dum1)*elProdCut;
}


double MollerBhabhaIonizationModel::SampleEnergyTransfer(const MaterialCuts *matcut, const double primekin,
                                                        const Geant::GeantTaskData* td) {
  const double tmin   = matcut->GetProductionCutsInEnergy()[1];
  const double tmax   = (fIsElectron) ? (0.5*primekin) : (primekin);
  const double xmin   = tmin/primekin;
  const double xmax   = tmax/primekin;
//  const double tau    = primekin/geant::kElectronMassC2;
  const double gamma  = primekin/geant::kElectronMassC2 + 1.0;
  const double gamma2 = gamma*gamma;
  const double beta2  = 1.-1./gamma2;
  //
  const double xminmax = xmin*xmax;
  //
  double dum;
  double deltaEkin = 0.;
  double *rndArray = td->fDblArray;
  if (fIsElectron) { //Moller (e-e-) scattering
    const double gg = (2.0*gamma-1.0)/gamma2;
    const double  y = 1.-xmax;
    const double gf = 1.0-gg*xmax+xmax*xmax*(1.0-gg+(1.0-gg*y)/(y*y));
    do {
      td->fRndm->uniform_array(2, rndArray);
      deltaEkin = xminmax/(xmin*(1.0 -rndArray[0])+xmax*rndArray[0]);
      const double xx = 1.0-deltaEkin;
      dum       = 1.0 - gg*deltaEkin + deltaEkin*deltaEkin*(1.0-gg+(1.0-gg*xx)/(xx*xx));
    } while (gf*rndArray[1]>dum);
  } else  {          //Bhabha (e+e-) scattering
    const double y     = 1.0/(1.0+gamma);
    const double y2    = y*y;
    const double y12   = 1.0-2.0*y;
    const double b1    = 2.0-y2;
    const double b2    = y12*(3.0+y2);
    const double y122  = y12*y12;
    const double b4    = y122*y12;
    const double b3    = b4+y122;
    const double xmax2 = xmax*xmax;
    const double gf    = 1.0 + (xmax2*b4 - xmin*xmin*xmin*b3 + xmax2*b2 - xmin*b1)*beta2;
    do {
      td->fRndm->uniform_array(2, rndArray);
      deltaEkin = xminmax/(xmin*(1.0 -rndArray[0])+xmax*rndArray[0]);
      const double xx = deltaEkin*deltaEkin;
      dum       = 1.0 + (xx*xx*b4 - deltaEkin*xx*b3 + xx*b2 - deltaEkin*b1)*beta2;
    } while (gf*rndArray[1]>dum);
  }
  deltaEkin *= primekin;
  return deltaEkin;
}


void MollerBhabhaIonizationModel::ClearSamplingTables() {
  size_t numST = fSamplingTables.size();
  for (size_t i=0; i<numST; ++i) {
    AliasDataMaterialCuts* st = fSamplingTables[i];
    if (st) {
      size_t numAT = st->fAliasData.size();
      for (size_t j=0; j<numAT; ++j) {
        LinAlias* la = st->fAliasData[j];
        if (la) {
          la->fXdata.clear();
          la->fYdata.clear();
          la->fAliasW.clear();
          la->fAliasIndx.clear();
          delete la;
        }
      }
      st->fAliasData.clear();
      delete st;
    }
  }
  fSamplingTables.clear();
}


void MollerBhabhaIonizationModel::InitSamplingTables() {
  // clear all sampling tables (if any)
  ClearSamplingTables();
  // determine global-to-local matcut indices:
  // - get number of different electron cuts (pdf do not depend on material nor on Z)
  // - allocate space and fill the global to local material-cut index map
  const std::vector<MaterialCuts*> &theMaterialCutsTable = MaterialCuts::GetTheMaterialCutsTable();
  int numMaterialCuts      = theMaterialCutsTable.size();
  int numDifferentMatECuts = 0;
  fGlobalMatECutIndxToLocal.resize(numMaterialCuts,-2);
  for (int i=0; i<numMaterialCuts; ++i) {
    // if the current MaterialCuts does not belong to the current active regions
    if (!IsActiveRegion(theMaterialCutsTable[i]->GetRegionIndex())) {
      continue;
    }
    bool isnew = true;
    int j = 0;
    for (; j<numDifferentMatECuts; ++j) {
      if (theMaterialCutsTable[i]->GetProductionCutsInEnergy()[1]==theMaterialCutsTable[j]->GetProductionCutsInEnergy()[1]) {
        isnew = false;
        break;
      }
    }
    if (isnew) {
     fGlobalMatECutIndxToLocal[i] = numDifferentMatECuts;
     ++numDifferentMatECuts;
    } else {
      fGlobalMatECutIndxToLocal[i] = fGlobalMatECutIndxToLocal[j];
    }
  }
  fSamplingTables.resize(numDifferentMatECuts,nullptr);
  // create an AliasTable object
  if (fAliasSampler) {
    delete fAliasSampler;
  }
  fAliasSampler = new AliasTable();
  // build tables per material-cuts
  for (int i=0; i<numMaterialCuts; ++i) {
    const MaterialCuts *matCut = theMaterialCutsTable[i];
    int indxLocal = fGlobalMatECutIndxToLocal[i];
    if (indxLocal>-1 && !(fSamplingTables[indxLocal])) {
      BuildSamplingTableForMaterialCut(matCut, indxLocal);
    }
  }
}


void MollerBhabhaIonizationModel::BuildSamplingTableForMaterialCut(const MaterialCuts *matcut, int indxlocal) {
  const double ecut   = (matcut->GetProductionCutsInEnergy())[1];
  const double minPrimEnergy = std::max(MinimumPrimaryEnergy(matcut,nullptr),GetLowEnergyUsageLimit());
  const double maxPrimEnergy = GetHighEnergyUsageLimit();
  if (minPrimEnergy>=maxPrimEnergy) {
    return;
  }
  double edge = ecut;
  if (fIsElectron) {
    edge += edge;
  }
  //
  // compute number of e-/e+ kinetic energy grid
  int numPrimEnergies  = fSTNumPrimaryEnergyPerDecade*std::lrint(std::log10(maxPrimEnergy/minPrimEnergy))+1;
  numPrimEnergies      = std::max(numPrimEnergies,3);
  double logEmin       = std::log(minPrimEnergy);
  double delta         = std::log(maxPrimEnergy/minPrimEnergy)/(numPrimEnergies-1.0);
  AliasDataMaterialCuts *dataMatCut = new AliasDataMaterialCuts(numPrimEnergies, logEmin, 1./delta);
  fSamplingTables[indxlocal] = dataMatCut;
  for (int ipe=0; ipe<numPrimEnergies; ++ipe) {
    double pekin = std::exp(logEmin+ipe*delta);
    if (ipe==0 && minPrimEnergy==edge) {
      pekin = minPrimEnergy+1.*geant::eV; // would be zero otherwise
    }
    if (ipe==numPrimEnergies-1) {
      pekin = maxPrimEnergy;
    }
    // create the alias data struct
    LinAlias *als = new LinAlias(fSTNumSamplingElecEnergies);
    const double adum = 1.0/(fSTNumSamplingElecEnergies-1.0);
    for (int i=0; i<fSTNumSamplingElecEnergies; ++i) {
      double xi = i*adum;
      if (i==0) {
        xi = 0.0;
      } else if (i==fSTNumSamplingElecEnergies-1) {
        xi = 1.0;
      }
      als->fXdata[i] = xi;
      if (fIsElectron) {
        als->fYdata[i] = ComputeMollerPDF(xi, ecut, pekin);
      } else {
        als->fYdata[i] = ComputeBhabhaPDF(xi, ecut, pekin);
      }
    }
    //
    fAliasSampler->PreparLinearTable(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]), &(als->fAliasIndx[0]),
                                     fSTNumSamplingElecEnergies);
    fSamplingTables[indxlocal]->fAliasData[ipe] = als;
  }
}


/**
 *  The (restricted) ionization cross section per electron is computed for Moller \f$[e^-+e^-\to e^-+e^-]\f$ and Bhabha
 *  \f$[e^++e^-\to e^++e^-]\f$ scattering.
 *  Need to multiply by atomic number(Z) of the target atom to get atomic cross section (in internal [lenght^2] units)
 *  and by the electron density of the material to get macroscopic cross section (in units of [1/length]).
 *  Integration of the Moller and Bhabha cross sections between a minimum(\f$T_{kin}^{min}\f$) and
 *  maximum(\f$T_{kin}^{max}\f$) kinetic energy for incident particle kinetic energy
 *  \f$ E_{kin} \in \{T_{kin}^{min},T_{kin}^{max}\}\f$ result in (restricted cross section per \f$e^-\f$):
 *
 *  - for Moller scattering:
 *  \f[
 *   \int_{T_{kin}^{min}=T_{pcut}^{e^-}}^{T_{kin}^{max}=0.5E_{kin}} \frac{\mathrm{d}\sigma^{Moller}(E_{kin})}{\mathrm{d}T}\mathrm{d}T =
 *      \frac{2\pi r_e^2 m_0c^2}{\beta^2 E_{kin}} \left\{
 *          [\varepsilon_2-\varepsilon_1] \left(1-C + \frac{1}{\varepsilon_2\varepsilon_1}
 *               +\frac{1}{(1-\varepsilon_1)(1-\varepsilon_2)}
 *               -C\ln\frac{\varepsilon_2(1-\varepsilon_1)}{\varepsilon_1(1-\varepsilon_2)}
 *          \right)
 *       \right\}
 *  \f]
 *  - for Bhabha scattering:
 *  \f[
 *    \int_{T_{kin}^{min}=T_{pcut}^{e^-}}^{T_{kin}^{max}=E_{kin}} \frac{\mathrm{d}\sigma^{Bhabha}(E_{kin})}{\mathrm{d}T}\mathrm{d}T =
 *      \frac{2\pi r_e^2 m_0c^2}{E_{kin}} \left\{
 *          [\varepsilon_2-\varepsilon_1] \left[
 *            \frac{1}{\beta^2\varepsilon_2\varepsilon_1} + B_2 -\frac{B_3}{2}(\varepsilon_1+\varepsilon_2)
 *            +\frac{B_4}{3}\left[(\varepsilon_1+\varepsilon_2)^2-\varepsilon_1\varepsilon_2 \right]
              \right] - B_1\ln\frac{\varepsilon_2}{\varepsilon_1}
 *       \right\}
 *  \f]
 *  where
 *  \f[
 *  \begin{array}{lcl}
 *  T_{pcut}^{e^-} &\to& \textrm{electron production cut (kinetic) energy} \\
 *  E_{kin}   &\to&  \textrm{kinetic energy of the incident particle i.e. e-(Moller) or e+(Bhabha)} \\
 *  T         &\to&  \textrm{kinetic energy of the scattered e-}\\
 *  r_e       &\to&  \textrm{classical electron radius} \\
 *  m_0c^2    &\to&  \textrm{electron rest mass energy} \\
 *  \gamma    &\to&  \textrm{Lorentz factor} =E_t/(m_0c^2) \textrm{ i.e. total energy of the incident particle (e-/e+)
 *                    in rest mass energy unit} \\
 *  \beta     &\to&  \textrm{ratio of the speed of the incident particle to } c;\quad = P_t/E_t \textrm{ with } P_t=pc
 *                   \textrm{ i.e. total momentum of the incident particle (e-/e+) in rest mass energy unit} \\
 *  C         &\to&  (2\gamma-1)/\gamma^2 \\
 *  \varepsilon_1 &\to& e^- \textrm{production cut (kinetic) energy in incident particle kinetic energy unit} \\
 *  \varepsilon_2 &\to& \textrm{maximum (kinetic)energy transfer in incident particle kinetic energy unit i.e.}\\
 *             &&    \quad = \begin{cases}
 *                     1           & \quad \textrm{Bhabha scattering i.e. incident particle is e+} \\
 *                     \frac{1}{2} & \quad \textrm{Moller scattering i.e. incident particle is e-}
 *                   \end{cases}\\
 *  y         &\to&  1/(1+\gamma)\\
 *  B_1       &\to&  2-y^2\\
 *  B_2       &\to&  (1-2y)(3+y^2)\\
 *  B_4       &\to&  (1-2y)^2\\
 *  B_3       &\to&  B_4+(1-2y)^2
 *  \end{array}
 *  \f]
 *
 * Note, that due to the indistinguishability of the incident and secondary electron in Moller scattering, the electron
 * with the higher energy after the interaction is considered to be the primary so the restricted cross section in
 * Moller scattering is obtained by the integration of the differential Moller scattering cross section between
 * \f$ T_{pcut}^{e^-} \f$ and \f$E_{kin}/2\f$. Therefore, the kinetic energy above a discrete Moller event can accour is
 * \f$ E_{kin}^{TH-Moller} = 2T_{pcut}^{e^-}\f$ (the post interaction primary electron kinetic energy, i.e. the one with
 * higher kinetic energy, will be alway above this threshold kinetic energy). In case of Bhabha scattering, the
 * restricted cross section is obtained by the integration of Bhabha differential cross section between
 * \f$ T_{pcut}^{e^-} \f$ and \f$E_{kin}\f$ i.e. the post interaction positron kinetic energy can be lower than
 * \f$ T_{pcut}^{e^-} \f$.
 */
double MollerBhabhaIonizationModel::ComputeXSectionPerElectron(const double pcutenergy, const double primekin) {
  constexpr double xsecFactor = geant::kTwoPi*geant::kClassicElectronRadius*geant::kClassicElectronRadius
                               *geant::kElectronMassC2;
  // secondary e- produced only above production cut energy T_c:
  //   - for Moller scattering: kinetic energy of incident e- should be higher than 2T_c
  //   - for Bhabha scattering: kinetic energy of incident e+ should be higher than T_c
  // otherwise the discrete restricted cross section is zero.
  double maxEkin = primekin;
  if (fIsElectron) {
    maxEkin *= 0.5;
  }
  double xsec = 0.0;
  if (maxEkin>pcutenergy) {
    //set min/max energies in incoming particle kinetic energy unit
    const double epsmin = pcutenergy/primekin;
    const double epsmax = maxEkin/primekin;
    // set other parameters
    const double tau       = primekin/geant::kElectronMassC2; // i.e. E_kin in (mc^2) units
    const double gamma     = tau + 1.0;            // = E_t/(mc^2)  i.e. E_t in (mc^2) units
    const double gamma2    = gamma*gamma;          // \gamma^2 = [E_t/(mc^2)]^2
    const double beta2     = tau*(tau+2.0)/gamma2; // \beta2 i.e. [P_t/E_t]^2
    if (fIsElectron) { // Moller scattering i.e. e- + e- -> e- + e-
      const double parC = (2.0*gamma-1.0)/gamma2;
      xsec  = (epsmax-epsmin)*(1.0-parC + 1.0/(epsmin*epsmax) + 1.0/((1.0-epsmin)*(1.0-epsmax)))
              - parC*std::log((epsmax*(1.0-epsmin))/(epsmin*(1.0-epsmax)));
      xsec /= beta2;
    } else {           // Bhabha scattering i.e. e+ + e- -> e+ + e-
      const double y     = 1.0/(1.0+gamma);
      const double y2    = y*y;
      const double ydum  = 1.0-2.0*y;
      const double ydum2 = ydum*ydum;
      const double b1    = 2.0-y2;
      const double b2    = ydum*(3.0+y2);
      const double b4    = ydum*ydum2;
      const double b3    = b4+ydum2;
      const double e1e2  = epsmin*epsmax;
      const double e1pe2 = epsmin+epsmax;
      xsec  = (epsmax-epsmin)*(1.0/(beta2*e1e2) + b2 - 0.5*b3*e1pe2 + b4*(e1pe2*e1pe2-e1e2)/3.0)
              - b1*std::log(epsmax/epsmin);
    }
  }
  xsec *= xsecFactor/primekin;
  return std::max(xsec,0.0);
}


/**
  *  The differential(in fractional energy transfer) atomic cross section for Moller scattering \cite moller1932theorie
  *  \f[
  *     \frac{\mathrm{d}\sigma}{\mathrm{d}\varepsilon} = \frac{2\pi r_e^2 Z}{\beta^2 (\gamma-1)}\left[
  *      C_1 + \frac{1}{\varepsilon}\left(\frac{1}{\varepsilon}-C_2\right)
  *          + \frac{1}{\varepsilon'}\left(\frac{1}{\varepsilon'}-C_2\right)
  *     \right]
  *  \f]
  *  where
  *  \f[
  *  \begin{array}{lcl}
  *  E_{kin}     &\to&  \textrm{kinetic energy of the incident e-} \\
  *  T           &\to&  \textrm{kinetic energy of the scattered e- } \in [T_{cut}^{e-},0.5E_{kin}] \textrm{ where }
  *                     T_{cut} \textrm{ is the electron kinetic energy production threshold}\\
  *  \varepsilon &\to&  \textrm{fractional energy transfer i.e. } \varepsilon=T/E_{kin} \in [T_{cut}^{e-}/E_{kin},0.5]\\
  *  \varepsilon' &\to&  \textrm{fraction of remaining kinetic energy i.e. } \varepsilon'=1-\varepsilon \\
  *  r_e         &\to&  \textrm{classical electron radius} \\
  *  Z           &\to&  \textrm{atomic number of the target atom} \\
  *  \gamma      &\to&  \textrm{Lorentz factor} =E_t/(m_0c^2) \textrm{ i.e. total energy of the incident e- in rest mass
  *                    energy unit} \\
  *  \beta       &\to&  \textrm{ratio of the speed of the incident e- to } c;\quad = P_t/E_t \textrm{ with } P_t=pc \\
  *  C_1         &\to&  \left(\frac{\gamma-1}{\gamma}\right)^2 \\
  *  C_2         &\to&  \frac{2\gamma-1}{\gamma^2}
  *  \end{array}
  *  \f]
  *
  *  The following variable transformations are applied:
  *  - first \f$\varepsilon \to u = \ln(\varepsilon)\f$ so \f$\varepsilon=e^u,\; \mathrm{d}\varepsilon/\mathrm{d}u=e^u=
  *    \varepsilon \f$ which leads to \f$p(u) \propto \varepsilon \left[
  *      C_1 + \frac{1}{\varepsilon}\left(\frac{1}{\varepsilon}-C_2\right)
  *          + \frac{1}{\varepsilon'}\left(\frac{1}{\varepsilon'}-C_2\right)
  *     \right] \f$
  *  - then \f$ u \to \xi = [u-\ln(T_{cut}^{e-}/E_{kin})]/[\ln(0.5E_{kin}/T_{cut}^{e-})] \in [0,1]\f$ so
  *    \f$ u = \xi\ln(0.5E_{kin}/T_{cut}^{e-})+\ln(T_{cut}^{e-}/E_{kin}),\;
  *    \mathrm{d}u/\mathrm{d}\xi = \ln(0.5E_{kin}/T_{cut}^{e-})\f$ which leads to \f$ p(\xi) \propto \varepsilon \left[
  *      C_1 + \frac{1}{\varepsilon}\left(\frac{1}{\varepsilon}-C_2\right)
  *          + \frac{1}{\varepsilon'}\left(\frac{1}{\varepsilon'}-C_2\right)
  *     \right] \ln(0.5E_{kin}/T_{cut}^{e-}) \f$ where the last factor is just a constant.
  *
  *  So the transformed distribution
  *  \f[
  *      p(\xi) \propto \left[
  *      \varepsilon C_1 + \left(\frac{1}{\varepsilon}-C_2\right)
  *          + \frac{\varepsilon }{\varepsilon'}\left(\frac{1}{\varepsilon'}-C_2\right)
  *     \right]
  *  \f]
  *  where the transformed variable in terms of \f$T\f$ kinetic energy transfer is
  *  \f$\xi = \ln(T/T_{cut}^{e-})/\ln(0.5E_{kin}/T_{cut}^{e-}) \f$ so after the sampling of
  *  \f$\xi\f$ the kinetic energy transfer can be obtained as \f$T=T_{cut}^{e-}e^{\xi\ln(0.5E_{kin}/T_{cut}^{e-})}\f$.
  *
  */
double MollerBhabhaIonizationModel::ComputeMollerPDF(const double xi, const double pcutenergy, const double primekin) {
  const double tau       = primekin/geant::kElectronMassC2; // i.e. E_kin in (mc^2) units
  const double gamma     = tau + 1.0;            // = E_t/(mc^2)  i.e. E_t in (mc^2) units
  const double gamma2    = gamma*gamma;          // \gamma^2 = [E_t/(mc^2)]^2
  //double beta2     = tau*(tau+2.0)/gamma2; // \beta2 i.e. [P_t/E_t]^2
  double C1              = (gamma-1.0)/gamma;
  C1 *=C1;
  const double C2        = (2.0*gamma-1.0)/gamma2;

  const double dum0      = pcutenergy/primekin;
  const double dum1      = std::log(0.5/dum0);
  const double a         = std::exp(xi*dum1)*dum0; // this is eps =  exp(xi*ln(0.5*T_0/T_cut))*T_cut/T_0
  const double b         = 1.0-a;                  // eps'
  return ((1.0/a-C2)+a*C1+a/b*(1.0/b-C2)) *dum0; // xdum0 is just scaling; this is the shape
}


/**
  *  The differential(in fractional energy transfer) atomic cross section for Bhabha scattering
  *  \cite bhabha1936scattering
  *  \f[
  *     \frac{\mathrm{d}\sigma}{\mathrm{d}\varepsilon} = \frac{2\pi r_e^2 Z}{(\gamma-1)}\left[
  *      \frac{1}{\beta^2\varepsilon^2} - \frac{B_1}{\varepsilon} +  B_2 - \varepsilon B_3 + \varepsilon^2 B_4
  *     \right]
  *  \f]
  *  where
  *  \f[
  *  \begin{array}{lcl}
  *  E_{kin}     &\to&  \textrm{kinetic energy of the incident e+} \\
  *  T           &\to&  \textrm{kinetic energy of the scattered e- } \in [T_{cut}^{e-},E_{kin}] \textrm{ where }
  *                     T_{cut} \textrm{ is the electron kinetic energy production threshold}\\
  *  \varepsilon &\to&  \textrm{fractional energy transfer i.e. } \varepsilon=T/E_{kin} \in [T_{cut}^{e-}/E_{kin},1]\\
  *  r_e         &\to&  \textrm{classical electron radius} \\
  *  Z           &\to&  \textrm{atomic number of the target atom} \\
  *  \gamma      &\to&  \textrm{Lorentz factor} =E_t/(m_0c^2) \textrm{ i.e. total energy of the incident e- in rest mass
  *                    energy unit} \\
  *  \beta       &\to&  \textrm{ratio of the speed of the incident e- to } c;\quad = P_t/E_t \textrm{ with } P_t=pc \\
  *  y           &\to&  1/(1+\gamma)\\
  *  B_1         &\to&  2-y^2\\
  *  B_2         &\to&  (1-2y)(3+y^2)\\
  *  B_4         &\to&  (1-2y)^2\\
  *  B_3         &\to&  B_4+(1-2y)^2
  *  \end{array}
  *  \f]
  *
  *  The following variable transformations are applied:
  *  - first \f$\varepsilon \to u = \ln(\varepsilon)\f$ so \f$\varepsilon=e^u,\; \mathrm{d}\varepsilon/\mathrm{d}u=e^u=
  *    \varepsilon \f$ which leads to \f$p(u) \propto \varepsilon \left[
  *      \frac{1}{\beta^2\varepsilon^2} - \frac{B_1}{\varepsilon} +  B_2 - \varepsilon B_3 + \varepsilon^2 B_4
  *     \right] \f$
  *  - then \f$ u \to \xi = [u-\ln(T_{cut}^{e-}/E_{kin})]/[\ln(E_{kin}/T_{cut}^{e-})] \in [0,1]\f$ so
  *    \f$ u = \xi\ln(E_{kin}/T_{cut}^{e-})+\ln(T_{cut}^{e-}/E_{kin}),\;
  *    \mathrm{d}u/\mathrm{d}\xi = \ln(E_{kin}/T_{cut}^{e-})\f$ which leads to \f$ p(\xi) \propto \varepsilon \left[
  *      \frac{1}{\beta^2\varepsilon^2} - \frac{B_1}{\varepsilon} +  B_2 - \varepsilon B_3 + \varepsilon^2 B_4
  *     \right] \ln(E_{kin}/T_{cut}^{e-}) \f$ where the last factor is just a constant.
  *
  *  So the transformed distribution
  *  \f[
  *      p(\xi) \propto \left[
  *      \frac{1}{\beta^2\varepsilon} - B_1 +  \varepsilon ( B_2 - \varepsilon ( B_3 + \varepsilon B_4))
  *     \right]
  *  \f]
  *  where the transformed variable in terms of \f$T\f$ kinetic energy transfer is
  *  \f$\xi = \ln(T/T_{cut}^{e-})/\ln(E_{kin}/T_{cut}^{e-}) \f$ so after the sampling of
  *  \f$\xi\f$ the kinetic energy transfer can be obtained as \f$T=T_{cut}^{e-}e^{\xi\ln(E_{kin}/T_{cut}^{e-})}\f$.
  *
  */
double MollerBhabhaIonizationModel::ComputeBhabhaPDF(const double xi, const double pcutenergy, const double primekin) {
  const double tau       = primekin/geant::kElectronMassC2; // i.e. E_kin in (mc^2) units
  const double gamma     = tau + 1.0;            // = E_t/(mc^2)  i.e. E_t in (mc^2) units
  const double gamma2    = gamma*gamma;          // \gamma^2 = [E_t/(mc^2)]^2
  const double beta2     = tau*(tau+2.0)/gamma2; // \beta2 i.e. [P_t/E_t]^2
  const double y         = 1.0/(1.0+gamma);
  const double y2        = y*y;
  const double ydum      = 1.0-2.0*y;
  const double ydum2     = ydum*ydum;
  const double b1        = 2.0-y2;
  const double b2        = ydum*(3.0+y2);
  const double b4        = ydum*ydum2;
  const double b3        = b4+ydum2;

  const double dum0      = pcutenergy/primekin;
  const double dum1      = std::log(1.0/dum0);
  const double a         = std::exp(xi*dum1)*dum0; // this is eps =  = exp(xi*ln(T_0/T_cut))*T_cut/T_0

  return ((1.0/(a*beta2)-b1) + a*(b2+a*(a*b4-b3)));//
}


}  // namespace geantphysics
