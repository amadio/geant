
#include "Geant/SeltzerBergerBremsModel.h"
// from material
#include "Geant/Types.h"

#include "Geant/PhysicalConstants.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/Material.h"
#include "Geant/Element.h"
#include "Geant/MaterialProperties.h"

#include "Geant/MaterialCuts.h"

#include "Geant/Spline.h"
#include "Geant/GLIntegral.h"
#include "Geant/AliasTable.h"

#include "Geant/PhysicsParameters.h"

#include "Geant/Gamma.h"

#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"

// from geantV
#include "Geant/TaskData.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <Geant/AliasTableAlternative.h>
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

const double SeltzerBergerBremsModel::gMigdalConst = 4.0 * geant::units::kPi * geant::units::kClassicElectronRadius *
                                                     geant::units::kRedElectronComptonWLenght *
                                                     geant::units::kRedElectronComptonWLenght;
std::vector<SeltzerBergerBremsModel::XsecDataZet *> SeltzerBergerBremsModel::fXsecDataPerZet;
std::vector<double> SeltzerBergerBremsModel::fXsecLimits;
std::vector<double> SeltzerBergerBremsModel::fLoadDCSElectronEnergyGrid;
std::vector<double> SeltzerBergerBremsModel::fLoadDCSReducedPhotonEnergyGrid;

SeltzerBergerBremsModel::SeltzerBergerBremsModel(bool iselectron, const std::string &modelname)
    : EMModel(modelname), fIsElectron(iselectron)
{
  fNGL                             = 64;
  fSecondaryInternalCode           = -1;  // will be set at init
  fDCSMaxZet                       = 0;   // will be set at init
  fLoadDCSNumElectronEnergies      = 0;   // will be set at init
  fLoadDCSNumReducedPhotonEnergies = 0;   // will be set at init
  fSTNumElectronEnergyPerDecade    = 7;   // ST=>sampling table
  fSTNumSamplingPhotEnergies       = 64;  // ST=>sampling table
  fLogLoadDCSMinElecEnergy         = -1.; // will be set at init
  fInvLogLoadDCSDeltaEnergy        = -1.; // will be set at init

  fAliasSampler = nullptr;
  fGL           = nullptr;
}

SeltzerBergerBremsModel::~SeltzerBergerBremsModel()
{
  ClearLoadDCSData();
  // clear sampling tables (if any)
  if (GetUseSamplingTables()) {
    ClearSamplingTables();
  }
  fGlobalMatGCutIndxToLocal.clear();
  if (fAliasSampler) {
    delete fAliasSampler;
  }
  if (fGL) {
    delete fGL;
  }
}

void SeltzerBergerBremsModel::Initialize()
{
  EMModel::Initialize();
  LoadDCSData();
  if (!fGL) {
    fGL = new GLIntegral(fNGL, 0.0, 1.0);
  }
  fSecondaryInternalCode = Gamma::Definition()->GetInternalCode();
  if (GetUseSamplingTables()) { // if sampling tables were requested
    InitSamplingTables();
    int NCutTables = (int)fSamplingTables.size();
    for (int i = 0; i < NCutTables; ++i) {
      // push empty data for equal spacing
      fAliasData.fILDelta.emplace_back();
      fAliasData.fLogEmin.emplace_back();
      fAliasData.fNData.emplace_back();
      fAliasData.fTablesPerMatCut.push_back(std::unique_ptr<AliasDataForMatCut>(nullptr));
      if (fSamplingTables[i] == nullptr) continue;

      AliasDataMaterialCuts *oldTable = fSamplingTables[i];
      fAliasData.fILDelta[i]          = oldTable->fILDelta;
      fAliasData.fLogEmin[i]          = oldTable->fLogEmin;
      fAliasData.fNData[i]            = oldTable->fNData;
      fAliasData.fTablesPerMatCut[i]  = std::unique_ptr<AliasDataForMatCut>(
          new AliasDataForMatCut((int)oldTable->fAliasData.size(), fAliasData.fLogEmin[i], fAliasData.fILDelta[i]));
      auto &newTable = *fAliasData.fTablesPerMatCut[i];

      for (int j = 0; j < (int)oldTable->fAliasData.size(); ++j) {
        auto &oldETable = oldTable->fAliasData[j];
        LinAliasCached newETable((int)oldETable->fAliasIndx.size());

        for (int k = 0; k < (int)oldETable->fAliasIndx.size() - 1; ++k) {
          auto &transAlias       = newETable.fLinAliasData[k];
          transAlias.fAliasIndx  = oldETable->fAliasIndx[k];
          transAlias.fAliasW     = oldETable->fAliasW[k];
          transAlias.fX          = oldETable->fXdata[k];
          transAlias.fYdata      = oldETable->fYdata[k];
          transAlias.fXdelta     = oldETable->fXdata[k + 1] - oldETable->fXdata[k];
          transAlias.fYdataDelta = (oldETable->fYdata[k + 1] - oldETable->fYdata[k]) / oldETable->fYdata[k];
          transAlias.fXdivYdelta = transAlias.fXdelta / transAlias.fYdataDelta;
        }

        newTable.fAliasData.push_back(newETable);
      }
    }
  } else { // rejection
    // we need element selectors per MaterialCuts
    InitialiseElementSelectors(this, nullptr, false);
  }
}

double SeltzerBergerBremsModel::ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle * /*particle*/,
                                            bool istotal)
{
  const Material *mat = matcut->GetMaterial();
  const double *cuts  = matcut->GetProductionCutsInEnergy();
  double gammacut     = cuts[0];
  if (istotal) {
    // for the total stopping power we just need a gamma production cut >=kinenergy
    gammacut = 1.01 * kinenergy;
  }
  return ComputeDEDXPerVolume(mat, gammacut, kinenergy);
}

double SeltzerBergerBremsModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy,
                                                           const Particle * /*particle*/)
{
  double xsec = 0.0;
  if (kinenergy < GetLowEnergyUsageLimit() || kinenergy > GetHighEnergyUsageLimit()) {
    return xsec;
  }
  const Material *mat = matcut->GetMaterial();
  const double *cuts  = matcut->GetProductionCutsInEnergy();
  double gammacut     = cuts[0];
  xsec                = ComputeXSectionPerVolume(mat, gammacut, kinenergy);
  return xsec;
}

double SeltzerBergerBremsModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut,
                                                       double kinenergy, const Particle *)
{
  double xsec = 0.0;
  if (kinenergy < GetLowEnergyUsageLimit() || kinenergy > GetHighEnergyUsageLimit()) {
    return xsec;
  }
  const Material *mat = matcut->GetMaterial();
  const double *cuts  = matcut->GetProductionCutsInEnergy();
  double gammacut     = cuts[0];
  xsec                = ComputeXSectionPerAtom(elem, mat, gammacut, kinenergy);
  return xsec;
}

int SeltzerBergerBremsModel::SampleSecondaries(LightTrack &track, geant::TaskData *td)
{
  int numSecondaries         = 0;
  double ekin                = track.GetKinE();
  const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(track.GetMaterialCutCoupleIndex());
  const Material *mat        = matCut->GetMaterial();
  const double gammacut      = (matCut->GetProductionCutsInEnergy())[0]; // gamma cut
  if (ekin < GetLowEnergyUsageLimit() || ekin > GetHighEnergyUsageLimit() || ekin <= gammacut) {
    return numSecondaries;
  }
  // sample gamma energy
  double gammaEnergy = 0.;
  if (GetUseSamplingTables()) {
    double *rndArray = td->fDblArray;
    td->fRndm->uniform_array(3, rndArray);
    gammaEnergy = SamplePhotonEnergy(matCut, ekin, rndArray[0], rndArray[1], rndArray[2]);
  } else {
    // sample target element
    const Vector_t<Element *> &theElements = mat->GetElementVector();
    double targetElemIndx                  = 0;
    if (theElements.size() > 1) {
      targetElemIndx = SampleTargetElementIndex(matCut, ekin, td->fRndm->uniform());
    }
    const double zet = theElements[targetElemIndx]->GetZ();
    // gamma energy
    gammaEnergy = SamplePhotonEnergy(ekin, gammacut, zet, mat, td);
  }
  // sample gamma scattering angle in the scattering frame i.e. which z-dir points to the orginal e-/e+ direction
  double *rndArray = td->fDblArray;
  td->fRndm->uniform_array(2, rndArray);
  double cosTheta = 1.0;
  double sinTheta = 0.0;
  SamplePhotonDirection(ekin, sinTheta, cosTheta, rndArray[0]);
  const double phi = geant::units::kTwoPi * (rndArray[1]);
  // gamma direction in the scattering frame
  double gamDirX = sinTheta * Math::Cos(phi);
  double gamDirY = sinTheta * Math::Sin(phi);
  double gamDirZ = cosTheta;
  // rotate gamma direction to the lab frame:
  Math::RotateToLabFrame(gamDirX, gamDirY, gamDirZ, track.GetDirX(), track.GetDirY(), track.GetDirZ());
  // create the secondary partcile i.e. the gamma
  numSecondaries         = 1;
  LightTrack &gammaTrack = td->fPhysicsData->InsertSecondary();
  gammaTrack.SetDirX(gamDirX);
  gammaTrack.SetDirY(gamDirY);
  gammaTrack.SetDirZ(gamDirZ);
  gammaTrack.SetKinE(gammaEnergy);
  gammaTrack.SetGVcode(fSecondaryInternalCode); // gamma GV code
  gammaTrack.SetMass(0.0);
  gammaTrack.SetTrackIndex(track.GetTrackIndex()); // parent Track index
  //
  // compute the primary e-/e+ post interaction direction: from momentum vector conservation
  const double elInitTotalMomentum = std::sqrt(ekin * (ekin + 2.0 * geant::units::kElectronMassC2));
  // final momentum of the e-/e+ in the lab frame
  double elDirX = elInitTotalMomentum * track.GetDirX() - gammaEnergy * gamDirX;
  double elDirY = elInitTotalMomentum * track.GetDirY() - gammaEnergy * gamDirY;
  double elDirZ = elInitTotalMomentum * track.GetDirZ() - gammaEnergy * gamDirZ;
  // normalisation
  const double norm = 1.0 / std::sqrt(elDirX * elDirX + elDirY * elDirY + elDirZ * elDirZ);
  // update primary track direction
  track.SetDirX(elDirX * norm);
  track.SetDirY(elDirY * norm);
  track.SetDirZ(elDirZ * norm);
  // update primary track kinetic energy
  track.SetKinE(ekin - gammaEnergy);
  // return with number of secondaries i.e. 1 gamma
  return numSecondaries;
}

double SeltzerBergerBremsModel::MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle *) const
{
  double mine = (matcut->GetProductionCutsInEnergy())[0]; // gamma production cut in the given material-cuts
  return std::max(mine, GetLowEnergyUsageLimit());
}

void SeltzerBergerBremsModel::ClearLoadDCSData()
{
  size_t numElems = fXsecDataPerZet.size();
  for (size_t iz = 0; iz < numElems; ++iz) {
    XsecDataZet *xs = fXsecDataPerZet[iz];
    if (xs) {
      size_t nprime = xs->fXsecDataVect.size();
      for (size_t ie = 0; ie < nprime; ++ie) {
        xs->fXsecDataVect[ie].clear();
        if (xs->fSplineVect[ie]) {
          delete xs->fSplineVect[ie];
        }
      }
      xs->fXsecDataVect.clear();
      xs->fSplineVect.clear();
      delete xs;
    }
  }
  fXsecDataPerZet.clear();
}

void SeltzerBergerBremsModel::LoadDCSData()
{
  // get the path to the main physics data directory
  char *path = std::getenv("GEANT_PHYSICS_DATA");
  if (!path) {
    std::cerr << "******   ERROR in SeltzerBergerBremsModel::LoadDCSData() \n"
              << "         GEANT_PHYSICS_DATA is not defined! Set the GEANT_PHYSICS_DATA\n"
              << "         environmental variable to the location of Geant data directory!\n"
              << std::endl;
    exit(1);
  }
  char baseFilename[512];
  sprintf(baseFilename, "%s/brems/NIST_BREM2/nist_brems_", path);
  //
  FILE *f = nullptr;
  char filename[512];
  sprintf(filename, "%sgrid", baseFilename);
  f = fopen(filename, "r");
  if (!f) {
    std::cerr << "******   ERROR in SeltzerBergerBremsModel::LoadDCSData() \n"
              << "         " << filename << " could not be found!\n"
              << std::endl;
    exit(1);
  }
  //
  int nmatches = fscanf(f, "%d%d%d", &fDCSMaxZet, &fLoadDCSNumElectronEnergies, &fLoadDCSNumReducedPhotonEnergies);
  (void)nmatches;
  //
  fLoadDCSElectronEnergyGrid.resize(fLoadDCSNumElectronEnergies, 0.);
  fLoadDCSReducedPhotonEnergyGrid.resize(fLoadDCSNumReducedPhotonEnergies, 0.);
  for (int i = 0; i < fLoadDCSNumElectronEnergies; ++i) {
    nmatches = fscanf(f, "%lf", &(fLoadDCSElectronEnergyGrid[i]));
    fLoadDCSElectronEnergyGrid[i] *= geant::units::MeV; // change to internal energy units
  }
  fLogLoadDCSMinElecEnergy  = Math::Log(fLoadDCSElectronEnergyGrid[0]);
  fInvLogLoadDCSDeltaEnergy = 1. / (Math::Log(fLoadDCSElectronEnergyGrid[1] / fLoadDCSElectronEnergyGrid[0]));
  for (int i = 0; i < fLoadDCSNumReducedPhotonEnergies; ++i)
    nmatches = fscanf(f, "%lf", &(fLoadDCSReducedPhotonEnergyGrid[i]));
  fclose(f);
  //
  ClearLoadDCSData();
  fXsecDataPerZet.resize(fDCSMaxZet, nullptr);
  fXsecLimits.resize(fDCSMaxZet, 0.);
  //
  const Vector_t<Element *> theElements = Element::GetTheElementTable();
  // std::cout<<theElements;
  size_t numElements = theElements.size();
  for (size_t i = 0; i < numElements; ++i) {
    int zet = std::lrint(theElements[i]->GetZ());
    zet     = std::min(zet, fDCSMaxZet);
    sprintf(filename, "%s%d", baseFilename, zet);
    f = fopen(filename, "r");
    if (!f) {
      std::cerr << "******   ERROR in SeltzerBergerBremsModel::LoadDCSData() \n"
                << "         " << filename << " could not be found!\n"
                << std::endl;
      exit(1);
    }
    // allocate space for this elemental DCS
    fXsecDataPerZet[zet - 1] = new XsecDataZet(fLoadDCSNumElectronEnergies, fLoadDCSNumReducedPhotonEnergies);
    for (int ie = 0; ie < fLoadDCSNumElectronEnergies; ++ie) {
      for (int ik = 0; ik < fLoadDCSNumReducedPhotonEnergies; ++ik) {
        double dum                                      = 0.;
        nmatches                                        = fscanf(f, "%lf", &dum);
        fXsecDataPerZet[zet - 1]->fXsecDataVect[ie][ik] = dum * geant::units::millibarn; // change to internal units
      }
      // set up Spline for this elektron energy
      fXsecDataPerZet[zet - 1]->fSplineVect[ie] =
          new Spline(&(fLoadDCSReducedPhotonEnergyGrid[0]), &(fXsecDataPerZet[zet - 1]->fXsecDataVect[ie][0]),
                     fLoadDCSNumReducedPhotonEnergies);
    }
    // determine ylimit i.e. value at kappa=0.97 and log(E) = 4*log(10MeV)
    fXsecLimits[zet - 1] = GetDXSECValue(zet, Math::Exp(4. * Math::Log(10. * geant::units::MeV)), 0.97);
    fclose(f);
  }
}

double SeltzerBergerBremsModel::GetEkinIndex(double &ekin, int &ie)
{
  double eresid = 0.;
  if (ekin < fLoadDCSElectronEnergyGrid[0]) {
    ekin   = fLoadDCSElectronEnergyGrid[0]; // y1 index
    eresid = 0.;                            // (y2-y1)*resid + y1 ==> y1
    ie     = 0;
  } else if (ekin >= fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies - 1]) {
    ekin   = fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies - 1];
    ie     = fLoadDCSNumElectronEnergies - 2; // y1 index
    eresid = 1.0;                             // (y2-y1)*resid + y1 ==> y2
  } else {
    double lekin = Math::Log(ekin);
    eresid       = (lekin - fLogLoadDCSMinElecEnergy) * fInvLogLoadDCSDeltaEnergy;
    ie           = (int)eresid; // y1 index
    eresid -= ie;               // (y2-y1)*resid + y1
  }
  return eresid;
}

double SeltzerBergerBremsModel::GetDXSECValue(int zet, int ie, double eresid, double kappa)
{
  double res = 0.;
  int ik     = 0;
  //  double kresid =  0.;
  if (kappa < fLoadDCSReducedPhotonEnergyGrid[0]) {
    //    kappa  = fLoadDCSReducedPhotonEnergyGrid[0];
    //    kresid =  0.;
    ik        = 0;
    double y1 = fXsecDataPerZet[zet - 1]->fXsecDataVect[ie][ik];
    double y2 = fXsecDataPerZet[zet - 1]->fXsecDataVect[ie + 1][ik];
    res       = (y2 - y1) * eresid + y1;
  } else if (kappa >= fLoadDCSReducedPhotonEnergyGrid[fLoadDCSNumReducedPhotonEnergies - 1]) {
    //    kappa  = fLoadDCSReducedPhotonEnergyGrid[fLoadDCSNumReducedPhotonEnergies-1];
    ik = fLoadDCSNumReducedPhotonEnergies - 1;
    //    kresid =  1.0;
    double y1 = fXsecDataPerZet[zet - 1]->fXsecDataVect[ie][ik];
    double y2 = fXsecDataPerZet[zet - 1]->fXsecDataVect[ie + 1][ik];
    res       = (y2 - y1) * eresid + y1;
  } else {
    ik = std::lower_bound(fLoadDCSReducedPhotonEnergyGrid.begin(), fLoadDCSReducedPhotonEnergyGrid.end(), kappa) -
         fLoadDCSReducedPhotonEnergyGrid.begin() - 1;
    //    kresid =
    //    (kappa-fLoadDCSReducedPhotonEnergyGrid[ik])/(fLoadDCSReducedPhotonEnergyGrid[ik+1]-fLoadDCSReducedPhotonEnergyGrid[ik]);
    double y1 = fXsecDataPerZet[zet - 1]->fSplineVect[ie]->GetValueAt(kappa, ik);
    double y2 = fXsecDataPerZet[zet - 1]->fSplineVect[ie + 1]->GetValueAt(kappa, ik);
    res       = (y2 - y1) * eresid + y1;
  }
  // spline on kappa
  return res;
}

double SeltzerBergerBremsModel::GetDXSECValue(int zet, double eprim, double kappa)
{
  double res = 0.;
  // linear interp on log e- energy
  int ie        = 0;
  double ekin   = eprim;
  double eresid = GetEkinIndex(ekin, ie);
  // determine kappa index
  int ik = 0;
  //  double kresid =  0.;
  if (kappa < fLoadDCSReducedPhotonEnergyGrid[0]) {
    //    kappa  = fLoadDCSReducedPhotonEnergyGrid[0];
    //    kresid =  0.;
    ik        = 0;
    double y1 = fXsecDataPerZet[zet - 1]->fXsecDataVect[ie][ik];
    double y2 = fXsecDataPerZet[zet - 1]->fXsecDataVect[ie + 1][ik];
    res       = (y2 - y1) * eresid + y1;
  } else if (kappa >= fLoadDCSReducedPhotonEnergyGrid[fLoadDCSNumReducedPhotonEnergies - 1]) {
    //    kappa  = fLoadDCSReducedPhotonEnergyGrid[fLoadDCSNumReducedPhotonEnergies-1];
    ik = fLoadDCSNumReducedPhotonEnergies - 1;
    //    kresid =  1.0;
    double y1 = fXsecDataPerZet[zet - 1]->fXsecDataVect[ie][ik];
    double y2 = fXsecDataPerZet[zet - 1]->fXsecDataVect[ie + 1][ik];
    res       = (y2 - y1) * eresid + y1;
  } else {
    ik = std::lower_bound(fLoadDCSReducedPhotonEnergyGrid.begin(), fLoadDCSReducedPhotonEnergyGrid.end(), kappa) -
         fLoadDCSReducedPhotonEnergyGrid.begin() - 1;
    //    kresid =
    //    (kappa-fLoadDCSReducedPhotonEnergyGrid[ik])/(fLoadDCSReducedPhotonEnergyGrid[ik+1]-fLoadDCSReducedPhotonEnergyGrid[ik]);
    double y1 = fXsecDataPerZet[zet - 1]->fSplineVect[ie]->GetValueAt(kappa, ik);
    double y2 = fXsecDataPerZet[zet - 1]->fSplineVect[ie + 1]->GetValueAt(kappa, ik);
    res       = (y2 - y1) * eresid + y1;
  }
  // spline on kappa
  return res;
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
double SeltzerBergerBremsModel::ComputeDEDXPerVolume(const Material *mat, double gammaprodcutenergy,
                                                     double electronekin)
{
  double dedx = 0.0;
  //
  if (electronekin < fLoadDCSElectronEnergyGrid[0] ||
      electronekin > fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies - 1])
    return dedx;

  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol() * gMigdalConst;

  double ptot2      = electronekin * (electronekin + 2.0 * geant::units::kElectronMassC2);
  double etot2      = ptot2 + geant::units::kElectronMassC2 * geant::units::kElectronMassC2;
  double ibeta2     = etot2 / ptot2;
  double densityCor = densityFactor * etot2; //*electronekin*electronekin;  // this is k_p^2

  // find electron energy grid index ie such that E_ie <= ekin <_E_ie+1
  int ie        = 0;
  double eresid = GetEkinIndex(electronekin, ie);
  // we will need the element composition of this material
  const Vector_t<Element *> theElements   = mat->GetElementVector();
  const double *theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int numElems                            = theElements.size();
  // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
  const std::vector<double> &glW = fGL->GetWeights();
  const std::vector<double> &glX = fGL->GetAbscissas();
  // integrate
  double integral   = 0.0;
  double upperlimit = gammaprodcutenergy;
  if (upperlimit > electronekin) {
    upperlimit = electronekin;
  }
  double kappacr = upperlimit / electronekin;
  for (int i = 0; i < fNGL; ++i) {
    double t      = glX[i] * kappacr; // kappa
    double egamma = glX[i] * upperlimit;
    double sum    = 0.0;
    for (int ielem = 0; ielem < numElems; ++ielem) {
      double zet = theElements[ielem]->GetZ();
      int izet   = std::lrint(zet);
      izet       = std::min(izet, fDCSMaxZet);
      double val = GetDXSECValue(izet, ie, eresid, t); // sp[ielem]->GetValueAt(t);
      if (!fIsElectron) {
        val *= PositronCorrection(electronekin, ibeta2, egamma / electronekin, zet);
      }
      sum += theAtomicNumDensityVector[ielem] * zet * zet * val;
    }
    integral += glW[i] * sum / (1. + densityCor / (egamma * egamma));
    // x 1/(1+k_p^2/k^2) i.e. density effect correction
  }
  dedx = upperlimit * ibeta2 * integral;
  return dedx;
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
double SeltzerBergerBremsModel::ComputeXSectionPerVolume(const Material *mat, double gammaprodcutenergy,
                                                         double electronekin)
{
  double xsec = 0.0;
  if (electronekin <= gammaprodcutenergy || electronekin < fLoadDCSElectronEnergyGrid[0] ||
      electronekin > fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies - 1])
    return xsec;
  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol() * gMigdalConst;
  double ptot2         = electronekin * (electronekin + 2.0 * geant::units::kElectronMassC2);
  double etot2         = ptot2 + geant::units::kElectronMassC2 * geant::units::kElectronMassC2;
  double ibeta2        = etot2 / ptot2;
  double densityCor    = densityFactor * etot2; // electronekin*electronekin;  // this is k_p^2
  double kappacr       = gammaprodcutenergy / electronekin;
  double logikappacr   = Math::Log(1. / kappacr);
  // find electron energy grid index ie such that E_ie <= ekin <_E_ie+1
  int ie        = 0;
  double eresid = GetEkinIndex(electronekin, ie);
  // we will need the element composition of this material
  const Vector_t<Element *> theElements   = mat->GetElementVector();
  const double *theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int numElems                            = theElements.size();
  // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
  const std::vector<double> &glW = fGL->GetWeights();
  const std::vector<double> &glX = fGL->GetAbscissas();
  // integrate
  double integral = 0.0;
  for (int i = 0; i < fNGL; ++i) {
    double dumx = 1.0 - Math::Exp(glX[i] * logikappacr) * kappacr;
    // double x = Math::Log(dumx+1.0e-12);
    double egamma = (1.0 - dumx) * electronekin; // 1-dumx => kappa
    double sum    = 0.0;
    for (int ielem = 0; ielem < numElems; ++ielem) {
      double zet = theElements[ielem]->GetZ();
      int izet   = std::lrint(zet);
      izet       = std::min(izet, fDCSMaxZet);
      double val = GetDXSECValue(izet, ie, eresid, 1.0 - dumx); // sp[ielem]->GetValueAt(x);
      if (!fIsElectron) {
        val *= PositronCorrection(electronekin, ibeta2, (1.0 - dumx), zet);
      }
      sum += theAtomicNumDensityVector[ielem] * zet * zet * val;
    }
    integral += glW[i] * sum / (1. + densityCor / (egamma * egamma));
  }
  xsec = logikappacr * ibeta2 * integral;
  return xsec;
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
 * evaluated (interpolated) at \f$\phi =\ln[1-\kappa_c e^{\xi\ln(1/\kappa_c)}+10^{-12}]\f$ at a given value of
 * \f$\xi\f$.
 */
double SeltzerBergerBremsModel::ComputeXSectionPerAtom(const Element *elem, const Material *mat,
                                                       double gammaprodcutenergy, double electronekin)
{
  double xsec = 0.0;
  if (electronekin <= gammaprodcutenergy || electronekin < fLoadDCSElectronEnergyGrid[0] ||
      electronekin > fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies - 1])
    return xsec;

  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  double densityFactor = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol() * gMigdalConst;

  double zet         = elem->GetZ();
  int izet           = std::lrint(zet);
  izet               = std::min(izet, fDCSMaxZet);
  double ptot2       = electronekin * (electronekin + 2.0 * geant::units::kElectronMassC2);
  double etot2       = ptot2 + geant::units::kElectronMassC2 * geant::units::kElectronMassC2;
  double ibeta2      = etot2 / ptot2;
  double densityCor  = densityFactor * etot2; // electronekin*electronekin;  // this is k_p^2
  double kappacr     = gammaprodcutenergy / electronekin;
  double logikappacr = Math::Log(1. / kappacr);
  // find electron energy grid index ie such that E_ie <= ekin <_E_ie+1
  int ie        = 0;
  double eresid = GetEkinIndex(electronekin, ie);
  // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
  const std::vector<double> &glW = fGL->GetWeights();
  const std::vector<double> &glX = fGL->GetAbscissas();
  // integrate
  double integral = 0.0;
  for (int i = 0; i < fNGL; ++i) {
    double dumx = 1.0 - Math::Exp(glX[i] * logikappacr) * kappacr; // ez a felso x -nek az exp(x)-e
    // double x = Math::Log(dumx+1.0e-12); // 1-dumx = kappa
    double egamma = (1.0 - dumx) * electronekin;
    double poscor = 1.0;
    if (!fIsElectron) {
      poscor *= PositronCorrection(electronekin, ibeta2, 1.0 - dumx, zet);
    }
    double val = GetDXSECValue(izet, ie, eresid, 1.0 - dumx);
    integral += glW[i] * val / (1. + densityCor / (egamma * egamma));
  }
  xsec = logikappacr * zet * zet * ibeta2 * integral;
  return xsec;
}

// correction for positrons : DCS must be multiplied by this for positrons
// ephoton is the reduced photon energy
double SeltzerBergerBremsModel::PositronCorrection(double ekinelectron, double ibeta2electron, double ephoton, double z)
{
  constexpr double dum1 = geant::units::kTwoPi * geant::units::kFineStructConst;
  double poscor         = 0.0;
  double ibeta1         = std::sqrt(ibeta2electron);
  double e2             = ekinelectron * (1.0 - ephoton);
  if (e2 > 0.0) {
    double ibeta2 = (e2 + geant::units::kElectronMassC2) / std::sqrt(e2 * (e2 + 2.0 * geant::units::kElectronMassC2));
    double dum0   = dum1 * z * (ibeta1 - ibeta2);
    if (dum0 < -12.0) {
      poscor = 0.0;
    } else {
      poscor = Math::Exp(dum0);
    }
  } else {
    poscor = 0.0;
  }
  return poscor;
}

// ephoton is the reduced photon energy
double SeltzerBergerBremsModel::PositronCorrection1(double ekinelectron, double ephoton, double gcutener, double z)
{
  constexpr double dum1 = geant::units::kTwoPi * geant::units::kFineStructConst;
  double poscor         = 0.0;
  double e1             = ekinelectron - gcutener; // here is the dif.
  double ibeta1 = (e1 + geant::units::kElectronMassC2) / std::sqrt(e1 * (e1 + 2.0 * geant::units::kElectronMassC2));
  double e2     = ekinelectron * (1.0 - ephoton);
  double ibeta2 = (e2 + geant::units::kElectronMassC2) / std::sqrt(e2 * (e2 + 2.0 * geant::units::kElectronMassC2));
  double ddum   = dum1 * z * (ibeta1 - ibeta2);
  if (ddum < -12.0) {
    poscor = 0.0;
  } else {
    poscor = Math::Exp(ddum);
  }
  return poscor;
}

double SeltzerBergerBremsModel::SamplePhotonEnergy(const MaterialCuts *matcut, double eekin, double r1, double r2,
                                                   double r3)
{
  double egamma = 0.;
  //
  const double gcut       = (matcut->GetProductionCutsInEnergy())[0];
  const double etot       = eekin + geant::units::kElectronMassC2;
  const double densityCor = matcut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol() *
                            gMigdalConst * (etot * etot); // k_p^2
  const int mcIndxLocal = fGlobalMatGCutIndxToLocal[matcut->GetIndex()];
  // determine electron energy lower grid point
  const double leekin = Math::Log(eekin);
  //
  int indxEekin = fSamplingTables[mcIndxLocal]->fAliasData.size() - 1;
  if (eekin < GetHighEnergyUsageLimit()) {
    const double val       = (leekin - fSamplingTables[mcIndxLocal]->fLogEmin) * fSamplingTables[mcIndxLocal]->fILDelta;
    indxEekin              = (int)val; // lower electron energy bin index
    const double pIndxHigh = val - indxEekin;
    if (r1 < pIndxHigh) ++indxEekin;
  }
  // sample the transformed variable
  const LinAlias *als = fSamplingTables[mcIndxLocal]->fAliasData[indxEekin];
  egamma = fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]), &(als->fAliasIndx[0]),
                                       fSTNumSamplingPhotEnergies, r2, r3);
  // transform back to gamma energy
  const double dum1 = gcut * gcut + densityCor;
  const double dum2 = (eekin * eekin + densityCor) / dum1;

  return std::sqrt(dum1 * Math::Exp(egamma * Math::Log(dum2)) - densityCor);
}

double SeltzerBergerBremsModel::SamplePhotonEnergy(double eekin, double gcut, double zet, const Material *mat,
                                                   geant::TaskData *td)
{
  double egamma = 0.;
  //
  const double etot = eekin + geant::units::kElectronMassC2;
  const double densityCor =
      mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol() * gMigdalConst * (etot * etot); // k_p^2
  const double kappac = gcut / eekin;
  int izet            = std::lrint(zet);
  izet                = std::min(izet, fDCSMaxZet);
  // find electron energy grid index ie such that E_ie <= ekin <_E_ie+1
  int ie        = 0;
  double eresid = GetEkinIndex(eekin, ie);
  // get max value for rejection at (kappa=kappa_c, eekin)
  double vmax = 1.02 * GetDXSECValue(izet, ie, eresid, kappac);
  //
  // majoranta corrected vmax for e-
  constexpr double epeaklimit = 300.0 * geant::units::MeV;
  constexpr double elowlimit  = 20.0 * geant::units::keV;
  if (fIsElectron && kappac < 0.97 && ((eekin > epeaklimit) || (eekin < elowlimit))) {
    vmax = std::max(vmax, std::min(fXsecLimits[izet - 1], 1.1 * GetDXSECValue(izet, ie, eresid, 0.97)));
  }
  if (kappac < 0.05) {
    vmax *= 1.2;
  }
  // comp. min/max of the tranformed variable xi (kappa) = ln[kappa^2 E_kin^2 +k_p^2]
  //    - kappa_min = kappa_c => xi(kappa_c) = ln[k_c^2 + k_p^2]
  //    - kappa_max = 1       => xi(1)       = ln[E_kin^2 + k_p^2]
  //    => xi \in [xi_min, xi_max] => [ ln[k_c^2 + k_p^2], ln[E_kin^2 + k_p^2]]
  const double minXi = Math::Log(gcut * gcut + densityCor);
  //  double maxXi     = Math::Log(eekin*eekin + densityCor);
  const double delXi = Math::Log(eekin * eekin + densityCor) - minXi;
  double val         = 0.;
  double *rndArray   = td->fDblArray;
  do {
    td->fRndm->uniform_array(2, rndArray);
    egamma             = std::sqrt(std::max(Math::Exp(minXi + rndArray[0] * delXi) - densityCor, 0.));
    const double kappa = egamma / eekin;
    // val    =  GetDXSECValue(izet, eekin, kappa);
    val = GetDXSECValue(izet, ie, eresid, kappa);
    // positron cor.
    if (!fIsElectron) {
      val *= PositronCorrection1(eekin, kappa, gcut, zet);
    }
  } while (val < vmax * rndArray[1]);
  return egamma;
}

// the simple DipBustgenerator
void SeltzerBergerBremsModel::SamplePhotonDirection(double elenergy, double &sinTheta, double &cosTheta, double rndm)
{
  const double c = 4. - 8. * rndm;
  double a       = c;
  double signc   = 1.;
  if (c < 0.) {
    signc = -1.;
    a     = -c;
  }
  const double delta = 0.5 * (std::sqrt(a * a + 4.) + a);
  //  delta += a;
  //  delta *= 0.5;

  const double cofA = -signc * Math::Exp(Math::Log(delta) / 3.0);
  cosTheta          = cofA - 1. / cofA;

  const double tau  = elenergy / geant::units::kElectronMassC2;
  const double beta = std::sqrt(tau * (tau + 2.)) / (tau + 1.);

  cosTheta = (cosTheta + beta) / (1. + cosTheta * beta);
  // check cosTheta limit
  cosTheta = std::min(1.0, cosTheta);
  // if (cosTheta>1.0) {
  //  cosTheta = 1.0;
  //}
  sinTheta = std::sqrt((1. - cosTheta) * (1. + cosTheta));
}

// clear all sampling tables (if any)
void SeltzerBergerBremsModel::ClearSamplingTables()
{
  size_t numST = fSamplingTables.size();
  for (size_t i = 0; i < numST; ++i) {
    AliasDataMaterialCuts *st = fSamplingTables[i];
    if (st) {
      size_t numAT = st->fAliasData.size();
      for (size_t j = 0; j < numAT; ++j) {
        LinAlias *la = st->fAliasData[j];
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

void SeltzerBergerBremsModel::InitSamplingTables()
{
  // clear all sampling tables (if any)
  ClearSamplingTables();
  // determine global-to-local matcut indices:
  // - get number of different material-gammacut pairs
  // - allocate space and fill the global to local material-cut index map
  const std::vector<MaterialCuts *> theMaterialCutsTable = MaterialCuts::GetTheMaterialCutsTable();
  int numMaterialCuts                                    = theMaterialCutsTable.size();
  int numDifferentMatGCuts                               = 0;
  fGlobalMatGCutIndxToLocal.resize(numMaterialCuts, -2);
  for (int i = 0; i < numMaterialCuts; ++i) {
    // if the current MaterialCuts does not belong to the current active regions
    if (!IsActiveRegion(theMaterialCutsTable[i]->GetRegionIndex())) {
      continue;
    }
    bool isnew = true;
    int j      = 0;
    for (; j < numDifferentMatGCuts; ++j) {
      if (theMaterialCutsTable[i]->GetMaterial()->GetIndex() == theMaterialCutsTable[j]->GetMaterial()->GetIndex() &&
          theMaterialCutsTable[i]->GetProductionCutsInEnergy()[0] ==
              theMaterialCutsTable[j]->GetProductionCutsInEnergy()[0]) {
        isnew = false;
        break;
      }
    }
    if (isnew) {
      fGlobalMatGCutIndxToLocal[i] = numDifferentMatGCuts;
      ++numDifferentMatGCuts;
    } else {
      fGlobalMatGCutIndxToLocal[i] = fGlobalMatGCutIndxToLocal[j];
    }
  }
  fSamplingTables.resize(numDifferentMatGCuts, nullptr);
  // create an AliasTable object
  if (fAliasSampler) {
    delete fAliasSampler;
  }
  fAliasSampler = new AliasTable();
  for (int i = 0; i < numMaterialCuts; ++i) {
    const MaterialCuts *matCut = theMaterialCutsTable[i];
    int indxLocal              = fGlobalMatGCutIndxToLocal[i];
    if (indxLocal > -1 && !(fSamplingTables[indxLocal])) {
      BuildSamplingTableForMaterialCut(matCut, indxLocal);
    }
  }
}

void SeltzerBergerBremsModel::BuildSamplingTableForMaterialCut(const MaterialCuts *matcut, int indxlocal)
{
  const Material *mat        = matcut->GetMaterial();
  const double gcut          = (matcut->GetProductionCutsInEnergy())[0];
  const double minElecEnergy = MinimumPrimaryEnergy(matcut, nullptr);
  const double maxElecEnergy = GetHighEnergyUsageLimit();
  if (minElecEnergy >= maxElecEnergy) {
    return;
  }
  //
  const double densityFactor = gMigdalConst * mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol();
  //
  const Vector_t<Element *> &theElements  = mat->GetElementVector();
  const double *theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  const int numElems                      = theElements.size();
  //
  // compute number of e-/e+ kinetic energy grid
  int numElecEnergies = fSTNumElectronEnergyPerDecade * std::lrint(Math::Log10(maxElecEnergy / minElecEnergy)) + 1;
  numElecEnergies     = std::max(numElecEnergies, 3);
  double logEmin      = Math::Log(minElecEnergy);
  double delta        = Math::Log(maxElecEnergy / minElecEnergy) / (numElecEnergies - 1.0);
  AliasDataMaterialCuts *dataMatCut = new AliasDataMaterialCuts(numElecEnergies, logEmin, 1. / delta);
  fSamplingTables[indxlocal]        = dataMatCut;
  for (int ie = 0; ie < numElecEnergies; ++ie) {
    double eekin = Math::Exp(dataMatCut->fLogEmin + ie * delta); // E_kin_i
    if (ie == 0 && minElecEnergy == gcut) {
      eekin = minElecEnergy + 1. * geant::units::eV; // would be zero otherwise
    }
    if (ie == numElecEnergies - 1) {
      eekin = maxElecEnergy;
    }
    const double etot       = eekin + geant::units::kElectronMassC2;
    const double densityCor = densityFactor * etot * etot; // this is k_p^2
    // xi(kappa) = log(kappa^2*E^2 + k_p^2)  ==> kappa =
    const double kappac = gcut / eekin;
    const double dum0   = gcut * gcut + densityCor; // xi(kappa=kappac)
    //    const double xiMax      = Math::Log(eekin*eekin+densityCor); // xi(kappa=1)
    const double xiNorm = Math::Log((eekin * eekin + densityCor) / dum0);
    // determine energy interpolation: low energy bin index and residual for the interpolation
    int ielow     = 0;
    double eresid = GetEkinIndex(eekin, ielow);
    // create the alias data struct
    LinAlias *als = new LinAlias(fSTNumSamplingPhotEnergies);
    // insert the first point: interpolated value at kappa_c
    int numData    = 1;
    als->fXdata[0] = 0.; // r(xi(kappa)) \in [0,1]
    for (int ikappa = 0; ikappa < fLoadDCSNumReducedPhotonEnergies; ++ikappa) {
      if (fLoadDCSReducedPhotonEnergyGrid[ikappa] > kappac &&
          std::abs(1. - fLoadDCSReducedPhotonEnergyGrid[ikappa] / kappac) > 1.e-8) {
        const double kappa   = fLoadDCSReducedPhotonEnergyGrid[ikappa];
        const double r       = Math::Log((kappa * kappa * eekin * eekin + densityCor) / dum0) / xiNorm;
        als->fXdata[numData] = r;
        ++numData;
      }
    }
    // get the corresponding scalled DCS values
    for (int id = 0; id < numData; ++id) {
      const double r     = als->fXdata[id];
      const double kappa = std::sqrt(Math::Exp(r * xiNorm) * dum0 - densityCor) / eekin;
      double dcs         = 0.;
      for (int ielem = 0; ielem < numElems; ++ielem) {
        const double zet = theElements[ielem]->GetZ();
        int izet         = std::lrint(zet);
        izet             = std::min(izet, fDCSMaxZet);
        double val       = GetDXSECValue(izet, ielow, eresid, kappa) * theAtomicNumDensityVector[ielem] * zet * zet;
        if (!fIsElectron) {
          val *= PositronCorrection1(eekin, std::min(kappa, 1. - 1e-12), gcut, zet);
        }
        dcs += val;
      }
      als->fYdata[id] = dcs;
    }
    // fill remaining points by minimizing the liear interpollation error
    while (numData < fSTNumSamplingPhotEnergies) {
      // find the lower index of the bin, where we have the biggest linear interp. error compared to spline
      double maxerr  = 0.0; // value of the current maximum error
      double thexval = 0.0;
      double theyval = 0.0;
      int maxerrindx = 0; // the lower index of the corresponding bin
      for (int i = 0; i < numData - 1; ++i) {
        const double xx    = 0.5 * (als->fXdata[i] + als->fXdata[i + 1]); // mid point
        const double yy    = 0.5 * (als->fYdata[i] + als->fYdata[i + 1]); // lin func val at the mid point
        const double kappa = std::sqrt(dum0 * Math::Exp(xx * xiNorm) - densityCor) / eekin;
        double spval       = 0.; // spline intp. val. at mid point
        for (int ielem = 0; ielem < numElems; ++ielem) {
          const double zet = theElements[ielem]->GetZ();
          int izet         = std::lrint(zet);
          izet             = std::min(izet, fDCSMaxZet);
          double val       = GetDXSECValue(izet, ielow, eresid, kappa) * theAtomicNumDensityVector[ielem] * zet * zet;
          if (!fIsElectron) {
            val *= PositronCorrection1(eekin, std::min(kappa, 1. - 1e-12), gcut, zet);
          }
          spval += val;
        }
        double err = std::fabs(yy - spval); // works better than fabs(1-yy/spval) might try constrained spline?
        if (err > maxerr) {
          maxerr     = err;
          maxerrindx = i;
          thexval    = xx;
          theyval    = spval;
        }
      }
      // extend x,y data by puting a spline interp.ted value at the mid point of the highest error bin
      // first shift all values to the right
      for (int j = numData; j > maxerrindx + 1; --j) {
        als->fXdata[j] = als->fXdata[j - 1];
        als->fYdata[j] = als->fYdata[j - 1];
      }
      // fill x mid point
      als->fXdata[maxerrindx + 1] = thexval;
      als->fYdata[maxerrindx + 1] = theyval;
      // increase number of data
      ++numData;
    }
    // prepare the sampling table at this eekin
    fAliasSampler->PreparLinearTable(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]), &(als->fAliasIndx[0]),
                                     fSTNumSamplingPhotEnergies);
    fSamplingTables[indxlocal]->fAliasData[ie] = als;
  }
}

void SeltzerBergerBremsModel::SampleSecondaries(LightTrack_v &tracks, geant::TaskData *td)
{
  const int N               = tracks.GetNtracks();
  double *gammaeEnergyArray = td->fPhysicsData->fPhysicsScratchpad.fEps;
  double *gammaCutArr       = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr;
  int *IZet                 = td->fPhysicsData->fPhysicsScratchpad.fIzet;
  double *Zet               = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr2;
  double *densityCorrArr    = td->fPhysicsData->fPhysicsScratchpad.fR0;

  for (int i = 0; i < N; i += kVecLenD) {

    Double_v primEkin = tracks.GetKinEVec(i);
    Double_v gammaCut;
    IndexD_v mcIndxLocal; // Used by alias method
    Double_v densityCor;
    for (int l = 0; l < kVecLenD; ++l) {
      const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(tracks.GetMaterialCutCoupleIndex(i + l));
      Set(gammaCut, l, matCut->GetProductionCutsInEnergy()[0]);
      Set(densityCor, l, matCut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol());
      if (GetUseSamplingTables()) {
        Set(mcIndxLocal, l, fGlobalMatGCutIndxToLocal[matCut->GetIndex()]);
      } else {
        // sample target element
        const Vector_t<Element *> &theElements = matCut->GetMaterial()->GetElementVector();
        double targetElemIndx                  = 0;
        if (theElements.size() > 1) {
          targetElemIndx = SampleTargetElementIndex(matCut, Get(primEkin, l), td->fRndm->uniform());
        }
        const double zet = theElements[targetElemIndx]->GetZ();
        IZet[i + l]      = std::min((int)std::lrint(zet), fDCSMaxZet);
        Zet[i + l]       = zet;
      }
    }
    Double_v totalEn = primEkin + geant::units::kElectronMassC2;
    densityCor *= gMigdalConst * (totalEn * totalEn);
    assert((primEkin >= gammaCut).isFull()); // Cut filtering should be applied up the call chain.
    if (GetUseSamplingTables()) {
      Double_v r1          = td->fRndm->uniformV();
      Double_v r2          = td->fRndm->uniformV();
      Double_v r3          = td->fRndm->uniformV();
      Double_v primEkin    = tracks.GetKinEVec(i);
      Double_v gammaEnergy = SampleEnergyTransfer(gammaCut, densityCor, mcIndxLocal, fAliasData.fLogEmin.data(),
                                                  fAliasData.fILDelta.data(), primEkin, r1, r2, r3);
      vecCore::Store(gammaEnergy, gammaeEnergyArray + i);
    } else {
      vecCore::Store(gammaCut, gammaCutArr + i);
      vecCore::Store(densityCor, densityCorrArr + i);
    }
  }

  if (!GetUseSamplingTables()) {
    tracks.GetKinEArr()[N] = tracks.GetKinEArr()[N - 1];
    gammaCutArr[N]         = gammaCutArr[N - 1];
    IZet[N]                = IZet[N - 1];
    Zet[N]                 = Zet[N - 1];
    densityCorrArr[N]      = densityCorrArr[N - 1];

    SampleEnergyTransfer(tracks.GetKinEArr(), gammaCutArr, IZet, Zet, gammaeEnergyArray, densityCorrArr, N, td);
  }

  for (int i = 0; i < N; i += kVecLenD) {
    Double_v ekin = tracks.GetKinEVec(i);
    Double_v gammaEnergy;
    vecCore::Load(gammaEnergy, gammaeEnergyArray + i);

    Double_v cosTheta = 1.0;
    Double_v sinTheta = 0.0;
    Double_v rnd0     = td->fRndm->uniformV();
    SamplePhotonDirection(ekin, sinTheta, cosTheta, rnd0);

    Double_v rnd1      = td->fRndm->uniformV();
    const Double_v phi = geant::units::kTwoPi * (rnd1);
    Double_v sinPhi, cosPhi;
    Math::SinCos(phi, sinPhi, cosPhi);
    // gamma direction in the scattering frame
    Double_v gamDirX = sinTheta * cosPhi;
    Double_v gamDirY = sinTheta * sinPhi;
    Double_v gamDirZ = cosTheta;
    // rotate gamma direction to the lab frame:
    Math::RotateToLabFrame(gamDirX, gamDirY, gamDirZ, tracks.GetDirXVec(i), tracks.GetDirYVec(i), tracks.GetDirZVec(i));

    LightTrack_v &secondaries = td->fPhysicsData->GetSecondarySOA();
    for (int l = 0; l < kVecLenD; ++l) {
      int idx = secondaries.InsertTrack();
      secondaries.SetKinE(Get(gammaEnergy, l), idx);
      secondaries.SetDirX(Get(gamDirX, l), idx);
      secondaries.SetDirY(Get(gamDirY, l), idx);
      secondaries.SetDirZ(Get(gamDirZ, l), idx);
      secondaries.SetGVcode(fSecondaryInternalCode, idx);
      secondaries.SetMass(0.0, idx);
      secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx);
    }

    // compute the primary e-/e+ post interaction direction: from momentum vector conservation
    const Double_v elInitTotalMomentum = Math::Sqrt(ekin * (ekin + 2.0 * geant::units::kElectronMassC2));
    // final momentum of the e-/e+ in the lab frame
    Double_v elDirX = elInitTotalMomentum * tracks.GetDirXVec(i) - gammaEnergy * gamDirX;
    Double_v elDirY = elInitTotalMomentum * tracks.GetDirYVec(i) - gammaEnergy * gamDirY;
    Double_v elDirZ = elInitTotalMomentum * tracks.GetDirZVec(i) - gammaEnergy * gamDirZ;

    Double_v norm = 1.0 / Math::Sqrt(elDirX * elDirX + elDirY * elDirY + elDirZ * elDirZ);

    // update primary track direction
    for (int l = 0; l < kVecLenD; ++l) {
      tracks.SetDirX(Get(elDirX * norm, l), i + l);
      tracks.SetDirY(Get(elDirY * norm, l), i + l);
      tracks.SetDirZ(Get(elDirZ * norm, l), i + l);
      tracks.SetKinE(Get(ekin - gammaEnergy, l), i + l);
    }
  }
}

Double_v SeltzerBergerBremsModel::SampleEnergyTransfer(Double_v gammaCut, Double_v densityCor, IndexD_v mcLocalIdx,
                                                          double *tableEmin, double *tableILDeta, Double_v primekin,
                                                          Double_v r1, Double_v r2, Double_v r3)
{
  Double_v lPrimEkin    = Math::Log(primekin);
  Double_v logEmin      = vecCore::Gather<Double_v>(tableEmin, mcLocalIdx);
  Double_v ilDelta      = vecCore::Gather<Double_v>(tableILDeta, mcLocalIdx);
  Double_v val          = (lPrimEkin - logEmin) * ilDelta;
  IndexD_v indxPrimEkin = (IndexD_v)val; // lower electron energy bin index
  Double_v pIndxHigh    = val - indxPrimEkin;
  MaskD_v mask          = r1 < pIndxHigh;
  if (!MaskEmpty(mask)) {
    vecCore::MaskedAssign(indxPrimEkin, mask, indxPrimEkin + 1);
  }

  Double_v egammaV;
  for (int l = 0; l < kVecLenD; ++l) {
    assert(indxPrimEkin[l] >= 0);
    assert(indxPrimEkin[l] <= fSamplingTables[mcLocalIdx[l]]->fAliasData.size() - 1);
    //    LinAliasCached& als = fAliasData.fTablesPerMatCut[mcLocalIdx[l]]->fAliasData[indxPrimEkin[l]];
    //    double xi = AliasTableAlternative::SampleLinear(als,fSTNumSamplingElecEnergies,r2[l],r3[l]);

    const LinAlias *als = fSamplingTables[Get(mcLocalIdx, l)]->fAliasData[Get(indxPrimEkin, l)];
    const double egamma =
        fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]), &(als->fAliasIndx[0]),
                                    fSTNumSamplingPhotEnergies, Get(r2, l), Get(r3, l));
    Set(egammaV, l, egamma);
  }

  const Double_v dum1 = gammaCut * gammaCut + densityCor;
  const Double_v dum2 = (primekin * primekin + densityCor) / dum1;

  return Math::Sqrt(dum1 * Math::Exp(egammaV * Math::Log(dum2)) - densityCor);
}

void SeltzerBergerBremsModel::SampleEnergyTransfer(const double *eEkin, const double *gammaCut, const int *IZet,
                                                      const double *zet, double *gammaEn, const double *densityCorArr,
                                                      int N, const geant::TaskData *td)
{

  //  // assert(N>=kVecLenD)
  int currN = 0;
  MaskD_v lanesDone;
  IndexD_v idx;
  for (int l = 0; l < kVecLenD; ++l) {
    AssignMaskLane(lanesDone, l, false);
    Set(idx, l, currN++);
  }

  //
  while (currN < N || !MaskFull(lanesDone)) {
    Double_v eekin       = vecCore::Gather<Double_v>(eEkin, idx);
    Double_v gcut        = vecCore::Gather<Double_v>(gammaCut, idx);
    Double_v densityCorr = vecCore::Gather<Double_v>(densityCorArr, idx);

    const Double_v kappac = gcut / eekin;

    Double_v lekin  = Math::Log(eekin);
    Double_v eresid = (lekin - fLogLoadDCSMinElecEnergy) * fInvLogLoadDCSDeltaEnergy;
    IndexD_v ie     = (IndexD_v)eresid; // y1 index
    eresid -= (Double_v)ie;             // (y2-y1)*resid + y1
    assert((eekin < fLoadDCSElectronEnergyGrid[0]).isEmpty());
    assert((eekin >= fLoadDCSElectronEnergyGrid[fLoadDCSNumElectronEnergies - 1]).isEmpty());

    Double_v vmax;
    for (int l = 0; l < kVecLenD; ++l) {
      Set(vmax, l, 1.02 * GetDXSECValue(IZet[Get(idx, l)], Get(ie, l), Get(eresid, l), Get(kappac, l)));
      //
      // majoranta corrected vmax for e-
      constexpr double epeaklimit = 300.0 * geant::units::MeV;
      constexpr double elowlimit  = 20.0 * geant::units::keV;
      if (fIsElectron && Get(kappac, l) < 0.97 && ((Get(eekin, l) > epeaklimit) || (Get(eekin, l) < elowlimit))) {
        vmax = std::max(vmax, std::min(fXsecLimits[IZet[Get(idx, l)] - 1],
                                       1.1 * GetDXSECValue(IZet[Get(idx, l)], Get(ie, l), Get(eresid, l), 0.97)));
      }
    }
    MaskD_v tmp = kappac < 0.05;
    vecCore::MaskedAssign(vmax, tmp, vmax * 1.2);

    const Double_v minXi = Math::Log(gcut * gcut + densityCorr);
    const Double_v delXi = Math::Log(eekin * eekin + densityCorr) - minXi;

    Double_v rnd0 = td->fRndm->uniformV();

    Double_v egamma      = Math::Sqrt(Math::Max(Math::Exp(minXi + rnd0 * delXi) - densityCorr, (Double_v)0.));
    const Double_v kappa = egamma / eekin;
    Double_v val         = 0.0;
    for (int l = 0; l < kVecLenD; ++l) {
      Set(val, l, GetDXSECValue(IZet[Get(idx, l)], Get(ie, l), Get(eresid, l), Get(kappa, l)));
    }

    if (!fIsElectron) {
      Double_v zetV = vecCore::Gather<Double_v>(zet, idx);
      val *= PositronCorrection1(eekin, kappa, gcut, zetV);
    }

    Double_v rnd1    = td->fRndm->uniformV();
    MaskD_v accepted = val > vmax * rnd1;

    if (!MaskEmpty(accepted)) {
      vecCore::Scatter(egamma, gammaEn, idx);
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

void SeltzerBergerBremsModel::SamplePhotonDirection(Double_v elenergy, Double_v &sinTheta, Double_v &cosTheta,
                                                       Double_v rndm)
{
  const Double_v c = 4. - 8. * rndm;
  Double_v a       = c;
  Double_v signc   = 1.;
  MaskD_v tmp      = c < 0.0;
  vecCore::MaskedAssign(signc, tmp, (Double_v)-1.0);
  vecCore::MaskedAssign(a, tmp, -c);

  const Double_v delta = 0.5 * (Math::Sqrt(a * a + 4.) + a);
  //  delta += a;
  //  delta *= 0.5;

  const Double_v cofA = -signc * Math::Exp(Math::Log(delta) / 3.0);
  cosTheta            = cofA - 1. / cofA;

  const Double_v tau  = elenergy * geant::units::kInvElectronMassC2;
  const Double_v beta = Math::Sqrt(tau * (tau + 2.)) / (tau + 1.);

  cosTheta = (cosTheta + beta) / (1. + cosTheta * beta);
  // check cosTheta limit
  cosTheta = Math::Min((Double_v)1.0, cosTheta);
  // if (cosTheta>1.0) {
  //  cosTheta = 1.0;
  //}
  sinTheta = Math::Sqrt((1. - cosTheta) * (1. + cosTheta));
}

Double_v SeltzerBergerBremsModel::PositronCorrection1(Double_v ekinelectron, Double_v ephoton, Double_v gcutener,
                                                         Double_v z)
{
  const Double_v dum1 = geant::units::kTwoPi * geant::units::kFineStructConst;
  Double_v poscor     = 0.0;
  Double_v e1         = ekinelectron - gcutener; // here is the dif.
  Double_v ibeta1 = (e1 + geant::units::kElectronMassC2) / Math::Sqrt(e1 * (e1 + 2.0 * geant::units::kElectronMassC2));
  Double_v e2     = ekinelectron * (1.0 - ephoton);
  Double_v ibeta2 = (e2 + geant::units::kElectronMassC2) / Math::Sqrt(e2 * (e2 + 2.0 * geant::units::kElectronMassC2));
  Double_v ddum   = dum1 * z * (ibeta1 - ibeta2);
  MaskD_v tmp     = ddum < -12.0;
  vecCore::MaskedAssign(poscor, tmp, (Double_v)-12.0);
  vecCore::MaskedAssign(poscor, !tmp, Math::Exp(ddum));
  return poscor;
}

bool SeltzerBergerBremsModel::IsModelUsable(const MaterialCuts *matCut, double ekin)
{
  const double gammaCut = matCut->GetProductionCutsInEnergy()[0];
  return ekin < GetHighEnergyUsageLimit() && ekin > GetLowEnergyUsageLimit() && ekin > gammaCut;
}

} // namespace geantphysics
