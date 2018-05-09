
#include "Geant/RelativisticBremsModel.h"

// from amterial
#include "Geant/Types.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/Material.h"
#include "Geant/MaterialProperties.h"
#include "Geant/Element.h"
#include "Geant/ElementProperties.h"

#include "Geant/MaterialCuts.h"

#include "Geant/GLIntegral.h"
#include "Geant/AliasTable.h"

#include "Geant/PhysicsParameters.h"

#include "Geant/Gamma.h"
#include "Geant/Electron.h"
#include "Geant/Positron.h"

#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"

// from geantV
#include "Geant/TaskData.h"

#include <cmath>
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

// use these elastic and inelatic form factors for light elements instead of TFM
// under the complete screening approximation
// Tsai Table.B2.
const double RelativisticBremsModel::gFelLowZet[]   = {0.0, 5.310, 4.790, 4.740, 4.710, 4.680, 4.620, 4.570};
const double RelativisticBremsModel::gFinelLowZet[] = {0.0, 6.144, 5.621, 5.805, 5.924, 6.012, 5.891, 5.788};

const double RelativisticBremsModel::gLPMFactor = 0.25 * geant::units::kFineStructConst *
                                                  geant::units::kElectronMassC2 * geant::units::kElectronMassC2 /
                                                  (geant::units::kPi * geant::units::kHBarPlanckCLight);
// this is k_p^2 / E_{t-electron}^2
const double RelativisticBremsModel::gDensityFactor = 4.0 * geant::units::kPi * geant::units::kClassicElectronRadius *
                                                      geant::units::kRedElectronComptonWLenght *
                                                      geant::units::kRedElectronComptonWLenght;

RelativisticBremsModel::LPMFuncs RelativisticBremsModel::gLPMFuncs;
std::vector<RelativisticBremsModel::ElementData *> RelativisticBremsModel::gElementData(gMaxZet, nullptr);

RelativisticBremsModel::RelativisticBremsModel(const std::string &modelname) : EMModel(modelname)
{
  fIsUseLPM              = true; // use LPM suppression by default
  fNGL                   = 64;
  fSecondaryInternalCode = -1; // internal gamma particle index (set at init.)
  // these are used only if sampling tables were requested
  fSTNumElectronEnergyPerDecade = 8; // number of sampling tables per decade (e-/e+ kinetic energy)
  fSTNumSamplingPhotEnergies    = 100;
  fAliasSampler                 = nullptr; // will be set at initialisation
  fGL                           = nullptr;
}

RelativisticBremsModel::~RelativisticBremsModel()
{
  // clear ElementData
  for (size_t i = 0; i < gElementData.size(); ++i) {
    if (gElementData[i]) {
      delete gElementData[i];
    }
  }
  gElementData.clear();
  //
  // clear LPMFunctions (if any)
  if (fIsUseLPM) {
    gLPMFuncs.fLPMFuncG.clear();
    gLPMFuncs.fLPMFuncPhi.clear();
  }
  //
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

void RelativisticBremsModel::Initialize()
{
  EMModel::Initialize();
  InitElementData();
  fSecondaryInternalCode = Gamma::Definition()->GetInternalCode();
  if (fIsUseLPM) {
    InitLPMFunctions(); // build table for fast LPM G(s), Phi(s) evaluation
  }
  if (!fGL) {
    fGL = new GLIntegral(fNGL, 0., 1.);
  }
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

double RelativisticBremsModel::ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle *, bool istotal)
{
  const Material *mat = matcut->GetMaterial();
  double gammacut     = (matcut->GetProductionCutsInEnergy())[0];
  if (istotal) {
    // for the total stopping power we just need a gamma production cut >=kinenergy
    gammacut = 1.01 * kinenergy;
  }
  return ComputeDEDXPerVolume(mat, gammacut, kinenergy);
}

double RelativisticBremsModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy,
                                                          const Particle *)
{
  double xsec = 0.0;
  if (kinenergy < GetLowEnergyUsageLimit() || kinenergy > GetHighEnergyUsageLimit()) {
    return xsec;
  }
  const Material *mat   = matcut->GetMaterial();
  const double gammacut = (matcut->GetProductionCutsInEnergy())[0];
  xsec                  = ComputeXSectionPerVolume(mat, gammacut, kinenergy);
  return xsec;
}

double RelativisticBremsModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut, double kinenergy,
                                                      const Particle *)
{
  double xsec = 0.0;
  if (kinenergy < GetLowEnergyUsageLimit() || kinenergy > GetHighEnergyUsageLimit()) {
    return xsec;
  }
  const Material *mat   = matcut->GetMaterial();
  const double gammacut = (matcut->GetProductionCutsInEnergy())[0];
  xsec                  = ComputeXSectionPerAtom(elem, mat, gammacut, kinenergy);
  return xsec;
}

template <typename R>
void RelativisticBremsModel::GetLPMFunctions(R &lpmGs, R &lpmPhis, const R s)
{
  vecCore::Mask<R> tmp = s < gLPMFuncs.fSLimit;

  R val                  = s / gLPMFuncs.fSDelta;
  vecCore::Index<R> ilow = vecCore::Index<R>(val);
  val -= ilow;
  for (int l = 0; l < (int)vecCore::VectorSize<R>(); ++l) {
    if (!Get(tmp, l)) Set(ilow, l, gLPMFuncs.fLPMFuncPhi.size() - 2); // above limit
  }
  R lpmG_lowp1   = vecCore::Gather<R>(gLPMFuncs.fLPMFuncG.data(), ilow + (vecCore::Index<R>)1);
  R lpmG_low     = vecCore::Gather<R>(gLPMFuncs.fLPMFuncG.data(), ilow);
  R lpmPhi_lowp1 = vecCore::Gather<R>(gLPMFuncs.fLPMFuncPhi.data(), ilow + (vecCore::Index<R>)1);
  R lpmPhi_low   = vecCore::Gather<R>(gLPMFuncs.fLPMFuncPhi.data(), ilow);
  lpmGs          = (lpmG_lowp1 - lpmG_low) * val + lpmG_low;
  lpmPhis        = (lpmPhi_lowp1 - lpmPhi_low) * val + lpmPhi_low;

  if (!MaskEmpty(!tmp)) {
    R ss = s * s;
    ss *= ss;
    vecCore::MaskedAssign(lpmPhis, !tmp, 1.0 - 0.01190476 / ss);
    vecCore::MaskedAssign(lpmGs, !tmp, 1.0 - 0.0230655 / ss);
  }
}

template <>
void RelativisticBremsModel::GetLPMFunctions(double &lpmGs, double &lpmPhis, const double s)
{
  if (s < gLPMFuncs.fSLimit) {
    double val     = s / gLPMFuncs.fSDelta;
    const int ilow = int(val);
    val -= ilow;
    lpmGs   = (gLPMFuncs.fLPMFuncG[ilow + 1] - gLPMFuncs.fLPMFuncG[ilow]) * val + gLPMFuncs.fLPMFuncG[ilow];
    lpmPhis = (gLPMFuncs.fLPMFuncPhi[ilow + 1] - gLPMFuncs.fLPMFuncPhi[ilow]) * val + gLPMFuncs.fLPMFuncPhi[ilow];
  } else {
    double ss = s * s;
    ss *= ss;
    lpmPhis = 1.0 - 0.01190476 / ss;
    lpmGs   = 1.0 - 0.0230655 / ss;
  }
}

int RelativisticBremsModel::SampleSecondaries(LightTrack &track, geant::TaskData *td)
{
  int numSecondaries         = 0;
  const double ekin          = track.GetKinE();
  const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(track.GetMaterialCutCoupleIndex());
  const double gammacut      = (matCut->GetProductionCutsInEnergy())[0];
  // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
  // check if kinetic energy is above fHighEnergyUsageLimit andd o nothing if yes
  // check if kinetic energy is below gamma production cut and do nothing if yes
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
    gammaEnergy = SamplePhotonEnergy(matCut, ekin, td);
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

double RelativisticBremsModel::MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle *) const
{
  double mine = (matcut->GetProductionCutsInEnergy())[0]; // gamma production cut in the given material-cuts
  return std::max(mine, GetLowEnergyUsageLimit());
}

// private methods

void RelativisticBremsModel::InitElementData()
{
  size_t numMatCuts = MaterialCuts::GetTheMaterialCutsTable().size();
  // get list of active region
  std::vector<bool> &isActiveInRegion = GetListActiveRegions();
  for (size_t imc = 0; imc < numMatCuts; ++imc) {
    const MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[imc];
    // if this MaterialCuts belongs to a region where this model is active:
    if (isActiveInRegion[matCut->GetRegionIndex()]) {
      // get the list of elements
      const Vector_t<Element *> &theElemVect = matCut->GetMaterial()->GetElementVector();
      size_t numElems                        = theElemVect.size();
      for (size_t ie = 0; ie < numElems; ++ie) {
        const Element *elem = theElemVect[ie];
        double zet          = elem->GetZ();
        int izet            = std::min(std::lrint(zet), gMaxZet - 1);
        if (!gElementData[izet]) {
          ElementData *elemData = new ElementData();
          const double fc       = elem->GetElementProperties()->GetCoulombCorrection();
          double Fel            = 1.;
          double Finel          = 1.;
          elemData->fLogZ       = Math::Log(zet);
          elemData->fFz         = elemData->fLogZ / 3. + fc;
          if (izet < 5) {
            Fel   = gFelLowZet[izet];
            Finel = gFinelLowZet[izet];
          } else {
            Fel   = Math::Log(184.15) - elemData->fLogZ / 3.;
            Finel = Math::Log(1194.) - 2. * elemData->fLogZ / 3.;
          }
          elemData->fZFactor1      = (Fel - fc) + Finel / zet;
          elemData->fZFactor2      = (1. + 1. / zet) / 12.;
          elemData->fVarS1         = Math::Pow(zet, 2. / 3.) / (184.15 * 184.15);
          elemData->fILVarS1Cond   = 1. / (Math::Log(std::sqrt(2.0) * elemData->fVarS1));
          elemData->fILVarS1       = 1. / Math::Log(elemData->fVarS1);
          elemData->fGammaFactor   = 100.0 * geant::units::kElectronMassC2 / Math::Pow(zet, 1. / 3.);
          elemData->fEpsilonFactor = 100.0 * geant::units::kElectronMassC2 / Math::Pow(zet, 2. / 3.);
          gElementData[izet]       = elemData;
        }
      }
    }
  }
}

/**
 *   The restricted atomic cross section for bremsstrahlung photon emeission for target element with atomic number
 *   \f$Z\f$, gamma photon production energy threshold \f$k_c\f$ is computed for e-/e+ kinetic energy \f$E\f$ according
 * to
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
template <typename R>
R RelativisticBremsModel::ComputeDXSecPerAtom(R egamma, R etotal, R zet)
{
  // constexpr R factor =
  // 16.*geant::units::kFineStructConst*geant::units::kClassicElectronRadius*geant::units::kClassicElectronRadius/3.;
  R dcs         = 0.;
  const R y     = egamma / etotal;
  const R onemy = 1. - y;
  R zFactor1, zFactor2, fz, logZ, gammaFactor, epsilonFactor;
  vecCore::Index<R> izetV;
  for (int l = 0; l < (int)vecCore::VectorSize<R>(); ++l) {
    int izet = std::lrint(Get(zet, l));
    Set(zFactor1, l, gElementData[izet]->fZFactor1);
    Set(zFactor2, l, gElementData[izet]->fZFactor2);
    Set(fz, l, gElementData[izet]->fFz);
    Set(logZ, l, gElementData[izet]->fLogZ);
    Set(gammaFactor, l, gElementData[izet]->fGammaFactor);
    Set(epsilonFactor, l, gElementData[izet]->fEpsilonFactor);
    Set(izetV, l, izet);
  }
  vecCore::Mask<R> smallZ = izetV < 5;
  if (!MaskEmpty(smallZ)) {
    vecCore::MaskedAssign(dcs, smallZ, (onemy + 0.75 * y * y) * zFactor1 + onemy * zFactor2);
  }
  if (!MaskEmpty(!smallZ)) {
    // Tsai: screening from Thomas-Fermi model of atom; Tsai Eq.(3.82)
    // variables gamma and epsilon from Tsai Eq.(3.30) and Eq.(3.31)
    const R invZ    = 1. / zet;
    const R dum0    = y / (etotal - egamma);
    const R gamma   = dum0 * gammaFactor;
    const R epsilon = dum0 * epsilonFactor;
    R phi1, phi1m2, xsi1, xsi1m2;
    ComputeScreeningFunctions(phi1, phi1m2, xsi1, xsi1m2, gamma, epsilon);
    vecCore::MaskedAssign(dcs, !smallZ,
                          (onemy + 0.75 * y * y) * ((0.25 * phi1 - fz) + (0.25 * xsi1 - 2. * logZ / 3.) * invZ) +
                              0.125 * onemy * (phi1m2 + xsi1m2 * invZ));
    // dcs *= factor*zet*zet/egamma;
  }

  return Math::Max(dcs, (R)0.0);
}

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
 * \f$k\f$ is the emitted photon energy, \f$y=k/E_t\f$ is the emitted photon energy in pre-interaction \f$e^-/e^+\f$
 * total
 *  energy (\f$ E_t\f$) units, \f$Z\f$ is the target atomic number, \f$f\f$ is the Coulomb correction \cite
 * davies1954theory
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
 *  24s^3\sum_{j=1}^{\infty} \frac{1}{(j+s)^2+s^2} = 12 s^2 i \left\{\psi^{(0)}[(1-i)s+1] - \psi^{(0)}[(1+i)s+1]\right\}
 * \\
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
 *   3\left\{1-\exp \left[-4s  - \frac{8s^2}{1+3.936s+4.97s^2-0.05s^3+7.5s^4} \right] \right\} -2\phi(s)  & \quad
 * \text{if}\; 0.1 < s < 2 \\
 *   1 - 0.0230655/s^4  & \quad \text{if}\; s \geq 2
 *  \end{cases}
 * \f]
 * where instead of \cite stanev1982development we used the series expansion around \f$s\f$ zero and
 * infinity respectively
 * \f[
 *  48s^3\sum_{j=0}^{\infty} \frac{1}{(j+s+1/2)^2+s^2} = 24 s^2 i \left\{ \psi^{0}((1-i)s+0.5) - \psi^{0}((1+i)s+0.5)
 * \right\} \approx \\
 *  \begin{cases}
 *   \approx 24\pi^2 s^3 + 48 \psi^{2}(0.5) s^4 & \quad \mathrm{around}\; s=0 \\
 *   \approx 12\pi s^2 -1 + 0.0230655/s^4   & \quad \mathrm{around}\; s=\infty
 *  \end{cases}
 * \f]
 * where \f$\psi^{2}(x)\f$ is the second derivative of the digamma function at 0.5 and \f$ 48 \psi^{2}(0.5) \approx
 * -807.782238923247359\f$.
 *
 * The variable \f$s\f$, (can be obtained from the LPM condition i.e. the relation of average photon emission angle to
 * the
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
double RelativisticBremsModel::ComputeURelDXSecPerAtom(double egamma, double etotal, double lpmenergy,
                                                       double densitycor, int izet)
{
  const double y     = egamma / etotal;
  const double onemy = 1. - y;
  const double dum0  = 0.25 * y * y;
  double funcGS, funcPhiS, funcXiS;
  ComputeLPMfunctions(funcXiS, funcGS, funcPhiS, lpmenergy, egamma, etotal, densitycor, izet);
  double dcs = funcXiS * (dum0 * funcGS + (onemy + 2.0 * dum0) * funcPhiS) * gElementData[izet]->fZFactor1 +
               onemy * gElementData[izet]->fZFactor2;
  return std::max(dcs, 0.0);
}

Double_v RelativisticBremsModel::ComputeURelDXSecPerAtom(Double_v egamma, Double_v etotal, Double_v lpmenergy,
                                                         Double_v densitycor, std::array<int, kVecLenD> izet)
{
  const Double_v y     = egamma / etotal;
  const Double_v onemy = 1. - y;
  const Double_v dum0  = 0.25 * y * y;
  Double_v funcGS, funcPhiS, funcXiS;
  ComputeLPMfunctions(funcXiS, funcGS, funcPhiS, lpmenergy, egamma, etotal, densitycor, izet);
  Double_v dcs;
  Double_v ZFactor1, ZFactor2;
  for (int l = 0; l < kVecLenD; ++l) {
    Set(ZFactor1, l, gElementData[izet[l]]->fZFactor1);
    Set(ZFactor2, l, gElementData[izet[l]]->fZFactor2);
  }
  dcs = funcXiS * (dum0 * funcGS + (onemy + 2.0 * dum0) * funcPhiS) * ZFactor1 + onemy * ZFactor2;
  return Math::Max(dcs, (Double_v)0.0);
}

template <typename R>
void RelativisticBremsModel::ComputeScreeningFunctions(R &phi1, R &phi1m2, R &xsi1, R &xsi1m2, const R gamma,
                                                       const R epsilon)
{
  const R gamma2 = gamma * gamma;
  phi1 =
      16.863 - 2.0 * Math::Log(1.0 + 0.311877 * gamma2) + 2.4 * Math::Exp(-0.9 * gamma) + 1.6 * Math::Exp(-1.5 * gamma);
  phi1m2           = 2.0 / (3.0 + 19.5 * gamma + 18.0 * gamma2); // phi1-phi2
  const R epsilon2 = epsilon * epsilon;
  xsi1             = 24.34 - 2.0 * Math::Log(1.0 + 13.111641 * epsilon2) + 2.8 * Math::Exp(-8.0 * epsilon) +
         1.2 * Math::Exp(-29.2 * epsilon);
  xsi1m2 = 2.0 / (3.0 + 120.0 * epsilon + 1200.0 * epsilon2); // xsi1-xsi2
}

void RelativisticBremsModel::ComputeLPMfunctions(double &funcXiS, double &funcGS, double &funcPhiS,
                                                 const double lpmenergy, const double egamma, const double etot,
                                                 const double densitycor, const int izet)
{
  static const double sqrt2 = std::sqrt(2.);
  const double redegamma    = egamma / etot;
  // const double varSprime = std::sqrt(0.125*redegamma/(1.0-redegamma)*lpmenergy/etot);
  const double varSprime = std::sqrt(0.125 * redegamma * lpmenergy / ((1.0 - redegamma) * etot));
  const double varS1     = gElementData[izet]->fVarS1;
  const double condition = sqrt2 * varS1;
  double funcXiSprime    = 2.0;
  if (varSprime > 1.0) {
    funcXiSprime = 1.0;
  } else if (varSprime > condition) {
    const double funcHSprime = Math::Log(varSprime) * gElementData[izet]->fILVarS1Cond;
    funcXiSprime             = 1.0 + funcHSprime -
                   0.08 * (1.0 - funcHSprime) * funcHSprime * (2.0 - funcHSprime) * gElementData[izet]->fILVarS1Cond;
  }
  const double varS = varSprime / std::sqrt(funcXiSprime);
  // - include dielectric suppression effect into s according to Migdal
  const double varShat = varS * (1.0 + densitycor / (egamma * egamma));
  funcXiS              = 2.0;
  if (varShat > 1.0) {
    funcXiS = 1.0;
  } else if (varShat > varS1) {
    funcXiS = 1.0 + Math::Log(varShat) * gElementData[izet]->fILVarS1;
  }
  GetLPMFunctions(funcGS, funcPhiS, varShat);
  // ComputeLPMGsPhis(funcGS, funcPhiS, varShat);
  //
  // MAKE SURE SUPPRESSION IS SMALLER THAN 1: due to Migdal's approximation on xi
  if (funcXiS * funcPhiS > 1. || varShat > 0.57) {
    funcXiS = 1. / funcPhiS;
  }
}

void RelativisticBremsModel::ComputeLPMfunctions(Double_v &funcXiS, Double_v &funcGS, Double_v &funcPhiS,
                                                 const Double_v lpmenergy, const Double_v egamma, const Double_v etot,
                                                 const Double_v densitycor, const std::array<int, kVecLenD> izet)
{
  static const Double_v sqrt2 = Math::Sqrt(2.);
  const Double_v redegamma    = egamma / etot;
  // const Double_v varSprime = std::sqrt(0.125*redegamma/(1.0-redegamma)*lpmenergy/etot);
  const Double_v varSprime = Math::Sqrt(0.125 * redegamma * lpmenergy / ((1.0 - redegamma) * etot));
  Double_v varS1, ILVarS1Cond, ILVarS1;
  for (int l = 0; l < kVecLenD; ++l) {
    Set(varS1, l, gElementData[izet[l]]->fVarS1);
    Set(ILVarS1Cond, l, gElementData[izet[l]]->fILVarS1Cond);
    Set(ILVarS1, l, gElementData[izet[l]]->fILVarS1);
  }
  const Double_v condition = sqrt2 * varS1;
  Double_v funcXiSprime    = 2.0;

  MaskD_v tmp1, tmp2;
  tmp1 = varSprime > 1.0;
  tmp2 = varSprime > condition;
  vecCore::MaskedAssign(funcXiSprime, tmp1, (Double_v)1.0);
  const Double_v funcHSprime = Math::Log(varSprime) * ILVarS1Cond;
  Double_v tmpFuncXiSprime =
      1.0 + funcHSprime - 0.08 * (1.0 - funcHSprime) * funcHSprime * (2.0 - funcHSprime) * ILVarS1Cond;
  vecCore::MaskedAssign(funcXiSprime, !tmp1 && tmp2, tmpFuncXiSprime);

  const Double_v varS = varSprime / Math::Sqrt(funcXiSprime);
  // - include dielectric suppression effect into s according to Migdal
  const Double_v varShat = varS * (1.0 + densitycor / (egamma * egamma));
  funcXiS                = 2.0;
  tmp1                   = varShat > 1.0;
  tmp2                   = varShat > varS1;
  vecCore::MaskedAssign(funcXiS, tmp1, (Double_v)1.0);
  vecCore::MaskedAssign(funcXiS, !tmp1 && tmp2, 1.0 + Math::Log(varShat) * ILVarS1);
  GetLPMFunctions(funcGS, funcPhiS, varShat);
  // ComputeLPMGsPhis(funcGS, funcPhiS, varShat);
  //
  // MAKE SURE SUPPRESSION IS SMALLER THAN 1: due to Migdal's approximation on xi
  tmp1 = funcXiS * funcPhiS > 1.0 || varShat > 0.57;
  vecCore::MaskedAssign(funcXiS, tmp1, 1.0 / funcPhiS);
}

void RelativisticBremsModel::ComputeLPMGsPhis(double &funcGS, double &funcPhiS, const double varShat)
{
  if (varShat < 0.01) {
    funcPhiS = 6.0 * varShat * (1.0 - geant::units::kPi * varShat);
    funcGS   = 12.0 * varShat - 2.0 * funcPhiS;
  } else {
    double varShat2 = varShat * varShat;
    double varShat3 = varShat * varShat2;
    double varShat4 = varShat2 * varShat2;
    if (varShat < 0.415827397755) { // use Stanev approximation: for \psi(s) and compute G(s)
      funcPhiS = 1.0 - Math::Exp(-6.0 * varShat * (1.0 + varShat * (3.0 - geant::units::kPi)) +
                                 varShat3 / (0.623 + 0.796 * varShat + 0.658 * varShat2));
      // 1-\exp \left\{  -4s  - \frac{8s^2}{1+3.936s+4.97s^2-0.05s^3+7.5s^4} \right\}
      const double funcPsiS =
          1.0 -
          Math::Exp(-4.0 * varShat -
                    8.0 * varShat2 / (1.0 + 3.936 * varShat + 4.97 * varShat2 - 0.05 * varShat3 + 7.5 * varShat4));
      // G(s) = 3 \psi(s) - 2 \phi(s)
      funcGS = 3.0 * funcPsiS - 2.0 * funcPhiS;
    } else if (varShat < 1.55) {
      funcPhiS = 1.0 - Math::Exp(-6.0 * varShat * (1.0 + varShat * (3.0 - geant::units::kPi)) +
                                 varShat3 / (0.623 + 0.796 * varShat + 0.658 * varShat2));
      const double dum0 = -0.16072300849123999 + 3.7550300067531581 * varShat - 1.7981383069010097 * varShat2 +
                          0.67282686077812381 * varShat3 - 0.1207722909879257 * varShat4;
      funcGS = std::tanh(dum0);
    } else {
      funcPhiS = 1.0 - 0.01190476 / varShat4;
      if (varShat < 1.9156) {
        const double dum0 = -0.16072300849123999 + 3.7550300067531581 * varShat - 1.7981383069010097 * varShat2 +
                            0.67282686077812381 * varShat3 - 0.1207722909879257 * varShat4;
        funcGS = std::tanh(dum0);
      } else {
        funcGS = 1.0 - 0.0230655 / varShat4;
      }
    }
  }
}

void RelativisticBremsModel::InitLPMFunctions()
{
  if (!gLPMFuncs.fIsInitialized) {
    //    gLPMFuncs.fSLimit = 2.;
    //    gLPMFuncs.fSDelta = 0.01;
    const int num = gLPMFuncs.fSLimit / gLPMFuncs.fSDelta + 1;
    gLPMFuncs.fLPMFuncG.resize(num);
    gLPMFuncs.fLPMFuncPhi.resize(num);
    for (int i = 0; i < num; ++i) {
      const double s = i * gLPMFuncs.fSDelta;
      ComputeLPMGsPhis(gLPMFuncs.fLPMFuncG[i], gLPMFuncs.fLPMFuncPhi[i], s);
    }
    gLPMFuncs.fIsInitialized = true;
  }
}

// clear all sampling tables (if any)
void RelativisticBremsModel::ClearSamplingTables()
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

void RelativisticBremsModel::InitSamplingTables()
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

void RelativisticBremsModel::BuildSamplingTableForMaterialCut(const MaterialCuts *matcut, int indxlocal)
{
  const Material *mat        = matcut->GetMaterial();
  const double gcut          = (matcut->GetProductionCutsInEnergy())[0];
  const double minElecEnergy = MinimumPrimaryEnergy(matcut, nullptr);
  const double maxElecEnergy = GetHighEnergyUsageLimit();
  if (minElecEnergy >= maxElecEnergy) {
    return;
  }
  //
  //
  const double densityFactor = gDensityFactor * mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol();
  const double lpmEnergy     = gLPMFactor * mat->GetMaterialProperties()->GetRadiationLength();
  // material dependent constant: the e-/e+ energy at which the corresponding k_LPM=k_p
  const double energyThLPM = std::sqrt(densityFactor) * lpmEnergy;
  //
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
    double eekin = Math::Exp(dataMatCut->fLogEmin + ie * delta);
    if (ie == 0 && minElecEnergy == gcut) {
      eekin = minElecEnergy + 1. * geant::units::eV; // would be zero otherwise
    }
    if (ie == numElecEnergies - 1) {
      eekin = maxElecEnergy;
    }
    const double etot       = eekin + geant::units::kElectronMassC2;
    const double densityCor = densityFactor * etot * etot; // this is k_p^2
    const bool isLPM        = (fIsUseLPM && etot > energyThLPM);
    //
    const double dumc0 = gcut * gcut + densityCor;
    const double dumc1 = Math::Log((eekin * eekin + densityCor) / dumc0);
    // create the alias data struct
    LinAlias *als = new LinAlias(fSTNumSamplingPhotEnergies);
    // fill 3 values at 0,0.8,1 of the transformed variable
    int numdata    = 3;
    als->fXdata[0] = 0.0;
    als->fXdata[1] = 0.8;
    als->fXdata[2] = 1.0;
    for (int i = 0; i < numdata; ++i) {
      double egamma = gcut;
      if (i == 0) {
        egamma = gcut;
      } else if (i == 2) {
        egamma = eekin;
      } else {
        egamma = std::sqrt(dumc0 * Math::Exp(als->fXdata[i] * dumc1) - densityCor);
      }
      double dcross = 0.0;
      for (int ielem = 0; ielem < numElems; ++ielem) {
        double zet = theElements[ielem]->GetZ();
        double val = 0.0;
        if (isLPM) {
          int izet = std::min(std::lrint(zet), gMaxZet - 1);
          val      = ComputeURelDXSecPerAtom(egamma, etot, lpmEnergy, densityCor, izet);
        } else {
          val = ComputeDXSecPerAtom(egamma, etot, zet);
        }
        dcross += theAtomicNumDensityVector[ielem] * zet * zet * val;
      }
      als->fYdata[i] = dcross;
    }
    // fill in up to fSTNumSamplingPhotEnergies
    while (numdata < fSTNumSamplingPhotEnergies) {
      // find the lower index of the bin, where we have the biggest linear interp. error compared to computed value
      double maxerr  = 0.0;
      double thexval = 0.0;
      double theyval = 0.0;
      int maxerrindx = 0; // the lower index of the corresponding bin
      for (int i = 0; i < numdata - 1; ++i) {
        const double xx     = 0.5 * (als->fXdata[i] + als->fXdata[i + 1]); // mid point
        const double yy     = 0.5 * (als->fYdata[i] + als->fYdata[i + 1]); // lin func val at the mid point
        const double egamma = std::sqrt(dumc0 * Math::Exp(xx * dumc1) - densityCor);
        double dcross       = 0.0;
        for (int ielem = 0; ielem < numElems; ++ielem) {
          double zet = theElements[ielem]->GetZ();
          double val = 0.0;
          if (isLPM) {
            int izet = std::min(std::lrint(zet), gMaxZet - 1);
            val      = ComputeURelDXSecPerAtom(egamma, etot, lpmEnergy, densityCor, izet);
          } else {
            val = ComputeDXSecPerAtom(egamma, etot, zet);
          }
          dcross += theAtomicNumDensityVector[ielem] * zet * zet * val;
        }
        double err = std::fabs((yy - dcross) * (1.0 - yy / dcross));
        if (err > maxerr) {
          maxerr     = err;
          maxerrindx = i;
          thexval    = xx;
          theyval    = dcross;
        }
      }
      // extend x,y data by puting a value at the mid point of the highest error bin
      // first shift all values to the right then insert x mid point
      for (int j = numdata; j > maxerrindx + 1; --j) {
        als->fXdata[j] = als->fXdata[j - 1];
        als->fYdata[j] = als->fYdata[j - 1];
      }
      als->fXdata[maxerrindx + 1] = thexval;
      als->fYdata[maxerrindx + 1] = theyval;
      // increase number of data
      ++numdata;
    }
    //
    fAliasSampler->PreparLinearTable(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]), &(als->fAliasIndx[0]),
                                     fSTNumSamplingPhotEnergies);
    fSamplingTables[indxlocal]->fAliasData[ie] = als;
  }
}

double RelativisticBremsModel::SamplePhotonEnergy(const MaterialCuts *matcut, double eekin, double r1, double r2,
                                                  double r3)
{
  const double gcut = (matcut->GetProductionCutsInEnergy())[0];
  const double etot = eekin + geant::units::kElectronMassC2;
  const double dum0 = gDensityFactor * matcut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol();
  const double densityCor = etot * etot * dum0; // this is k_p^2
  const int mcIndxLocal   = fGlobalMatGCutIndxToLocal[matcut->GetIndex()];
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
  double gammae       = fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]),
                                              &(als->fAliasIndx[0]), fSTNumSamplingPhotEnergies, r2, r3);
  // transform back to gamma energy
  const double dum1 = gcut * gcut + densityCor;
  const double dum2 = (eekin * eekin + densityCor) / dum1;

  return std::sqrt(dum1 * Math::Exp(gammae * Math::Log(dum2)) - densityCor);
}

Double_v RelativisticBremsModel::SampleEnergyTransfer(Double_v gammaCut, Double_v densityCor, IndexD_v mcLocalIdx,
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
    auto &als     = fAliasData.fTablesPerMatCut[Get(mcLocalIdx, l)]->fAliasData[Get(indxPrimEkin, l)];
    double egamma = als.Sample(Get(r2, l), Get(r3, l));

    //    Standard version
    //    const LinAlias *als = fSamplingTables[Get(mcLocalIdx, l)]->fAliasData[Get(indxPrimEkin, l)];
    //    const double egamma =
    //        fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]),
    //        &(als->fAliasIndx[0]),
    //                                    fSTNumSamplingPhotEnergies, Get(r2, l), Get(r3, l));
    Set(egammaV, l, egamma);
  }

  const Double_v dum1 = gammaCut * gammaCut + densityCor;
  const Double_v dum2 = (primekin * primekin + densityCor) / dum1;

  return Math::Sqrt(dum1 * Math::Exp(egammaV * Math::Log(dum2)) - densityCor);
}

double RelativisticBremsModel::SamplePhotonEnergy(const MaterialCuts *matcut, double eekin, geant::TaskData *td)
{
  const double gcut       = (matcut->GetProductionCutsInEnergy())[0];
  const Material *mat     = matcut->GetMaterial();
  const double etot       = eekin + geant::units::kElectronMassC2;
  const double dumc0      = gDensityFactor * mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol();
  const double densityCor = etot * etot * dumc0; // this is k_p^2
  const double lpmEnergy  = mat->GetMaterialProperties()->GetRadiationLength() * gLPMFactor;
  // material dependent constant: the e-/e+ energy at which the corresponding k_LPM=k_p
  const double energyThLPM = std::sqrt(dumc0) * lpmEnergy;
  const bool isLPM         = (fIsUseLPM && etot > energyThLPM);

  // sample target element
  const Vector_t<Element *> &theElements = mat->GetElementVector();
  double targetElemIndx                  = 0;
  if (theElements.size() > 1) {
    targetElemIndx = SampleTargetElementIndex(matcut, eekin, td->fRndm->uniform());
  }
  // sample the transformed variable: u(k) = ln(k^2+k_p^2) that is in [ln(k_c^2+k_p^2), ln(E_k^2+k_p^2)]
  const double minVal   = Math::Log(gcut * gcut + densityCor);
  const double valRange = Math::Log(eekin * eekin + densityCor) - minVal;
  const double zet      = theElements[targetElemIndx]->GetZ();
  const int izet        = std::min(std::lrint(zet), gMaxZet - 1);
  const double funcMax  = gElementData[izet]->fZFactor1 + gElementData[izet]->fZFactor2;
  double egamma, funcVal;

  double *rndArray = td->fDblArray;
  do {
    td->fRndm->uniform_array(2, rndArray);
    egamma  = std::sqrt(std::max(Math::Exp(minVal + rndArray[0] * valRange) - densityCor, 0.));
    funcVal = (isLPM) ? ComputeURelDXSecPerAtom(egamma, etot, lpmEnergy, densityCor, izet)
                      : ComputeDXSecPerAtom(egamma, etot, zet);
  } while (funcVal < rndArray[1] * funcMax);

  return egamma;
}

void RelativisticBremsModel::SampleEnergyTransfer(const double *eEkin, const double *gammaCut, const double *zetArr,
                                                  const double *densityCorConstArr, const double *lpmEnergyArr,
                                                  double *gammaEn, int N, const geant::TaskData *td)
{

  // assert(N>=kVecLenD)
  // Init variables for lane refiling bookkeeping
  int currN = 0;
  MaskD_v lanesDone;
  IndexD_v idx;
  for (int l = 0; l < kVecLenD; ++l) {
    AssignMaskLane(lanesDone, l, false);
    Set(idx, l, currN++);
  }

  while (currN < N || !MaskFull(lanesDone)) {
    // Gather data for one round of sampling
    Double_v eekin            = vecCore::Gather<Double_v>(eEkin, idx);
    Double_v gcut             = vecCore::Gather<Double_v>(gammaCut, idx);
    Double_v densityCorrConst = vecCore::Gather<Double_v>(densityCorConstArr, idx);
    Double_v lpmEnergy        = vecCore::Gather<Double_v>(lpmEnergyArr, idx);
    Double_v zet              = vecCore::Gather<Double_v>(zetArr, idx);

    Double_v etot        = eekin + geant::units::kElectronMassC2;
    Double_v densityCorr = densityCorrConst * etot * etot;

    Double_v energyThLPM = Math::Sqrt(densityCorrConst) * lpmEnergy;
    MaskD_v isLPM        = MaskD_v(fIsUseLPM) && (etot > energyThLPM);

    Double_v minVal   = Math::Log(gcut * gcut + densityCorr);
    Double_v valRange = Math::Log(eekin * eekin + densityCorr) - minVal;
    Double_v funcMax;
    std::array<int, kVecLenD> izetV;
    for (int l = 0; l < kVecLenD; ++l) {
      const int izet = std::min(std::lrint(Get(zet, l)), gMaxZet - 1);
      Set(funcMax, l, gElementData[izet]->fZFactor1 + gElementData[izet]->fZFactor2);
      izetV[l] = izet;
    }

    Double_v egamma =
        Math::Sqrt(Math::Max(Math::Exp(minVal + td->fRndm->uniformV() * valRange) - densityCorr, (Double_v)0.));
    Double_v funcVal;
    if (!MaskEmpty(isLPM)) {
      vecCore::MaskedAssign(funcVal, isLPM, ComputeURelDXSecPerAtom(egamma, etot, lpmEnergy, densityCorr, izetV));
    }
    if (!MaskEmpty(!isLPM)) {
      vecCore::MaskedAssign(funcVal, !isLPM, ComputeDXSecPerAtom(egamma, etot, zet));
    }

    MaskD_v accepted = funcVal > td->fRndm->uniformV() * funcMax;

    // Scatter data sampled in this round
    if (!MaskEmpty(accepted)) {
      vecCore::Scatter(egamma, gammaEn, idx);
    }

    // Refill index lanes with next value if there is some particles to sample
    // if there is no work left - fill index with sentinel value(all input arrays should contain sentinel values after
    // end)
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

// the simple DipBustgenerator
template <typename R>
void RelativisticBremsModel::SamplePhotonDirection(R elenergy, R &sinTheta, R &cosTheta, R rndm)
{
  const R c            = 4. - 8. * rndm;
  R a                  = c;
  R signc              = 1.;
  vecCore::Mask<R> tmp = c < 0.0;
  vecCore::MaskedAssign(signc, tmp, (R)-1.0);
  vecCore::MaskedAssign(a, tmp, -c);

  const R delta = 0.5 * (Math::Sqrt(a * a + 4.) + a);
  //  delta += a;
  //  delta *= 0.5;

  const R cofA = -signc * Math::Exp(Math::Log(delta) / 3.0);
  cosTheta     = cofA - 1. / cofA;

  const R tau  = elenergy * geant::units::kInvElectronMassC2;
  const R beta = Math::Sqrt(tau * (tau + 2.)) / (tau + 1.);

  cosTheta = (cosTheta + beta) / (1. + cosTheta * beta);
  // check cosTheta limit
  cosTheta = Math::Min((R)1.0, cosTheta);
  // if (cosTheta>1.0) {
  //  cosTheta = 1.0;
  //}
  sinTheta = Math::Sqrt((1. - cosTheta) * (1. + cosTheta));
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
double RelativisticBremsModel::ComputeDEDXPerVolume(const Material *mat, double gammacut, double eekin)
{
  const double xsecFactor = 16.0 * geant::units::kFineStructConst * geant::units::kClassicElectronRadius *
                            geant::units::kClassicElectronRadius / 3.0;
  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  const double etot       = eekin + geant::units::kElectronMassC2; // it's always e-/e+
  const double dumc0      = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol() * gDensityFactor;
  const double densityCor = dumc0 * etot * etot;
  const double lpmEnergy  = mat->GetMaterialProperties()->GetRadiationLength() * gLPMFactor;
  // material dependent constant: the e-/e+ energy at which the corresponding k_LPM=k_p
  const double energyThLPM = std::sqrt(dumc0) * lpmEnergy;
  const bool isLPM         = (fIsUseLPM && etot > energyThLPM);

  // upper limit of the integral
  const double upperlimit   = std::min(gammacut, eekin);
  const double kappacr      = upperlimit / etot;
  const double log1mKappaCr = Math::Log(1.0 - kappacr);
  //
  double dedx = 0.0;
  // we will need the element composition of this material
  const Vector_t<Element *> &theElements  = mat->GetElementVector();
  const double *theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int numElems                            = theElements.size();
  // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
  const std::vector<double> &glW = fGL->GetWeights();
  const std::vector<double> &glX = fGL->GetAbscissas();
  // do the integration on reduced photon energy transformd to log(kappa) that transformed to integral between 0-1
  double integral = 0.0;
  for (int i = 0; i < fNGL; ++i) {
    const double dumx   = 1.0 - Math::Exp(glX[i] * log1mKappaCr);
    const double egamma = etot * dumx;
    double sum          = 0.0;
    for (int ielem = 0; ielem < numElems; ++ielem) {
      double zet = theElements[ielem]->GetZ();
      double val = 0.0;
      if (isLPM) {
        int izet = std::min(std::lrint(zet), gMaxZet - 1);
        val      = ComputeURelDXSecPerAtom(egamma, etot, lpmEnergy, densityCor, izet);
      } else {
        val = ComputeDXSecPerAtom(egamma, etot, zet);
      }
      sum += theAtomicNumDensityVector[ielem] * zet * zet * val;
    }
    integral += glW[i] * sum * (egamma - etot) / (1. + densityCor / (egamma * egamma));
  }
  dedx = log1mKappaCr * xsecFactor * integral;
  return dedx;
}

/**
 *   The restricted atomic cross section for bremsstrahlung photon emeission for target element with atomic number
 *   \f$Z\f$, gamma photon production energy threshold \f$k_c\f$ is computed for e-/e+ kinetic energy \f$E\f$ according
 * to
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
double RelativisticBremsModel::ComputeXSectionPerAtom(const Element *elem, const Material *mat, double gcut,
                                                      double eekin)
{
  double xsec = 0.0;
  if (eekin <= gcut) {
    return xsec;
  }
  const double xsecFactor = 16.0 * geant::units::kFineStructConst * geant::units::kClassicElectronRadius *
                            geant::units::kClassicElectronRadius / 3.0;
  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  const double etot       = eekin + geant::units::kElectronMassC2; // it's always e-/e+
  const double dumc0      = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol() * gDensityFactor;
  const double densityCor = dumc0 * etot * etot;
  const double lpmEnergy  = mat->GetMaterialProperties()->GetRadiationLength() * gLPMFactor;
  // material dependent constant: the e-/e+ energy at which the corresponding k_LPM=k_p
  const double energyThLPM = std::sqrt(dumc0) * lpmEnergy;
  const bool isLPM         = (fIsUseLPM && etot > energyThLPM);

  const double lKappaPrimePerCr = Math::Log(eekin / gcut);
  //
  // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
  const std::vector<double> &glW = fGL->GetWeights();
  const std::vector<double> &glX = fGL->GetAbscissas();
  // do the integration on reduced photon energy transformd to log(kappa) that transformed to integral between 0-1
  const double zet = elem->GetZ();
  const int izet   = std::min(std::lrint(zet), gMaxZet - 1);
  double integral  = 0.0;
  for (int i = 0; i < fNGL; ++i) {
    const double egamma = Math::Exp(glX[i] * lKappaPrimePerCr) * gcut;
    double val          = 0.0;
    if (isLPM) {
      val = ComputeURelDXSecPerAtom(egamma, etot, lpmEnergy, densityCor, izet);
    } else {
      val = ComputeDXSecPerAtom(egamma, etot, zet);
    }
    integral += glW[i] * val / (1. + densityCor / (egamma * egamma));
  }
  return lKappaPrimePerCr * xsecFactor * zet * zet * integral;
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
double RelativisticBremsModel::ComputeXSectionPerVolume(const Material *mat, double gcut, double eekin)
{
  double xsec = 0.0;
  if (eekin <= gcut) {
    return xsec;
  }
  const double xsecFactor = 16.0 * geant::units::kFineStructConst * geant::units::kClassicElectronRadius *
                            geant::units::kClassicElectronRadius / 3.0;
  // this is k_p^2 / E_{t-electron}^2 that is independent from E_{t-electron}^2
  const double etot       = eekin + geant::units::kElectronMassC2; // it's always e-/e+
  const double dumc0      = mat->GetMaterialProperties()->GetTotalNumOfElectronsPerVol() * gDensityFactor;
  const double densityCor = dumc0 * etot * etot;
  const double lpmEnergy  = mat->GetMaterialProperties()->GetRadiationLength() * gLPMFactor;
  // material dependent constant: the e-/e+ energy at which the corresponding k_LPM=k_p
  const double energyThLPM = std::sqrt(dumc0) * lpmEnergy;
  const bool isLPM         = (fIsUseLPM && etot > energyThLPM);

  const double lKappaPrimePerCr = Math::Log(eekin / gcut);
  // we will need the element composition of this material
  const Vector_t<Element *> &theElements  = mat->GetElementVector();
  const double *theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  const int numElems                      = theElements.size();
  // we need the abscissas and weights for the numerical integral 64-points GL between 0-1
  const std::vector<double> &glW = fGL->GetWeights();
  const std::vector<double> &glX = fGL->GetAbscissas();
  // do the integration on reduced photon energy transformd to log(kappa) that transformed to integral between 0-1
  double integral = 0.0;
  for (int i = 0; i < fNGL; ++i) {
    const double egamma = Math::Exp(glX[i] * lKappaPrimePerCr) * gcut;
    double sum          = 0.0;
    for (int ielem = 0; ielem < numElems; ++ielem) {
      const double zet = theElements[ielem]->GetZ();
      double val       = 0.0;
      if (isLPM) {
        int izet = std::min(std::lrint(zet), gMaxZet - 1);
        val      = ComputeURelDXSecPerAtom(egamma, etot, lpmEnergy, densityCor, izet);
      } else {
        val = ComputeDXSecPerAtom(egamma, etot, zet);
      }
      sum += theAtomicNumDensityVector[ielem] * zet * zet * val;
    }
    integral += glW[i] * sum / (1. + densityCor / (egamma * egamma));
  }
  return lKappaPrimePerCr * xsecFactor * integral;
}

void RelativisticBremsModel::SampleSecondaries(LightTrack_v &tracks, geant::TaskData *td)
{
  // Prepare temporary arrays for SIMD processing
  const int N                 = tracks.GetNtracks();
  double *gammaeEnergyArray   = td->fPhysicsData->fPhysicsScratchpad.fEps;
  double *gammaCutArr         = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr;
  double *zet                 = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr2;
  double *densityCorrConstArr = td->fPhysicsData->fPhysicsScratchpad.fR0;
  double *lpmEnergyArray      = td->fPhysicsData->fPhysicsScratchpad.fR1;

  // Gather data for sampling into array form
  for (int i = 0; i < N; i += kVecLenD) {

    Double_v primEkin = tracks.GetKinEVec(i);
    Double_v gammaCut;
    IndexD_v mcIndxLocal; // Used by alias method
    Double_v densityCorConst;
    Double_v lpmEnergy;
    for (int l = 0; l < kVecLenD; ++l) {
      const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(tracks.GetMaterialCutCoupleIndex(i + l));
      Set(gammaCut, l, matCut->GetProductionCutsInEnergy()[0]);
      Set(densityCorConst, l,
          matCut->GetMaterial()->GetMaterialProperties()->GetTotalNumOfElectronsPerVol() * gDensityFactor);
      Set(lpmEnergy, l, matCut->GetMaterial()->GetMaterialProperties()->GetRadiationLength() * gLPMFactor);
      if (GetUseSamplingTables()) {
        Set(mcIndxLocal, l, fGlobalMatGCutIndxToLocal[matCut->GetIndex()]);
      } else {
        // sample target element
        const Vector_t<Element *> &theElements = matCut->GetMaterial()->GetElementVector();
        double targetElemIndx                  = 0;
        if (theElements.size() > 1) {
          targetElemIndx = SampleTargetElementIndex(matCut, Get(primEkin, l), td->fRndm->uniform());
        }
        zet[i + l] = theElements[targetElemIndx]->GetZ();
      }
    }
    Double_v totalEn    = primEkin + geant::units::kElectronMassC2;
    Double_v densityCor = densityCorConst * (totalEn * totalEn);
    assert((primEkin >= gammaCut).isFull()); // Cut filtering should be applied up the call chain.
    if (GetUseSamplingTables()) {
      // Sample alias table here
      Double_v r1          = td->fRndm->uniformV();
      Double_v r2          = td->fRndm->uniformV();
      Double_v r3          = td->fRndm->uniformV();
      Double_v primEkin    = tracks.GetKinEVec(i);
      Double_v gammaEnergy = SampleEnergyTransfer(gammaCut, densityCor, mcIndxLocal, fAliasData.fLogEmin.data(),
                                                  fAliasData.fILDelta.data(), primEkin, r1, r2, r3);
      vecCore::Store(gammaEnergy, gammaeEnergyArray + i);
    } else {
      vecCore::Store(gammaCut, gammaCutArr + i);
      vecCore::Store(densityCorConst, densityCorrConstArr + i);
      vecCore::Store(lpmEnergy, lpmEnergyArray + i);
    }
  }

  if (!GetUseSamplingTables()) {
    // Always create fake particle at the end of the input arrays to vector rejection sampling method
    tracks.GetKinEArr()[N] = tracks.GetKinEArr()[N - 1];
    gammaCutArr[N]         = gammaCutArr[N - 1];
    zet[N]                 = zet[N - 1];
    densityCorrConstArr[N] = densityCorrConstArr[N - 1];
    lpmEnergyArray[N]      = lpmEnergyArray[N - 1];

    SampleEnergyTransfer(tracks.GetKinEArr(), gammaCutArr, zet, densityCorrConstArr, lpmEnergyArray, gammaeEnergyArray,
                         N, td);
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

bool RelativisticBremsModel::IsModelUsable(const MaterialCuts *matCut, double ekin)
{
  const double gammaCut = matCut->GetProductionCutsInEnergy()[0];
  return ekin < GetHighEnergyUsageLimit() && ekin > GetLowEnergyUsageLimit() && ekin > gammaCut;
}

} // namespace geantphysics
