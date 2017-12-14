
#include "KleinNishinaComptonModel.h"

#include "PhysicalConstants.h"
#include "Material.h"
#include "Element.h"
#include "MaterialProperties.h"

#include "MaterialCuts.h"
#include "AliasTable.h"

#include "Gamma.h"
#include "Electron.h"


#include "LightTrack.h"
#include "PhysicsData.h"

// from geantV
#include "GeantTaskData.h"

namespace geantphysics {


KleinNishinaComptonModel::KleinNishinaComptonModel(const std::string &modelname) : EMModel(modelname) {
  fSecondaryInternalCode              = -1;    // will be set at init

  fSTNumPhotonEnergiesPerDecade       = 10;    // ST=>SamplingTables
  fSTNumDiscreteEnergyTransferVals    = 60;    // ST=>SamplingTables
  fSTNumPhotonEnergies                = -1;    // ST=>SamplingTables: will be set at init

  fSTLogMinPhotonEnergy               = -1.;   // ST=>SamplingTables: will be set at init
  fSTILDeltaPhotonEnergy              = -1.;   // ST=>SamplingTables: will be set at init

  fAliasSampler                       = nullptr;

  SetLowestSecondaryEnergy(100.0*geant::eV);    // zero by default in the base class
}


KleinNishinaComptonModel::~KleinNishinaComptonModel() {
  if (GetUseSamplingTables()) {
    ClearSamplingTables();
  }
  if (fAliasSampler) {
    delete fAliasSampler;
  }
}


void KleinNishinaComptonModel::Initialize() {
  EMModel::Initialize();
  fSecondaryInternalCode = Electron::Definition()->GetInternalCode();
  if (GetUseSamplingTables()) {
    InitSamplingTables();
  }
}


double KleinNishinaComptonModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy,
                                                            const Particle* /*particle*/) {
  double xsec = 0.0;
  // this equal just to get the G4 behaviour at the lowest edge lambda table
  if (kinenergy<=GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
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
  return xsec;
}


double KleinNishinaComptonModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts* /*matcut*/,
                                                        double kinenergy, const Particle* /*particle*/) {
   double xsec  = 0.0;
   // this equal just to get the G4 behaviour at the lowest edge lambda table
   if (kinenergy<=GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
     return xsec;
   }
   // compute the parametrized atomic cross section: depends only on target Z and gamma energy.
   xsec = ComputeAtomicCrossSection(elem->GetZ(), kinenergy);
   return xsec;
}


int KleinNishinaComptonModel::SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td) {
  int    numSecondaries = 0;
  const double ekin     = track.GetKinE();
  // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
  // check if kinetic energy is above fHighEnergyUsageLimit andd o nothing if yes
  if (ekin<GetLowEnergyUsageLimit() || ekin>GetHighEnergyUsageLimit()) {
    return numSecondaries;
  }
  // sample post interaction gamma energy
  // here we need 3 random number + later one more for sampling phi
  double *rndArray  = td->fDblArray;
  td->fRndm->uniform_array(4, rndArray);
  double eps, oneMinusCost, sint2;
  if (GetUseSamplingTables()) {
    eps = SampleReducedPhotonEnergy(ekin, rndArray[0], rndArray[1], rndArray[2]);
    const double kappa = ekin/geant::kElectronMassC2;
    oneMinusCost = (1./eps-1.)/kappa;
    sint2        = oneMinusCost*(2.-oneMinusCost);
  } else {
    eps = SampleReducedPhotonEnergy(ekin, oneMinusCost, sint2, td);
  }
  // compute gamma scattering angle (realtive to the origininal dir i.e. Z)
  sint2 = std::max(0.,sint2);
  const double cost  = 1.0-oneMinusCost;
  const double sint  = std::sqrt(sint2);
  const double phi   = geant::kTwoPi*(rndArray[3]);
  // direction of the scattered gamma in the scattering frame
  double dirX  = sint*std::cos(phi);
  double dirY  = sint*std::sin(phi);
  double dirZ  = cost;
  // rotate back to lab frame
  RotateToLabFrame(dirX, dirY, dirZ, track.GetDirX(), track.GetDirY(), track.GetDirZ());
  //
  // keep org. gamma dir in lab frame: will be updated but will be needed later
  const double orgGamDirX = track.GetDirX();
  const double orgGamDirY = track.GetDirY();
  const double orgGamDirZ = track.GetDirZ();
  // Update primary gamma track properties i.e. the scattered gamma
  double eDeposit = 0.0;
  const double postGammaE = ekin*eps;
  if (postGammaE>GetLowestSecondaryEnergy()) {
    // update primary track kinetic energy
    track.SetKinE(postGammaE);
    // update primary track direction
    track.SetDirX(dirX);
    track.SetDirY(dirY);
    track.SetDirZ(dirZ);
  } else {
    eDeposit += postGammaE;
    track.SetKinE(0.0);
    track.SetTrackStatus(LTrackStatus::kKill);
  }
  //
  // Compute secondary e- properties: first the enrgy to check
  const double elEnergy = ekin-postGammaE; // E_el = E_0-E_1
  if (elEnergy>GetLowestSecondaryEnergy()) {
    // compute the secondary e- direction: from momentum vector conservation
    // final momentum of the secondary e- in the lab frame: = P_1-P_0 (binary col.)
    double elDirX = ekin*orgGamDirX - postGammaE*dirX;
    double elDirY = ekin*orgGamDirY - postGammaE*dirY;
    double elDirZ = ekin*orgGamDirZ - postGammaE*dirZ;
    // normalisation factor
    const double norm  = 1.0/std::sqrt(elDirX*elDirX + elDirY*elDirY + elDirZ*elDirZ);
    elDirX *= norm;
    elDirY *= norm;
    elDirZ *= norm;
    // create the secondary partcile i.e. the e-
    numSecondaries = 1;
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
    // this is known since it is a secondary track
    //  sectracks[secIndx].SetTrackStatus(LTrackStatus::kNew);
    std::vector<LightTrack>& sectracks = td->fPhysicsData->GetListOfSecondaries();
    sectracks[secIndx].SetDirX(elDirX);
    sectracks[secIndx].SetDirY(elDirY);
    sectracks[secIndx].SetDirZ(elDirZ);
    sectracks[secIndx].SetKinE(elEnergy);
    sectracks[secIndx].SetGVcode(fSecondaryInternalCode);  // e- GV code
    sectracks[secIndx].SetMass(geant::kElectronMassC2);
    sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); // parent GeantTrack index
  } else {
    eDeposit += elEnergy;
  }
  //
  // set (possible) energy deposit
  track.SetEnergyDeposit(eDeposit);
  //
  // return with number of secondaries i.e. 1 or 0
  return numSecondaries;
}



/**
 * @internal
 *
 * The method computes atomic Compton scattering cross sections based on the Geant4 \cite cirrone2010validation
 * \cite agostinelli2003geant4 parametrization. Note, that the Klein-Nishina model \cite klein1929streuung describes the
 * interaction of photons with free electrons at rest i.e. atomic binding and momentum distribution (that would result
 * in Doppler broadening) of the target atomic electrons are not included. However, the Geant4 parametrization of the
 * atomic cross section is based on numerical cross section tables \cite storm1970photon \cite hubbell1980pair that
 * include an incoherent scattering function, \f$S(\bar{q},Z)\f$ (~ probablity that the atom will be raised to an
 * excited or ionized state when the photon imparts a recoil momentum \f$ \bar{q} \f$ to one of the atomic electrons)
 * i.e. account the elecrton binding effect that will result in more accurate cross sections at lower photon energies
 * compared to the pure Klein-Nishina cross sections.
 *
 * According to the Geant4 documentation \cite g4physref, the numerical atomic cross sections \f$ \sigma(Z,E_0) \f$ as
 * a function of the target atomic numberg \f$ Z \f$ and photon energy \f$ E_0 \f$ given in \cite storm1970photon
 * \cite hubbell1980pair where approximated by the following function
 * \f[
 *    \sigma(Z,E_0) \approx \tilde{\sigma}(Z,E_0) \equiv P_1(Z)\frac{\ln[1+2\kappa]}{\kappa}
 *           + \frac{ P_2(Z) + \kappa P_3(Z) + \kappa^2P_4(Z)}{ 1 + a\kappa +b\kappa^2 +c\kappa^3}
 * \f]
 * where \f$ \kappa=E_0/(m_ec^2) \f$ i.e. photon energy in electron rest mass units and \f$P_i(Z)=d_iZ+e_iZ^2+f_iZ^3\f$
 * and the parameter values
 * \f[
 *   \begin{array}{lcrr}
 *    a   & = &  20.0 &\\
 *    b   & = & 230.0 &\\
 *    c   & = & 440.0 &\\
 *    &&&\\
 *    d_1 & = & +2.7965e-1  &  \text{[barn]} \\
 *    d_2 & = & -1.8300e-1  &  \text{[barn]} \\
 *    d_3 & = & +6.7527e+0  &  \text{[barn]} \\
 *    d_4 & = & -1.9798e+1  &  \text{[barn]} \\
 *    &&&\\
 *    e_1 & = & +1.9756e-5  &  \text{[barn]} \\
 *    e_2 & = & -1.0205e-2  &  \text{[barn]} \\
 *    e_3 & = & -7.3913e-2  &  \text{[barn]} \\
 *    e_4 & = & +2.7079e-2  &  \text{[barn]} \\
 *    &&&\\
 *    f_1 & = & -3.9178e-7  &  \text{[barn]} \\
 *    f_2 & = & +6.8241e-5  &  \text{[barn]} \\
 *    f_3 & = & +6.0480e-5  &  \text{[barn]} \\
 *    f_4 & = & +3.0274e-4  &  \text{[barn]} \\
 *   \end{array}
 * \f]
 * where determined by fitting the approximate formula to 511 data points from the tables \cite storm1970photon
 * \cite hubbell1980pair in the intervals \f$ Z \in [1,100], E_0 \in [10 \text{ keV}, 100 \text{ GeV}] \f$. The accuracy
 * of the approximated cross sections \f$ \tilde{\sigma}(Z,E_0) \f$ is estimated to be \cite cirrone2010validation
 * \cite g4physref :
 * \f[
 *  \begin{cases}
 *    \approx 10  \%  & \quad \text{ for } & E_0 \simeq 10 \text{ keV } - 20 \text{ keV } \\
 *    \leq    5-6 \%  & \quad \text{ for } & E_0 >      20 \text{ keV }  \\
 *  \end{cases}
 * \f]
 * @endinternal
 */
double KleinNishinaComptonModel::ComputeAtomicCrossSection(double z, double egamma) {
  double xsec = 0.0;
  //
  constexpr double a  =  20.0;
  constexpr double b  = 230.0;
  constexpr double c  = 440.0;
  //
  constexpr double d1 =   2.7965e-1*geant::barn;
  constexpr double d2 =  -1.8300e-1*geant::barn;
  constexpr double d3 =   6.7527   *geant::barn;
  constexpr double d4 =  -1.9798e+1*geant::barn;
  //
  constexpr double e1 =   1.9756e-5*geant::barn;
  constexpr double e2 =  -1.0205e-2*geant::barn;
  constexpr double e3 =  -7.3913e-2*geant::barn;
  constexpr double e4 =   2.7079e-2*geant::barn;
  //
  constexpr double f1 =  -3.9178e-7*geant::barn;
  constexpr double f2 =   6.8241e-5*geant::barn;
  constexpr double f3 =   6.0480e-5*geant::barn;
  constexpr double f4 =   3.0274e-4*geant::barn;
  //
  const double z2  = z*z;
  const double p1Z = z*(d1 + e1*z + f1*z2);
  const double p2Z = z*(d2 + e2*z + f2*z2);
  const double p3Z = z*(d3 + e3*z + f3*z2);
  const double p4Z = z*(d4 + e4*z + f4*z2);
  //
  double t0  = 15.0*geant::keV;
  if (z<1.5) {
    t0 = 40.0*geant::keV;
  }
  //
  double kappa  = std::max(egamma,t0)/geant::kElectronMassC2;
  double kappa2 = kappa*kappa;
  double kappa3 = kappa2*kappa;
  xsec  = p1Z*std::log(1. + 2.*kappa)/kappa + (p2Z + p3Z*kappa + p4Z*kappa2)/(1. + a*kappa + b*kappa2 + c*kappa3);
  // low energy correction:
  if (egamma<t0) {
    constexpr double dt0 = 1.*geant::keV;
    kappa   = (t0+dt0)/geant::kElectronMassC2;
    kappa2  = kappa*kappa;
    kappa3  = kappa2*kappa;
    const double sigma = p1Z*std::log(1. + 2*kappa)/kappa + (p2Z + p3Z*kappa + p4Z*kappa2)/(1. + a*kappa + b*kappa2 + c*kappa3);
    const double    c1 = -t0*(sigma-xsec)/(xsec*dt0);
    double c2 = 0.15;
    if (z>1.5) {
      c2 = 0.375-0.0556*std::log(z);
    }
    const double y = std::log(egamma/t0);
    xsec *= std::exp(-y*(c1+c2*y));
  }
  return xsec;
}


/**
 * @internal
 *
 * This method is called from the SampleSecondaries() method to sample the post interaction reduced photon energy
 * \f$ \epsilon = E_1/E_0 \f$ for given initial photon energy \f$ E_1 \f$. The sampling is done by using sampling tables
 * built over an initial photon energy grid at the initialisation of the model. It is ensured in the SampleSecondaries()
 * method that the primary photon energy is within the min/max values of this photon energy grid. After the initial
 * photon energy bin is determined, one sampling table is selected (i.e. the lower or higher bin edge table) by using
 * statistical interpolation on log photon energy scale. After the sampling table is selected, the reduced photon
 * energy related transformed variable \f$ \xi(\epsilon) \f$ is sampled from the table (by using discrete alias sampling
 * and linear approximation of the transformed p.d.f.) and transformed back by using the relation
 * \f$ \epsilon(\xi)=\exp(\ln(1+2\kappa)(\xi-1))\f$ with \f$ \kappa \equiv E_0/(m_ec^2) \f$ i.e. the initial photon
 * energy in electron rest mass units [see more about the variable transformation at the
 * ComputeDXSecPerAtom() method].
 *
 * @endinternal
 */
double KleinNishinaComptonModel::SampleReducedPhotonEnergy(const double egamma, const double r1, const double r2,
                                                            const double r3) {
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
  const LinAlias *als = fSamplingTables[indxEgamma];
  // sample the transformed variable xi=[\alpha-ln(ep)]/\alpha (where \alpha=ln(1/(1+2\kappa)))
  // that is in [0,1] when ep is in [ep_min=1/(1+2\kappa),ep_max=1] (that limits comes from energy and momentum
  // conservation in case of scattering on free electron at rest).
  // where ep = E_1/E_0 and kappa = E_0/(mc^2)
  const double  xi    = fAliasSampler->SampleLinear(&(als->fXdata[0]), &(als->fYdata[0]), &(als->fAliasW[0]),
                                                    &(als->fAliasIndx[0]), fSTNumDiscreteEnergyTransferVals, r2, r3);
  // transform it back to eps = E_1/E_0
  // \epsion(\xi) = \exp[ \alpha(1-\xi) ] = \exp [\ln(1+2\kappa)(\xi-1)]
  const double kappa  = egamma/geant::kElectronMassC2;
  return std::exp(std::log(1.+2.*kappa)*(xi-1.)); // eps = E_1/E_0
}


double KleinNishinaComptonModel::SampleReducedPhotonEnergy(const double egamma, double &onemcost, double &sint2,
                                                          const Geant::GeantTaskData *td) {
  const double kappa = egamma/geant::kElectronMassC2;
  const double eps0  = 1./(1.+2.*kappa);
  const double eps02 = eps0*eps0;
  const double al1   = -std::log(eps0);
  const double cond  = al1/(al1+0.5*(1.-eps02));
  //
  double eps, eps2, gf;
  double *rndArray = td->fDblArray;
  do {
    td->fRndm->uniform_array(3, rndArray);
    if (cond>rndArray[0]) {
      eps  = std::exp(-al1*rndArray[1]);
      eps2 = eps*eps;
    } else {
      eps2 = eps02+(1.-eps02)*rndArray[1];
      eps  = std::sqrt(eps2);
    }
    onemcost = (1.-eps)/(eps*kappa);
    sint2    = onemcost*(2.-onemcost);
    gf       = 1.-eps*sint2/(1.+eps2);
  } while (gf<rndArray[2]);
  return eps;
}


/**
  * @internal
  *
  * The Klein-Nishina \cite klein1929streuung differential cross section for incoherent scattering of photons with
  * energy \f$ E_0\f$ on free electron at rest
  * \f[
  *   \frac{\mathrm{d}\sigma}{\mathrm{d}\epsilon} = \pi r_0^2 Z \frac{m_ec^2}{E_0}
  *      \left[ \frac{1}{\epsilon}+\epsilon \right]
  *      \left[ 1- \frac{\epsilon \sin^2(\theta)}{1+\epsilon^2} \right]
  * \f]
  * where \f$ r_0\f$ is the calssical electron radius, \f$ Z\f$ is the target atomic number, \f$ m_e\f$ is the electron
  * rest mass, \f$ c \f$ is the speed of the light, \f$ \theta \f$ is the scattering angle  and
  * \f$ \epsilon \equiv E_1/E_0\f$ with \f$ E_1 \f$ being the post interaction energy of the photon that is (due to
  * energy and momentum conservation in case of scattering on free electron at rest)
  * \f[
  *    E_1=E_0 \frac{m_ec^2}{m_ec^2 + E_0[1-\cos(\theta)]} = \frac{E_0}{1+\kappa[1-\cos(\theta)]}
  * \f]
  * with \f$ \kappa \equiv E_0/(m_ec^2) \f$) i.e. the initial photon energy in electron rest mass units. The above
  * equation gives the kinematical limits of \f$ \epsilon \f$ by dividing both side with \f$ E_0\f$:
  * \f[
  *  \begin{matrix}
  *    \epsilon_{min} & = & \epsilon(\theta=\pi)   & =  & 1/(1+2\kappa) \\
  *    \epsilon_{max} & = & \epsilon(\theta=0)     & =  & 1
  *  \end{matrix}
  * \f]
  * Eliminating \f$ \sin^2(\theta)\f$, the differential cross section can be written as
  * \f[
  *   \frac{\mathrm{d}\sigma}{\mathrm{d}\epsilon} = \pi r_0^2 Z \frac{m_ec^2}{E_0}
  *      \left[ \frac{1}{\epsilon}+\epsilon \right]
  *      \left[ 1- \frac{\epsilon \beta[\beta+2]]}{1+\epsilon^2} \right]
  * \f]
  * where \f$ \beta \equiv (1-1/\epsilon)/\kappa\f$.
  * This cross section, differential in \f$ \epsilon \f$ can be transformed to
  * \f[
  *   \frac{\mathrm{d}\sigma}{\mathrm{d}\xi} = \pi r_0^2 Z \frac{m_ec^2}{E_0}(-\alpha)
  *      \left[ \frac{1}{\epsilon(\xi)}+\epsilon(\xi) \right]
  *      \left[ 1- \frac{\epsilon(\xi) \beta[\beta+2]]}{1+\epsilon(\xi)^2} \right]\epsilon(\xi)
  * \f]
  * where \f$ \alpha = \ln[1/(1+2\kappa)]\f$, \f$ \xi(\epsilon) = [\alpha-\ln(\epsilon)]/\alpha \f$ and
  * \f$ \epsilon(\xi) = \exp\{ \alpha[1-\xi] \}\f$. Note, that \f$ \xi(\epsilon_{min}) = 0 \f$ and
  * \f$ \xi(\epsilon_{max}) = 1 \f$ so \f$ \xi(\epsilon) \in [0,1] \f$ when \f$ \epsilon \in [1/(1+2\kappa),1] \f$.
  * This cross section, differential in \f$ \xi \f$ can be written as
  * \f[
  *   \frac{\mathrm{d}\sigma}{\mathrm{d}\xi} = \eta  \left( \frac{\mathrm{d}\sigma}{\mathrm{d}\xi} \right)^*
  * \f]
  * with the \f$ \xi \f$ independent, constant
  * \f[
  *    \eta = - \pi r_0^2 Z \frac{\alpha}{\kappa}
  * \f]
  * and the \f$ \xi \f$ dependent
  * \f[
  *   \left( \frac{\mathrm{d}\sigma}{\mathrm{d}\xi} \right)^* =
  *      \left[ \frac{1}{\epsilon(\xi)}+\epsilon(\xi) \right]
  *      \left[ 1- \frac{\epsilon(\xi) \beta[\beta+2]]}{1+\epsilon(\xi)^2} \right]\epsilon(\xi)
  * \f]
  * parts.
  *
  * This method computes the \f$ \xi \f$ dependent \f$ \left( \frac{\mathrm{d}\sigma}{\mathrm{d}\xi} \right)^* \f$ part
  * of the transformed differential cross section for the given \f$ \xi \in [0,1] \text{ and } \kappa \f$ input values.
  * It's used to compute the distribution of the transformed variable \f$ \xi \f$ in the BuildOneLinAlias()
  * method to build alias sampling tables at discrete \f$ E_0 \f$ photon energies. These sampling tables are used for
  * for run-time sampling of the transformed \f$ \xi \f$ variable in the SampleReducedPhotonEnergy() method and the
  * sampled value is transformed back to the corresponding \f$ \epsilon \f$ value by using the expression for
  * \f$ \epsilon(\xi) = \exp\{ \alpha[1-\xi] \}\f$ given above.
  *
  * @endinternal
  */
double KleinNishinaComptonModel::ComputeDXSecPerAtom(double xi, double kappa) {
  const double inv2Kappa  = 1./(1.+2.*kappa);
  const double linv2Kappa = std::log(inv2Kappa);
  const double eps        = std::exp(linv2Kappa*(1.-xi));
  const double invEps     = 1./eps;
  const double beta       = (1.-invEps)/kappa;

  return eps*(eps+invEps)*(1.+eps*beta*(2.+beta)/(1.+eps*eps));
}



void KleinNishinaComptonModel::InitSamplingTables() {
  ClearSamplingTables();
  // set number of primary gamma energy grid points
  const double minEprim = GetLowEnergyUsageLimit();
  const double maxEprim = GetHighEnergyUsageLimit();
  fSTNumPhotonEnergies = fSTNumPhotonEnergiesPerDecade*std::lrint(std::log10(maxEprim/minEprim))+1;
  fSTNumPhotonEnergies = std::max(fSTNumPhotonEnergies,3);
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
  // 4. set up the container that stores sampling tables for all the materials (init to nullptr-s)
  fSamplingTables.resize(fSTNumPhotonEnergies,nullptr);
  // 5. prepare sampling tables one-by-one
  for (int i=0; i<fSTNumPhotonEnergies; ++i) {
    const double egamma = primEVect[i];
    const double kappa  = egamma/geant::kElectronMassC2;
    BuildOneLinAlias(i,kappa);
  }
  primEVect.clear();
}


void KleinNishinaComptonModel::BuildOneLinAlias(int indx, double kappa) {
  LinAlias *tb          = new LinAlias(fSTNumDiscreteEnergyTransferVals);
  fSamplingTables[indx] = tb;
  // note: the transformd variable (xi) is in [0,1] when eps=E_1/E_0 in [eps_min, eps_max]
  // so fill 3 initial values of xi:
  //  -  xi_0 = x_min = 0
  //  -  xi_1 = (x_max-x_min)/2 = 0.5
  //  -  xi_2 = x_max = 1
  // and the corresponding y(i.e.~PDF) values
  tb->fXdata[0] = 0.0;
  tb->fXdata[1] = 0.5;
  tb->fXdata[2] = 1.0;
  tb->fYdata[0] = ComputeDXSecPerAtom(tb->fXdata[0],kappa);
  tb->fYdata[1] = ComputeDXSecPerAtom(tb->fXdata[1],kappa);
  tb->fYdata[2] = ComputeDXSecPerAtom(tb->fXdata[2],kappa);
  int curNumData = 3;
  // expand the data up to numdata points
  while (curNumData<fSTNumDiscreteEnergyTransferVals) {
    // find the lower index of the bin, where we have the biggest linear interp. error compared to spline
    double maxerr     = 0.0; // value of the current maximum error
    double thexval    = 0.0;
    double theyval    = 0.0;
    int    maxerrindx = 0;   // the lower index of the corresponding bin
    for (int i=0; i<curNumData-1; ++i) {
      const double xx    = 0.5*(tb->fXdata[i]+tb->fXdata[i+1]);    // mid x point
      const double yy    = 0.5*(tb->fYdata[i]+tb->fYdata[i+1]);    // lin. interpolated pdf value at the mid point
      const double val   = ComputeDXSecPerAtom(xx,kappa); // real pdf value at the mid point
      const double err   = std::abs(1.-(yy/val));
      if (err>maxerr) {
        maxerr     = err;
        maxerrindx = i;
        thexval    = xx;
        theyval    = val;
      }
    }
    // extend x,y data by puting a new real value at the mid point of the highest error bin
    // first shift all values to the right
    for (int j=curNumData; j>maxerrindx+1; --j) {
      tb->fXdata[j] = tb->fXdata[j-1];
      tb->fYdata[j] = tb->fYdata[j-1];
    }
    // fill x mid point
    tb->fXdata[maxerrindx+1] = thexval;
    tb->fYdata[maxerrindx+1] = theyval;
    // increase number of data
    ++curNumData;
  } // end while
  // prepare the alias data for this PDF(x,y)
  fAliasSampler->PreparLinearTable(&(tb->fXdata[0]), &(tb->fYdata[0]), &(tb->fAliasW[0]), &(tb->fAliasIndx[0]),
                                   fSTNumDiscreteEnergyTransferVals);
}


void KleinNishinaComptonModel::ClearSamplingTables() {
  size_t num = fSamplingTables.size();
  for (size_t itb=0; itb<num; ++itb) {
    LinAlias *tb = fSamplingTables[itb];
    if (tb) {
      tb->fXdata.clear();
      tb->fYdata.clear();
      tb->fAliasW.clear();
      tb->fAliasIndx.clear();
    }
    delete tb;
  }
  fSamplingTables.clear();
}

}  // namespace geantphysics
