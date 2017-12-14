
#ifndef KLEINNISHINACOMPTONMODEL_H
#define KLEINNISHINACOMPTONMODEL_H

#include "EMModel.h"

// from geantV
#include "Geant/Config.h"
namespace Geant {
  inline namespace GEANT_IMPL_NAMESPACE {
  class GeantTaskData;
}
}

namespace geantphysics {
  inline namespace GEANT_IMPL_NAMESPACE {
    class Material;
    class Element;
  }
}

#include <string>
#include <vector>

namespace geantphysics {

/**
 * @brief   Compton (incoherent) scattering model of photons on atomic electrons.
 * @class   KleinNishinaComptonModel
 * @author  M Novak
 * @date    april 2017
 *
 * Model for describing Compton scattering of photons on atomic electrons. The model is based on parametrized atomic
 * cross sections (Geant4 \cite cirrone2010validation \cite agostinelli2003geant4) and the Klein-Nishina
 * \cite klein1929streuung differential cross section for sampling the final state photon energy. Note, that the
 * Klein-Nishina model describes the interaction of photons with free electrons at rest i.e. atomic binding and momentum
 * distribution (that results in Doppler broadening) of the target atomic electrons are not included. However, the
 * parametrization of the atomic cross section is based on numerical cross section tables \cite storm1970photon
 * \cite hubbell1980pair that include an incoherent scattering function, \f$S(\bar{q},Z)\f$ (~ probablity that the atom
 * will be raised to an excited or ionized state when the photon imparts a recoil momentum \f$ \bar{q} \f$ to one of the
 * atomic electrons) i.e. account the elecrton binding effect that will result in more accurate cross sections at lower
 * photon energies compared to the pure Klein-Nishina cross sections.
 *
 */

//class Material;
class MaterialCuts;
//class Element;
class AliasTable;
class Particle;
class LightTrack;

class KleinNishinaComptonModel : public EMModel {
public:
/**
* @name Constructor, destructor:
*/
//@{
  /**
  * @brief Constructor.
  *
  * @param[in] modelname   Name of the model.
  */
  KleinNishinaComptonModel(const std::string &modelname="ComptonKleinNishina");
  /** @brief Destructor. */
  virtual ~KleinNishinaComptonModel();
//@}

/**
* @name Implemented EMModel base class methods:
*/
//@{
  virtual void   Initialize();
  virtual double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *particle);
  virtual double ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut, double kinenergy, const Particle *particle);
  virtual int    SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td);
//@}



private:
  /** @brief Copy constructor  (deleted) */
  KleinNishinaComptonModel(const KleinNishinaComptonModel&) = delete;
  /** @brief Operator=  (deleted) */
  KleinNishinaComptonModel &operator=(const KleinNishinaComptonModel&) = delete;


/**
* @name Model specific private methods.
*/
//@{

  /**
   * @brief Internal method to compute parametrized atomic Compton scattering cross section for a given target atom,
   *        gamma particle kinetic energy.
   *
   * @param[in] z        Atomic number of the target atom.
   * @param[in] egamma   Kinetic energy of the gamma particle in internal [energy] units.
   * @return    The computed atomic Compton scattering cross section in internal [lenght^2] units.
   */
  double ComputeAtomicCrossSection(double z, double egamma);


  /**
    * @brief Internal method to sample post interaction reduced photon energy from the prepared sampling tables.
    *
    *  @param[in] egamma     Kinetic energy of the primary gamma particle \f$ E_0 \f$.
    *  @param[in] r1         Random number distributed uniformly in [0,1].
    *  @param[in] r2         Random number distributed uniformly in [0,1].
    *  @param[in] r3         Random number distributed uniformly in [0,1].
    *  @return    Sampled post interaction reduced photon energy \f$ \epsilon = E_1/E_0\f$.
    */
  double SampleReducedPhotonEnergy(const double egamma, const double r1, const double r2, const double r3);


  /**
    * @brief Internal method to sample post interaction reduced photon energy.
    *
    *  @param[in]     egamma     Kinetic energy of the primary gamma particle \f$ E_0 \f$.
    *  @param[in,out] onemcost   One minus cos theta.
    *  @param[in,out] sint2      Sin theta square.
    *  @param[in]     td         Pointer to the GeantV thread local data object (used to get random numbers).
    *  @return    Sampled post interaction reduced photon energy \f$ \epsilon = E_1/E_0\f$.
    */
  double SampleReducedPhotonEnergy(const double egamma, double &onemcost, double &sint2, const Geant::GeantTaskData *td);

  /**
   * @brief Internal method to compute distribution of reduced (post interaction) photon energy related transformed
   *        variable.
   *
   * @param[in] xi       The transformed variable value (\f$ \xi \in [0,1] \f$).
   * @param[in] kappa    The constant \f$ \kappa=E_0/(m_ec^2)\f$ i.e. initial photon energy in electron rest mass units.
   * @return    The computed \f$ \xi \f$ dependent \f$ \left( \frac{\mathrm{d}\sigma}{\mathrm{d}\xi} \right)^* \f$ part
   *            of the transformed differential cross section.
   */
  double ComputeDXSecPerAtom(double xi, double kappa);


  /** @brief Internal method to build reduced (post interaction) photon energy related sampling tables.*/
  void   InitSamplingTables();


  /** @brief Internal method to build one reduced (post interaction) photon energy sampling tables for a given initial
   *         gamma energy.
   *
   *  The initial gamma energy \f$ E_0 \f$ will determine the \f$ \kappa=E_0/(m_ec^2)\f$ input variable.
   *  This method is used by the InitSamplingTables() method to build sampling table (using Walker's discrete alias
   *  method combined with <em>linear approximation of the p.d.f.</em>) at a given initial gamma energy. The method will
   *  call the ComputeDXSecPerAtom() internal method to compute the post interaction photon energy related transformed
   *  variable distribution.
   *
   *  @param[in]  indx   Index of the alias table data structure to build in the container.
   *  @param[in]  kappa  Initial photon energy \f$ E_0 \f$ dependent input variable \f$ \kappa=E_0/(m_ec^2)\f$.
   */
  void   BuildOneLinAlias(int indx, double kappa);

  void   ClearSamplingTables();
//@}

  /** @brief Internal data structure to store data for sampling the post interaction gaamma energy related transformd
    *        variable.
    *
    *  This data structure is set up at initialisation for each initial gamma kinetic energy grid point to sample the
    *  post interaction gamma energy related transformed variable using a combination  of Walker's alias sampling and
    *  liner approximation. We will have as many data structure as discrete primary gamma energy grid points.
    */
  struct LinAlias{
    LinAlias(int num) { fXdata.resize(num); fYdata.resize(num); fAliasW.resize(num); fAliasIndx.resize(num); }
    /** @brief Post interaction gamma energy related transformed variable values. */
    std::vector<double> fXdata;
    /** @brief The pdf values (not necessarily normalised) over the energy transfer related variable values. */
    std::vector<double> fYdata;
    /** @brief The alias probabilities (not necessarily normalised) over the energy transfer related variables. */
    std::vector<double> fAliasW;
    /** @brief The alias indices over the energy transfer related transformed variable values. */
    std::vector<int>    fAliasIndx;
  };


// data members
private:
  /** @brief  Internal code of the secondary partcile (e-). */
  int    fSecondaryInternalCode;

  int    fSTNumPhotonEnergiesPerDecade;    // ST=>SamplingTables
  int    fSTNumDiscreteEnergyTransferVals; // ST=>SamplingTables
  int    fSTNumPhotonEnergies;             // ST=>SamplingTables

  double fSTLogMinPhotonEnergy;            // ST=>SamplingTables
  double fSTILDeltaPhotonEnergy;           // ST=>SamplingTables

  /** @brief Container to store pointers to LinAlias data structures.*/
  std::vector<LinAlias*>   fSamplingTables;
  /** @brief An alias sampler used at run-time sampling of the post interaction gamma kinetic energy related transfered
    *        variable from a LinAlias data structure (prepared at initialisation).
    */
  AliasTable *fAliasSampler;
};

}        // namespace geantphysics

#endif   // KLEINNISHINACOMPTONMODEL_H
