
#ifndef RELATIVISTICBREMSMODEL_H
#define RELATIVISTICBREMSMODEL_H

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

//class Material;
//class Element;
class AliasTable;
class MaterialCuts;
class Particle;
class LightTrack;

/**
 * @brief   High energy Bremsstrahlung models for electron/positron.
 * @class   RelativisticBremsModel
 * @author  M Novak, A Ribon
 * @date    march 2016
 *
 * Relativistic and ultra relativistic Bremsstrahlung models for electron/positron.
 * Models are implemented based on \cite tsai1974pair. Landau-Pomeranchuk-Migdal
 * effect is included in the ultra relativistic model under the complete
 * screening approximation. Dielectric suppression is taken into account in both
 * cases. The electron/positron kinetic energy range of the model is between 1[keV] - 100[TeV]
 * but typically used between 1-2[GeV] - 100[TeV] since at lower energies (where LPM suppression is
 * not active) the Seltzer-Berger model (SeltzerBergerBremsModel) is more accurate.
 */

class RelativisticBremsModel: public EMModel {
public:

  virtual void Initialize(); // from EMModel
  virtual double ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle* particle, bool istotal=false);
  virtual double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *particle);
  virtual double ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut, double kinenergy, const Particle *particle);
  virtual int    SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td);

  virtual double MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle *part) const;

/**
* @name Constructor, destructor:
*/
//@{
   /**
    * @brief Constructor.
    *
    * @param[in] isuselpm  Flag to indicate if LPM suppression should be included.
    * @param[in] modelname Name of the model.
    */
  RelativisticBremsModel(bool isuselpm=true, const std::string &modelname="eRelativisticBrems");
   /** @brief Destructor. */
 ~RelativisticBremsModel();
//@}

  /**
   * @brief Public method to initilise the model.
   *
   * During the initialisation, internal sampling tables for run time sampling of the emitted photon energy will be
   * created. This method always need to be invoked between the construction and usage of the model.
   *
   */
  void Initialise();

  /**
   * @brief Public method to obtain (restricted) atomic cross sections.
   *
   * Since the model includes some material dependent corrections, for consistency reasons one also needs to provide
   * the material that the element, that the atomic cross section is requested, belongs to.
   *
   * @param[in] elem               Pointer to the element object that the atomic cross section is required.
   * @param[in] mat                Pointer to the material that the element blongs to.
   * @param[in] gammaprodcutenergy Kinetic energy threshold for gamma particle production.
   * @param[in] particleekin       Kinetic energy of the incident particle i.e. e-/e+.
   * @return                       (Restricted) Bremsstrahlung atomic cross section for the specified configuration
   *                                in internal [lenght^2] units.
   */
  double ComputeXSectionPerAtom(const Element *elem, const Material *mat, double gammaprodcutenergy, double particleekin);

  /**
   * @brief Public method to obtain (restricted) macroscopic cross sections.
   *
   * @param[in] mat                Pointer to the material that the cross section is requested.
   * @param[in] gammaprodcutenergy Kinetic energy threshold for gamma particle production.
   * @param[in] particleekin       Kinetic energy of the incident particle i.e. e-/e+.
   * @return                       (Restricted) Macroscopic bremsstrahlung cross section for the specified configuration
   *                                in internal [1/lenght] units.
   */
  double ComputeXSectionPerVolume(const Material *mat, double gammaprodcutenergy, double particleekin);

  /**
   * @brief Public method to obtain (restricted) radiative stopping power.
   *
   * @param[in] mat                Pointer to the material that the cross section is requested.
   * @param[in] gammaprodcutenergy Kinetic energy threshold for gamma particle production.
   * @param[in] particleekin       Kinetic energy of the incident particle i.e. e-/e+.
   * @return                       (Restricted) Radiative stopping power for the specified configuration in internal
   *                               [energy/lenght] units.
   */
  double ComputeDEDXPerVolume(const Material *mat, double gammaprodcutenergy, double particleekin);
//  double ComputeDEDXPerVolume1(Material *mat, double gammaprodcutenergy, double particleekin);

  /**
   * @brief Public method to sample the emitted (restricted) bremsstrahlung photon energy.
   *
   * The emitted bremsstrahlung photon energy is always higher than the gamma particle kinetic energy production
   * threshold in the specified material-cut pair and not higher than the incident particle kinetic energy.
   *
   * @param[in] matcut    Pointer to the material-cut pair in which the interaction takes place.
   * @param[in] eekin     Kinetic energy of the incident particle i.e. e-/e+.
   * @param[in] r1        Random number distributed uniformly on the 0 1 interval.
   * @param[in] r2        Random number distributed uniformly on the 0 1 interval.
   * @param[in] r3        Random number distributed uniformly on the 0 1 interval.
   * @return              An emitted bremsstrahlung photon energy sampled from the distribution specified by
   *                      given configuration and the model in internal [energy] units if the particle (e-/e+)
   *                      kinetic energy is higher than the gamma particle production energy threshold. Zero otherwise.
   */
  double SamplePhotonEnergy(const MaterialCuts *matcut, double eekin, double r1, double r2, double r3);

  //
  void   SamplePhotonDirection(double elenergy, double &sinTheta, double &cosTheta, double rndm);

/**
 * @name Field related setters/getters:
 */
//@{
  /** @brief Get the flag that indicates if LPM suppression is active (true by default). */
  bool               GetLPMFlag()  {return fIsUseLPM;}
  /** @brief Set the flag that indicates if LPM suppression is active (true by default).
    *
    * Must be set before initialisation if not the default value required.
    */
  void               SetLPMFlag(bool islpm=true) {fIsUseLPM = islpm;}
//@}


private:
  /** @brief Internal method to build emitted photon energy sampling tables under <em>linear approximation of the p.d.f.</em>.
   *
   *  Used at initialisation of the model to prepare emitted photon energy related transformed variable
   *  sampling tables for all different material-gamma production cut pairs over an e-/e+ kinetic energy grid.
   *  These tables are prepared for Walker's alias method combined with <em>linear approximation of the p.d.f.</em>
   *  within the bins.
   */
  void InitSamplingTables();

  /** @brief Internal method to build one emitted photon energy sampling tables under <em>linear approximation of the p.d.f.</em>.
   *
   *  This method is used by RelativisticBremsModel::InitSamplingTables() to build emitted photon energy related
   *  sampling tables for Walker's alias method combined with <em>linear approximation of the p.d.f.</em>. Grid points
   *  are adaptively inserted such that the error between the real and linearly interpolated p.d.f. is minimised.
   *  Number of grid points will be RelativisticBremsModel::LinAlias::fNumdata which is always the maximum i.e.
   *  RelativisticBremsModel::fNumSamplingPhotEnergies .
   */
  void BuildOneLinAlias(int ialias, const Material *mat, double gcut);

  /** @brief Internal method to build emitted photon energy sampling tables based on <em>rational interpolation of the c.d.f.</em>.
   *
   *  Used at initialisation of the model to prepare emitted photon energy related transformed variable
   *  sampling tables for all different material-gamma production cut pairs over an e-/e+ kinetic energy grid.
   *  These tables are prepared for Walker's alias method combined with <em>rational interpolation based inverse transform
   *  of the c.d.f.</em>.
   *
   */
  void InitSamplingTables1();

  /** @brief Internal method to build one emitted photon energy sampling tables based on <em>rational interpolation of the c.d.f.</em>.
   *
   *  This method is used by RelativisticBremsModel::InitSamplingTables1() to build emitted photon energy related
   *  sampling tables for Walker's alias method combined with <em>rational interpolation based inverse transform
   *  of the c.d.f.</em>. Grid points are adaptively inserted such that the error between the real p.d.f. and that obtaind
   *  from the rational approximation is minimised. Number of grid points will be RelativisticBremsModel::LinAlias::fNumdata
   *  which is at most RelativisticBremsModel::fNumSamplingPhotEnergies .
   */
  void BuildOneRatinAlias1(int ialias, const Material *mat, double gcut);

  /** @brief Internal method to compute (part of the) relativistic differential cross section for Bremsstrahlung photon
    *        emission when Landau-Pomeranchuk-Migdal effect is inactive.
    *
    * @param[in] gammaenergy Energy of the emitted gamma photon (\f$k\f$) at which the differential cross section is
    *                        required.
    * @param[in] totalenergy Total energy (\f$E_t\f$) of the pre-interaction projectile (e-/e+) for which the
    *                        differential cross section is required.
    * @param[in] z           Atomic number of the target (\f$Z\f$).
    * @return    Differential cross section, or more exactly the
    *            \f$\left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \f$ part of
    *            \f$ \frac{\mathrm{d}\sigma}{\mathrm{d}k} = \frac{16 \alpha r_e^2 Z^2}{3k}
    *            \left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \f$. Note, that the \f$ 1/\Gamma \f$ factor i.e.
    *            the leading factor of the dielectric suppression will be also applied in this case later (see more
    *            details on \f$ \Gamma \f$ at RelativisticBremsModel::ComputeURelDXSecPerAtom() and e.g.
    *            RelativisticBremsModel::ComputeXSectionPerVolume()).
    */
  double ComputeDXSecPerAtom(double gammaenergy, double totalenergy, double z);

  /** @brief Internal method to compute (part of the) relativistic differential cross section for Bremsstrahlung photon
    *        emission when Landau-Pomeranchuk-Migdal effect is active.
    *
    * @param[in] gammaenergy Energy of the emitted gamma photon (\f$k\f$) at which the differential cross section is
    *                        required.
    * @param[in] totalenergy Total energy (\f$E_t\f$) of the pre-interaction projectile (e-/e+) for which the
    *                        differential cross section is required.
    * @param[in] z           Atomic number of the target (\f$Z\f$).
    * @param[in] lpmenergy   LPM energy i.e. material dependent constant (\f$E_{LPM}\f$ see the more detailed
    *                        description).
    * @param[in] densitycor  Characteristic photon energy square i.e. (\f$k_p^2\f$ see the more detailed description).
    * @return    Differential cross section, or more exactly the
    *            \f$ \left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \f$ part of
    *            \f$ \frac{\mathrm{d}\sigma}{\mathrm{d}k} = \frac{16 \alpha r_e^2 Z^2}{3k}
    *            \left(\frac{\mathrm{d}\sigma}{\mathrm{d}k}\right)^* \f$
    */
  double ComputeURelDXSecPerAtom(double gammaenergy, double totalenergy, double z, double lpmenergy, double densitycor);

// data members
private:
  // use these elastic and inelatic form factors for light elements instead of TFM
  // under the complete screening approximation
  // Tsai Table.B2.
  /** @brief Elastic form factors for elements Z<8 from \cite tsai1974pair Table B2. */
  static const double gFelLowZet[8];
  /** @brief Inelastic form factors for elements Z<8 from \cite tsai1974pair Table B2. */
  static const double gFinelLowZet[8];

  // secondary related data
  int      fSecondaryInternalCode; // gamma GV code set at initialization

  /** @brief Flag to indicate if LPM effect should be included (true by default).
    *
    * Must be set before initialisation if not the default value is required by using the RelativisticBremsModel::SetLPMFlag(bool) method.
    */
  bool        fIsUseLPM;                 // flag to indicate if LPM suppression should be used (must be set before initialisation)
  /** @brief Flag to indicate if linear approximation based sampling tables are required to use (default true)
    *
    * If this flag is true (default value true), then the emitted photon energy is sampled by using a combination of Walker's alias
    * method with linear approximation of the probability density function within the bins. Otherwise, a rational interpolation based
    * inverse transform of the cumulative distribution function is used instead the linear approximation of the probability density
    * function within the bins.
    */
  bool        fIsUseLinearSamplingTable;

  // data members that describe internal data structures of sampling tables built at initialisation:
  // At initialisation, we build sampling tables over an e-/e+ kinetic energy grid (fNumSamplingElecEnergies grid
  // point with logarithmic spaceing between fMinElecEnergy and fMaxElecEnergy) for each different
  // material-gamma production cut pairs. At each kinetic energy grid point we build a sampling table with
  // fNumSamplingPhotEnergies photon energy point transformed to some other variable on the [0,1] interval.

  /** @brief Number of e-/e+ kinetic energy grid points in [RelativisticBremsModel::fMinElecEnergy, RelativisticBremsModel::fMaxElecEnergy] (default 41). */
  int     fNumSamplingElecEnergies;  // number of e-/e+ kinetic energy bins
  /** @brief Number of transformed emitted photon energy related variable in [0,1] (default 100). */
  int     fNumSamplingPhotEnergies;  // number of photon energy related bins
  /** @brief Minimum of the e-/e+ kinetic energy grid (default 1.0 [GeV]) */
  double  fMinElecEnergy;            // minimum e-/e+ kinetic energy
  /** @brief Maximum of the e-/e+ kinetic energy grid (default 100.0 [TeV]) */
  double  fMaxElecEnergy;            // maximum e-/e+ kinetic energy
  /** @brief Logarithm of RelativisticBremsModel::fMinElecEnergy i.e. ln(fMinElecEnergy) . */
  double  fElEnLMin;     // log min electron energy
  /** @brief Inverse of the e-/e+ kinetic energy grid delta i.e. ln[fMaxElecEnergy/fMinElecEnergy]/(fNumSamplingElecEnergies-1) */
  double  fElEnILDelta;  // 1/log(delta-electron energy)
  /** @brief The logarithmically spaced e-/e+ kinetic energy grid.
    *
    * RelativisticBremsModel::fNumSamplingElecEnergies points in [RelativisticBremsModel::fMinElecEnergy, RelativisticBremsModel::fMaxElecEnergy])
    */
  double *fSamplingElecEnergies;     // the e-/e+ kinetic energy grid with size of fNumSamplingElecEnergies
  /** @brief The logarithm of RelativisticBremsModel::fSamplingElecEnergies array.
    *
    * RelativisticBremsModel::fNumSamplingElecEnergies points in [ln(RelativisticBremsModel::fMinElecEnergy), ln(RelativisticBremsModel::fMaxElecEnergy)])
    */
  double *fLSamplingElecEnergies;    // log of e-/e+ kinetic energy grid with size of fNumSamplingElecEnergies

  // data to map all material-gamma production cut pair indices to local indices including only the subset of all
  // material-gamma production cut that are different. These data used only internally by the model.
  /** @brief Number of different material-gamma production cut pairs. */
  int     fNumDifferentMaterialGCuts;  // number of different matrial-gammacut pairs
  /** @brief Map from global to local material-gamma production cut indices. The size of the array is
    *        RelativisticBremsModel::numMaterialCuts
    */
  int    *fGlobalMatGCutIndxToLocal;  // maps the global mat.-cut indices to local indices that are used here

  /** @brief Internal data structure to store data for sampling the emitted photon energy when
    *        RelativisticBremsModel::fIsUseLinearSamplingTable is true (default).
    *
    *  This data structure is set up at initialisation for each different material-gamma production cut pairs
    *  over the e-/e+ kinetic energy grid to sample the emitted photon energy distribution using a combination
    *  of Walker's alias sampling and liner approximation. At most we will have as many data structure as
    *  RelativisticBremsModel::fNumDifferentMaterialGCuts times RelativisticBremsModel::fNumSamplingElecEnergies
    *  and these data structure pointers are stored in the RelativisticBremsModel::fAliasData linear array.
    *  However, data structures are created only for the possible e-/e+ kinetic energy - emitted photon energy
    *  combinations (i.e. for those e-/e+ kinetic energy grid points in RelativisticBremsModel::fSamplingElecEnergies
    *  that are above the gamma production cut value) and the other elements of RelativisticBremsModel::fAliasData
    *  linear array are left to be nullptr.
    */
  struct LinAlias{
    /** @brief Number of data points i.e. size of the arrays. */
    int     fNumdata;
    /** @brief Reduced photon energy related transformed variable values. */
    double *fXdata;
    /** @brief The probability density function values (not necessarily normalised) over the reduced photon energy
      *        related transformed variable values.
      */
    double *fYdata;
    /** @brief The alias probabilities (not necessarily normalised) over the reduced photon energy related transformed
      *        variable values.
      */
    double *fAliasW;
    /** @brief The alias indices over the reduced photon energy related transformed variable values. */
    int    *fAliasIndx;
  };
  /** @brief Linear array to store pointers to LinAlias data structures.
    *
    * The size is RelativisticBremsModel::fNumDifferentMaterialGCuts times RelativisticBremsModel::fNumSamplingElecEnergies.
    * Some of the stored pointers are left to be nullptr (that correspond to e-/e+ kinetic energy grid points that are below
    * the gamma production cut value).
    */
  LinAlias   **fAliasData;


  /** @brief Internal data structure to store data for sampling the emitted photon energy when RelativisticBremsModel::fIsUseLinearSamplingTable is false.
  *
  *  This data structure is set up at initialisation for each different material-gamma production cut pairs
  *  over the e-/e+ kinetic energy grid to sample the emitted photon energy distribution using a combination
  *  of Walker's alias sampling and rational interpolation based inverse transform. At most we will have as many
  *  data structure as RelativisticBremsModel::fNumDifferentMaterialGCuts times RelativisticBremsModel::fNumSamplingElecEnergies
  *  and these data structure pointers are stored in the RelativisticBremsModel::fRatinAliasData linear array.
  *  However, data structures are created only for the possible e-/e+ kinetic energy - emitted photon energy
  *  combinations (i.e. for those e-/e+ kinetic energy grid points in RelativisticBremsModel::fSamplingElecEnergies
  *  that are above the gamma production cut value) and the other elements of RelativisticBremsModel::fRatinAliasData
  *  linear array are left to be nullptr.
  */
  struct RatinAlias{
    /** @brief Number of data points i.e. size of the arrays. */
    int     fNumdata;
    /** @brief Reduced photon energy related transformed variable values. */
    double *fXdata;
    /** @brief The alias probabilities over the reduced photon energy related transformed variable values. */
    double *fAliasW;
    /** @brief The cumulative distribution function over the reduced photon energy related transformed variable values. */
    double *fC;
    /** @brief Parameters for rational interpolation based inverse transform. */
    double *fA;
    /** @brief Parameters for rational interpolation based inverse transform. */
    double *fB;
    /** @brief The alias indices over the reduced photon energy related transformed variable values. */
    int    *fAliasIndx;
  };
  RatinAlias   **fRatinAliasData; //alias data for each matrial-gammacut pairs

  /** @brief A general(linera/rational) alias sampler used at run-time sampling of the emitted photon energy. */
  AliasTable  *fAliasSampler;

};

} // namespace geantphysics

#endif // RELATIVISTICBREMSMODEL_H
