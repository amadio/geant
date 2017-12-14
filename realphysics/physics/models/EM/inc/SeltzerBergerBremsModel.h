
#ifndef SELTZERBERGERBREMSMODEL_H
#define SELTZERBERGERBREMSMODEL_H

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

#include <vector>
#include <string>

namespace geantphysics {

//class Material;
//class Element;
class MaterialCuts;
class AliasTable;
class Particle;
class LightTrack;
class GLIntegral;

/**
 * @brief   Low energy Bremsstrahlung models for electron/positron.
 * @class   SeltzerBergerBremsModel
 * @author  M Novak, A Ribon
 * @date    march 2016
 *
 * Bremsstrahlung model for electron/positron based on Seltzer-Berger numerical differential cross sections
 * \cite seltzer1985bremsstrahlung \cite seltzer1986bremsstrahlung.
 * Dielectric suppression is taken into account(see more at RelativisticBremsModel::ComputeURelDXSecPerAtom).
 * The electron/positron kinetic energy range of the model is between 1[keV] - 10[GeV] but typically used between
 * 1[keV] - 1 [GeV] since at higher energies LPM suppression can be important.
 */

class SeltzerBergerBremsModel : public EMModel {
public:

  virtual void   Initialize(); // from EMModel
  virtual double ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle* particle,bool istotal=false);
  virtual double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *particle);
  virtual double ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut, double kinenergy, const Particle *particle);
  virtual int    SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td);


  virtual double MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle *part) const;



/**
* @name Constructor, destructor:
*/
//@{
   /**
    * @brief Constructor to build a model based on the numerical differential cross sections stored in files.
    *
    * There are 3 different representation of the Seltzer-Berger differential cross sections:
    * - datafileindx = 0 indicates the interpolated NIST data base which is the same as in Geant4.
    * - datafileindx = 1 indicates the NIST data base
    * - datafileindx = 2 indicates that the NRC data base is used: the electron-nuclear part of the differential
    *                    cross section is taken from the NIST data base but the electron-electron part is from
    *                    \cite tessier2008calculation (unlike the NIST free electron-electron interaction this
    *                    later computation includes the effect of binding; differencies expected at lower < 1MeV
    *                    energies).
    *
    * @param[in] iselectron   Flag to indicate that the model is for electron(true) or for psitron(false).
    * @param[in] datafileindx Index to specify the data set.
    * @param[in] modelname    Name of the model.
    */
  SeltzerBergerBremsModel(bool iselectron, int datafileindx = 0, const std::string &modelname = "eSeltzerBergerBrems");

  /** @brief Destructor. */
  virtual ~SeltzerBergerBremsModel();
//@}


  /**
   * @brief Public method to initilise the model.
   *
   * During the initialisation, Seltzer-Berger numerical differential cross section for bremsstrahlung photon emission
   * will be loaded from files, internal sampling tables for run time sampling of the emitted photon energy will be
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
   * @param[in] electronekin       Kinetic energy of the incident particle i.e. e-/e+.
   * @return                       (Restricted) Bremsstrahlung atomic cross section for the specified configuration
   *                               in internal [lenght^2] units.
   */
  double ComputeXSectionPerAtom(const Element *elem, const Material *mat, double gammaprodcutenergy, double electronekin);

  /**
   * @brief Public method to obtain (restricted) macroscopic cross sections.
   *
   * @param[in] mat                Pointer to the material that the cross section is requested.
   * @param[in] gammaprodcutenergy Kinetic energy threshold for gamma particle production.
   * @param[in] electronekin       Kinetic energy of the incident particle i.e. e-/e+.
   * @return                       (Restricted) Macroscopic bremsstrahlung cross section for the specified configuration
   *                               in internal [1/lenght] units.
   */
  double ComputeXSectionPerVolume(const Material *mat, double gammaprodcutenergy, double electronekin);

  /**
   * @brief Public method to obtain (restricted) radiative stopping power.
   *
   * @param[in] mat                Pointer to the material that the cross section is requested.
   * @param[in] gammaprodcutenergy Kinetic energy threshold for gamma particle production.
   * @param[in] electronekin       Kinetic energy of the incident particle i.e. e-/e+.
   * @return                       (Restricted) Radiative stopping power for the specified configuration in internal
   *                               [energy/lenght] units.
   */
  double ComputeDEDXPerVolume(const Material *mat, double gammaprodcutenergy, double electronekin);

  /**
   * @brief Public method to sample the emitted (restricted) bremsstrahlung photon energy.
   *
   * The emitted bremsstrahlung photon energy is always higher than the gamma particle kinetic energy production
   * threshold in the specified material-cut pair and not higher than the incident particle kinetic energy.
   *
   * @param[in] matcut     Pointer to the material-cut pair in which the interaction takes place.
   * @param[in] eekin      Kinetic energy of the incident particle i.e. e-/e+.
   * @param[in] r1         Random number distributed uniformly on the 0 1 interval.
   * @param[in] r2         Random number distributed uniformly on the 0 1 interval.
   * @param[in] r3         Random number distributed uniformly on the 0 1 interval.
   * @return               Emitted bremsstrahlung photon energy sampled from the distribution specified by the
   *                       given configuration and the model in internal [energy] units if the particle (e-/e+)
   *                       kinetic energy is higher than the gamma particle production energy threshold. Zero otherwise.
   */
  double SamplePhotonEnergy(const MaterialCuts *matcut, double eekin, double r1, double r2, double r3);

  //
  void   SamplePhotonDirection(double elenergy, double &sinTheta, double &cosTheta, double rndm);


private:
  /** @brief Copy constructor  (deleted) */
  SeltzerBergerBremsModel(const SeltzerBergerBremsModel&) = delete;
  /** @brief Operator=  (deleted) */
  SeltzerBergerBremsModel &operator=(const SeltzerBergerBremsModel&) = delete;

  /**
   * @brief Internal method to load Seltzer-Berger atomic differential cross sections for bremsstrahlung photon emission
   *        from file.
   *
   *        Used at the initialisation of the model.
   */
  void   LoadDCSData();

  /** @brief Internal method to build emitted photon energy sampling tables under <em>linear approximation of
   *         the p.d.f.</em>.
   *
   *  Used at initialisation of the model to prepare emitted photon energy related transformed variable
   *  sampling tables for all different material-gamma production cut pairs over an e-/e+ kinetic energy grid.
   *  These tables are prepared for Walker's alias method combined with <em>linear approximation of the p.d.f.</em>
   *  within the bins.
   */
  void   InitSamplingTables();

  /** @brief Internal method to build one emitted photon energy sampling tables under <em>linear approximation of
   *         the p.d.f.</em>.
   *
   *  This method is used by SeltzerBergerBremsModel::InitSamplingTables() to build emitted photon energy related
   *  sampling tables for Walker's alias method combined with <em>linear approximation of the p.d.f.</em>. Grid points
   *  are adaptively inserted such that the error between the spline and linearly interpolated p.d.f. is minimised.
   *  Number of grid points will be SeltzerBergerBremsModel::LinAlias::fNumdata which is always the maximum i.e.
   *  SeltzerBergerBremsModel::fNumSamplingPhotEnergies .
   */
  void   BuildOneLinAlias(int indxlalias, const Material *mat, double gcut);

  /**
   * @brief Correction for accounting some differencies between positrons and electrons.
   */
  double PositronCorrection(double ekinelectron, double ibeta2electron, double ephoton, double z);

  /**
   * @brief Correction for accounting some differencies between positrons and electrons.
   */
  double PositronCorrection1(double ekinelectron, double ephoton, double gcutener, double z);

private:
 /** @brief Flag to indicate if the model is for e-(true) or for e+(false). Must be set before initialisation. */
 bool     fIsElectron;                      // flag to indicate if the model is for electron(true) or for positron(flase)
 /** @brief Index of the data file where the numerical DCS are taken from. Must be set before initialisation.
   *        (see more at the description of the constructor).
   */
 int      fDataFileIndx;                    // flag to indicate which SB data should be used

 int      fNGL;
 // secondary related data
 int      fSecondaryInternalCode; // gamma GV code set at initialization

 //
 // these are for the DCS data loaded from the file
 /** @brief Maximum atomic number that numerical DCS are available in the file. */
 int      fDCSMaxZet;                       // max Z that we have DCS data in files
 /** @brief Number of electron energies that numerical DCS are available in the file. */
 int      fLoadDCSNumElectronEnergies;      // the elelectron energy grid dimension of the DCS loaded
 /** @brief Number of reduced photon energies at each electron energy that numerical DCS are available in the file. */
 int      fLoadDCSNumReducedPhotonEnergies; // the reduced photon energy grid dimension of the DCS loaded
 /** @brief The electron energy grid that numerical DCS data are given in the file. Loaded from file.
   *        Size of the array is fLoadDCSNumElectronEnergies.
   */
 double  *fLoadDCSElectronEnergyGrid;       // the electron energy grid of the DCS
 /** @brief The reduced photon energy grid that numerical DCS data are given at each electron energy. Loaded from file.
   *         Size of the array is fLoadDCSNumReducedPhotonEnergies.
   */
 double  *fLoadDCSReducedPhotonEnergyGrid;  // the reduced photon energy grid of the DCS
 /** @brief The numerical DCS data loaded from file for each element that used in the detector. Size of the array is
   *         fDCSMaxZet. Each element of the array is a pointer to a double array with size of
   *         fLoadDCSNumElectronEnergies times fLoadDCSNumReducedPhotonEnergies.
   */
 double **fLoadDCSForElements;              // container to store DCS data for the used elements

 ///
 // these are to describe and define the common electron kinetic energy grid above we build sampling tables for run-time
 // samling of emitted photon energy.
 /** @brief Number of e-/e+ kinetic energy grid points in [SeltzerBergerBremsModel::fMinElecEnergy,
   *        SeltzerBergerBremsModel::fMaxElecEnergy] (default 71).
   */
 int     fNumSamplingElecEnergies;
 /** @brief Number of transformed emitted photon energy related variable in [0,1] (default 54). */
 int     fNumSamplingPhotEnergies;
 /** @brief Minimum of the e-/e+ kinetic energy grid (default 1.0 [keV] that is the minimum available.) */
 double  fMinElecEnergy;
 /** @brief Maximum of the e-/e+ kinetic energy grid (default 10.0 [GeV] that is the maximum available.) */
 double  fMaxElecEnergy;
 /** @brief Logarithm of SeltzerBergerBremsModel::fMinElecEnergy i.e. ln(fMinElecEnergy) . */
 double  fElEnLMin;                         // log min electron energy
 /** @brief Inverse of the e-/e+ kinetic energy grid delta i.e.
   *        ln[fMaxElecEnergy/fMinElecEnergy]/(fNumSamplingElecEnergies-1)
   */
 double  fElEnILDelta;                      // 1 log delta electron energy of the electron energy grid
 /** @brief The logarithmically spaced e-/e+ kinetic energy grid.
   *
   *        Size of the array is SeltzerBergerBremsModel::fNumSamplingElecEnergies points in the
   *        [SeltzerBergerBremsModel::fMinElecEnergy, SeltzerBergerBremsModel::fMaxElecEnergy] interval.
   */
 double *fSamplingElecEnergies;             // the common electron energy grid which we build sampling tables above
 /** @brief The logarithm of SeltzerBergerBremsModel::fSamplingElecEnergies grid.
   *
   *        Size of the array is SeltzerBergerBremsModel::fNumSamplingElecEnergies points in the
   *        [ln(SeltzerBergerBremsModel::fMinElecEnergy), ln(SeltzerBergerBremsModel::fMaxElecEnergy)] interval.
   */
 double *fLSamplingElecEnergies;            // log of sampling electron/positron energies

 //
 // data to map all material-gamma production cut pair indices to local indices including only the subset of all
 // material-gamma production cut that are different. These data used only internally by the model.
 /** @brief Number of different material-gamma production cut pairs. */
 int     fNumDifferentMaterialGCuts;        // number of different matrial-gammacut pairs
 /** @brief Map from global to local material-gamma production cut indices. The size of the array is
   *        SeltzerBergerBremsModel::numMaterialCuts.
   */
 int    *fGlobalMatGCutIndxToLocal;         // maps the global mat.-cut indices to local indices that are used here

 /** @brief Internal data structure to store data for sampling the emitted photon energy.
   *
   *  This data structure is set up at initialisation for each different material-gamma production cut pairs
   *  over the e-/e+ kinetic energy grid to sample the emitted photon energy distribution using a combination
   *  of Walker's alias sampling and liner approximation. At most we will have as many data structure as
   *  SeltzerBergerBremsModel::fNumDifferentMaterialGCuts times SeltzerBergerBremsModel::fNumSamplingElecEnergies
   *  and these data structure pointers are stored in the SeltzerBergerBremsModel::fAliasData linear array.
   *  However, data structures are created only for the possible e-/e+ kinetic energy - emitted photon energy
   *  combinations (i.e. for those e-/e+ kinetic energy grid points in SeltzerBergerBremsModel::fSamplingElecEnergies
   *  that are above the gamma production cut value) and the other elements of SeltzerBergerBremsModel::fAliasData
   *  linear array are left to be nullptr.
   */
 struct LinAlias{
   /** @brief Number of data points i.e. size of the arrays = SeltzerBergerBremsModel::fNumSamplingPhotEnergies. */
   int     fNumdata;                        // size of the arrays
   /** @brief Reduced photon energy related transformed variable values. */
   double *fXdata;                          // reduced photon energies
   /** @brief The probability density function values (not necessarily normalised) over the reduced photon energy
     *        related transformed variable values.
     */
   double *fYdata;                          // p.d.f (not necessarily norm.)
   /** @brief The alias probabilities (not necessarily normalised) over the reduced photon energy related transformed
     *        variable values.
     */
   double *fAliasW;                         // alias probs (not necessarily norm.)
   /** @brief The alias indices over the reduced photon energy related transformed variable values. */
   int    *fAliasIndx; // alias indices
 };
 /** @brief Linear array to store pointers to LinAlias data structures.
   *
   * The size is SeltzerBergerBremsModel::fNumDifferentMaterialGCuts times
   * SeltzerBergerBremsModel::fNumSamplingElecEnergies. Some of the stored pointers are left to be nullptr (that
   * correspond to e-/e+ kinetic energy grid points that are below the gamma production cut value).
   */
 LinAlias   **fAliasData;                   //alias data structure for all different matrial-gammacut pairs
 /** @brief An alias sampler used at run-time sampling of the emitted photon energy. */
 AliasTable  *fAliasSampler;

 GLIntegral  *fGL;

};

}      // namespace geantphysics

#endif // SELTZERBERGERBREMSMODEL_H
