
#ifndef MollerBhabhaIONIZATIONMODEL_H
#define MollerBhabhaIONIZATIONMODEL_H

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

namespace geantphysics {

/**
 * @brief   Ionization model for electron/positron.
 * @class   MollerBhabhaIonizationModel
 * @author  M Novak, A Ribon
 * @date    may 2016
 *
 * Model for Moller \cite moller1932theorie \cite crawford1970electron \f$[e^-+e^-\to e^-+e^-]\f$ and Bhabha
 * \cite bhabha1936scattering \cite crawford1970electron \f$[e^++e^-\to e^++e^-]\f$ scattering.
 */

//class Material;
class MaterialCuts;
//class Element;
class AliasTable;
class Particle;
class LightTrack;

class MollerBhabhaIonizationModel : public EMModel {
public:

  virtual void Initialize(); // from EMModel
  virtual double ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle* particle, bool istotal=false);
  virtual double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *particle);
  virtual double ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut, double kinenergy, const Particle *particle);
  virtual int    SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td);

  virtual double MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle *part) const;

public:
/**
* @name Constructor, destructor:
*/
//@{
   /**
    * @brief Constructor.
    *
    * @param[in] iselectron  Flag to indicate that the model is for electron or for psitron.
    * @param[in] modelname   Name of the model.
    */
  MollerBhabhaIonizationModel(bool iselectron, const std::string &modelname="eIonization");
  /** @brief Destructor. */
 ~MollerBhabhaIonizationModel();
//@}

  /**
   * @brief Public method to initilise the model.
   *
   * Internal sampling tables for run time sampling of the kinetic energy transfered to the e- will be created during
   * the initialisation. This method always need to be invoked between the construction and usage of the model.
   *
   */
  void   Initialise();

  /**
   * @brief Public method to obtain (restricted) collision stopping power for a given material, particle kinetic energy
   *        and e- production cut energy.
   *
   * @param[in] mat            Pointer to specify the material.
   * @param[in] prodcutenergy  Production cut for e- in internal energy units.
   * @param[in] particleekin   Kinetic energy of the particle (e-/e+) internal energy units.
   * @return    The computed (restricted) collision stopping power value in internal [energy/lenght] units.
   */
  double ComputeDEDXPerVolume(const Material *mat, double prodcutenergy, double particleekin);

  /**
   * @brief Public method to obtain (restricted) atomic ionization cross section for a given element of a material,
   *        particle kinetic energy and e- production cut energy.
   *
   * There is no material dependent correction in this case so one could call this method for simple element by passing
   * nullptr for the material.
   *
   * @param[in] elem           Pointer to specify the target atom(element).
   * @param[in] mat            Pointer to specify the material in whihc the target element is.
   * @param[in] prodcutenergy  Production cut (kinetic) energy for e- in internal energy units.
   * @param[in] particleekin   Kinetic energy of the particle (e-/e+) internal energy units.
   * @return    The computed (restricted) atomic ionization scross section in internal [lenght^2] units.
   */
  double ComputeXSectionPerAtom(const Element *elem, const Material *mat, double prodcutenergy, double particleekin);

  /**
   * @brief Public method to obtain (restricted) macroscopic ionization cross section for a given material, particle
   *        kinetic energy and e- production cut energy.
   *
   * @param[in] mat            Pointer to specify the material.
   * @param[in] prodcutenergy  Production cut (kinetic) energy for e- in internal energy units.
   * @param[in] particleekin   Kinetic energy of the particle (e-/e+) internal energy units.
   * @return    The computed (restricted) macroscopic ionization cross section in internal [1/lenght] units.
   */
  double ComputeXSectionPerVolume(const Material *mat, double prodcutenergy, double particleekin);

  /**
    * @brief Public method to sample (restricted) kinetic energy transfered to the electron in Moller/Bhabha scattering.
    *
    *  @param[in] matcut     Pointer to the material-production cut object in which the interaction happens.
    *  @param[in] primekin   Kinetic energy of the primary particle i.e. e-/e+.
    *  @param[in] r1         Random number distributed uniformly in [0,1].
    *  @param[in] r2         Random number distributed uniformly in [0,1].
    *  @param[in] r3         Random number distributed uniformly in [0,1].
    *  @return    Sample of kinetic energy transfered to the electron (in internal [energy] units) in Moller/Bhabha
    *             interaction, if it is a kinematically allowed combination of the primary particle kinetic energy -
    *             current electron production threshold. Zero otherwise.
    */
  double SampleEnergyTransfer(const MaterialCuts *matcut, double primekin, double r1, double r2, double r3);


private:
  /** @brief Internal method to build energy transfer sampling tables under <em>linear approximation of
   *         the p.d.f.</em>.
   *
   *  Used at initialisation of the model to prepare energy transfer related transformed variable
   *  sampling tables for all different material-electron production cut pairs over an e-/e+ kinetic energy grid.
   *  These tables are prepared for Walker's alias method combined with <em>linear approximation of the p.d.f.</em>
   *  within the e-/e+ kinetic energy bins.
   */
  void   InitSamplingTables();

  /** @brief Internal method to build energy transfer sampling tables for one given material-electron production kinetic
   *         energy threshold under <em>linear approximation of the p.d.f.</em>.
   *
   *  This method is used by MollerBhabhaIonizationModel::InitSamplingTables() to build energy transfer related
   *  sampling tables for Walker's alias method combined with <em>linear approximation of the p.d.f.</em>. The
   *  kinematically allowed energy transfer range is transformed to the [0,1] inteval and there are a fixed number of
   *  grid points. The transformed variable grid is determined by the value of
   *  MollerBhabhaIonizationModel::fNumSamplingElecEnergies member.
   */
  void   BuildOneLinAlias(int indxlalias, double elecprodcut);

  /** @brief Internal method to compute (restricted) ionization atomic cross section per electron.
    *
    * This method is used by MollerBhabhaIonizationModel::ComputeXSectionPerAtom,
    * MollerBhabhaIonizationModel::ComputeXSectionPerVolume to compute (restricted) atomic and macroscopic ionization
    * cross sections for Moller and Bhabha scattering.
    */
  double ComputeXSectionPerElectron(double prodcutenergy, double particleekin);

  /** @brief Internal method to compute (unnormalized) probaility density function of kinetic energy transfer related
    *        transformed variable in Moller scattering.
    *
    * This method is used at initialisation by MollerBhabhaIonizationModel::BuildOneLinAlias to prepare energy transfer
    * related sampling tables for run-time sampling of energy transfer in Moller scattering.
    *
    * @param[in] xi             Kinetic energy transfer related transformed variable.
    * @param[in] prodcutenergy Electron kinetic energy production threshold.
    * @param[in] particleekin  Kinetic energy of the primary electron.
    * @return    (Unnormalized) probability of the given kinetic energy transfer related transformed variable.
    */
  double ComputeMollerPDF(double xi, double prodcutenergy, double particleekin);

  /** @brief Internal method to compute (unnormalized) probaility density function of kinetic energy transfer related
    *        transformed variable in Bhabha scattering.
    *
    * This method is used at initialisation by MollerBhabhaIonizationModel::BuildOneLinAlias to prepare energy transfer
    * related sampling tables for run-time sampling of energy transfer in Bhabha scattering.
    *
    * @param[in] xi            Kinetic energy transfer related transformed variable.
    * @param[in] prodcutenergy Electron kinetic energy production threshold.
    * @param[in] particleekin  Kinetic energy of the primary positron.
    * @return    (Unnormalized) probability of the given kinetic energy transfer related transformed variable.
    */
  double ComputeBhabhaPDF(double xi, double prodcutenergy, double particleekin);


// data members
private:
  /** @brief Flag to indicate if the model is for electron or positron interaction. */
  bool        fIsElectron;

  int fSecondaryInternalCode; // internal code i.e. GV code of the secondary partcile i.e. e-

  ///
  // these are to describe and define the common e-/e+ kinetic energy grid above we build sampling tables for run-time
  // samling of emitted photon energy.
  /** @brief Number of e-/e+ kinetic energy grid points in [MollerBhabhaIonizationModel::fMinPrimEnergy,
    *        MollerBhabhaIonizationModel::fMaxPrimEnergy].
    */
  int fNumSamplingPrimEnergies;
  /** @brief Number of energy transfer related transformed variable in [0,1]. */
  int fNumSamplingElecEnergies;
  /** @brief Minimum of the e-/e+ kinetic energy grid. */
  double  fMinPrimEnergy;
  /** @brief Maximum of the e-/e+ kinetic energy grid. */
  double  fMaxPrimEnergy;
  /** @brief Logarithm of MollerBhabhaIonizationModel::fMinPrimEnergy i.e. ln(fMinPrimEnergy) . */
  double  fPrimEnLMin;
  /** @brief Inverse of the e-/e+ kinetic energy grid delta i.e.
    *        ln[fMaxPrimEnergy/fMinPrimEnergy]/(fNumSamplingElecEnergies-1)
    */
  double  fPrimEnILDelta;
  /** @brief The logarithmically spaced e-/e+ kinetic energy grid.
    *
    *        Size of the array is MollerBhabhaIonizationModel::fNumSamplingElecEnergies points in the
    *        [MollerBhabhaIonizationModel::fMinPrimEnergy, MollerBhabhaIonizationModel::fMaxPrimEnergy] interval.
    */
  double *fSamplingPrimEnergies;
  /** @brief The logarithm of MollerBhabhaIonizationModel::fSamplingElecEnergies grid.
    *
    *        Size of the array is MollerBhabhaIonizationModel::fNumSamplingElecEnergies points in the
    *        [ln(MollerBhabhaIonizationModel::fMinPrimEnergy), ln(MollerBhabhaIonizationModel::fMaxPrimEnergy)] interval.
    */
  double *fLSamplingPrimEnergies;

  //
  // data to map all material-production cut pair indices to local indices including only the subset of all
  // electron production cuts that are different. These data used only internally by the model.
  /** @brief Number of different electron production cut pairs. */
  int     fNumDifferentElecCuts;
  /** @brief Map from global to local material-production cut indices. The size of the array is
    *        MollerBhabhaIonizationModel::numMaterialCuts.
    */
  int    *fGlobalMatCutIndxToLocal;


  /** @brief Internal data structure to store data for sampling the energy transfer in Moller/Bhabha scattering.
    *
    *  This data structure is set up at initialisation for each different material-electron production cut pairs
    *  over the e-/e+ kinetic energy grid to sample the energy transfered to the electron using a combination
    *  of Walker's alias sampling and liner approximation. At most we will have as many data structure as
    *  MollerBhabhaIonizationModel::fNumDifferentElecCuts times MollerBhabhaIonizationModel::fNumSamplingPrimEnergies
    *  and these data structure pointers are stored in the MollerBhabhaIonizationModel::fAliasData linear array.
    *  However, data structures are created only for the possible e-/e+ kinetic energy - energy transfer
    *  combinations (i.e. for those e-/e+ kinetic energy grid points in
    *  MollerBhabhaIonizationModel::fSamplingPrimEnergies that results in kinematically allowed combinations) and the
    *  other elements of MollerBhabhaIonizationModel::fAliasData linear array are left to be nullptr.
    */
  struct LinAlias{
    /** @brief Number of data points i.e. size of the arrays = SeltzerBergerBremsModel::fNumSamplingElecEnergies. */
    int     fNumdata;
    /** @brief Energy transfer related transformed variable values. */
    double *fXdata;
    /** @brief The probability density function values (not necessarily normalised) over the energy transfer related
      *        transformed variable values.
      */
    double *fYdata;
    /** @brief The alias probabilities (not necessarily normalised) over the energy transfer related transformed variable
      *        values.
      */
    double *fAliasW;
    /** @brief The alias indices over the energy transfer related transformed variable values. */
    int    *fAliasIndx;
  };
  /** @brief Linear array to store pointers to LinAlias data structures.
    *
    * The size is MollerBhabhaIonizationModel::fNumDifferentElecCuts times
    * MollerBhabhaIonizationModel::fNumSamplingPrimEnergies. Some of the stored pointers are left to be nullptr (that
    * correspond to e-/e+ kinetic energy grid points that results in a kinematically not allowed e-/e+ kinetic energy
    * electron production cut combination).
    */
  LinAlias   **fAliasData;
  /** @brief An alias sampler used at run-time sampling of the kinetic energy transfered to the electron. */
  AliasTable  *fAliasSampler;
};

} // namespace geantphysics

#endif // MollerBhabhaIONIZATIONMODEL_H
