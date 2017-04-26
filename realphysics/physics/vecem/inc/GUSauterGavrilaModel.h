#ifndef GUSAUTERGABRILAMODEL_H
#define GUSAUTERGABRILAMODEL_H

#include "EMModel.h"
#include "base/VecPhys.h"

//vecphys
#include "PhotoElectronSauterGavrila.h"

#include <string>

namespace geantphysics {

  using namespace vecphys;

/**
 * @brief   PhotoElectroic Effect model with the SauterGavrila angular distributionfor gamma.
 * @class   GUSauterGavrilaModel
 * @date    march 2017
 *
 * Model for KleinNishina
 */

class Material;
class MaterialCuts;
class Element;
class AliasTable;
class Particle;
class LightTrack;

class GUSauterGavrilaModel : public EMModel
{

 public:

  void Initialize() override final; // from EMModel

  double ComputeXSectionPerAtom(const Element *elem,
                                const MaterialCuts *matcut,
                                double kinenergy,
                                const Particle *particle) override final;
  
  double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *particle)
     override final;
  
  int    SampleSecondaries(LightTrack &track, std::vector<LightTrack> &sectracks, Geant::GeantTaskData *td)
     override final;

  double MinimumPrimaryEnergy(const MaterialCuts * /*matcut*/, const Particle * /*part*/) const override final;

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
  GUSauterGavrilaModel(bool iselectron, const std::string &modelname="gPhotoElectric");
  /** @brief Destructor. */
 ~GUSauterGavrilaModel();
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
// data members
private:
  vecphys::PhotoElectronSauterGavrila *fVectorModel;

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
  /** @brief Number of all material-production cut pairs. */
  int     fNumMaterialCuts;
  /** @brief Number of different electron production cut pairs. */
  int    *fGlobalMatCutIndxToLocal;

};

} // namespace geantphysics

#endif // MollerBhabhaIONIZATIONMODEL_H
