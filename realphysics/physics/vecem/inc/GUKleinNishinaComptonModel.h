#ifndef GUKLENINNISHINACOMPTONMODEL_H
#define GUKLENINNISHINACOMPTONMODEL_H

#include "EMModel.h"
#include "base/VecPhys.h"

//vecphys
#include "ComptonKleinNishina.h"

#include <string>

namespace geantphysics {

  using namespace vecphys;

/**
 * @brief   KleinNishina model for gamma.
 * @class   GUKleinNishinaComptonModel
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

class GUKleinNishinaComptonModel : public EMModel {
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
    * @param[in] modelname   Name of the model.
    */
  GUKleinNishinaComptonModel(const std::string &modelname="gCompton");
  /** @brief Destructor. */
 ~GUKleinNishinaComptonModel();
//@}

  /**
   * @brief Public method to initialise the model.
   *
   * Creates internal sampling tables for run time sampling of the kinetic energy transfered to the e-
   * Must be invoked between construction and using of the model.
   */
  void   Initialise();

  /**
   * @brief Public method to obtain (restricted) macroscopic ionization cross section for a given
   *        material, particle, kinetic energy and e- production cut energy.
   *
   * @param[in] mat            Pointer to specify the material.
   * @param[in] prodcutenergy  Production cut (kinetic) energy for e- in internal energy units.
   * @param[in] particleekin   Kinetic energy of the particle (e-/e+) internal energy units.
   * @return    The computed (restricted) macroscopic ionization cross section in internal [1/lenght] units.
   */
  double ComputeXSectionPerVolume(const Material *mat, double prodcutenergy, double particleekin);


private:

  void   InitSamplingTables();

// data members
private:
  vecphys::ComptonKleinNishina *fVectorModel;

  int fSecondaryInternalCode; // internal code i.e. GV code of the secondary partcile i.e. e-

  ///
  // The numbers below are for the common kinetic energy grid  (??)
  //   we build sampling tables for run-time samling of emitted photon energy.
  /** @brief Number of e-/e+ kinetic energy grid points in [fMinPrimEnergy,
    *        fMaxPrimEnergy].
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
