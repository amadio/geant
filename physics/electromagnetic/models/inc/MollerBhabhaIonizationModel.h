
#ifndef MOLLERBHABHAIONIZATIONMODELc1_H
#define MOLLERBHABHAIONIZATIONMODELc1_H

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
  MollerBhabhaIonizationModel(bool iselectron, const std::string &modelname="ieioniMollerBhabha");
  /** @brief Destructor. */
  virtual ~MollerBhabhaIonizationModel();
//@}

/**
* @name Implemented EMModel base class methods:
*/
//@{
    /** @brief Interface method to initilize the model. */
    virtual void   Initialize();
    virtual double ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle* particle,
                               bool istotal=false);
    virtual double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *particle);
    virtual double ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut, double kinenergy,
                                          const Particle *particle);
    virtual int    SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td);

    virtual double MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle *part) const;
//
//@}


private:

  /** @brief Copy constructor  (deleted) */
  MollerBhabhaIonizationModel(const MollerBhabhaIonizationModel&) = delete;
  /** @brief Operator=  (deleted) */
  MollerBhabhaIonizationModel &operator=(const MollerBhabhaIonizationModel&) = delete;


  /**
   * @brief Public method to obtain (restricted) collision stopping power for a given material, particle kinetic energy
   *        and e- production cut energy.
   *
   * @param[in] mat        Pointer to specify the material.
   * @param[in] pcutenergy Production cut for e- in internal energy units.
   * @param[in] primekin   Kinetic energy of the particle (e-/e+) internal energy units.
   * @return    The computed (restricted) collision stopping power value in internal [energy/lenght] units.
   */
  double ComputeDEDXPerVolume(const Material *mat, const double pcutenergy, const double primkin);

  /**
   * @brief Public method to obtain (restricted) atomic ionization cross section for a given element of a material,
   *        particle kinetic energy and e- production cut energy.
   *
   * There is no material dependent correction in this case so one could call this method for simple element by passing
   * nullptr for the material.
   *
   * @param[in] elem       Pointer to specify the target atom(element).
   * @param[in] mat        Pointer to specify the material in whihc the target element is.
   * @param[in] pcutenergy Production cut (kinetic) energy for e- in internal energy units.
   * @param[in] primekin   Kinetic energy of the particle (e-/e+) internal energy units.
   * @return    The computed (restricted) atomic ionization scross section in internal [lenght^2] units.
   */
  double ComputeXSectionPerAtom(const Element *elem, const Material *mat, const double pcutenergy, const double primekin);

  /**
   * @brief Public method to obtain (restricted) macroscopic ionization cross section for a given material, particle
   *        kinetic energy and e- production cut energy.
   *
   * @param[in] mat        Pointer to specify the material.
   * @param[in] pcutenergy Production cut (kinetic) energy for e- in internal energy units.
   * @param[in] primekin   Kinetic energy of the particle (e-/e+) internal energy units.
   * @return    The computed (restricted) macroscopic ionization cross section in internal [1/lenght] units.
   */
  double ComputeXSectionPerVolume(const Material *mat, const double pcutenergy, const double primekin);

  /**
    * @brief Internal method to sample (restricted) kinetic energy transfered to the electron in Moller/Bhabha
    *        scattering (from sampling tables).
    *
    *  @param[in] matcut     Pointer to the material-production cut object in which the interaction happens.
    *  @param[in] primekin   Kinetic energy of the primary particle i.e. e-/e+.
    *  @param[in] r1         Random number distributed uniformly in [0,1].
    *  @param[in] r2         Random number distributed uniformly in [0,1].
    *  @param[in] r3         Random number distributed uniformly in [0,1].
    *  @return    Sample of kinetic energy transfered to the electron (in internal [energy] units) in Moller/Bhabha
    *             interaction.
    */
  double SampleEnergyTransfer(const MaterialCuts *matcut, const double primekin, const double r1, const double r2,
                              const double r3);

  /**
    * @brief Internal method to sample (restricted) kinetic energy transfered to the electron in Moller/Bhabha
    *        scattering (with rejection).
    *
    *  @param[in] matcut     Pointer to the material-production cut object in which the interaction happens.
    *  @param[in] primekin   Kinetic energy of the primary particle i.e. e-/e+.
    *  @param[in] td         Pointer to the (thread local) Geant task data (used to generate random numbers).
    *  @return    Sample of kinetic energy transfered to the electron (in internal [energy] units) in Moller/Bhabha
    *             interaction.
    */
  double SampleEnergyTransfer(const MaterialCuts *matcut, const double primekin, const Geant::GeantTaskData* td);

  /** @brief Internal method to build energy transfer sampling tables under <em>linear approximation of
   *         the p.d.f.</em>.
   *
   *  Used at initialisation of the model to prepare energy transfer related transformed variable
   *  sampling tables for all different material-electron production cut pairs over an e-/e+ kinetic energy grid.
   *  These tables are prepared for Walker's alias method combined with <em>linear approximation of the p.d.f.</em>
   *  within the e-/e+ kinetic energy bins.
   */
  void   InitSamplingTables();

  void   BuildSamplingTableForMaterialCut(const MaterialCuts *matcut, int indxlocal);

  /** @brief Internal method to compute (restricted) ionization atomic cross section per electron.
    *
    * @param[in] pcutenergy Electron kinetic energy production threshold.
    * @param[in] primekin   Kinetic energy of the primary electron.
    * @return               Cross section per electron.
    */
  double ComputeXSectionPerElectron(const double pcutenergy, const double primekin);

  /** @brief Internal method to compute (unnormalized) probaility density function of kinetic energy transfer related
    *        transformed variable in Moller scattering.
    *
    * This method is used at initialisation by #BuildOneLinAlias to prepare energy transfer
    * related sampling tables for run-time sampling of energy transfer in Moller scattering.
    *
    * @param[in] xi         Kinetic energy transfer related transformed variable.
    * @param[in] pcutenergy Electron kinetic energy production threshold.
    * @param[in] primekin   Kinetic energy of the primary electron.
    * @return    (Unnormalized) probability of the given kinetic energy transfer related transformed variable.
    */
  double ComputeMollerPDF(const double xi, const double pcutenergy, const double primekin);

  /** @brief Internal method to compute (unnormalized) probaility density function of kinetic energy transfer related
    *        transformed variable in Bhabha scattering.
    *
    * This method is used at initialisation by #BuildOneLinAlias to prepare energy transfer
    * related sampling tables for run-time sampling of energy transfer in Bhabha scattering.
    *
    * @param[in] xi            Kinetic energy transfer related transformed variable.
    * @param[in] prodcutenergy Electron kinetic energy production threshold.
    * @param[in] particleekin  Kinetic energy of the primary positron.
    * @return    (Unnormalized) probability of the given kinetic energy transfer related transformed variable.
    */
  double ComputeBhabhaPDF(const double xi, const double pcutenergy, const double primekin);


  void   ClearSamplingTables();

  struct LinAlias {
    LinAlias(int num) {fXdata.resize(num); fYdata.resize(num); fAliasW.resize(num); fAliasIndx.resize(num);}
    /** @brief Reduced photon energy related transformed variable values. */
    std::vector<double> fXdata;
    /** @brief The probability density function values (not necessarily normalised) over the reduced photon energy
      *        related transformed variable values.
      */
    std::vector<double> fYdata;
    /** @brief The alias probabilities (not necessarily normalised) over the reduced photon energy related transformed
      *        variable values.
      */
    std::vector<double> fAliasW;
    /** @brief The alias indices over the reduced photon energy related transformed variable values. */
    std::vector<int>    fAliasIndx;
  };

  struct AliasDataMaterialCuts {
    AliasDataMaterialCuts(int ntables, double lemin, double ildel) : fNData(ntables), fLogEmin(lemin), fILDelta(ildel) {
      fAliasData.resize(ntables,nullptr);
    }
    int    fNData;
    double fLogEmin;
    double fILDelta;
    std::vector<LinAlias*>  fAliasData;
  };

private:
  /** @brief Flag to indicate if the model is for electron or positron interaction. */
  bool    fIsElectron;

  int     fSecondaryInternalCode;

  int     fSTNumPrimaryEnergyPerDecade;   // ST=> sampling tables
  int     fSTNumSamplingElecEnergies;     // ST=> sampling tables

  AliasTable*                          fAliasSampler;
  std::vector<int>                     fGlobalMatECutIndxToLocal;
  std::vector<AliasDataMaterialCuts*>  fSamplingTables;
};

} // namespace geantphysics

#endif // MOLLERBHABHAIONIZATIONMODELc1_H
