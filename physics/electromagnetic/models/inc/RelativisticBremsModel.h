
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
class GLIntegral;

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
/**
* @name Constructor, destructor:
*/
//@{
     /**
      * @brief Constructor.
      *
      * @param[in] modelname Name of the model.
      */
    RelativisticBremsModel(const std::string &modelname="eRelativisticBrems");
     /** @brief Destructor. */
    virtual ~RelativisticBremsModel();
//@}

/**
* @name Implemented EMModel base class methods:
*/
//@{
    /** @brief Interface method to initilize the model. */
    virtual void   Initialize();

    virtual double ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle* particle, bool istotal=false);
    virtual double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *particle);
    virtual double ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut, double kinenergy, const Particle *particle);
    virtual int    SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td);

    virtual double MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle *part) const;
//
//@}

/**
 * @name Field related setters/getters:
 */
//@{
  /** @brief Get the flag that indicates if LPM suppression is active (true by default). */
  bool  GetLPMFlag()  {return fIsUseLPM;}
  /** @brief Set the flag that indicates if LPM suppression is active (true by default).
    *
    * Must be set before initialisation if not the default value required.
    */
  void  SetLPMFlag(bool islpm) {fIsUseLPM = islpm;}
//@}


private:

  /** @brief Copy constructor  (deleted) */
  RelativisticBremsModel(const RelativisticBremsModel&) = delete;
  /** @brief Operator=  (deleted) */
  RelativisticBremsModel &operator=(const RelativisticBremsModel&) = delete;

  /**
   * @brief Internal method to obtain (restricted) radiative stopping power.
   *
   * @param[in] mat         Pointer to the material that the cross section is requested.
   * @param[in] gammacut    Kinetic energy threshold for gamma particle production.
   * @param[in] eekin       Kinetic energy of the incident particle i.e. e-/e+.
   * @return                (Restricted) Radiative stopping power for the specified configuration in internal
   *                        [energy/lenght] units.
   */
  double ComputeDEDXPerVolume(const Material *mat, double gammacut, double eekin);

  /**
   * @brief Internal method to obtain (restricted) macroscopic cross sections.
   *
   * @param[in] mat         Pointer to the material that the cross section is requested.
   * @param[in] gcut        Kinetic energy threshold for gamma particle production.
   * @param[in] eekin       Kinetic energy of the incident particle i.e. e-/e+.
   * @return                (Restricted) Macroscopic bremsstrahlung cross section for the specified configuration
   *                        in internal [1/lenght] units.
   */
  double ComputeXSectionPerVolume(const Material *mat, double gcut, double eekin);

  /**
   * @brief Internal method to obtain (restricted) atomic cross sections.
   *
   * Since the model includes some material dependent corrections, for consistency reasons one also needs to provide
   * the material that the element, that the atomic cross section is requested, belongs to.
   *
   * @param[in] elem        Pointer to the element object that the atomic cross section is required.
   * @param[in] mat         Pointer to the material that the element blongs to.
   * @param[in] gcut        Kinetic energy threshold for gamma particle production.
   * @param[in] eekin       Kinetic energy of the incident particle i.e. e-/e+.
   * @return                (Restricted) Bremsstrahlung atomic cross section for the specified configuration
   *                        in internal [lenght^2] units.
   */
  double ComputeXSectionPerAtom(const Element *elem, const Material *mat, double gcut, double eekin);

  /**
   * @brief Internal method to sample the emitted (restricted) bremsstrahlung photon energy (with sampling table).
   *
   * The emitted bremsstrahlung photon energy is always higher than the gamma particle kinetic energy production
   * threshold in the specified material-cut pair and not higher than the incident particle kinetic energy.
   *
   * @param[in] matcut    Pointer to the material-cut pair in which the interaction takes place.
   * @param[in] eekin     Kinetic energy of the incident particle i.e. e-/e+.
   * @param[in] r1        Random number distributed uniformly on the 0 1 interval.
   * @param[in] r2        Random number distributed uniformly on the 0 1 interval.
   * @param[in] r3        Random number distributed uniformly on the 0 1 interval.
   * @return              An emitted bremsstrahlung photon energy (sampled from the distribution specified by
   *                      given configuration and the model) in internal [energy] units.
   */
  double SamplePhotonEnergy(const MaterialCuts *matcut, double eekin, double r1, double r2, double r3);

  /**
   * @brief Internal method to sample the emitted (restricted) bremsstrahlung photon energy (with rejection).
   *
   * The emitted bremsstrahlung photon energy is always higher than the gamma particle kinetic energy production
   * threshold in the specified material-cut pair and not higher than the incident particle kinetic energy.
   *
   * @param[in] matcut    Pointer to the material-cut pair in which the interaction takes place.
   * @param[in] eekin     Kinetic energy of the incident particle i.e. e-/e+.
   * @param[in] td        Pointer to the (thread local) Geant task data (used to generate random numbers).
   * @return              An emitted bremsstrahlung photon energy (sampled from the distribution specified by
   *                      given configuration and the model) in internal [energy] units.
   */
  double SamplePhotonEnergy(const MaterialCuts *matcut, double eekin, Geant::GeantTaskData* td);

  /**
   * @brief Internal method to sample the emitted (restricted) bremsstrahlung photon energy (with rejection).
   *
   * The emitted bremsstrahlung photon energy is always higher than the gamma particle kinetic energy production
   * threshold in the specified material-cut pair and not higher than the incident particle kinetic energy.
   *
   * @param[in]     eekin     Kinetic energy of the incident particle i.e. e-/e+.
   * @param[in,out] sinTheta  Result of the emitted photon direction sampling.
   * @param[in,out] cosTheta  Result of the emitted photon direction sampling.
   * @param[in]     rndm      A random number distributed uniformly on [0,1).
   */
  void   SamplePhotonDirection(double eekin, double &sinTheta, double &cosTheta, double rndm);


  void   InitElementData();
  double ComputeDXSecPerAtom(double egamma, double etotal, double zet);
  double ComputeURelDXSecPerAtom(double egamma, double etotal, double lpmenergy, double densitycor, int izet);

  void   ComputeScreeningFunctions(double &phi1, double &phi1m2, double &xsi1, double &xsi1m2, const double gamma,
                                   const double epsilon);
  //
  void   ComputeLPMfunctions(double &funcXiS, double &funcGS, double &funcPhiS, const double lpmenergy,
                             const double egamma, const double etot, const double densitycor, const int izet);
  void   ComputeLPMGsPhis(double &funcGS, double &funcPhiS, const double varShat);
  void   InitLPMFunctions();
  void   GetLPMFunctions(double &lpmGs, double &lpmPhis, const double s);
  //
  void   ClearSamplingTables();
  void   InitSamplingTables();
  void   BuildSamplingTableForMaterialCut(const MaterialCuts *matCut, int indxlocal);



  struct ElementData {
    /** @brief \f$ \ln(Z) \f$  */
    double  fLogZ;
    /** @brief \f$ \ln(Z)/3 + f_c \f$  */
    double  fFz;
    /** @brief \f$ ((Fel-fc)+Finel*invZ)\f$  */
    double  fZFactor1;
    /** @brief \f$ (1.0+invZ)/12  \f$  */
    double  fZFactor2;
    // LPM variables
    double  fVarS1;
    double  fILVarS1;
    double  fILVarS1Cond;
    // constant factors to the screening function evaluations
    double  fGammaFactor;
    double  fEpsilonFactor;
  };

  struct LPMFuncs {
    LPMFuncs() : fIsInitialized(false), fSDelta(0.01), fSLimit(2.) {}
    bool   fIsInitialized;
    double fSDelta;
    double fSLimit;
    std::vector<double>  fLPMFuncG;
    std::vector<double>  fLPMFuncPhi;
  };

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

// data members
private:

  static const long   gMaxZet = 121;
  static const double gFelLowZet[8];
  static const double gFinelLowZet[8];
  static const double gLPMFactor;
  static const double gDensityFactor;
  static LPMFuncs     gLPMFuncs;
  static std::vector<ElementData*>  gElementData;

  bool    fIsUseLPM;
  int     fNGL;
  int     fSecondaryInternalCode;
  int     fSTNumElectronEnergyPerDecade;   // ST=> sampling tables
  int     fSTNumSamplingPhotEnergies;      // ST=> sampling tables

  GLIntegral *fGL;

  AliasTable*                          fAliasSampler;
  std::vector<int>                     fGlobalMatGCutIndxToLocal;
  std::vector<AliasDataMaterialCuts*>  fSamplingTables;
};


}       // namespace geantphysics

#endif  // RELATIVISTICBREMSMODEL_H
