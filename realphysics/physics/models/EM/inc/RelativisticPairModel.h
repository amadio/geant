
#ifndef RELATIVISTICPAIRMODEL_H
#define RELATIVISTICPAIRMODEL_H

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
 * @brief   Model for conversion of gamma particles into e-/e+ pair including LPM suppresion effects.
 * @class   RelativisticPairModel
 * @author  F Hariri, M Novak
 * @date    November 2017
 *
 * The model is based on the Bethe-Heitler \cite bethe1934stopping differential cross section (DCS) corrected for
 * various effects like screening, Coulomb correction, conversion in the field of atomic electrons (see more at the
 * #ComputeDXsectionPerAtom() method). However, triplet production is not generated explicitely. The Landau-Pomeranchuk-
 * Migdal (LPM) \cite migdal1956bremsstrahlung  suppression effect is included in the high energy DCS (see more at the
 * #ComputeLPMDXsectionPerAtom() method). The atomic and macroscopic cross sections are obtained by the direct
 * integration of the DCS and the final state energies (total energy transfered to one of the e-/e+ pair) are computed
 * based on the DCS (with or without LPM effect depending on the photon energy). The model can be used from 50 [MeV].
 * However, it is recommended to use the alternative BetheHeitlerPairModel at lower energies and this model above 80
 * [GeV].
 */



//class Material;
class MaterialCuts;
//class Element;
class AliasTable;
class Particle;
class LightTrack;


class RelativisticPairModel : public EMModel {
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

  RelativisticPairModel(const std::string &modelname="PairRelativisticLPM");
  /** @brief Destructor. */
  ~RelativisticPairModel();
//@}

/**
* @name Implemented EMModel base class methods:
*/
//@{
  /** @brief Interface method to initilize the model. */
  virtual void   Initialize();

  /**
    * @brief Interface method to obtain macroscopic cross sections.
    *
    * @param[in] matcut      Pointer to the MaterialCuts object in which the macroscopic cross section must be computed.
    * @param[in] kinenergy   Kinetic energy of the gamma particle at which the macroscopic cross section must be computed.
    * @return    Macroscopic pair-production cross section in internal [1/length] units.
    */
  virtual double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *);
  /**
    * @brief Interface method to obtain atomic cross sections.
    *
    * @param[in] elem        Pointer to the Element object for which the atomic cross section must be computed.
    * @param[in] kinenergy   Kinetic energy of the gamma particle at which the atomic cross section must be computed.
    * @return    Atomic pair-production cross section in internal [length^2] units.
    */
  virtual double ComputeXSectionPerAtom(const Element *elem, const MaterialCuts*, double kinenergy, const Particle*);
  /**
    * @brief Interface method to generate final state of the interaction.
    *
    * @param[in,out] track     Primary track. At input, it stores the pre-interaction primary particle properties and
    *                          some information about the current material-cut couple. It is updated by the method and
    *                          it stores the post-interaction primary track properties at output.
    * @param[in,out] td        Pointer to a Geant thread local data object. At output, its fPhysicsData object will store
    *                          the seconadry tracks generated in the interaction.
    * @return                  Number of secondary tracks generated in the interaction.
    */
  virtual int    SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td);
//
//@}

/**
* @name Model specific public methods for customize sampling tables (if they are used):
*/
//@{
  /** @brief Public method to set the number of primary photon energy grid points per decade for sampling tables.
    *        Will be used only if sampling tables were required and must be set before initialisation.
    */
  int   GetNumberOfPhotonEnergiesPerDecade()   const { return fNumSamplingPrimEnergiesPerDecade; }
  /** @brief Public method to get the number of primary photon energy grid points per decade for sampling tables.*/
  void  SetNumberOfPhotonEnergiesPerDecade(int val)  { fNumSamplingPrimEnergiesPerDecade = val;  }
  /** @brief Public method to set the number of discrete samples used by the sampling tables.
    *        Will be used only if sampling tables were required and must be set before initialisation.
    */
  int   GetNumberOfDiscretePDFSamples()        const { return fNumSamplingEnergies;              }
  /** @brief Public method to get the number of discrete samples used by the sampling tables.*/
  void  SetNumberOfDiscretePDFSamples(int val)       { fNumSamplingEnergies = val;               }
//@}



private:
  /** @brief Internal method to initilise the model.*/
  void   InitialiseModel();

  /** @brief Internal method to build data collection of some frequently used target atom dependent variables. */
  void   InitialiseElementData();

  /**
   * @brief Internal method to build sampling table data structures for fast run-time sampling of the reduced total
   *        energy transfered to one of the e-/e+ pair for all possible target atom.
   *
   * Sampling tables are built at initialization for all materials that appear in regions in which the model is active
   * if sampling tables were requested. First the common, discrete, primary photon energy grid is generated between
   * #fMinPrimEnergy and #fMaxPrimEnergy using #fNumSamplingPrimEnergies discrete points on log-linear scale. Then
   * sampling tables are built for each material at all dicrete primary photon energies that are stored in separate
   * #RatinAliasDataPerMaterial data structures for the different target Materials. Such data structures for all
   * possible target materials will be stored in the #fRatinAliasDataForAllMaterials data structure.
   */
  void   InitSamplingTables();

  /**
   * @brief Internal method to sample reduced total energy transfered to one of the e-/e+ pair using sampling tables
   *        prepared at initialization}.
   *
   * @param[in] primekin Gamma photon energy in internal [energy] units.
   * @param[in] matindx  Index of the target material.
   * @param[in] r1       Random number distributed uniformly in [0,1].
   * @param[in] r2       Random number distributed uniformly in [0,1].
   * @param[in] r3       Random number distributed uniformly in [0,1].
   * @return             The sampled reduced total energy transfered to one of the e-/e+ pair.
   */
  double SampleTotalEnergyTransfer(double primekin, int matindx, double r1, double r2, double r3);

  /**
   * @brief Internal method to sample reduced total energy transfered to one of the e-/e+ pair using rejection.
   *
   * @param[in]  epsmin      Minimum value of the reduced total energy transfer
   *                         \f$ \epsilon_{\text{min}}=\text{max}[\epsilon_1(Z_{min},\epsilon_0] \f$.
   * @param[in]  eps0        Kinematical minimum of the reduced total energy transfer \f$ m_e c^2/E_{gamma}\f$
   *                         \f$ \epsilon_0 \equiv m_ec^2/E_{\gamma} \f$.
   * @param[in]  fz          The Coulomb correction factor \f$ F(Z) \f$ (see #ComputeDXsectionPerAtom()).
   * @param[in]  z23         Target atom dependent constant \f$ Z^{2/3}\f$.
   * @param[in]  egamma      Primary gamma particle energy.
   * @param[in]  lpmenergy   LPM energy (see more at #ComputeLPMDXsectionPerAtom()).
   * @param[in]  deltafactor Target atom dependent constant \f$ 136 Z^{-1/3}\f$.
   * @param[in]  td          Pointer to the GeantV thread local data object (used to get random numbers).
   * @return                 The sampled reduced total energy transfered to one of the e-/e+ pair.
   */
  double SampleTotalEnergyTransfer(double epsmin, double eps0, double fz, double z23, double egamma, double lpmenergy,
                                   double deltafactor, Geant::GeantTaskData *td);
  /**
   * @brief Internal helper method to integrate the DCS in order to get the atomic cross scection.
   *
   * @param[in]  elem      Target element on which atomic cross section is required.
   * @param[in]  mat       Material in which the atomic cross section is required for the given target element.
   * @param[in]  egamma    Primary gamma particle energy.
   * @return               The computed atomic cross section in internal [lenght^2] units.
   */
  double ComputeAtomicCrossSection(const Element *elem, const Material *mat, double egamma);

  // the 2 DCS (with and without LPM) plus a third method that computes these two transformed DCS for the sampling table
  double ComputeDXsectionPerAtom(double epsmin, double eps0, double df, double fz, double xi, bool istsai=false);
  double ComputeLPMDXsectionPerAtom(double epsmin, double df, double fz, double lpmenergy, double z23, double egamma, double xi, bool istsai=false);
  double ComputeDXsection(const Material *mat, double egamma, double epsmin, double xi, bool istsai); // for the sampling table

  //
  void ComputeScreeningFunctions(double &phi1, double &phi2, double delta, bool istsai=false);
  void ComputeLPMFunctions(double &funcPhiS, double &funcGS, double &funcXiS, double z23, double lpmenergy, double eps, double egamma);

  // these 2 are used n the rejection
  double ScreenFunction1(double delta, bool istsai);
  double ScreenFunction2(double delta, bool istsai);



// builds sampling tables for a given material over the discrete photon energy grid
void BuildSamplingTablesForMaterial(const Material *mat);
// builds one sampling table for a given material ata a given photon energy 
void BuildOneRatinAlias(double egamma, const Material *mat, double *pdfarray, int egammaindx, int ilowestz);

// data members
private:
  /** @brief Size of some containers that store data per elements (\f$ Z_{\text{max}} = gMaxZet-1)\f$. */
  static const int gMaxZet = 120; // max Z+1

  /** @brief Elastic form factors for elements Z<8 from \cite tsai1974pair Table B2. */
  static const double gFelLowZet[8];
  /** @brief Inelastic form factors for elements Z<8 from \cite tsai1974pair Table B2. */
  static const double gFinelLowZet[8];

  bool fIsUseTsaisScreening; // setter/getter!

  bool fIsUseLPM; // setter/getter !

  /** @brief  Internal code of the secondary e-. */
  int fElectronInternalCode;
  /** @brief  Internal code of the secondary e+. */
  int fPositronInternalCode;

  double fLPMEnergyLimit;

  /** Data collection that stores some frequently used target atom specific constants for one atom. */
  struct ElementData {
    /** @brief \f$ 136*Z^{-1/3} \f$ */
    double  fDeltaFactor;
    /** @brief Coulomb correction \f$ f_c \f$ as in \cite davies1954theory [Eqs(36-38)] */
    double  fCoulombCor;
    /** @brief \f$ 8\ln(Z)/3 + f_c \f$ */
    double  fFz;
    /** @brief \f$ \exp[(42.24-(8\ln(Z)/3+8f_c))/8.368]-0.952 \f$ */
    double  fDeltaMax;
    /** @brief \f$   1.36\sqrt{\exp(0.5*20.863-2.-0.25*(8\ln(Z)/3+8f_c))-1.}/0.55846 \f$ */
    double  fDeltaMaxTsai;
    /** @brief \f$ Z(Z+\eta(Z)) \f$ with \f$ \eta(Z) = L_{inel}/[L_{el}-f_c]\f$*/
    double  fEtaValue;
  };

  /** @brief Container to store target atom specific data collections (ElementData ) for all target atoms which
    *        the model needs to respond.
    *
    * The size of the container will be equal to #gMaxZet. After initialisation (i.e. calling the InitialiseElementData()
    * method) an element with index \f$Z\f$ will contain data collection for target atom with atomic number \f$ Z \f$.
    * Only those elements will be non-nullptr that corresponds to an atom with atomic number \f$ Z \f$ that the model
    * needs to provide response(final state): i.e. that appears in materials that belongs to regions inwhich the model
    * is active.
    */
  ElementData  **fElementData;


  /**
   * @name Members to describe the common discrete photon energy grid for sampling tables:
   *
   * These variables describe and define the primary gamma kinetic energy grid above we build sampling tables in the
   * InitSamplingTables() method for run-time samling of the total energy transferd to one particle. The min of the table
   * is set to #fMinPrimEnergy and the maximum is #fMaxPrimEnergy (min/max energy usage of model). The number of discrete
   * gamma energy points will be determined by the value of #fNumSamplingPrimEnergiesPerDecade variable. The default value
   * is 10 and it can be changed by the user through the SetNumberOfPhotonEnergiesPerDecade() public method (must be set
   * before the the initialisation of the model!).
   */
  //@{
    /** @brief Number of primary gamma kinetic energy grid points in [#fMinPrimEnergy,#fMaxPrimEnergy].*/
    int fNumSamplingPrimEnergies;
    /** @brief Number of primary gamma kinetic energy grid points per decade. */
    int fNumSamplingPrimEnergiesPerDecade;
    /** @brief Number of "total energy transfered to one of the secondaries" related discrete transformed variable PDF
      *        points in [0,1].
      */
    int fNumSamplingEnergies;
    /** @brief Minimum of the primary gamma kinetic energy grid. */
    double  fMinPrimEnergy;
    /** @brief Maximum of the primary gamma kinetic energy grid. */
    double  fMaxPrimEnergy;
    /** @brief Logarithm of #fMinPrimEnergy ie ln(#fMinPrimEnergy) . */
    double  fPrimEnLMin;
    /** @brief Inverse of the primary gamma kinetic energy grid delta ie
      *        ln[#fMaxPrimEnergy/#fMinPrimEnergy]/(#fNumSamplingPrimEnergies-1)
      */
    double  fPrimEnILDelta;
    /** @brief The logarithmically spaced primary gamma kinetic energy grid.
      *
      * Size of the array is #fNumSamplingPrimEnergies points in the [#fMinPrimEnergy, #fMaxPrimEnergy] interval.
      */
    double *fSamplingPrimEnergies;
    /** @brief The logarithm of #fSamplingPrimEnergies grid.
      *
      *  Size of the array is #fNumSamplingPrimEnergies points in the [ln(#fMinPrimEnergy), ln(#fMaxPrimEnergy)] interval.
      */
    double *fLSamplingPrimEnergies;
  //@}

  /** @brief Internal data structure to store data that are required to sample the energy (transfered to one of the e-/e+
    *        pair) related transformed variable by using the combination of Walker's alias method and rational
    *        interpolation based numerical inversion of the cumulative function (by using AliasTable::SampleRatin() method).
    */
  struct RatinAliasData {
    /** @brief Number of dicrete data points ie size of the arrays = #fNumSamplingEnergies. */
    int     fNumdata;
    /** @brief Total energy (transfered to one of the secondaries) related transformed variable values. */
    double *fXdata;
    /** @brief The cumulative distribution function values over the energy transfer related transformed variable values.*/
    double *fCumulative;
    /** @brief Interpolation parameters over the energy transfer related transformed variable values.*/
    double *fParaA;
    /** @brief Interpolation parameters over the energy transfer related transformed variable values.*/
    double *fParaB;
    /** @brief The alias probabilities over the energy transfer related transformed variable values.*/
    double *fAliasW;
    /** @brief The alias indices over the energy transfer related transformed variable values. */
    int    *fAliasIndx;
  };

  /** @brief Internal data structure to store RatinAliasData for a given target material. */
  struct RatinAliasDataPerMaterial {
    int  fILowestZ; // lowest Z in the corresponding material (used for the variable transform)
    /** @brief Container that stores #fNumSamplingEnergies pointers to RatinAliasData structure over the primary gamma
      *         energy gird #fSamplingPrimEnergies for a specific material built at initialisation if sampling tables
      *         were requested.
      *  The indices of the RatinAliasData pointers correspond to the primary gamma energies in #fSamplingPrimEnergies.
      */
    RatinAliasData **fRatinAliasDataForOneMaterial;
  };

  /** @brief Container to store pointers to RatinAliasDataPerMaterial data structures for all target materials that the
    *        model needs to respond.
    *  Size of the container is equal to number of materials in the detector and indexed by the material index. The
    *  corresponding RatinAliasDataPerElement data structures are built at initialization if sampling tables were required.
    *  Non-nullptr only at indices that corresponds to materials (with the same index) that belongs to regions in
    *  which the model is active.
    */
  std::vector<RatinAliasDataPerMaterial*> fRatinAliasDataForAllMaterials;

  /** @brief Pointer to an AliasTable uitility object that is used both to prepare(at initialization) and to provide
    *        samples(at run-time) from the sampling tables. (Used only if sampling tables were required).
    */
  AliasTable *fAliasSampler;


};


}         // namespace geantphysics

#endif    // RELATIVISTICPAIRMODEL_H
