
#ifndef BETHEHEITLERPAIRMODEL_H
#define BETHEHEITLERPAIRMODEL_H

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
 * @brief   Model for conversion of gamma particles into e-/e+ pair in the the field of nucleus and atomic electrons.
 * @class   BetheHeitlerPairModel
 * @author  F Hariri, M Novak
 * @date    May 2017
 *
 * The model is based on atomic cross sections computed by using the Geant4 parametrization
 * \cite cirrone2010validation \cite agostinelli2003geant4 of the numerical cross section values in
 * \cite hubbell1980pair. These numerical cross sections include pair production both in the nuclear (with screening,
 * Coulomb and radiative corrections) and electron fields (with approximate screening, radiative, exchange and
 * retardation effects). The parametrization gives accurate results up to 70-90 [GeV]. See more at the
 * #ComputeAtomicCrossSection() method.
 * The final state energies (total energy transfered to one of the e-/e+ pair) are computed based on the Bethe-Heitler
 * \cite bethe1934stopping differential cross section (DCS) corrected for various effects like screening, Coulomb
 * correction, conversion in the field of atomic electrons (see more at the #ComputeDXSection() method). However,
 * triplet production is not generated and the asymmetry between the e-/e+ at low photon energies are not accounted.
 */

//class Material;
class MaterialCuts;
//class Element;
class AliasTable;
class Particle;
class LightTrack;


class BetheHeitlerPairModel : public EMModel {
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
  BetheHeitlerPairModel(const std::string &modelname="PairBetheHeitler");
  /** @brief Destructor. */
  virtual ~BetheHeitlerPairModel();
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

  /**
   * @brief Method to obtain the minimum primary gamma energy at which the interaction can happen.
   *
   * It's used by the target element selector to build partial (per element) lambda tables that are used at run-time to
   * select target atom for the interaction:
   * \f$ E_{\gamma}^{\text{min}} = \text{max}\{2m_ec^2,\text{min-energy-usage} \} \f$
   *
   * @return \f$ E_{\gamma}^{\text{min}} \f$
   */
  virtual double MinimumPrimaryEnergy(const MaterialCuts * /*matcut*/, const Particle * /*part*/) const {
    return fMinimumPrimaryEnergy;
  }
//
//@}



private:
  /** @brief Copy constructor  (deleted) */
  BetheHeitlerPairModel(const BetheHeitlerPairModel&) = delete;
  /** @brief Operator=  (deleted) */
  BetheHeitlerPairModel &operator=(const BetheHeitlerPairModel&) = delete;

/**
* @name Model specific private methods.
*/
//@{
  /**
   * @brief Internal method to compute parametrized atomic cross section of gammas conversion into e-/e+ pair
   *        for a given target atom, gamma particle kinetic energy.
   *
   * @param[in] z        Atomic number of the target atom.
   * @param[in] egamma   Kinetic energy of the gamma particle in internal [energy] units.
   * @return    The computed atomic cross section in internal [lenght^2] units.
   */
  double ComputeAtomicCrossSection(double z, double egamma);

  /** @brief Internal method to build data collection of some frequently used target atom dependent variables. */
  void InitialiseElementData();

  /**
   * @brief Internal method to build sampling table data structures for fast run-time sampling of the reduced total
   *        energy transfered to one of the e-/e+ pair for all possible target atom.
   *
   * Sampling tables are built at initialization for all target elements that belongs to materials that appear in
   * regions in which the model is active if sampling tables were requested. First the common, discrete, primary photon
   * energy grid is generated between. Then sampling tables are built for each element at all dicrete primary photon
   * energies that are stored in separate #RatinAliasDataPerElement data structures for the different target elements.
   * Such data structures for all possible target elements will be stored in the #fSamplingTables data structure.
   */
   void InitSamplingTables();

   /**
    * @brief Internal method to sample reduced total energy transfered to one of the e-/e+ pair using sampling tables
    *        prepared at initialization}.
    *
    * @param[in] egmma    Gamma photon energy in internal [energy] units.
    * @param[in] izet     Index of the target element sampling table collection in #fRatinAliasDataForAllElements (Z).
    * @param[in] r1       Random number distributed uniformly in [0,1].
    * @param[in] r2       Random number distributed uniformly in [0,1].
    * @param[in] r3       Random number distributed uniformly in [0,1].
    * @return             The sampled reduced total energy transfered to one of the e-/e+ pair.
    */
   double SampleTotalEnergyTransfer(const double egamma, const int izet, const double r1, const double r2, const double r3);

   /**
    * @brief Internal method to sample reduced total energy transfered to one of the e-/e+ pair using rejection.
    *
    * @param[in] egmma    Gamma photon energy in internal [energy] units.
    * @param[in] izet     Index of the target element sampling table collection in #fRatinAliasDataForAllElements (Z).
    * @param[in]  td      Pointer to the GeantV thread local data object (used to get random numbers).
    * @return             The sampled reduced total energy transfered to one of the e-/e+ pair.
    */
   double SampleTotalEnergyTransfer(const double epsmin, const int izet, const Geant::GeantTaskData *td);


   /**
    * @brief Internal method to compute the screening functions.
    *
    * @param[in,out]  phi1 Screening function value \f$ \phi(\delta)_1 \f$ .
    * @param[in,out]  phi2 Screening function value \f$ \phi(\delta)_2 \f$ .
    * @param[in]      delta  Screening variable \f$ \delta(\epsilon) = 136 Z^{-1/3} \epsilon_0/(\epsilon(1-\epsilon)) \f$.
    * @param[in]      istsai Flag to indicate if Tsai's screening function approximation needs to be used.
    */
   void   ComputeScreeningFunctions(double &phi1, double &phi2, const double delta, const bool istsai);

   /**
    * @brief Same as ScreenFunction1() and ScreenFunction2() at once (used in the rejection algorithm).
    */
   void   ScreenFunction12(double &val1, double &val2, const double delta, const bool istsai);

   /**
    * @brief Internal method used in the rejection algorithm to compute the screening funtion related part.
    *
    * @param[in] delta  Screening variable \f$ \delta(\epsilon) = 136 Z^{-1/3} \epsilon_0/(\epsilon(1-\epsilon)) \f$.
    * @param[in] istsai Flag to indicate if Tsai's screening function approximation needs to be used.
    * @return    \f$ 3 \Phi_1(\delta) - \Phi_2(\delta) \f$ (see more at #ComputeDXSection()).
    */
   double ScreenFunction1(const double delta, const bool istsai);
   /**
    * @brief Internal method used in the rejection algorithm to compute the screening funtion related part.
    *
    * @param[in] delta   Screening variable \f$ \delta(\epsilon) = 136 Z^{-1/3} \epsilon_0/(\epsilon(1-\epsilon)) \f$.
    * @param[in] istsai  Flag to indicate if Tsai's screening function approximation needs to be used.
    * @return    \f$ 1.5 \Phi_1(\delta) +0.5 \Phi_2(\delta) \f$ (see more at #ComputeDXSection()).
    */
   double ScreenFunction2(const double delta, const bool istsai);


   /**
     * @brief Internal method to build sampling tables for a given target element at all discrete photon energies.
     *
     * @param elem Pointer to the given target element object.
     */
   void BuildSamplingTablesForElement(const Element *elem, const std::vector<double> &primevect);

   /**
     * @brief Internal method to build one sampling table for a given photon energy and target element.
     *
     * @param egamma     Gamma photon energy in internal [energy] units.
     * @param elem       Pointer to the target element object.
     * @param thepdfdata Array to store some temporary pdf values (size is fNumSamplingEnergies).
     * @param egammaindx Index of the gamma photon energy in the gamma energy grid.
     */
   void BuildOneRatinAlias(const double egamma, const int izet, double *thepdfdata, const int egammaindx);

   /**
     * @brief Internal method to compute the total energy, transferd to one of the e-/e+ pair, related transformed
     *        variable dependent part of the DCS.
     *
     *  The method is used to prepare sampling tables at initialization for the fast sampling of the reduced total
     *  energy, transfered to one of the e-/e+ pair, at run-time.
     *
     *  @param[in] epsmin      Minimum value of the reduced total energy transfer \f$ \epsilon_{\text{min}} =
     *                         \text{max}[\epsilon_0,\epsilon'] \f$ .
     *  @param[in] eps0        Kinematical minimum of the reduced total energy transfer
     *                         \f$ \epsilon_0 \equiv m_ec^2/E_{\gamma} \f$ .
     *  @param[in] deltafactor Target atomic number dependent constant \f$ 136Z^{-1/3} \f$
     *  @param[in] fz          The target atomic number dependent Coulomb correction function
     *                          \f[
     *                             F(Z) =
     *                               \begin{cases}
     *                                 \frac{8}{3} \ln(Z) & \quad \text{if } E_{\gamma} < 50 \text{[MeV]} \\
     *                                 \frac{8}{3} \ln(Z) + 8 f_c(Z) & \quad \text{if } E_{\gamma} \geq 50 \text{[MeV]}
     *                               \end{cases}
     *                           \f]
     *  @param[in] xi           Total energy transfer related transformed variable \f$ \xi \in [0,1] \f$ when
     *                          \f$ \epsilon \in [\epsilon_{\text{min}},0.5] \f$  and  \f$ \xi(\epsilon) =
     *                               \ln[\epsilon/\epsilon_{\text{min}}]/\ln[0.5/\epsilon_{\text{min}}]\f$
     * @param[in] istsai        Flag to indicate if Tsai's screening function approximation needs to be used.
     *  @return                 The \f$\xi\f$ depends part of the transformed DCS \f$ \left( \mathrm{d}\sigma(Z,
     *                          \epsilon(\xi))/\mathrm{d}\xi \right)^* \f$ at the given input parameters.
     */
   double ComputeDXSection(double epsmin, double eps0, double deltafactor, double fz, double xi);

   void   ClearSamplingTables();
//@}


  /** Data collection that stores some frequently used target atom specific constants for one atom. */
  struct ElementData {
    /** @brief \f$ 136*Z^{-1/3} \f$ */
    double  fDeltaFactor;
    /** @brief Coulomb correction \f$ f_c \f$ as in \cite davies1954theory [Eqs(36-38)] */
    double  fCoulombCor;
    /** @brief \f$ 8\ln(Z)/3 \f$ */
    double  fFzLow;
    /** @brief \f$ 8\ln(Z)/3 + f_c \f$ */
    double  fFzHigh;
    /** @brief \f$ \exp[(42.24-8\ln(Z)/3)/8.368]-0.952 \f$ */
    double  fDeltaMaxLow;
    /** @brief \f$ \exp[(42.24-(8\ln(Z)/3+f_c))/8.368]-0.952 \f$ */
    double  fDeltaMaxHigh;
    /** @brief \f$ 1.36\sqrt{\exp(0.5*20.863-2.-0.25*(8\ln(Z)/3))-1.}/0.55846 \f$ */
    double  fDeltaMaxLowTsai;
    /** @brief \f$ 1.36\sqrt{\exp(0.5*20.863-2.-0.25*(8\ln(Z)/3+8f_c))-1.}/0.55846 \f$ */
    double  fDeltaMaxHighTsai;
  };

  /** @brief Internal data structure to store data that are required to sample the energy (transfered to one of the e-/e+
    *        pair) related transformed variable by using the combination of Walker's alias method and rational
    *        interpolation based numerical inversion of the cumulative function (by using AliasTable::SampleRatin() method).
    */
  struct RatinAliasData {
    RatinAliasData(size_t n) {
      fXdata.resize(n); fCumulative.resize(n); fParaA.resize(n); fParaB.resize(n);
      fAliasW.resize(n); fAliasIndx.resize(n);
    }
    /** @brief Total energy (transfered to one of the secondaries) related transformed variable values. */
    std::vector<double> fXdata;
    /** @brief The cumulative distribution function values over the energy transfer related transformed variable values.*/
    std::vector<double> fCumulative;
    /** @brief Interpolation parameters over the energy transfer related transformed variable values.*/
    std::vector<double> fParaA;
    /** @brief Interpolation parameters over the energy transfer related transformed variable values.*/
    std::vector<double> fParaB;
    /** @brief The alias probabilities over the energy transfer related transformed variable values.*/
    std::vector<double> fAliasW;
    /** @brief The alias indices over the energy transfer related transformed variable values. */
    std::vector<int>    fAliasIndx;
  };

  /** @brief Internal data structure to store RatinAliasData for a given target element. */
  struct RatinAliasDataPerElement {
    RatinAliasDataPerElement(int numeprims) {
      fRatinAliasData.resize(numeprims,nullptr);
    }
    /** @brief Container that stores #fSTNumPhotonEnergies pointers to RatinAliasData structure over the primary gamma
      *         energy gird for a specific target atom built at initialisation if sampling tables were requested.
      */
    std::vector<RatinAliasData*> fRatinAliasData;
  };


// data members
private:
  /** @brief Size of some containers that store data per elements (\f$ Z_{\text{max}} = gMaxZet-1)\f$. */
  static const long  gMaxZet = 121; // max Z+1
  static std::vector<ElementData*>  gElementData;

  bool   fIsUseTsaisScreening;

  /** @brief  Internal code of the secondary e-. */
  int    fElectronInternalCode;
  /** @brief  Internal code of the secondary e+. */
  int    fPositronInternalCode;

  int    fSTNumPhotonEnergiesPerDecade;    // ST=>SamplingTables
  int    fSTNumDiscreteEnergyTransferVals; // ST=>SamplingTables
  int    fSTNumPhotonEnergies;             // ST=>SamplingTables

  double fSTLogMinPhotonEnergy;            // ST=>SamplingTables
  double fSTILDeltaPhotonEnergy;           // ST=>SamplingTables

  /** @brief Minimum photon energy at which the discrete interaction can happen:
             \f$ E_{\gamma}^{\text{min}} = \text{max}\{2m_ec^2,\text{min-energy-usage} \} \f$
    */
  double fMinimumPrimaryEnergy;

  /** simplified sampling below this energy limit: 2 MeV */
  double fGammaEneregyLimit;

  /** @brief Container to store pointers to RatinAliasDataPerElement data structures for all target elements that the
    *        model needs to respond.
    *  Size of the container is equal to #gMaxZet and indexed by the atomic number of the target elements. The
    *  corresponding RatinAliasDataPerElement data structures are built at initialization if sampling tables were required.
    *  Non-nullptr only at indices that corresponds to target atomic number that appears in materials that belongs to
    *  regions in which the model is active.
    */
  std::vector<RatinAliasDataPerElement*> fSamplingTables;

  /** @brief Pointer to an AliasTable uitility object that is used both to prepare(at initialization) and to provide
    *        samples(at run-time) from the sampling tables. (Used only if sampling tables were required).
    */
  AliasTable *fAliasSampler;


};

}      // namespace geantphysics

#endif // BETHEHEITLERPAIRMODEL_H
