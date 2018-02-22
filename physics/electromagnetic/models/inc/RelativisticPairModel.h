
#ifndef RELATIVISTICPAIRMODEL_H
#define RELATIVISTICPAIRMODEL_H

#include "EMModel.h"

// from geantV
#include "Geant/Config.h"
namespace geant {
  inline namespace GEANT_IMPL_NAMESPACE {
  class TaskData;
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
class GLIntegral;

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
  virtual ~RelativisticPairModel();
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
  virtual int    SampleSecondaries(LightTrack &track, geant::TaskData *td);
//
//@}



private:
  /** @brief Copy constructor  (deleted) */
  RelativisticPairModel(const RelativisticPairModel&) = delete;
  /** @brief Operator=  (deleted) */
  RelativisticPairModel &operator=(const RelativisticPairModel&) = delete;

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
   * @param[in] egamma   Gamma photon energy in internal [energy] units.
   * @param[in] matindx  Index of the target material.
   * @param[in] r1       Random number distributed uniformly in [0,1].
   * @param[in] r2       Random number distributed uniformly in [0,1].
   * @param[in] r3       Random number distributed uniformly in [0,1].
   * @return             The sampled reduced total energy transfered to one of the e-/e+ pair.
   */
  double SampleTotalEnergyTransfer(const double egamma, const int matindx, const double r1, const double r2, const double r3);

  /**
   * @brief Internal method to sample reduced total energy transfered to one of the e-/e+ pair using rejection.
   *
   * @param[in]  egamma      Primary gamma particle energy.
   * @param[in]  lpmenergy   LPM energy (see more at #ComputeLPMDXsectionPerAtom()).
   * @param[in]  td          Pointer to the GeantV thread local data object (used to get random numbers).
   * @return                 The sampled reduced total energy transfered to one of the e-/e+ pair.
   */
  double SampleTotalEnergyTransfer(const double egamma, const double lpmenergy, const int izet,
                                   const geant::TaskData *td);
  /**
   * @brief Internal helper method to integrate the DCS in order to get the atomic cross scection.
   *
   * @param[in]  elem      Target element on which atomic cross section is required.
   * @param[in]  mat       Material in which the atomic cross section is required for the given target element.
   * @param[in]  egamma    Primary gamma particle energy.
   * @return               The computed atomic cross section in internal [lenght^2] units.
   */
  double ComputeAtomicCrossSection(const Element *elem, const Material *mat, const double egamma);


  // the 2 DCS (with and without LPM) plus a third method that computes these two transformed DCS for the sampling table
  double ComputeDXsectionPerAtom(const double epsmin, const double egamma, const double xi, const int izet,
                                 bool istsai=false);
  double ComputeLPMDXsectionPerAtom(const double epsmin, const double egamma, const double xi, const double lpmenergy,
                                    const int izet,  bool istsai=false);
  // for the sampling table
  double ComputeDXsection(const Material *mat, double egamma, double epsmin, double xi, bool istsai);


  //
  void   ComputeLPMfunctions(double &funcXiS, double &funcGS, double &funcPhiS, const double lpmenergy, const double eps,
                             const double egamma, const int izet);
  void   ComputeLPMGsPhis(double &funcGS, double &funcPhiS, const double varShat);
  void   InitLPMFunctions();
  void   GetLPMFunctions(double &lpmGs, double &lpmPhis, const double s);


  void   ComputeScreeningFunctions(double &phi1, double &phi2, const double delta, const bool istsai);
  // these 3 are used only in the rejection
  void   ScreenFunction12(double &val1, double &val2, const double delta, const bool istsai);
  double ScreenFunction1(const double delta, const bool istsai);
  double ScreenFunction2(const double delta, const bool istsai);


  void   ClearSamplingTables();

  // builds sampling tables for a given material over the discrete photon energy grid
  void   BuildSamplingTablesForMaterial(const Material *mat, const std::vector<double> &primevect);
  // builds one sampling table for a given material ata a given photon energy
  void   BuildOneRatinAlias(const double egamma, const Material *mat, double *pdfarray, const int egammaindx, const int ilowestz);


  /** Data collection that stores some frequently used target atom specific constants for one atom. */
  struct ElementData {
    /** @brief Coulomb correction \f$ f_c \f$ as in \cite davies1954theory [Eqs(36-38)] */
    double  fCoulombCor;
    /** @brief \f$ 136*Z^{-1/3} \f$ */
    double  fDeltaFactor;
    /** @brief \f$ 8\ln(Z)/3 + f_c \f$ */
    double  fFz;
    /** @brief \f$ \exp[(42.24-(8\ln(Z)/3+8f_c))/8.368]-0.952 \f$ */
    double  fDeltaMax;
    /** @brief \f$   1.36\sqrt{\exp(0.5*20.863-2.-0.25*(8\ln(Z)/3+8f_c))-1.}/0.55846 \f$ */
    double  fDeltaMaxTsai;
    /** @brief \f$ Z(Z+\eta(Z)) \f$ with \f$ \eta(Z) = L_{inel}/[L_{el}-f_c]\f$*/
    double  fEtaValue;
    /** @brief \f$ s_1 = \sqrt{2} Z^{2/3}/(184.15^2)*/
    double  fVarS1Cond;
    /** @brief \f$  1/\ln(s_1) */
    double  fILVarS1Cond;
  };

  struct LPMFuncs {
    LPMFuncs() : fIsInitialized(false), fSDelta(0.01), fSLimit(2.) {}
    bool   fIsInitialized;
    double fSDelta;
    double fSLimit;
    std::vector<double>  fLPMFuncG;
    std::vector<double>  fLPMFuncPhi;
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

  /** @brief Internal data structure to store RatinAliasData for a given target material. */
  struct RatinAliasDataPerMaterial {
    RatinAliasDataPerMaterial(int numeprims) {
      fILowestZ = 200;
      fRatinAliasData.resize(numeprims,nullptr);
    }
    int  fILowestZ; // lowest Z in the corresponding material (used for the variable transform)
    /** @brief Container that stores pointers to RatinAliasData structure over the primary gamma
      *         energy gird for a specific material built at initialisation if sampling tables were requested.
      */
    std::vector<RatinAliasData*> fRatinAliasData;
  };


// data members
private:
  /** @brief Size of some containers that store data per elements (\f$ Z_{\text{max}} = gMaxZet-1)\f$. */
  static const long                 gMaxZet = 121; // max Z+1
  /** @brief Elastic form factors for elements Z<8 from \cite tsai1974pair Table B2. */
  static const double               gFelLowZet[8];
  /** @brief Inelastic form factors for elements Z<8 from \cite tsai1974pair Table B2. */
  static const double               gFinelLowZet[8];
  static const double               gLPMFactor;
  static LPMFuncs                   gLPMFuncs;
  static std::vector<ElementData*>  gElementData;


  bool   fIsUseTsaisScreening;
  bool   fIsUseLPM;

  /** @brief  GL integral number of abscissas and weights. */
  int    fNGL;
  /** @brief  Internal code of the secondary e-. */
  int    fElectronInternalCode;
  /** @brief  Internal code of the secondary e+. */
  int    fPositronInternalCode;

  int    fSTNumPhotonEnergiesPerDecade;    // ST=>SamplingTables
  int    fSTNumDiscreteEnergyTransferVals; // ST=>SamplingTables
  int    fSTNumPhotonEnergies;             // ST=>SamplingTables

  double fLPMEnergyLimit;  // setter/getter

  double fSTLogMinPhotonEnergy;            // ST=>SamplingTables
  double fSTILDeltaPhotonEnergy;           // ST=>SamplingTables

  /** @brief Container to store pointers to RatinAliasDataPerMaterial data structures for all target materials that the
    *        model needs to respond.
    *  Size of the container is equal to number of materials in the detector and indexed by the material index. The
    *  corresponding RatinAliasDataPerElement data structures are built at initialization if sampling tables were required.
    *  Non-nullptr only at indices that corresponds to materials (with the same index) that belongs to regions in
    *  which the model is active.
    */
  std::vector<RatinAliasDataPerMaterial*> fSamplingTables;

  /** @brief Pointer to an AliasTable uitility object that is used both to prepare(at initialization) and to provide
    *        samples(at run-time) from the sampling tables. (Used only if sampling tables were required).
    */
  AliasTable *fAliasSampler;

  /** @brief A GL numerical integral for integrations. */
  GLIntegral *fGL;
};


}         // namespace geantphysics

#endif    // RELATIVISTICPAIRMODEL_H
