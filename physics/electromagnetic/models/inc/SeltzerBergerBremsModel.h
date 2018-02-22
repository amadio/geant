
#ifndef SELTZERBERGERBREMSMODEL_H
#define SELTZERBERGERBREMSMODEL_H

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
class Spline;

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
  /**
  * @name Constructor, destructor:
  */
  //@{
     /**
      * @brief Constructor to build a model based on the numerical differential cross sections stored in files.
      *
      * @param[in] iselectron   Flag to indicate that the model is for electron(true) or for psitron(false).
      * @param[in] modelname    Name of the model.
      */
    SeltzerBergerBremsModel(bool iselectron, const std::string &modelname = "eSeltzerBergerBrems");

    /** @brief Destructor. */
    virtual ~SeltzerBergerBremsModel();
  //@}

/**
* @name Implemented EMModel base class methods:
*/
//@{
  virtual void   Initialize(); // from EMModel

  virtual double ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle* particle,bool istotal=false);
  virtual double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *particle);
  virtual double ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut, double kinenergy, const Particle *particle);
  virtual int    SampleSecondaries(LightTrack &track, geant::TaskData *td);
  virtual double MinimumPrimaryEnergy(const MaterialCuts *matcut, const Particle *part) const;
//
//@}


private:
  /** @brief Copy constructor  (deleted) */
  SeltzerBergerBremsModel(const SeltzerBergerBremsModel&) = delete;
  /** @brief Operator=  (deleted) */
  SeltzerBergerBremsModel &operator=(const SeltzerBergerBremsModel&) = delete;

  double   ComputeDEDXPerVolume(const Material *mat, double gammaprodcutenergy, double electronekin);
  double   ComputeXSectionPerVolume(const Material *mat, double gammaprodcutenergy, double electronekin);
  double   ComputeXSectionPerAtom(const Element *elem, const Material *mat, double gammaprodcutenergy,double electronekin);
  double   SamplePhotonEnergy(const MaterialCuts *matcut, double eekin, double r1, double r2, double r3);
  double   SamplePhotonEnergy(double eekin, double gcut, double zet, const Material *mat, geant::TaskData *td);
  void     SamplePhotonDirection(double elenergy, double &sinTheta, double &cosTheta, double rndm);

  void     ClearLoadDCSData();
  void     LoadDCSData();
  void     InitSamplingTables();
  void     ClearSamplingTables();
  void     BuildSamplingTableForMaterialCut(const MaterialCuts *matcut, int indxlocal);
//  double   GetDXSECValue(double eprim, double kappa, int ie, int ik);
  double   GetEkinIndex(double &ekin, int &ie);
  double   GetDXSECValue(int zet, int ie, double eresid, double kappa);
  double   GetDXSECValue(int zet, double eprim, double kappa);

  /**
   * @brief Correction for accounting some differencies between positrons and electrons.
   */
  double PositronCorrection(double ekinelectron, double ibeta2electron, double ephoton, double z);

  /**
   * @brief Correction for accounting some differencies between positrons and electrons.
   */
  double PositronCorrection1(double ekinelectron, double ephoton, double gcutener, double z);


private:

  struct XsecDataZet {
    XsecDataZet(int nprime, int nphote) {
      fSplineVect.resize(nprime,nullptr);
      fXsecDataVect.resize(nprime);
      for (int i=0; i<nprime; ++i)
        fXsecDataVect[i].resize(nphote,0.);
    }
    std::vector< std::vector <double> > fXsecDataVect;  // LoadDCSNumElectronEnergies x fLoadDCSNumReducedPhotonEnergies
    std::vector<Spline*>                fSplineVect;    // LoadDCSNumElectronEnergies
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


private:

  static const double gMigdalConst;

  static std::vector<XsecDataZet*>  fXsecDataPerZet;
  static std::vector<double>        fXsecLimits;
  static std::vector<double>        fLoadDCSElectronEnergyGrid;
  static std::vector<double>        fLoadDCSReducedPhotonEnergyGrid;


  /** @brief Flag to indicate if the model is for e-(true) or for e+(false). Must be set before initialisation. */
  bool     fIsElectron;
  /** @brief Number of points used in the GL integrals. */
  int      fNGL;
  /** @brief Secondary particle i.e. gamma particle internal code. */
  int      fSecondaryInternalCode;
  //
  // these are for the DCS data loaded from the file
  /** @brief Maximum atomic number that numerical DCS are available in the file. */
  int      fDCSMaxZet;
  /** @brief Number of electron energies that numerical DCS are available in the file. */
  int      fLoadDCSNumElectronEnergies;
  /** @brief Number of reduced photon energies at each electron energy that numerical DCS are available in the file. */
  int      fLoadDCSNumReducedPhotonEnergies;

  int      fSTNumElectronEnergyPerDecade;   // ST=> sampling tables
  int      fSTNumSamplingPhotEnergies;      // ST=> sampling tables

  double   fLogLoadDCSMinElecEnergy;
  double   fInvLogLoadDCSDeltaEnergy;

  AliasTable*                          fAliasSampler;
  std::vector<int>                     fGlobalMatGCutIndxToLocal;
  std::vector<AliasDataMaterialCuts*>  fSamplingTables;

  GLIntegral  *fGL;

};

}      // namespace geantphysics


#endif // SELTZERBERGERBREMSMODEL_H
