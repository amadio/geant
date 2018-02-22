
#ifndef POSITRONTO2GAMMAMODEL_H
#define POSITRONTO2GAMMAMODEL_H

#include "EMModel.h"

// from geantV
#include "Geant/Config.h"
namespace geant {
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

/**
 * @brief   Positron annihilation into 2 gamma model.
 * @class   PositronTo2GammaModel
 * @author  M Novak
 * @date    january 2018
 *
 * The model is based on Heitler's differential cross section of two-photon positron-electron annihilation
 * \cite heitler1954quantum. One-photon annihilation is not considered since it's assumed that the target e- is free and
 * at rest i.e. there are no binding effects. Similarly, positron annihilation into three or more photons are also
 * ignored.
 *
 */


class PositronTo2GammaModel : public EMModel {
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
    PositronTo2GammaModel(const std::string &modelname="e2GammaAnnih");
     /** @brief Destructor. */
    virtual ~PositronTo2GammaModel();
//@}

/**
* @name Implemented EMModel base class methods:
*/
//@{
    /** @brief Interface method to initilize the model. */
    virtual void   Initialize();
    virtual double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *particle);
    virtual double ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *matcut, double kinenergy, const Particle *particle);
    virtual int    SampleSecondaries(LightTrack &track, geant::GeantTaskData *td);
//
//@}

private:

  /** @brief Copy constructor  (deleted) */
  PositronTo2GammaModel(const PositronTo2GammaModel&) = delete;
  /** @brief Operator=  (deleted) */
  PositronTo2GammaModel &operator=(const PositronTo2GammaModel&) = delete;

  double ComputeXsectionPerElectron(double pekin);
  double SampleEnergyTransfer(double pekin, double gamma, double r1, double r2, double r3);
  double SampleEnergyTransfer(double gamma, geant::GeantTaskData *td);


  /** @brief Internal method to build energy transfer (to one of the gammas) related sampling tables.*/
  void   InitSamplingTables();


  /** @brief Internal method to build one sampling tables.
   *
   *  @param[in]  indx   Index of the alias table data structure to build in the container.
   *  @param[in]  gamma  Initial e+ total energy in rest mass units \f$ \gamma=(E_k+m_ec^2)/(m_ec^2)\f$.
   */
  void   BuildOneLinAlias(int indx, double gamma);

  void   ClearSamplingTables();

  double ComputeTransfDXSec(double xi, double gamma, double mineps, double maxeps);
//@}

  /** @brief Internal data structure to store data for sampling tables. */
  struct LinAlias{
    LinAlias(int num) { fXdata.resize(num); fYdata.resize(num); fAliasW.resize(num); fAliasIndx.resize(num); }
    /** @brief Transformed variable values. */
    std::vector<double> fXdata;
    /** @brief The pdf values (not necessarily normalised) over the energy transfer related variable values. */
    std::vector<double> fYdata;
    /** @brief The alias probabilities (not necessarily normalised) over the energy transfer related variables. */
    std::vector<double> fAliasW;
    /** @brief The alias indices over the energy transfer related transformed variable values. */
    std::vector<int>    fAliasIndx;
  };

private:

  int    fSecondaryInternalCode;

  int    fSTNumPositronEnergiesPerDecade;  // ST=>SamplingTables
  int    fSTNumDiscreteEnergyTransferVals; // ST=>SamplingTables
  int    fSTNumPositronEnergies;           // ST=>SamplingTables

  double fSTLogMinPositronEnergy;          // ST=>SamplingTables
  double fSTILDeltaPositronEnergy;         // ST=>SamplingTables

  /** @brief Container to store pointers to LinAlias data structures.*/
  std::vector<LinAlias*>   fSamplingTables;
  /** @brief An alias sampler used at run-time sampling of the energy transfered to one of the gammas
    *        variable from a LinAlias data structure (prepared at initialisation).
    */
  AliasTable *fAliasSampler;


};

}      // namespace geantphysics

#endif // POSITRONTO2GAMMAMODEL_H
