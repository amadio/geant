#ifndef NeutronNudyInelasticModel_H
#define NeutronNudyInelasticModel_H

#include "Geant/TNudyEndfRecoPoint.h"
#include "Geant/HadronicFinalStateModel.h"

#include <string>

namespace geantphysics {
/**
 * @brief   Endf based elastic model
 * @class   NeutronNudyInelasticModel
 * @author  H Kumawat
 * @date    June 2018
 *
 *
 */

namespace geantphysics {
inline namespace GEANT_IMPL_NAMESPACE {
class Isotope;
}
} // namespace geantphysics

class LightTrack;

class NeutronNudyInelasticModel : public HadronicFinalStateModel {
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
  NeutronNudyInelasticModel(const std::string &modelname = "NeutronNudyInelastic");
  /** @brief Destructor. */
  ~NeutronNudyInelasticModel();
  //@}

  /**
   * @name Implemented HadronicFinalStateModel base class methods:
   */
  //@{
  /** @brief Interface method to initilize the model. */
  virtual void Initialize();

  /** @brief Interface method to generate final state of the interaction. */
  virtual int SampleFinalState(LightTrack &track, Isotope *targetisotope, geant::TaskData *td);

  // sample momentum transfer in the CMS system
  double SampleInvariantT(double mass, double plab, Isotope *targetisotope, geant::TaskData *td);

  //
  //@}

private:
  /**
   * @name Model specific private methods.
   */
  //@{

  //@}

  // data members
private:
  const char *fRENDF;
  std::string filename;
  NudyPhysics::TNudyEndfRecoPoint recopoint;
};

} // namespace geantphysics

#endif // DIFFUSEELASTICMODEL_H
