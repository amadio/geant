#ifndef NeutronNudyCaptureModel_H
#define NeutronNudyCaptureModel_H

#include "Geant/TNudyEndfRecoPoint.h"
#include "Geant/HadronicFinalStateModel.h"

#include <string>

namespace geantphysics {
/**
 * @brief   Endf based elastic model
 * @class   NeutronNudyCaptureModel
 * @author  H Kumawat
 * @date    June 2018
 *
 * 
 */

namespace geantphysics {
inline namespace GEANT_IMPL_NAMESPACE {
class Isotope;
}
}

class LightTrack;

class NeutronNudyCaptureModel : public HadronicFinalStateModel {
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
  NeutronNudyCaptureModel(const std::string &modelname = "NeutronNudyCapture");
  /** @brief Destructor. */
  ~NeutronNudyCaptureModel();
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
