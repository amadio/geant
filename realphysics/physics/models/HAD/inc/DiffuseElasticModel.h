#ifndef DIFFUSEELASTICMODEL_H
#define DIFFUSEELASTICMODEL_H

#include "HadronicFinalStateModel.h"

#include <string>

namespace geantphysics {
  /**
   * @brief   Optical diffuse elastic scattering with 4-momentum balance. Based on the original implementation by V Grishine.
   * @class   DiffuseElasticModel
   * @author  W Pokorski
   * @date    June 2017
   *
   * detailed description will be added later
   */

  namespace geantphysics {
    inline namespace GEANT_IMPL_NAMESPACE {
      class Isotope;
    }
  }

  class LightTrack;


  class DiffuseElasticModel : public HadronicFinalStateModel {
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
    DiffuseElasticModel(const std::string &modelname="DiffuseElastic");
    /** @brief Destructor. */
    ~DiffuseElasticModel();
    //@}

    /**
     * @name Implemented HadronicFinalStateModel base class methods:
     */
    //@{
    /** @brief Interface method to initilize the model. */
    virtual void   Initialize();
  
    /** @brief Interface method to generate final state of the interaction. */
    virtual int  SampleFinalState(LightTrack &track, Isotope* targetisotope, 
				  Geant::GeantTaskData *td);

    // sample momentum transfer in the CMS system 
    double SampleInvariantT(double mass, double plab, Isotope* targetisotope, Geant::GeantTaskData *td);  

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
  };

}      // namespace geantphysics

#endif // DIFFUSEELASTICMODEL_H
