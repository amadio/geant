
#ifndef GSMSCMODEL_H
#define GSMSCMODEL_H

#include <string>

#include "MSCModel.h"

// from geantV
#include "Geant/Config.h"
#include "GeantTrack.h"

namespace Geant {
  inline namespace GEANT_IMPL_NAMESPACE {
  class GeantTaskData;
}
}


namespace geantphysics {

/**
 * @brief   Model for multiple scattering of e-/e+ based on pre-computed angular distributions computed by means of
 *          Goudsmit-Saunderson theory based on screened Rutherford DCS.
 * @class   GSMSCModel
 * @author  M Novak
 * @date    November 2017
 */


class GSMSCTable;
class GSPWACorrections;
class Particle;

class GSMSCModel : public MSCModel {
public:
  GSMSCModel(bool iselectron = true, const std::string &name = "Goudsmit-Saunderson");
 ~GSMSCModel();

 // implemented base class method
  virtual void  Initialize();
 // implemented MSC base class model methods
  virtual void  StepLimit(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td);
  virtual bool  SampleScattering(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td);
  virtual void  ConvertTrueToGeometricLength(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td);
  virtual void  ConvertGeometricToTrueLength(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td);

 // model specifc method
 void SetOptionPWACorrection(bool opt)  { fIsUsePWACorrection = opt; }
 bool GetOptionPWACorrection() const    { return fIsUsePWACorrection; }

 void SetOptionMottCorrection(bool opt) { fIsUseMottCorrection = opt; }
 bool GetOptionMottCorrection() const   { return fIsUseMottCorrection; }


// make it public for testing
void   ComputeParameters(const MaterialCuts *matcut, double ekin, double &lambel, double &lambtr1,
                         double &scra, double &g1, double &mccor1, double &mccor2);

private:
  double RandomizeTrueStepLength(Geant::GeantTaskData *td, double tlimit);
  void   SampleMSC(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td);
  double GetTransportMeanFreePathOnly(const MaterialCuts *matcut, double ekin);

//  void   ComputeParameters(const MaterialCuts *matcut, double ekin, double &lambel, double &lambtr1,
//                           double &scra, double &g1);
// data members
private:
  bool   fIsElectron          = false;      // is the model for e- (e+ otherwise)
  bool   fIsUsePWACorrection  = true;       // use screening that gives back pwa first transport mean free path
  bool   fIsUseMottCorrection = false;      // use Mott correction
  bool   fIsUseAccurate       = true;       // use accurate step limits
  bool   fIsOptimizationOn    = true;       // use optimisation in the step limit: check current range and pre-safety

  double fCharge             = 0.0;

  double fTauSmall           = 1.e-16;
  double fTauLim             = 1.e-6;
  double fTLimitMinfix2      = 1.*geant::nm;
  double fDtrl               = 0.05;

  Particle* fParticle        = nullptr;    //e-/e+
  Geant::TrackToken fMSCdata;   // Handle for MSCData

  GSMSCTable*       fGSTable        = nullptr;
  GSPWACorrections* fPWACorrection  = nullptr;

};

}      // namespace geantphysics

#endif // GSMSCMODEL_H
