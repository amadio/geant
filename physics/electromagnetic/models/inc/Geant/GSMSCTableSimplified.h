
#ifndef GSMSCTABLESIMPLIFIED_H
#define GSMSCTABLESIMPLIFIED_H

// from geantV
#include "Geant/Config.h"
namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {
class TaskData;
}
}

#include <vector>

namespace geantphysics {

/**
 * @brief   Helper class for the GSMSCModel to handle the pre-computed angular distributions.
 * @class   GSMSCTableSimplified
 * @author  M Novak
 * @date    November 2017
 */

class GSMottCorrection;
class MaterialCuts;

class GSMSCTableSimplified {
public:
  GSMSCTableSimplified(bool iselectron);
  ~GSMSCTableSimplified();

  void Initialize(double lownergylimit, double highenergylimit, const std::vector<bool> &activeregionv);

  // structure to store one GS transformed angular distribution (for a given s/lambda_el,s/lambda_elG1)
  struct GSMSCAngularDtr {
    int fNumData; // # of data points
    double fDelta;
    struct Data {
      double fUValues; // array of transformed variables
      double fParamA;  // array of interpolation parameters a
      double fParamB;  // array of interpolation parameters b
    };
    std::vector<Data> fData;
    GSMSCAngularDtr() : fNumData(0), fData(0), fDelta(0.0) {}
  };

  bool Sampling(double lambdaval, double qval, double scra, double &cost, double &sint, GSMSCAngularDtr **gsDtr,
                double &transfPar, geant::TaskData *td, bool isfirst);

  void SampleTheta12(const double *lambdaval, const double *qval, const double *scra, double *cost1, double *cost2,
                     int N, geant::TaskData *td);

  double SampleCosTheta(double lambdaval, double qval, double scra, GSMSCAngularDtr **gsDtr, double &transfPar,
                        geant::TaskData *td, bool isfirst);

  double SampleGSSRCosTheta(const GSMSCAngularDtr *gsDrt, double transfpar, geant::TaskData *td);

  void SampleGSSRCosThetaVector(GSMSCAngularDtr **gsDrt, const double *transfpar, double *cost, int N,
                                geant::TaskData *td);

  GSMSCAngularDtr *GetGSAngularDtr(double scra, double lambdaval, double qval, double &transfpar, geant::TaskData *td);
  GSMSCAngularDtr *GetGSAngularDtr(double scra, double lambdaval, double lLambda, double qval, double &transfpar,
                                   geant::TaskData *td);

  // material dependent MSC parameters (computed at initialisation) regarding
  // Moliere's screening parameter
  double GetMoliereBc(const int matindx) { return gMoliere[matindx].Bc; }

  double GetMoliereXc2(const int matindx) { return gMoliere[matindx].Xc2; }

  struct MoliereData {
    double Bc, Xc2;
  };

private:
  void LoadMSCData();
  // initialisation of material dependent Moliere's MSC parameters
  void InitMoliereMSCParams();

private:
  double SampleWithSingle(double scra, geant::TaskData *td);
  double SampleMSCWithSingle(double expn, double lambdaval, double scra, double rndTheta, geant::TaskData *td);

  static bool gIsInitialised;                  // are the precomputed angular distributions already loaded in?
  static constexpr int gLAMBNUM    = 64;       // # L=s/lambda_el in [fLAMBMIN,fLAMBMAX]
  static constexpr int gQNUM1      = 15;       // # Q=s/lambda_el G1 in [fQMIN1,fQMAX1] in the 1-st Q grid
  static constexpr int gQNUM2      = 32;       // # Q=s/lambda_el G1 in [fQMIN2,fQMAX2] in the 2-st Q grid
  static constexpr double gLAMBMIN = 1.0;      // minimum s/lambda_el
  static constexpr double gLAMBMAX = 100000.0; // maximum s/lambda_el
  static constexpr double gQMIN1   = 0.001;    // minimum s/lambda_el G1 in the 1-st Q grid
  static constexpr double gQMAX1   = 0.99;     // maximum s/lambda_el G1 in the 1-st Q grid
  static constexpr double gQMIN2   = 0.99;     // minimum s/lambda_el G1 in the 1-st Q grid
  static constexpr double gQMAX2   = 7.99;     // maximum s/lambda_el G1 in the 1-st Q grid
  //
  bool fIsElectron; // GS-table for e- (for e+ otherwise)
  //
  double fLogLambda0;        // ln(gLAMBMIN)
  double fLogDeltaLambda;    // ln(gLAMBMAX/gLAMBMIN)/(gLAMBNUM-1)
  double fInvLogDeltaLambda; // 1/[ln(gLAMBMAX/gLAMBMIN)/(gLAMBNUM-1)]
  double fInvDeltaQ1;        // 1/[(gQMAX1-gQMIN1)/(gQNUM1-1)]
  double fDeltaQ2;           // [(gQMAX2-gQMIN2)/(gQNUM2-1)]
  double fInvDeltaQ2;        // 1/[(gQMAX2-gQMIN2)/(gQNUM2-1)]

  // vector to store all GS transformed angular distributions (cumputed based on the Screened-Rutherford DCS)
  static std::vector<GSMSCAngularDtr> gGSMSCAngularDistributions1;
  static std::vector<GSMSCAngularDtr> gGSMSCAngularDistributions2;

  /** Precomputed \f$ b_lambda_{c} $\f and \f$ \chi_c^{2} $\f material dependent
  *   Moliere parameters that can be used to compute the screening parameter,
  *   the elastic scattering cross section (or \f$ \lambda_{e} $\f) under the
  *   screened Rutherford cross section approximation. (These are used in
  *   GSMSCModel if gIsUsePWATotalXsecData is FALSE.)
  */
  static std::vector<MoliereData> gMoliere;
};

} // namespace geantphysics

#endif // GSMSCTABLE_H
