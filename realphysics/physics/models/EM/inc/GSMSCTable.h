

#ifndef GSMSCTABLE_H
#define GSMSCTABLE_H

// from geantV
#include "Geant/Config.h"
namespace Geant {
  inline namespace GEANT_IMPL_NAMESPACE {
  class GeantTaskData;
}
}

#include <vector>


namespace geantphysics {

class GSMSCTable {
public:
  GSMSCTable();
 ~GSMSCTable();

  void Initialize();

  // structure to store one GS transformed angular distribution (for a given s/lambda_el,s/lambda_elG1)
  struct GSMSCAngularDtr {
    int     fNumData;    // # of data points
    double  fQScale;
    double *fUValues;    // array of transformed variables
    double *fParamA;     // array of interpolation parameters a
    double *fParamB;     // array of interpolation parameters b
  };

  // only for testing
  GSMSCAngularDtr* GetOne(int indx) {return fGSMSCAngularDistributions1[indx];}

  void          Sampling(double lambdaval, double qval, double scra, double &cost, double &sint,  Geant::GeantTaskData *td);
  double SampleCosTheta (double lambdaval, double qval, double scra, double rndm1, double rndm2, double rndm);
  double SampleCosTheta1(double lambdaval, double qval, double scra, double rndm1, double rndm2, double rndm);
  double SampleCosTheta2(double lambdaval, double qval, double scra, double rndm1, double rndm2, double rndm);
  double GetScreeningParam(double G1);

  // material dependent MSC parameters (computed at initialisation) regarding
  // Moliere's screening parameter
  double GetMoliereBc (int matindx) { return fMoliereBc[matindx]; }
  double GetMoliereXc2(int matindx) { return fMoliereXc2[matindx]; }

private:
  void LoadMSCData();
  // initialisation of material dependent Moliere's MSC parameters
  void InitMoliereMSCParams();


private:
  static bool             gIsInitialised;       // are the precomputed angular distributions already loaded in?
  static constexpr int    gLAMBNUM = 64;        // # L=s/lambda_el in [fLAMBMIN,fLAMBMAX]
  static constexpr int    gQNUM1   = 15;        // # Q=s/lambda_el G1 in [fQMIN1,fQMAX1] in the 1-st Q grid
  static constexpr int    gQNUM2   = 32;        // # Q=s/lambda_el G1 in [fQMIN2,fQMAX2] in the 2-st Q grid
  static constexpr int    gNUMSCR1 = 201;       // # of screening parameters in the A(G1) function
  static constexpr int    gNUMSCR2 = 51;        // # of screening parameters in the A(G1) function
  static constexpr double gLAMBMIN = 1.0;       // minimum s/lambda_el
  static constexpr double gLAMBMAX = 100000.0;  // maximum s/lambda_el
  static constexpr double gQMIN1   = 0.001;     // minimum s/lambda_el G1 in the 1-st Q grid
  static constexpr double gQMAX1   = 0.99;      // maximum s/lambda_el G1 in the 1-st Q grid
  static constexpr double gQMIN2   = 0.99;      // minimum s/lambda_el G1 in the 1-st Q grid
  static constexpr double gQMAX2   = 7.99;      // maximum s/lambda_el G1 in the 1-st Q grid
  // precomputed A(G1) function with its interpolation parameters
  static constexpr double gSCRMIN1 = 1.93214991408357e-12;
  static constexpr double gSCRMAX1 = 2.42974344203683e-01;
  static constexpr double gSCRMAX2 = 5.50564555556202e+01;
  //
  static const double gG1Values1[];
  static const double gScrAValues1[];
  static const double gScrBValues1[];
  static const double gG1Values2[];
  static const double gScrAValues2[];
  static const double gScrBValues2[];

  double fLogLambda0;          // ln(gLAMBMIN)
  double fLogDeltaLambda;      // ln(gLAMBMAX/gLAMBMIN)/(gLAMBNUM-1)
  double fInvLogDeltaLambda;   // 1/[ln(gLAMBMAX/gLAMBMIN)/(gLAMBNUM-1)]
  double fInvDeltaQ1;          // 1/[(gQMAX1-gQMIN1)/(gQNUM1-1)]
  double fDeltaQ2;             // [(gQMAX2-gQMIN2)/(gQNUM2-1)]
  double fInvDeltaQ2;          // 1/[(gQMAX2-gQMIN2)/(gQNUM2-1)]
  // for the precumputed A(G1) function
  double fLogG1FuncMin1;
  double fInvLogDeltaG1Func1;
  double fLogG1FuncMin2;
  double fInvLogDeltaG1Func2;

   // vector to store all GS transformed angular distributions
   std::vector<GSMSCAngularDtr*> fGSMSCAngularDistributions1;
   std::vector<GSMSCAngularDtr*> fGSMSCAngularDistributions2;

   /** Precomputed \f$ b_lambda_{c} $\f and \f$ \chi_c^{2} $\f material dependent
   *   Moliere parameters that can be used to compute the screening parameter,
   *   the elastic scattering cross section (or \f$ \lambda_{e} $\f) under the
   *   screened Rutherford cross section approximation. (These are used in
   *   GSMSCModel if gIsUsePWATotalXsecData is FALSE.)
   */
   std::vector<double> fMoliereBc;
   std::vector<double> fMoliereXc2;
};

}      // namespace geantphysics

#endif // GSMSCTABLE_H
