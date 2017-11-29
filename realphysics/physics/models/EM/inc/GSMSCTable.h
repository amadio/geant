
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

/**
 * @brief   Helper class for the GSMSCModel to handle the pre-computed angular distributions.
 * @class   GSMSCTable
 * @author  M Novak
 * @date    November 2017
 */

class GSMottCorrection;
class MaterialCuts;

class GSMSCTable {
public:
  GSMSCTable(bool iselectron);
 ~GSMSCTable();

  void Initialize(double lownergylimit, double highenergylimit, const std::vector<bool>& activeregionv);

  // structure to store one GS transformed angular distribution (for a given s/lambda_el,s/lambda_elG1)
  struct GSMSCAngularDtr {
    int     fNumData;    // # of data points
    double *fUValues;    // array of transformed variables
    double *fParamA;     // array of interpolation parameters a
    double *fParamB;     // array of interpolation parameters b
  };

  bool   Sampling(double lambdaval, double qval, double scra, double &cost, double &sint, double lekin, double beta2,
                  int matindx, GSMSCAngularDtr **gsDtr, int &mcekini, int &mcdelti, double &transfPar,
                  Geant::GeantTaskData *td, bool isfirst);

  double SampleCosTheta(double lambdaval, double qval, double scra, double lekin, double beta2, int matindx,
                        GSMSCAngularDtr **gsDtr, int &mcekini, int &mcdelti, double &transfPar,
                        Geant::GeantTaskData *td, bool isfirst);

  double SampleGSSRCosTheta(const GSMSCAngularDtr* gsDrt, double transfpar, Geant::GeantTaskData *td);

  double SingleScattering(double lambdaval, double scra, double lekin, double beta2, int matindx,
                          Geant::GeantTaskData *td);

  GSMSCAngularDtr* GetGSAngularDtr(double scra, double &lambdaval, double &qval, double &transfpar,
                                   Geant::GeantTaskData *td);

  // material dependent MSC parameters (computed at initialisation) regarding
  // Moliere's screening parameter
  double GetMoliereBc(const int matindx)  { return gMoliereBc[matindx];  }

  double GetMoliereXc2(const int matindx) { return gMoliereXc2[matindx]; }

  void   GetMottCorrectionFactors(double logekin, double beta2, int matindx, double &mcToScr, double &mcToQ1,
                                  double &mcToG2PerG1);

  // set option to activate/inactivate Mott-correction
  void   SetOptionMottCorrection(bool val) { fIsMottCorrection = val; }
  // set option to activate/inactivate PWA-correction
  void   SetOptionPWACorrection(bool val)  { fIsPWACorrection = val;  }

  // this method returns with the scattering power correction (to avoid double counting of sub-threshold deflections)
  // interpolated from tables prepared at initialisation
  double ComputeScatteringPowerCorrection(const MaterialCuts *matcut, double ekin);

  void   InitSCPCorrection(const std::vector<bool>& activeregionv);

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
  //
  bool   fIsElectron;          // GS-table for e- (for e+ otherwise)
  bool   fIsMottCorrection;    // flag to indicate if Mott-correction was requested to be used
  bool   fIsPWACorrection;     // flag to indicate is PWA corrections were requested to be used
  //
  double fLogLambda0;          // ln(gLAMBMIN)
  double fLogDeltaLambda;      // ln(gLAMBMAX/gLAMBMIN)/(gLAMBNUM-1)
  double fInvLogDeltaLambda;   // 1/[ln(gLAMBMAX/gLAMBMIN)/(gLAMBNUM-1)]
  double fInvDeltaQ1;          // 1/[(gQMAX1-gQMIN1)/(gQNUM1-1)]
  double fDeltaQ2;             // [(gQMAX2-gQMIN2)/(gQNUM2-1)]
  double fInvDeltaQ2;          // 1/[(gQMAX2-gQMIN2)/(gQNUM2-1)]
  //
  double fLowEnergyLimit;
  double fHighEnergyLimit;
  //
  int      fNumSPCEbinPerDec;    // scattering power correction energy grid bins per decade
  struct SCPCorrection {
    bool   fIsUse;               //
    double fPrCut;               // sec. e- production cut energy
    double fLEmin;               // log min energy
    double fILDel;               // inverse log delta kinetic energy
    std::vector<double> fVSCPC;  // scattering power correction vector
  };
  std::vector<SCPCorrection*>  fSCPCPerMatCuts;
  std::vector<bool>            fActiveRegionsVector;

  // vector to store all GS transformed angular distributions (cumputed based on the Screened-Rutherford DCS)
  static std::vector<GSMSCAngularDtr*> gGSMSCAngularDistributions1;
  static std::vector<GSMSCAngularDtr*> gGSMSCAngularDistributions2;

   /** Precomputed \f$ b_lambda_{c} $\f and \f$ \chi_c^{2} $\f material dependent
   *   Moliere parameters that can be used to compute the screening parameter,
   *   the elastic scattering cross section (or \f$ \lambda_{e} $\f) under the
   *   screened Rutherford cross section approximation. (These are used in
   *   GSMSCModel if gIsUsePWATotalXsecData is FALSE.)
   */
   static std::vector<double> gMoliereBc;
   static std::vector<double> gMoliereXc2;
   //
   //
   GSMottCorrection   *fMottCorrection;
};

}      // namespace geantphysics

#endif // GSMSCTABLE_H
