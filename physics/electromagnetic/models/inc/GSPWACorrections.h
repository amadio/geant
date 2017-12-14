
#ifndef GSPWACORRECTIONS_H
#define GSPWACORRECTIONS_H

#include "SystemOfUnits.h"

#include <vector>
#include <string>

// from geantV
#include "Geant/Config.h"

namespace geantphysics {
  inline namespace GEANT_IMPL_NAMESPACE {
    class Material;
    class Element;
  }
}



namespace geantphysics {

/**
 * @class   GSPWACorrections
 * @author  M Novak
 * @date    November 2017
 *
 *  Class to describe and store correction factors to the integrated quantities of GSMscModel (screening parameter,
 *   first and second moments) derived by using accurate Dirac-PWA based integrated quantities.
 */

class GSPWACorrections {
public:
  GSPWACorrections(bool iselectron=true);

 ~GSPWACorrections();

  void     Initialise(const std::vector<bool>& activeregionv);

  void     GetPWACorrectionFactors(double logekin, double beta2, int matindx, double &corToScr, double &corToQ1,
                                   double &corToG2PerG1);
private:
  void     InitDataPerElement(const std::vector<bool>& activeregionv);

  void     InitDataPerMaterials(const std::vector<bool>& activeregionv);

  void     LoadDataElement(const Element*);

  void     InitDataMaterial(const Material*);

  void     ClearDataPerElement();

  void     ClearDataPerMaterial();

  // either per material or per Z
  struct DataPerMaterial {
    std::vector<double>   fCorScreening;    // correction factor to Moliere screening parameter
    std::vector<double>   fCorFirstMoment;  // correction factor to first moment
    std::vector<double>   fCorSecondMoment; // correction factor to second
  };


// data members
private:
  bool   fIsElectron;
  static constexpr int     gMaxZet    = 98;                 // max. Z for which correction data were computed (98)
  static constexpr int     gNumEkin   = 31;                 // number of kinetic energy grid points for Mott correction
  static constexpr int     gNumBeta2  = 16;                 // \beta^2 values between [fMinBeta2-fMaxBeta2]
  static constexpr double  gMinEkin   =   1.*geant::keV;    // minimum kinetic energy value
  static constexpr double  gMidEkin   = 100.*geant::keV;    // kinetic energy at the border of the E_{kin}-\beta^2 grids
  static constexpr double  gMaxBeta2  =   0.9999;           // maximum \beta^2 value
  //
  double                   fMaxEkin;        // from max fMaxBeta2 = 0.9999 (~50.5889 [MeV])
  double                   fLogMinEkin;     // \ln[fMinEkin]
  double                   fInvLogDelEkin;  // 1/[\ln(fMidEkin/fMinEkin)/(fNumEkin-fNumBeta2)]
  double                   fMinBeta2;       // <= E_{kin}=100 [keV] (~0.300546)
  double                   fInvDelBeta2;    // 1/[(fMaxBeta2-fMinBeta2)/(fNumBeta2-1)]
  //
  static const std::string   gElemSymbols[];
  //
  std::vector<DataPerMaterial*>  fDataPerElement;   // size will be gMaxZet+1; won't be null only at used Z indices
  std::vector<DataPerMaterial*>  fDataPerMaterial;  // size will #materials; won't be null only at used mat. indices

};

}      // namespace geantphysics

#endif // GSPWACORRECTIONS_H
