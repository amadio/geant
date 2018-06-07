//===-- Nudy/TNudyEndfEnergyAng.h - Instruction class definition -------*- C++ -*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \class TNudyEndfEnergyAng
/// \author H. Kumawat
/// \date July 2016
//===----------------------------------------------------------------------===//
#ifndef TNudyEndfEnergyAng_H
#define TNudyEndfEnergyAng_H

#include "Geant/RngWrapper.h"
#include "Geant/TNudyEndfSigma.h"

namespace NudyPhysics {
class TNudyEndfSigma;
}

namespace Nudy {
class TNudyEndfFile;
}
#define PI acos(-1.0)
typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
typedef std::vector<rowd> matrixd2;
typedef std::vector<std::vector<rowd>> matrixd3;
typedef std::vector<std::vector<std::vector<rowd>>> matrixd4;
#ifdef USE_ROOT
#include "Rtypes.h"
#endif

namespace NudyPhysics {
class TNudyEndfEnergyAng : public NudyPhysics::TNudyEndfSigma {

public:
  TNudyEndfEnergyAng();
  /// \brief Default constructure
  TNudyEndfEnergyAng(Nudy::TNudyEndfFile *file, double[]);
  /// \brief constructure to be called in TNudyEndfSigma
  virtual void ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, rowint &fint,rowd &x1, rowd &x2);
  /// \brief process tab1 entry
  virtual void ProcessTab2(Nudy::TNudyEndfTab2 *tab2, int &NR,int &NP,rowint &fnbt, rowint &fint); 
  /// \brief process tab2 entry 
  virtual void ModifyTab1(Nudy::TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x);
  /// \brief modify Tab1 record after linearization
  virtual void CreateTab2(Nudy::TNudyEndfTab2 *secTab2, int &NE);
  /// \brief create tab2 record in file5 (MF = 5) to make LF=1 structure
  virtual void CreateTab1(Nudy::TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, int mtf[], int law, double &x);
  /// \brief create tab1 record in file5 (MF = 5) to make LF=1 structure
  virtual ~TNudyEndfEnergyAng();
  
private:
  double RecursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2);
  /// \brief recursive linear for file 5 type probability distribution
  double RecursionLinearProb(double x1, double x2, double pdf1, double pdf2);
  /// \brief recursive linear for file 6 type probability distribution
  double RecursionLinearLeg1D(double x1, double x2, double pdf1, double pdf2);
  /// \brief linearization of angular distribution if it is given in terms if the legendre polynomial
  void LegendreToCosine(rowd &x1, rowd &x2);
  /// \brief converting legendre coefficient into cosine and pdf distribution 
  
  int fNR, fNP;
  /// \brief standard ENDF parameters for range and no. of point for interpolation
  rowd fCosFile4, fEnergyFile5;
  /// \brief cosine and energy
  rowd fCosPdfFile4, fEnergyPdfFile5;
  /// \brief pdf for cosine and energy
  rowd fCosCdfFile4, fEnergyCdfFile5;
  /// \brief cdf for cosine and energy
  rowd fEdes6, fF06, fR6, fA6;
  /// \brief standard file 6 parameters from ENDF manual
  rowd fE1, fP1, fE2, fP2, fE3, fP3;
  /// \brief standard file 6 parameters from ENDF manual
  rowint fNbt1, fInt1;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  rowd fEint, fYi;
  /// \brief incident energy and yield of the reaction
  rowint fNbt2, fInt2;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  int fNr2, fNp2;
  /// \brief standard ENDF parameters for no .of regions and points for interpolation
  rowint fNbt3, fInt3;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  int fNr3, fNp3;
  /// \brief standard ENDF parameters for no .of regions and points for interpolation
  rowd fEin, fLegendCoef, fCosIn, fCosInPdf;
  /// \brief temp. energy, pdf and cdf for energy distribution
  rowd fEoute, fPdfe;
  /// \brief temp. energy, pdf and cdf for energy distribution
  matrixd2 fCos2d, fCosinpdf2d;
  /// \brief temp. energy, pdf and cdf for energy distribution
  matrixd2 fEout2de, fPdf2de,fLegendCoef1;
  /// \brief temp. energy, pdf and cdf for energy distribution
  matrixd3 fEout3de, fPdf3de;
  /// \brief temp. energy, pdf and cdf for energy distribution
  
  geant::RngWrapper fRng;
#ifdef USE_ROOT
  ClassDef(TNudyEndfEnergyAng, 1) // class for an ENDF reconstruction
#endif
};

} // namespace
#endif
