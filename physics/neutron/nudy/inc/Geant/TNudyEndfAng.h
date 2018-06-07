//===-- Nudy/TNudyEndfAng.h - Instruction class definition -------*- C++ -*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \class TNudyEndfAng
/// \author H. Kumawat
/// \date July 2016
//===----------------------------------------------------------------------===//
#ifndef TNudyEndfAng_H
#define TNudyEndfAng_H

#include "Geant/TNudyEndfSigma.h"

namespace NudyPhysics {
class TNudyEndfSigma;
}

namespace Nudy {
class TNudyEndfFile;
class TNudyEndfSec;
class TNudyEndfTab1;
class TNudyEndfTab2;
}

typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
#define PI acos(-1.0)

#ifdef USE_ROOT
#include "Rtypes.h"
#endif

namespace NudyPhysics {
class TNudyEndfAng : public NudyPhysics::TNudyEndfSigma {

public:
  TNudyEndfAng();
  /// \brief default constructure
  TNudyEndfAng(Nudy::TNudyEndfFile *file);
  /// \brief constructure to be used
  virtual void ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, rowint &fint,rowd &x1, rowd &x2);
  /// \brief process tab1 entry
  virtual void ModifyTab1(Nudy::TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x);
  /// \brief modify Tab1 record after linearization
  virtual void CreateTab2(Nudy::TNudyEndfTab2 *secTab2, int &NE);
  /// \brief create tab2 record
  virtual void CreateTab1(Nudy::TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, int mtf[], int law, double &x);
  /// \brief create tab1 record 
  virtual ~TNudyEndfAng();

private:
  double RecursionLinearLeg1D(double x1, double x2, double pdf1, double pdf2);
  /// \brief linearization of angular distribution if it is given in terms if the legendre polynomial
  double RecursionLinearProb(double x1, double x2, double pdf1, double pdf2);
  /// \brief linearization of angular distribution if it is given in terms if the probability distribution
  void LegendreToCosine(TIter &iter, Nudy::TNudyEndfSec *sec);
  /// \brief converting legendre coefficient into cosine and pdf distribution for neutron with tab2+tab1 records
  void CosineToLinearCosine(TIter &iter);
  /// \brief converting consine into linearized probability distribution to be used in TNudyEndfReco class
  void LegendreToCosinePlusCosine(TIter &recIter, Nudy::TNudyEndfSec *sec);
  /// \brief converting legendre coefficient to probability distribution and 
  ///  consine into linearized probability distribution to be used in TNudyEndfReco class
  rowd fEin, fLegendCoef;
  /// \brief Temp. variables for energy and legendre coefficient for angular distribution
  rowd fCosFile4, fCosPdfFile4, fCosCdfFile4;
  /// \brief cosine, pdf and cdf for angular distribution
  rowint fNbt1, fInt1;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  int fNR, fNP;
  /// \brief standard ENDF parameters for no .of regions and points for interpolation
  int fLTT, fNLTT;
  /// \brief flag to identify legendre vs probability distribution and no. of energy point in the LTT=3 legendre
  #ifdef USE_ROOT
  ClassDef(TNudyEndfAng, 1) // class for an ENDF reconstruction
#endif
};

} // namespace
#endif
