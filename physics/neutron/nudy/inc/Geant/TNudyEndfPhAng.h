//===-- Nudy/TNudyEndfPhAng.h - Instruction class definition -------*- C++
//-*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \file This class processes neutron induced photon angular dist. it makes 
///  linearly interpolable if possible. In this data are given in lab frame only
/// \class TNudyEndfPhAng
/// \author H. Kumawat
/// \date April 2018
//===----------------------------------------------------------------------===//

#ifndef TNudyEndfPhAng_H
#define TNudyEndfPhAng_H

#include "Geant/TNudyEndfSigma.h"

namespace Nudy {
class TNudyEndfFile;
class TNudyEndfSec;
class TNudyEndfTab1;
class TNudyEndfTab2;
}

namespace NudyPhysics {
class TNudyEndfSigma;
}

typedef std::vector<double> rowd;
typedef std::vector<int> rowint;

#ifdef USE_ROOT
#include "Rtypes.h"
#endif

namespace NudyPhysics {

  class TNudyEndfPhAng : public TNudyEndfSigma{

public:
  TNudyEndfPhAng();
  /// \brief Default constructure
  TNudyEndfPhAng(Nudy::TNudyEndfFile *file);
  /// \brief constructure to be used
  virtual void CreateTab1(Nudy::TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, int mtf[], int law, double &x);
  /// \brief create tab1 record 
  virtual void ModifyTab1(Nudy::TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x);
  /// \brief modify Tab1 record after linearization
  virtual void CreateTab2(Nudy::TNudyEndfTab2 *secTab2, int &NE);
  /// \brief create tab2 record
  virtual ~TNudyEndfPhAng();

private:
  double RecursionLinearLeg1D(double x1, double x2, double pdf1, double pdf2);
  /// \brief linearization of angular distribution if it is given in terms if
  /// the legendre polynomial
  double RecursionLinearProb(double x1, double x2, double pdf1, double pdf2);
  /// \brief linearization of angular distribution if it is given in terms if
  /// the probability distribution
  void GenerateProb(Nudy::TNudyEndfTab1 *tab);
  /// \brief this function creates cosine probabilities and calls  
  /// linearization program
  rowd fEin, fLegendCoef;
  /// \brief Temp. variables for energy and legendre coefficient for
  /// angular distribution
  rowd fCosFile4, fCosPdfFile4;
  /// \brief cosine, pdf for angular distribution
  rowint fNbt1, fInt1;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  int fNr, fNp;
  /// \brief standard ENDF parameters for no .of regions and points for
  /// interpolation
  int fNK, fNI;
  /// \brief Total no. of discrete + continuous photon distributions, discrete distributions
  
#ifdef USE_ROOT
  ClassDef(TNudyEndfPhAng, 1) // class for an ENDF reconstruction
#endif
};

} // namespace
#endif
