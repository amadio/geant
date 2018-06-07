//===-- Nudy/TNudyEndfPhProd.h - Instruction class definition -------*- C++
//-*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \file This class processes neutron induced photon production cross-section 
/// from file 13 of the endf format. It makes linearly interpolable. The cross
/// sections are independent from file 3 (absolute values are given here)
/// \class TNudyEndfPhProd
/// \author H. Kumawat
/// \date April 2018
//===----------------------------------------------------------------------===//
#ifndef TNudyEndfPhProd_H
#define TNudyEndfPhProd_H

#include "Geant/TNudyEndfSigma.h"

namespace NudyPhysics {
class TNudyEndfSigma;
}

namespace Nudy {
class TNudyEndfFile;
class TNudyEndfTab1;
}

typedef std::vector<double> rowd;
typedef std::vector<int> rowi;

#ifdef USE_ROOT
#include "Rtypes.h"
#endif

namespace NudyPhysics {
class TNudyEndfPhProd : public NudyPhysics::TNudyEndfSigma {

public:
  TNudyEndfPhProd();
  /// \brief Default constructure
  TNudyEndfPhProd(Nudy::TNudyEndfFile *file);
  /// \brief constructure to be use
  virtual void ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, rowint &fint,rowd &x1, rowd &x2);
  /// \brief process tab1 entry
  virtual void ModifyTab1(Nudy::TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x);
  /// \brief modify Tab1 record after linearization
  virtual ~TNudyEndfPhProd();

private:
  double RecursionLinear(double x1, double x2, double pdf1, double pdf2);
  rowd fEIN, fSigPh;
  /// \brief Neutron energy and cross-section
  rowi fNBT, fINT;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  double fLP, fLF;
  /// \brief fLP = flag for primary, non-primary photon, photon energy if == 0,1 Egp = eg +awr/(1+awr) if ==2,
  /// fLF = energy distribution law, tabulated in file 15 if == 1, discrete if == 2 
  double fES, fEG;
  /// \brief Energy of emission level of photon, photon energy
  int fNR, fNP;
  /// \brief fNR = No. of interpolation regions, fNP = number of energy points

  ClassDef(TNudyEndfPhProd, 1)
};

} // namespace
#endif
