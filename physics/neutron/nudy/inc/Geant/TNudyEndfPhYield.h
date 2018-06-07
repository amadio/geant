//===-- Nudy/TNudyEndfPhYield.h - Instruction class definition -------*- C++
//-*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \file This class processes neutron induced photon production multiplicities 
/// and transition probability arrays from file 12 of the endf format. 
/// It makes linearly interpolable.
/// \class TNudyEndfPhYield
/// \author H. Kumawat
/// \date April 2018
//===----------------------------------------------------------------------===//
#ifndef TNudyEndfPhYield_H
#define TNudyEndfPhYield_H

#include "Geant/TNudyEndfSigma.h"

namespace NudyPhysics {
class TNudyEndfSigma;
}

namespace Nudy {
class TNudyEndfFile;
}

typedef std::vector<double> rowd;
typedef std::vector<int> rowi;
#ifdef USE_ROOT
#include "Rtypes.h"
#endif

namespace NudyPhysics {
class TNudyEndfPhYield : public TNudyEndfSigma {

public:
  TNudyEndfPhYield();
  /// \brief Default constructure
  TNudyEndfPhYield(Nudy::TNudyEndfFile *file);
  /// \brief constructure to be use
  virtual void ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, rowint &fint,rowd &x1, rowd &x2);
  /// \brief process tab1 entry
  virtual void ModifyTab1(Nudy::TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x);
  /// \brief modify Tab1 record after linearization
  virtual ~TNudyEndfPhYield();

private:
  double RecursionLinear(double x1, double x2, double pdf1, double pdf2);
  rowd fEIN, fMultiPh;
  /// \brief photon energy and multiplicity
  rowd fENS, fTP, fGP;                          
  /// \brief ENS = energy of NS level
  /// TP = direct photon transition probability, 
  /// GP = conditional photon transition probability
  rowi fNBT, fINT;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  double fES, fEG;
  /// \brief Energy of emission level of photon and photon energy
  int fNR, fNP;
  /// \brief fNR = No. of interpolation regions, fNP = number of energy points
  int fLP, fLF;
  /// \brief fLP = flag for primary, non-primary photon, photon energy if == 0,1 Egp = eg +awr/(1+awr) if ==2,
  /// fLF = energy distribution law, tabulated in file 15 if == 1, discrete if == 2 

  ClassDef(TNudyEndfPhYield, 1) 
};

} // namespace
#endif
