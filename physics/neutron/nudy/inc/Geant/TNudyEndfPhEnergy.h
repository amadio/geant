//===-- Nudy/TNudyEndfPhEnergy.h - Instruction class definition -------*- C++
//-*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \file This class processes neutron induced photon enerrgy dist. it makes 
///  linearly interpolable. data are given in LAB frame only
/// \class TNudyEndfPhEnergy
/// \author H. Kumawat
/// \date April 2018
//===----------------------------------------------------------------------===//
#ifndef TNudyEndfPhEnergy_H
#define TNudyEndfPhEnergy_H

#include "Geant/TNudyEndfSigma.h"

namespace Nudy {
class TNudyEndfFile;
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
class TNudyEndfPhEnergy : public NudyPhysics::TNudyEndfSigma {

public:
  TNudyEndfPhEnergy();
  /// \brief Default constructure
  TNudyEndfPhEnergy(Nudy::TNudyEndfFile *file);
  /// \brief constructor to be called for in Reco point
  virtual void ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, rowint &fint,rowd &x1, rowd &x2);
  /// \brief process tab1 entry
  virtual void ProcessTab2(Nudy::TNudyEndfTab2 *tab2, int &NR,int &NP,rowint &fnbt, rowint &fint); 
  /// \brief process tab2 entry 
  virtual void ModifyTab1(Nudy::TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x);
  /// \brief modify Tab1 record after linearization
  virtual ~TNudyEndfPhEnergy();

private:
  double RecursionLinearFile5Prob(double x1, double x2, double pdf1,
                                  double pdf2);
  /// \brief linearization function if probabilities are given
  int fNr1, fNp1;
  /// \brief standard ENDF parameters for range and interpolation
  rowd fE1, fP1, fE2, fP2, fE3, fP3;
  /// \brief standard file 5 parameters from ENDF manual
  rowint fNbt1, fInt1;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  rowint fNbt2, fInt2;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  int fNr2, fNp2;
  /// \brief standard ENDF parameters for no .of regions and points for
  /// interpolation
  rowint fNbt3, fInt3;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  int fNr3, fNp3;
  /// \brief standard ENDF parameters for no .of regions and points for
  /// interpolation
  rowd fEin;
  /// \brief Temp. energy for energy distribution for the file 15
  rowd fEnergyFile5, fEnergyPdfFile5;
  /// \brief energy, pdf for energy distribution for the file 15 with
  /// integrated angular distribution

  ClassDef(TNudyEndfPhEnergy, 1) // class for an ENDF energy reconstruction
};

} // namespace
#endif
