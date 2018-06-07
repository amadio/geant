//===-- Nudy/TNudyEndfEnergy.h - Instruction class definition -------*- C++ -*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \class TNudyEndfEnergy
/// \author H. Kumawat
/// \date July 2016
//===----------------------------------------------------------------------===//
#ifndef TNudyEndfEnergy_H
#define TNudyEndfEnergy_H

#include "Geant/TNudyEndfSigma.h"

namespace NudyPhysics {
class TNudyEndfSigma;
}

namespace Nudy {
class TNudyEndfFile;
class TNudyEndfTab1;
class TNudyEndfTab2;
}

#define PI acos(-1.0)
typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
typedef std::vector<rowd> matrixd2;
typedef std::vector<std::vector<rowd>> matrixd3;

#ifdef USE_ROOT
#include "Rtypes.h"
#endif

namespace NudyPhysics {

class TNudyEndfEnergy : public NudyPhysics::TNudyEndfSigma {

public:
  TNudyEndfEnergy();
  /// \brief Default constructure
  TNudyEndfEnergy(Nudy::TNudyEndfFile *file);
  /// \brief constructor to be called for in Reco point
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
  virtual ~TNudyEndfEnergy();

private:
  double RecursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2);
  /// \brief linearization function if probabilities are given
  double RecursionLinearFile5GenEva(double x1, double x2, double pdf1, double pdf2, double energy);
  /// \brief linearization function if General evaporation spectra is to be constructed
  double RecursionLinearFile5Maxwell(double x1, double x2, double pdf1, double pdf2, double energy);
  /// \brief linearization function if Maxwellian spectra is to be constructed
  double RecursionLinearFile5Watt(double x1, double x2, double pdf1, double pdf2, double energy);
  /// \brief linearization function if Watt spectra is to be constructed
  
  int fNR, fNP;
  /// \brief standard ENDF parameters for range and interpolation
  rowint fMtNumbers;
  /// \brief MT numbers temporary
  rowd fE1, fP1, fE2, fP2, fE3, fP3, INorm;
  /// \brief standard file 5 parameters from ENDF manual
  rowint fNbt1, fInt1;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  rowint fNbt2, fInt2;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  int fNr2, fNp2;
  /// \brief standard ENDF parameters for no .of regions and points for interpolation
  rowint fNbt3, fInt3;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  int fNr3, fNp3;
  /// \brief standard ENDF parameters for no .of regions and points for interpolation
  rowd fEnergyFile5, fEnergyPdfFile5;
  /// \brief energy and pdf for energy distribution for the file 5 with integrated angular distribution
  rowd fEin;
  /// \brief Temp. energy for energy distribution for the file 5
#ifdef USE_ROOT
  ClassDef(TNudyEndfEnergy, 1) // class for an ENDF energy reconstruction
#endif
};

} // namespace
#endif
