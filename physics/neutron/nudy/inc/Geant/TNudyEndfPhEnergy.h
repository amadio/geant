//===-- Nudy/TNudyEndfPhEnergy.h - Instruction class definition -------*- C++
//-*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \file This class processes neutron induced photon enerrgy dist. it makes 
///  linearly interpolable if possible and also makes in LAB frame if possible
/// \class TNudyEndfPhEnergy
/// \author H. Kumawat
/// \date March 2016
//===----------------------------------------------------------------------===//
#ifndef TNudyEndfPhEnergy_H
#define TNudyEndfPhEnergy_H

#include "Geant/TNudyEndfRecoPoint.h"

namespace Nudy {
class TNudyEndfFile;
}

namespace NudyPhysics {
class TNudyEndfRecoPoint;
}

#define PI acos(-1.0)
typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
typedef std::vector<rowd> matrixd2;
typedef std::vector<std::vector<rowd> > matrixd3;

#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
#endif

namespace NudyPhysics {
class TNudyEndfPhEnergy : public NudyPhysics::TNudyEndfRecoPoint {

public:
  TNudyEndfPhEnergy();
  /// \brief Default constructure
  TNudyEndfPhEnergy(Nudy::TNudyEndfFile *file);
  /// \brief constructor to be called for in Reco point
  virtual double GetEnergy5(int elemid, int mt, double energyK);
  /// \brief getter for the energy distribution based in the element ID, MT and
  /// energy of the neutron
  virtual double GetDelayedFraction(int ielemId, int mt, double energyK);
  /// \brief getter for delayed neutron fraction based in the element ID, MT and
  /// energy of the neutron
  virtual ~TNudyEndfPhEnergy();

private:
  double RecursionLinearFile5Prob(double x1, double x2, double pdf1,
                                  double pdf2);
  /// \brief linearization function if probabilities are given
  double RecursionLinearFile5GenEva(double x1, double x2, double pdf1,
                                    double pdf2, double energy);
  /// \brief linearization function if General evaporation spectra is to be
  /// constructed
  double RecursionLinearFile5Maxwell(double x1, double x2, double pdf1,
                                     double pdf2, double energy);
  /// \brief linearization function if Maxwellian spectra is to be constructed
  double RecursionLinearFile5Watt(double x1, double x2, double pdf1,
                                  double pdf2, double energy);
  /// \brief linearization function if Watt spectra is to be constructed
  void FillPdf1D();
  /// \brief filling 1 dimentional pdf for energy distribution
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
  /// \brief standard ENDF parameters for no .of regions and points for
  /// interpolation
  rowint fNbt3, fInt3;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  int fNr3, fNp3;
  /// \brief standard ENDF parameters for no .of regions and points for
  /// interpolation
  rowd fEnergyFile5, fEnergyPdfFile5, fEnergyCdfFile5;
  /// \brief energy, pdf and cdf for energy distribution for the file 5 with
  /// integrated angular distribution
  rowd fEin, fEneE, fCdf, fPdf;
  /// \brief Temp. energy, cdf and pdf for energy distribution for the file 5
  matrixd2 fEne2D, fFrac2D, fCdf2D, fPdf2D, fEin2D;
  /// \brief Temp. energy, cdf and pdf for energy distribution for the file 5
  matrixd3 fEne3D, fCdf3D, fPdf3D;
  /// \brief Temp. energy, cdf and pdf for energy distribution for the file 5
#ifdef USE_ROOT
  TRandom3 *fRnd;
#endif
  ClassDef(TNudyEndfPhEnergy, 1) // class for an ENDF energy reconstruction
};

} // namespace
#endif
