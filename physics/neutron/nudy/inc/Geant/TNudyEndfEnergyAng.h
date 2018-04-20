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

#include "Geant/TNudyEndfRecoPoint.h"

namespace NudyPhysics {
class TNudyEndfRecoPoint;
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
class TRandom3;
#endif

namespace NudyPhysics {
class TNudyEndfEnergyAng : public NudyPhysics::TNudyEndfRecoPoint {

public:
  TNudyEndfEnergyAng();
  /// \brief Default constructure
  TNudyEndfEnergyAng(Nudy::TNudyEndfFile *file, double[]);
  /// \brief constructure to be called in recopoint
  virtual double GetCos64(int elemid, int mt, double energyK);
  virtual double GetCos6(int elemid, int mt, double energyK);
  /// \brief getting cosine from file 6
  virtual double GetEnergy6(int elemid, int mt, double energyK);
  /// \brief getting energy from file 6
  virtual int GetLaw6(int ielemId, int mt);
  /// \brief getting rule to get the energy and angle from file 6
  virtual int GetZd6(int ielemId, int mt);
  /// \brief getting Z value after particle evaporation
  virtual int GetAd6(int ielemId, int mt);
  /// \brief getting A value after particle evaporation
  virtual int GetMt6Neutron(int ielemId, int mt);
  /// \brief getting flag if energy-angle distribution is given for the neutron or photon
  virtual ~TNudyEndfEnergyAng();

private:
  double RecursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2);
  /// \brief recursive linear for file 5 type probability distribution
  double RecursionLinearFile4(int i, double x1, double x2, double pdf1, double pdf2);
  /// \brief recursive linear for file 4 type probability distribution
  double RecursionLinearProb(double x1, double x2, double pdf1, double pdf2);
  /// \brief recursive linear for file 6 type probability distribution

  int fNR, fNP;
  /// \brief standard ENDF parameters for range and no. of point for interpolation
  int fNR2, fNE2;
  /// \brief standard ENDF parameters for range and no. of point for interpolation
  rowd fCosFile4, fEnergyFile5;
  /// \brief cosine and energy
  rowd fCosPdfFile4, fEnergyPdfFile5;
  /// \brief pdf for cosine and energy
  rowd fCosCdfFile4, fEnergyCdfFile5;
  /// \brief cdf for cosine and energy
  rowd fEdes6, fF06, fR6, fA6;
  /// \brief standard file 6 parameters from ENDF manual
  rowd fE1, fP1, fE2, fP2, fE3, fP3, INorm;
  /// \brief standard file 6 parameters from ENDF manual
  rowint fZd, fAd;
  /// \brief z and A
  rowint fLaw;
  /// \brief law6 numbers for endf file 6
  rowint fMtNumbers, fMtNumbers6, fMtNumbers4;
  /// \brief temp. MT numbers
  rowint fMtNumNeutron, fMtNumPhoton, fMtNumCharge;
  /// \brief MT numbers either for neutron, photon or charge particles
  rowint fMtLct;
  /// \brief temp. LCT numbers (flag for angular distribution in cm or lab system)
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
  rowd fEin, fCosc, fCdfc, fPdfc, fLegendCoef1, fCosIn, fCosInPdf, fCosInCdf;
  /// \brief temp. energy, pdf and cdf for energy distribution
  rowd fEoute, fCdfe, fPdfe;
  /// \brief temp. energy, pdf and cdf for energy distribution
  matrixint fMt6Values;
  /// \brief MT numbers in file 6
  matrixint fMt6Neutron, fMt6Photon, fMt6Charge;
  /// \brief MT numbers if it is for neutron photon or charge particle
  matrixd2 fCos2d, fCosinpdf2d, fCosincdf2d, fCos2dc, fPdf2dc, fCdf2dc, fLegendCoef, fEin2d, fEin2dc;
  /// \brief temp. energy, pdf and cdf for energy distribution
  matrixd3 fCos3d, fCosinpdf3d, fCosincdf3d, fCos3dc, fPdf3dc, fCdf3dc;
  /// \brief temp. energy, pdf and cdf for energy distribution
  matrixd2 fEout2de, fPdf2de, fCdf2de;
  /// \brief temp. energy, pdf and cdf for energy distribution
  matrixd3 fEout3de, fPdf3de, fCdf3de;
  /// \brief temp. energy, pdf and cdf for energy distribution
  matrixd4 fEout4de, fPdf4de, fCdf4de;
/// \brief energy, pdf and cdf for energy distribution
#ifdef USE_ROOT
  TRandom3 *fRnd;
#endif
#ifdef USE_ROOT
  ClassDef(TNudyEndfEnergyAng, 1) // class for an ENDF reconstruction
#endif
};

} // namespace
#endif
