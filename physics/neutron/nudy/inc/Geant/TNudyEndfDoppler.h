//===-- Nudy/TNudyEndfDoppler.h - Instruction class definition -------*- C++ -*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \brief The class reads doppler broadened cross-sections from the root file
///  and gives all the corss-sections/ secondary parameters forthe outgoing particles
/// \class TNudyEndfDoppler
/// \author H. Kumawat
/// \date March 2016
//===----------------------------------------------------------------------===//
#ifndef TNudyEndfDoppler_H
#define TNudyEndfDoppler_H

#include "Geant/TNudyEndfSigma.h"

namespace NudyPhysics {

class TNudyEndfDoppler : public TNudyEndfSigma {

public:
  TNudyEndfDoppler();
  /// \brief Default constructure
  TNudyEndfDoppler(double isigDiff, double aw, double t1, double t2, std::vector<double> &x1, std::vector<double> &x2);
  /// \brief constructor to be called for in TNudyEndfSigma class
  virtual ~TNudyEndfDoppler() {}
  std::vector<double> fSigma;
  /// \brief temp. cross-section variable
private:
  double fOVSQPI = 0.564189583547756279;
  /// \brief constant from Cullen's description of the doppler broadening with continuous cross-sections
  double fBoltz = 8.617385E-05;
  /// \brief boltzman constant
  double fONE = 1.0, fTWO = 2.0, fTHH = 3.0 / fTWO, fTW3 = fTWO / 3.0, fZERO = 0.0, fHALF = fONE / fTWO;
  /// \brief numerical constants
  double fZLIMI = 5;
  /// \brief limit on integration for the doppler
  double fTk;
  /// \brief temperature difference in kelvin
  double fALPHA;
  /// \brief constant
  double fAWRI;
  /// \brief mass in units of neutron mass
  double fZKT, fY2;
  /// \brief temp. variable
  double fZK2P, fZK22P, fEXPAP, fF0K2P, fF1K2P, fF2K2P, fF3K2P, fF4K2P;
  /// \brief temp. variables
  double fZK2, fZK22, fEXPA, fF0K2, fF1K2, fF2K2, fF3K2, fF4K2;
  /// \brief temp. variables
  double fZK1, fF0K1, fF1K1, fF2K1, fF3K1, fF4K1;
  /// \brief temp. variables
  double fFACT, fAK, fCK, fCKY, fCKY2;
  /// \brief temp. variables
  double fFTAIL1 = 0, fFTAIL2 = 0;
  /// \brief temp. variables
  double fE1, fE2, fS1, fS2;
  /// \brief temp. variables
  double fZK1P, fF0K1P, fF1K1P, fF2K1P, fF3K1P, fF4K1P;
  /// \brief temp. variables
  double fRATLOW = 1, fRATHIG = 2;
  /// \brief constants
  double fRATHLF;
  /// \brief temp. variable
  double fHTEST;
  /// \brief temp. variable
  double fRATMAX;
  /// \brief temp. variable
  double fY;
  /// \brief temp. variable
  double fXSUM, fSigDiff, fXss;
  /// \brief temp. variables
  int fNcrs, fIPP, fKPP, fSize, fMipp, fJLoop, fMLoop;
  /// \brief temp. variables
  double RecursionLinear1(std::vector<double> &x1, std::vector<double> &x2, double x, double y, double sig, double xd,
                          double yd, double sigd);
  /// \brief recursive linear function to linearize far from peak region
  double BroadMore(std::vector<double> &x1, std::vector<double> &x2, double xp);
/// \brief Broadening far from peak region
#ifdef USE_ROOT
  ClassDef(TNudyEndfDoppler, 1)
#endif
};

} // namespace
#endif
