//===-- Nudy/TNudyEndfNuPh.h - Instruction class definition -------*- C++ -*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \class TNudyEndfNuPh
/// \author H. Kumawat
/// \date July 2016
//===----------------------------------------------------------------------===//
#ifndef TNudyEndfNuPh_H
#define TNudyEndfNuPh_H

#include "Geant/TNudyEndfRecoPoint.h"

#define PI acos(-1.0)
typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
#ifdef USE_ROOT
#include "Rtypes.h"
#endif

namespace Nudy {
class TNudyEndfFile;
}

namespace NudyPhysics {
class TNudyEndfNuPh : public TNudyEndfRecoPoint {

public:
  TNudyEndfNuPh();
  /// \brief Default constructure
  TNudyEndfNuPh(Nudy::TNudyEndfFile *file);
  /// \brief constructure to be called in recopoint
  virtual double GetNuTotal(int elemid, double energyK);
  /// \brief getting fission neutron multiplicity (Total)
  virtual double GetNuPrompt(int elemid, double energyK);
  /// \brief getting fission neutron multiplicity (Prompt)
  virtual double GetNuDelayed(int elemid, double energyK);
  /// \brief getting fission neutron multiplicity (Delayed)
  virtual double GetFissHeat(int elemid, double energyK);
  /// \brief getting fission Heat
  virtual double GetLambdaD(int elemid, int time);
  /// \brief getting fission neutron Lambda for the decaying family
  virtual ~TNudyEndfNuPh();

private:
  double RecursionLinearNuPh(double x1, double x2, double sig1, double sig2, std::vector<double> x,
                             std::vector<double> sig);
  /// \brief recursive linear for cross-section
  int fNR, fNP;
  /// \brief standard ENDF parameters for range and interpolation
  rowd fEintFile1, fNutFile1, fEinFile1, fNuFile1;
  /// \brief energy, total fission neutron multiplicity, energy, prompt fission neutron multiplicity
  rowd fEindFile1, fNudFile1, fEinPhFile1, fPhFile1;
  /// \brief energy, delayed fission neutron multiplicity, energy, photon
  rowd fEinfFile1, fHeatFile1;
  /// \brief energy, fission heat
  rowd fCnc, fNui;
  /// \brief coefficients for getting neutron multiplicity \cite ENDF manual
  matrixd2 fEint, fNut;
  /// \brief incident energy and total nu,  all elements
  rowint fNbt1, fInt1;
  /// \brief endf interpolation parameter
  double fSigDiff;
  /// \brief precision/tolerance for cross-section reconstruction while linearization from true values
  ClassDef(TNudyEndfNuPh, 1) // class for an ENDF reconstruction
};

} // namespace
#endif
