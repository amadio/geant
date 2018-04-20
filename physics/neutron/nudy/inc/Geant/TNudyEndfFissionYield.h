//===-- Nudy/TNudyEndfFissionYield.h - Instruction class definition -------*- C++ -*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \class TNudyEndfFissionYield
/// \author H. Kumawat
/// \date July 2016
//===----------------------------------------------------------------------===//
#ifndef TNudyEndfFissionYield_H
#define TNudyEndfFissionYield_H

#include "Geant/TNudyEndfRecoPoint.h"

namespace Nudy {
class TNudyEndfFile;
}

typedef std::vector<double> rowd;
typedef std::vector<rowd> matrixd2;
#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
#endif

namespace NudyPhysics {
class TNudyEndfFissionYield : public NudyPhysics::TNudyEndfRecoPoint {

public:
  TNudyEndfFissionYield();
  /// \brief Default constructure
  TNudyEndfFissionYield(Nudy::TNudyEndfFile *file);
  /// \brief constructure to be called in recopoint
  virtual double GetFisYield(int elemid, double energyK);
  /// \brief getting fission yield for the element and neutron energy
  virtual ~TNudyEndfFissionYield();

private:
  rowd fEin, fEinc;
  /// \brief incident energy
  matrixd2 fZafp, fFps, fZafpc, fFpsc, fYi, fCyi, fDyi, fYc, fDyc;
  /// \brief charge, mass, yield (independent and cummulative and error)
  rowd fZafp1, fFps1, fZafpc1, fFpsc1, fYi1, fCyi1, fDyi1, fYc1, fDyc1;
/// \brief charge, mass, yield (independent and cummulative)
#ifdef USE_ROOT
  TRandom3 *fRnd;
#endif
  ClassDef(TNudyEndfFissionYield, 1) // class for an ENDF fission yield reconstruction
};

} // namespace
#endif
