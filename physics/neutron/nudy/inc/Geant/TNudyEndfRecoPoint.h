//===-- Nudy/TNudyEndfRecoPoint.h - Instruction class definition -------*- C++ -*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \brief The class reads doppler broadened cross-sections from the root file
///  and gives all the corss-sections/ secondary parameters forthe outgoing particles
/// \class TNudyEndfRecoPoint
/// \author H. Kumawat
/// \date March 2016
//===----------------------------------------------------------------------===//
#ifndef __TNudyEndfRecoPoint__
#define __TNudyEndfRecoPoint__

#include <vector>
#include <fstream>
namespace Nudy {
class TNudyEndfFile;
class TNudyEndfList;
}

class TList;

namespace NudyPhysics {
class TNudyEndfNuPh;
class TNudyEndfFissionYield;
class TNudyEndfEnergy;
class TNudyEndfEnergyAng;
class TNudyEndfAng;
class TNudyEndfPhYield;
class TNudyEndfPhProd;
class TNudyEndfPhAng;
class TNudyEndfPhEnergy;
}

#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
#endif

namespace NudyPhysics {

#define PI acos(-1.0)
typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
typedef std::vector<rowint> matrixint;
typedef std::vector<rowd> matrixd2;
typedef std::vector<std::vector<rowd>> matrixd3;
typedef std::vector<std::vector<std::vector<rowd>>> matrixd4;
typedef std::vector<std::vector<std::vector<std::vector<rowd>>>> matrixd5;

class TNudyEndfRecoPoint {

public:
  TNudyEndfRecoPoint();
  /// \brief Default constructure
  TNudyEndfRecoPoint(int ielemId, const char *irENDF);
  /// \brief constructure to be called for any interface class to get data
  virtual ~TNudyEndfRecoPoint();
  void GetData(int elemid, const char *irENDF);
  /// \brief main function to process the data
  double GetSigmaTotal(int elemid, double energyK);
  /// \brief getting total XSec.
  double GetSigmaPartial(int elemid, int i, double energyK);
  /// \brief getting partial cross-section
  virtual double GetCos4(int elemid, int mt, double energyK);
  /// \brief getting cosine from file 4
  virtual double GetCos64(int elemid, int mt, double energyK);
  /// \brief getting cosine from file 6 of type 4
  virtual int GetCos4Lct(int elemid, int mt);
  /// \brief getting CM vs. LAB flag
  virtual double GetEnergy5(int elemid, int mt, double energyK);
  /// \brief getting secondary neutron energy from file 5
  virtual double GetCos6(int elemid, int mt, double energyK);
  /// \brief getting cosine from file 6
  virtual double GetEnergy6(int elemid, int mt, double energyK);
  /// \brief getting energy from file 6
  virtual double GetQValue(int elemid, int mt);
  /// \brief getting q value for the reaction specially for the excited state emission
  virtual double GetMt4(int elemid, int mt);
  /// \brief getting MT values for file 4 for which angular distributions are given
  virtual double GetMt5(int elemid, int mt);
  /// \brief getting MT values for which energy distributions are given in file 5
  virtual double GetMt6(int elemid, int mt);
  /// \brief getting MT values for which angle-energy correlated distributions are given in file 6
  virtual int GetLaw6(int ielemId, int mt);
  /// \brief getting LAW of the energy distribution in file 6
  virtual int GetZd6(int ielemId, int mt);
  /// \brief getting Z value
  virtual int GetAd6(int ielemId, int mt);
  /// \brief getting A value
  virtual int GetMt6Neutron(int ielemId, int mt);
  /// \brief getting MT values for neutron only
  virtual double GetNuTotal(int elemid, double energyK);
  /// \brief getting total fission neutron multiplicity
  virtual double GetNuPrompt(int elemid, double energyK);
  /// \brief getting prompt fission neutron multiplicity
  virtual double GetNuDelayed(int elemid, double energyK);
  /// \brief getting delayed fission neutron multiplicity
  virtual double GetFissHeat(int elemid, double energyK);
  /// \brief getting total fission Heat
  virtual double GetFisYield(int elemid, double energyK);
  /// \brief getting fission yield
  virtual double GetLambdaD(int elemid, int time);
  /// \brief getting lambda of delayed neutron emitter family
  virtual double GetDelayedFraction(int ielemId, int mt, double energyK);
  /// \brief getting delayed fraction for each group of the emitter family
  std::fstream out, outtotal;
  std::string outstring, outstringTotal;
  matrixint MtValues;
  /// \brief MT values for which cross-section/ heating values are given

protected:
  int fElemId;
  /// \brief serial no of the element in the material list of the simulation
  const char *rENDF;
  /// \brief root file name
  matrixd2 fEneUni, fSigUniT;
  /// \brief unionization of energy and total cross-section
  matrixd3 fSigUniOfMt;
  /// \brief Xsec for each reaction after unionization of energy for all reactions
  matrixint fEnergyLocMtId;
  /// \brief MT wise starting energy locations for all cross-sections
  matrixd4 fCos4OfMts;
  /// \brief cosine from file 4 for each reaction
  matrixd4 fCosPdf4OfMts;
  /// \brief pdf from file 4 for each reaction
  matrixd4 fCosCdf4OfMts;
  /// \brief cdf from file 4 for each reaction
  matrixd3 fEnergy4OfMts;
  /// \brief incident energy in file 4 for each reaction
  matrixint fMt4Values;
  ///  \brief MT values for which angular distributions are given in file 4
  matrixint fMt4Lct;
  /// \brief CM and Lab flag for angular distributions as given in file 4
  matrixd4 fEnergyOut5OfMts;
  /// \brief energy from file 5 for each reaction
  matrixd4 fEnergyPdf5OfMts;
  /// \brief pdf from file 5 for each reaction
  matrixd4 fEnergyCdf5OfMts;
  /// \brief cdf from file 5 for each reaction
  matrixd4 fCos6OfMts;
  /// \brief cosine 6 for each reaction
  matrixd4 fCosin6Pdf, fCosin6Cdf;
  /// \brief pdf cdf file 6 for each reaction
  matrixd5 fEnergyOut6OfMts;
  /// \brief energy from file 6 for each reaction
  matrixd5 fEnergyPdf6OfMts;
  /// \brief pdf from file 6 for each reaction
  matrixd5 fEnergyCdf6OfMts;
  /// \brief cdf from file 6 for each reaction
  matrixd3 fEnergy5OfMts;
  /// \brief incident energy in file 5 for each reaction
  matrixd3 fFraction5OfMts;
  /// \brief fraction for incident energy in file 5 for each reaction
  matrixint fMt5Values;
  /// \brief MT values for which energy distributions are given in file 5
  matrixint fLaw6;
  /// \brief law 6 for angular-energy distributions are given in file 6
  matrixint fZD6, fAD6;
  /// \brief Z, A of law 6 distributions
  matrixd2 fEint, fNut;
  /// \brief total incident energy and nu
  matrixd2 fEinp, fNup;
  /// \brief prompt incident energy and nu
  matrixd2 fEind, fNud, fLambdaD;
  /// \brief delayed incident energy, nu and lambda
  matrixd2 fEinFissHeat, fFissHeat;
  /// \brief fission incident energy and heat
  matrixd2 fEinfId, fQvalue;
  /// \brief incident energy for fission yield and q-value
  matrixd3 fZafId, fPdfYieldId, fCdfYieldId;
  /// \brief za and yield fission
  double AWRI;
  /// \brief Mass in units of the neutron

private:
  void ReadFile2(Nudy::TNudyEndfFile *file);
  /// \brief function to read file 2
  void ReadFile3(Nudy::TNudyEndfFile *file);
  /// \brief function to read file 3
  void FixupTotal(rowd &x1);
  /// \brief making total cross-section

  int fFlagRead = -1;
  /// \brief flag for reading charge particle and photon production cross-sections
  double fQValue[999];
  /// \brief q-value for all the reactions
  int fNR, fNP;
  /// \brief standard ENDF parameters for range and interpolation
  matrixint fMt4, fMt5, fMt6;
  /// \brief MT values for which angular, energy/ angular-energy distributions are given in file 4, 5, 6
  matrixd2 fSigmaOfMts;
  /// \brief sigma for each reaction
  matrixd2 fSigmaUniOfMts;
  /// \brief sigma for each reaction after unionization of energy
  rowint fEnergyLocationMts;
  /// \breif MT wise starting energy for cross-section
  rowint fMtNumbers, fMtNum4, fMtNum5, fMtNum6;
  /// \breif MT numbers
  rowd fEnergyMts, fSigmaMts, fQvalueTemp;
  /// \brief MT numbers for sigma in file3
  rowd fELinearFile3, fXLinearFile3;
  /// \breif energy and Xsec
  rowd fEneTemp, fSigTemp;
  /// \breif temporary vectors to store energy and sigma
  NudyPhysics::TNudyEndfAng *recoAng;
  NudyPhysics::TNudyEndfEnergy *recoEnergy;
  NudyPhysics::TNudyEndfEnergyAng *recoEnergyAng;
  NudyPhysics::TNudyEndfNuPh *recoNuPh;
  NudyPhysics::TNudyEndfFissionYield *recoFissY;
  NudyPhysics::TNudyEndfPhYield *recoPhYield;
  NudyPhysics::TNudyEndfPhProd *recoPhProd;
  NudyPhysics::TNudyEndfPhAng *recoPhAng;
  NudyPhysics::TNudyEndfPhEnergy *recoPhEnergy;
#ifdef USE_ROOT
  TRandom3 *fRnd;
#endif

#ifdef USE_ROOT
  ClassDef(TNudyEndfRecoPoint, 1) // class for an ENDF reconstruction
#endif
};
}
#endif
