//-----------------------------------------------------------------
// @file TNudyInterface.h
// @brief prototype Nudy interface for GV
// @author Abhijit Bhattacharyya
//----------------------------------------------------------------
#ifndef TNUDY_INTERFACE_H
#define TNUDY_INTERFACE_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>

/*

//#include "Material.h"
//#include "Element.h"

//#include "SystemOfUnits.h"
//#include "PhysicalConstants.h"
// #include "HadronicProcess.h"
// #include "NudyProcess.h"
// #include "NudyCrossSection.h"

*/

#ifdef USE_ROOT
#include "Geant/TNudyDB.h"
#include "Geant/TNudyAlias.h"
#include "Geant/TNudyElementTable.h"
#include "Geant/TNudyENDF.h"
#include "Geant/TNudyEndfCont.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfINTG.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyEndfMat.h"
#include "Geant/TNudyEndfTape.h"
#include "Geant/TNudyEndfRecord.h"
#include "Geant/TNudyEndfSec.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyEndfTab2.h"
#include "Geant/TNudyLibrary.h"
#include "Geant/TNudyManager.h"
#include "Geant/TNudySubLibrary.h"
#include "Geant/TVNudyModel.h"
#include "Geant/TNudyEndfEnergy.h"
#include "Geant/TNudyEndfSigma.h"
#include "Geant/TNudyEndfRecoPoint.h"
#include "Geant/TNudyEndfTape.h"
#include "Geant/TNudyEndfAng.h"
//#include "NudyXSProcess.h"

#include "TSystem.h"
#include "TFile.h"
#include "TList.h"
#include "TKey.h"
class TList;
#endif

using Nudy::TNudyAlias;
using Nudy::TNudyCore;
using Nudy::TNudyDB;
using Nudy::TNudyElementTable;
using Nudy::TNudyEndfFile;
using Nudy::TNudyENDF;
using Nudy::TNudyEndfMat;
using Nudy::TNudyEndfList;
using Nudy::TNudyEndfSec;
using Nudy::TNudyEndfCont;
using Nudy::TNudyEndfRecord;
using Nudy::TNudyEndfTape;
using Nudy::TNudyEndfTab1;
using Nudy::TNudyEndfTab2;
using Nudy::TNudyLibrary;
using Nudy::TNudyManager;
using Nudy::TVNudyModel;
using Nudy::TNudyEndfINTG;

using NudyPhysics::TNudyEndfDoppler;
using NudyPhysics::TNudyEndfEnergyAng;
using NudyPhysics::TNudyEndfEnergy;
using NudyPhysics::TNudyEndfFissionYield;
using NudyPhysics::TNudyEndfNuPh;
using NudyPhysics::TNudyEndfPhAng;
using NudyPhysics::TNudyEndfPhEnergy;
using NudyPhysics::TNudyEndfPhProd;
using NudyPhysics::TNudyEndfPhYield;
using NudyPhysics::TNudyEndfRecoPoint;
using NudyPhysics::TNudyEndfSigma;
// using NudyPhysics::NudyXSProcess;

namespace geantphysics {
//  class LightTrack;
// class HadronicCrossSection;
// class HadronicProcess;
//  class HadronicProcessType;
//  class NudyCrossSection;
//  class NudyProcess;
enum class NudyProcessType {
  kNotDefined,
  kElastic,
  kFission,
  kInelastic,
  kCapture,
  kRadioActiveDecay,
  kQuasiElastic,
  kLeptoNuclear,
  kUserDefined
};

// class HadronicProcess;
//  inline namespace GEANT_IMPL_NAMESPACE {
//    class Isotope;
//    class Material;
//    class Element;
//  }
}

namespace NudyPhysics {

class TNudyInterface {
public:
  TNudyInterface();
  TNudyInterface(const int projCode, const double projKE, const double temp, const std::string isoName, const int tZ,
                 const int tA, geantphysics::NudyProcessType processType);
  virtual ~TNudyInterface();

public:
  double GetXS(int projCode, double projKE, double temp, std::string isoName, int tZ, int tA,
               geantphysics::NudyProcessType pType);
  std::string SetDataFileNameENDF(int projCode, std::string isoName); //, int tZ, int tA );
  std::string findENDFFileName(std::string ele);                      //, int tZ, int tA ) ;
  std::string GetDataFileName(std::string str1, std::string str2);
  std::string FixRootDataFile(std::string str1);           // ENDF filename without path and extension
  std::string SetDataFileNameROOT(std::string isoName);    //, int tZ, int tA );
  std::string SetDataFileNameENDFSUB(std::string isoName); //, int tZ, int tA );
  std::string GetCWD();
  bool GetFisCha(int inKey);
  void SetProjIDFn(int projCode, std::string style);
  double ComputeCrossSection();
  void SetMTValues(geantphysics::NudyProcessType pType);
  void setFileNames(std::string fIN, std::string fOUT, std::string fSubName);
  void DumpEndf2Root(std::string fIN, std::string fEndfSub, std::string fOUT, double temp);
  // void ConvertENDF2ROOT(std::string fENDF, std::string fEndfSub, std::string rENDF, double temp);
  void printE2RErr();

public:
  inline std::string GetIsotopeName();
  inline int GetProjectileCode();
  inline int GetZ();
  inline int GetA();
  inline double GetProjectileKE();
  inline double GetTemp();
  inline double GetCrossSection();
  inline bool GetIsFissKey();
  inline std::string GetProjID();
  // inline std::vector<double> GetXSTable();

  inline void SetIsotopeName(const std::string &isoName);
  inline void SetProjectileCode(const int projCode);
  inline void SetZ(const int tZValue);
  inline void SetA(const int tAvalue);
  inline void SetProjectileKE(const double projKE);
  inline void SetTemp(const double temp);
  inline void SetCrossSection(const double XSvalue);
  inline void SetEndfDataFileName(std::string fileName);
  inline void SetEndfSubDataFileName(std::string fileName);
  inline void SetRootFileName(std::string fileName);
  inline void SetIsFissKey(const bool theKey);
  inline void SetProjID(const std::string &theID);
  inline void AppendXS(const double xsvalue);
  inline void SetProcessType(const geantphysics::NudyProcessType ptype);

private:
  unsigned int fNumberOfReactionChannels = 895;
  geantphysics::NudyProcessType fProcType;
  std::string fIsoName;
  int fMTValue;
  int fProjCode;
  int ftZ;
  int ftA;
  int MTChargeFlag[10];
  double fProjKE;
  double fTemperature;
  double fXS;
  std::string fProjID;
  const char *fEndfFileN;
  const char *fEndfSubDataFileName;
  const char *fRootFileName;
  std::vector<int> fFissionFragmentsMass;
  std::vector<int> fFissionFragmentCharge;
  std::vector<double> fChannelXSArray;
  bool fIsFiss;
  std::vector<int> fChannelFiss{18, 19, 20, 21, 38, 151, 452, 454, 455, 456, 457, 458, 459};
};

//--------- GETTERS -------
inline std::string TNudyInterface::GetIsotopeName()
{
  return fIsoName;
}
inline int TNudyInterface::GetProjectileCode()
{
  return fProjCode;
}
inline int TNudyInterface::GetZ()
{
  return ftZ;
}
inline int TNudyInterface::GetA()
{
  return ftA;
}
inline double TNudyInterface::GetProjectileKE()
{
  return fProjKE;
}
inline double TNudyInterface::GetTemp()
{
  return fTemperature;
}
inline bool TNudyInterface::GetIsFissKey()
{
  return fIsFiss;
}
inline std::string TNudyInterface::GetProjID()
{
  return fProjID;
}
// inline std::vector<double> TNudyInterface::GetXSTable() { return fChannelXSArray; }

//--------- SETTERS ---------
inline void TNudyInterface::SetIsotopeName(const std::string &isoName)
{
  fIsoName = isoName;
}
inline void TNudyInterface::SetProjectileCode(const int projCode)
{
  fProjCode = projCode;
}
inline void TNudyInterface::SetZ(const int tZValue)
{
  ftZ = tZValue;
}
inline void TNudyInterface::SetA(const int tAvalue)
{
  ftA = tAvalue;
}
inline void TNudyInterface::SetProjectileKE(const double projKE)
{
  fProjKE = projKE;
}
inline void TNudyInterface::SetTemp(const double temp)
{
  fTemperature = temp;
}
inline void TNudyInterface::SetCrossSection(const double XSvalue)
{
  fXS = XSvalue;
}
inline void TNudyInterface::SetEndfDataFileName(std::string fileName)
{
  fEndfFileN = fileName.c_str();
}
inline void TNudyInterface::SetEndfSubDataFileName(std::string fileName)
{
  fEndfSubDataFileName = fileName.c_str();
}
inline void TNudyInterface::SetRootFileName(std::string fileName)
{
  fRootFileName = fileName.c_str();
}
inline void TNudyInterface::SetIsFissKey(const bool theKey)
{
  fIsFiss = theKey;
}
inline void TNudyInterface::SetProjID(const std::string &theID)
{
  fProjID = theID;
}
// inline void TNudyInterface::AppendXS ( const double xsvalue ) { fChannelFiss.push_back(xsvalue); }
inline void TNudyInterface::SetProcessType(const geantphysics::NudyProcessType ptype)
{
  fProcType = ptype;
}
} // namespace TNudyPhysics

#endif
