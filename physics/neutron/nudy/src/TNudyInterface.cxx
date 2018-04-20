/**********************************************
 * @file TNudyInterface.cxx
 * @author Abhijit Bhattacharyya
 * @brief This interface links to the GV through NUDYCrossSection of realphysics/physics/HAD/CrossSection
 **********************************************/
#include <iostream>
#include "Geant/TNudyInterface.h"

using namespace geantphysics;
using namespace Nudy;
using namespace NudyPhysics;

NudyPhysics::TNudyInterface::TNudyInterface() //:
    // fProjCode(2112), fProjKE(4.0e+6), fTemperature(293.60608), fIsoName("Fe"), ftZ(26), ftN(56),
    // fEndfDataFileName(""), fEndfSubDataFileName(""), fRootFileName("")
    {};

NudyPhysics::TNudyInterface::TNudyInterface(const int projCode, const double projKE, double temp,
                                            const std::string isoName, const int tZ, const int tA,
                                            geantphysics::NudyProcessType ProcessType)
{
  SetProjectileCode(projCode);
  SetProjectileKE(projKE);
  SetTemp(temp);
  SetIsotopeName(isoName);
  SetZ(tZ);
  SetA(tA);
  SetProcessType(ProcessType);
  SetMTValues(ProcessType);
};

NudyPhysics::TNudyInterface::~TNudyInterface()
{
}

void NudyPhysics::TNudyInterface::setFileNames(std::string fIN, std::string fOUT, std::string fSubName)
{
  TNudyInterface::fEndfFileN = fIN.c_str();
  if (!fOUT.length()) {
    fOUT = fIN + ".root";
  }
  fRootFileName        = fOUT.c_str();
  fEndfSubDataFileName = fSubName.c_str();
  std::cout << "endf: " << TNudyInterface::fEndfFileN << "     root: " << fRootFileName << std::endl;
}

////////////////////////
// void NudyPhysics::TNudyInterface::DumpEndf2Root(std::string fIN, std::string fOUT, std::string fSUBName, double temp)
// {
void NudyPhysics::TNudyInterface::DumpEndf2Root(std::string fIN, std::string fEndfSub, std::string fOUT, double temp)
{
  SetTemp(temp);
  SetIsFissKey(false);

  std::string fENDFD1;
  std::string fENDFD2;
  fENDFD1 = "n-" + fIN;

  fENDFD1 = GetDataFileName("neutrons", fENDFD1);
  fENDFD1 += ".endf";
  fEndfFileN = fENDFD1.c_str();

  fENDFD2 = "nfy-" + fEndfSub;
  fENDFD2 = GetDataFileName("nfy", fENDFD2);
  fENDFD2 += ".endf";
  fEndfSubDataFileName = fENDFD2.c_str();

  if (!fOUT.length()) fOUT = fIN + ".root";
  fRootFileName            = fOUT.c_str();

  Nudy::TNudyENDF *proc = new Nudy::TNudyENDF(fEndfFileN, fRootFileName, "recreate");
  proc->SetPreProcess(0);
  proc->SetLogLev(0);
  proc->Process();

  if (fEndfSub.find("ZZ") == std::string::npos) {
    bool LFIval = proc->GetLFI();
    SetIsFissKey(LFIval);
    proc->SetEndfSub(fEndfSubDataFileName);
    proc->Process();
  } else {
    // Nudy::TNudyENDF *tn = new Nudy::TNudyENDF(fEndfFileN, fRootFileName, "recreate");
    // tn->SetLogLev(2);
    // tn->Process();
    return;
  }

  double iSigDiff                   = 0.001; // documentation required
  NudyPhysics::TNudyEndfSigma *xsec = new NudyPhysics::TNudyEndfSigma(fRootFileName, iSigDiff);
  xsec->SetsigPrecision(iSigDiff);
  xsec->SetPreProcess(0);
  xsec->SetInitTempDop(0.0);
  xsec->SetOutTempDop(temp);
  xsec->GetData(fRootFileName, iSigDiff);
  // NudyPhysics::TNudyEndfRecoPoint *recoPoint = new NudyPhysics::TNudyEndfRecoPoint(0, fRootFileName);
}

////////////////////////////////////
double NudyPhysics::TNudyInterface::GetXS(int projCode, double projKE, double temp, std::string isoName, int tZ, int tA,
                                          geantphysics::NudyProcessType pType)
{
  SetProcessType(pType);
  SetMTValues(pType);
  SetA(tA);
  SetZ(tZ);
  SetProjectileKE(projKE);
  SetTemp(temp);
  SetIsotopeName(isoName);
  SetIsFissKey(false); // initializing to false first

  // Nudy::TNudyENDF *proc;

  //  Fix the name of the ENDF, ROOT and ENDFSUB filenames here
  std::string fileENDF1 = SetDataFileNameENDF(projCode, isoName);
  SetEndfDataFileName(fileENDF1);
  fileENDF1 += ".endf";
  std::string fileENDF2 = SetDataFileNameROOT(isoName);
  SetRootFileName(fileENDF2);
  std::string fileENDF3 = SetDataFileNameENDFSUB(isoName);
  SetEndfSubDataFileName(fileENDF3);
  fileENDF3 += ".endf";

  fEndfFileN           = fileENDF1.c_str();
  fRootFileName        = fileENDF2.c_str();
  fEndfSubDataFileName = fileENDF3.c_str();

  // Create and process with NUDY with keywords
  Nudy::TNudyENDF *procX = new Nudy::TNudyENDF(fEndfFileN, fRootFileName, "recreate");

  procX->SetPreProcess(0); // make comment if wrong

  procX->SetLogLev(0);
  procX->Process();

  bool LFIval = procX->GetLFI();
  SetIsFissKey(LFIval);
  procX->SetEndfSub(fEndfSubDataFileName);

  procX->Process();

  // SetProjIDFn(projCode, projKE,"");
  double XSvalue = ComputeCrossSection(); // call routines from Nudy
  return XSvalue;
}

void NudyPhysics::TNudyInterface::printE2RErr()
{
  std::cout << "SYNTAX: convertENDF2ROOT < ENDF file name > " << std::endl
            << "        convertENDF2ROOT < ENDF file name > < ROOT file name > " << std::endl
            << "......................... Exiting." << std::endl;
  exit(0);
}

std::string NudyPhysics::TNudyInterface::SetDataFileNameENDF(int projCode, std::string isoName)
{
  std::string fstyle = "ENDF";
  SetProjIDFn(projCode, fstyle);
  std::string DataFileNameString = findENDFFileName(isoName);
  std::string fileENDFName1      = GetDataFileName(fProjID, DataFileNameString);
  return fileENDFName1;
}

std::string NudyPhysics::TNudyInterface::SetDataFileNameROOT(std::string isoName)
{
  std::string DataFileNameString = findENDFFileName(isoName);
  std::string fileENDFName2      = FixRootDataFile(DataFileNameString);
  return fileENDFName2;
}

// This method is to be modified
void NudyPhysics::TNudyInterface::SetProjIDFn(int projCode, std::string fstyle)
{
  bool isChFiss = GetFisCha(fMTValue);
  //  if ( projKE < 0.03*geant::eV && projCode == 2112 ) SetProjID("thermal_scatt");
  if (projCode == 2112) SetProjID("neutrons");
  if (isChFiss && fstyle.compare("ENDFSUB") == 0) SetProjID("nfy");
}

std::string NudyPhysics::TNudyInterface::SetDataFileNameENDFSUB(std::string isoName)
{
  std::string fileENDFName3;
  std::string DataFileNameString;
  SetProjID("nfy");
  std::string prjId  = GetProjID();
  DataFileNameString = findENDFFileName(isoName);
  fileENDFName3      = GetDataFileName(prjId, DataFileNameString);
  return fileENDFName3;
}

// Actual Nudy CrossSection computation method
double NudyPhysics::TNudyInterface::ComputeCrossSection()
{
  int iElementID  = 0; //<------------- confusing testing by Abhijit 419 ?
  double xsvalue  = 0.0;
  double iSigDiff = 0.001; // trial value for test documentation reqd.

  NudyPhysics::TNudyEndfSigma *xsec = new TNudyEndfSigma(fRootFileName, iSigDiff);
  xsec->SetsigPrecision(iSigDiff);
  // iSigDiff = xsec->GetsigPrecision ();
  xsec->SetPreProcess(0);
  xsec->SetInitTempDop(0.0);
  xsec->SetOutTempDop(293.6);
  xsec->GetData(fRootFileName, iSigDiff);

  NudyPhysics::TNudyEndfRecoPoint *recoPoint = new TNudyEndfRecoPoint(iElementID, fRootFileName);
  // This is  under testing to check Harphool code for interfacing to GV :: Abhijit

  for (unsigned int crsp = 0; crsp < recoPoint->MtValues[iElementID].size(); crsp++) {
    int mtNow = recoPoint->MtValues[iElementID][crsp];
    if (mtNow == fMTValue) {
      xsvalue = recoPoint->GetSigmaPartial(iElementID, crsp, fProjKE);
      break;
    }
  }
  return xsvalue;
}

// selects the data file name for ENDF data
std::string NudyPhysics::TNudyInterface::GetDataFileName(std::string str1, std::string str2)
{
  std::string EndfDataPath = "";
  if (std::getenv("ENDFDATADIR") != NULL) {
    EndfDataPath               = std::getenv("ENDFDATADIR");
    std::string ENDFDataString = EndfDataPath + "/" + str1 + "/" + str2; // + ".endf";
    return ENDFDataString;
  } else {
    std::cout << " Please set environment ENDFDATADIR pointing to root of ENDF-B-VII data directory ." << std::endl;
    exit(99);
  }
  return EndfDataPath;
}

// selects name of root data file in the current working directory
std::string NudyPhysics::TNudyInterface::FixRootDataFile(std::string str1)
{
  std::string cwdPath      = GetCWD();
  std::string rootENDFFile = cwdPath + "/" + str1 + ".root";
  return rootENDFFile;
}

// store root data file in current working directory
// This is a bad technique. Using for the testing as quick solution.
std::string NudyPhysics::TNudyInterface::GetCWD()
{
  char *tempDIRarray  = new char[1024];
  std::string cwdPath = (getcwd(tempDIRarray, 1024)) ? std::string(tempDIRarray) : std::string("");
  delete[] tempDIRarray;
  return cwdPath;
}

std::string NudyPhysics::TNudyInterface::findENDFFileName(std::string elementName)
{
  std::stringstream ss;
  std::string fName = "";
  if (fProjID == "thermal_scatt") {
    fName = "tsl-";
  } else if (fProjID == "neutrons") {
    fName = "n-";
  } else {
    fName = fProjID + "-";
  }

  ss << ftZ;
  std::string stZ = ss.str();
  ss.str("");
  ss << ftA;
  std::string stA = ss.str();
  ss.str("");

  switch (stZ.length()) {
  case 0:
    return "";
  case 1:
    stZ = "00" + stZ;
    break;
  case 2:
    stZ = "0" + stZ;
    break;
  case 3:
    stZ = stZ;
  }

  if (ftZ == 12) {
    stA = "000";
  } else {
    switch (stA.length()) {
    case 0:
      return "";
    case 1:
      stA = "00" + stA;
      break;
    case 2:
      stA = "0" + stA;
      break;
    case 3:
      stA = stA;
    }
  }
  fName = fName + stZ + "_" + elementName + "_" + stA;
  return fName;
}

// Depending on the process type set the MT value for cross section
void NudyPhysics::TNudyInterface::SetMTValues(geantphysics::NudyProcessType pType)
{
  SetProcessType(pType);
  switch (pType) {
  case geantphysics::NudyProcessType::kElastic:
    fMTValue = 2;
    break;
  case geantphysics::NudyProcessType::kInelastic:
    fMTValue = 3;
    break;
  case geantphysics::NudyProcessType::kFission:
    fMTValue = 18;
    break;
  case geantphysics::NudyProcessType::kCapture:
    fMTValue = 102;
    break;
  case geantphysics::NudyProcessType::kRadioActiveDecay:
    fMTValue = 457;
    break;
  case geantphysics::NudyProcessType::kNotDefined:
  case geantphysics::NudyProcessType::kQuasiElastic:
  case geantphysics::NudyProcessType::kUserDefined:
  case geantphysics::NudyProcessType::kLeptoNuclear:
    fMTValue = 2; // to be changed later
  }
}

bool NudyPhysics::TNudyInterface::GetFisCha(int inKey)
{
  auto p = std::find(fChannelFiss.begin(), fChannelFiss.end(), inKey);
  return (p != fChannelFiss.end());
}
