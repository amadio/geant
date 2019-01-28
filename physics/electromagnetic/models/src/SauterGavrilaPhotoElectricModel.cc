#include "Geant/SauterGavrilaPhotoElectricModel.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/Material.h"
#include "Geant/Element.h"
#include "Geant/MaterialProperties.h"

#include "Geant/MaterialCuts.h"

#include "Geant/Spline.h"
#include "Geant/GLIntegral.h"
#include "Geant/AliasTable.h"
#include "Geant/XSectionsVector.h"

#include "Geant/PhysicsParameters.h"
#include "Geant/Gamma.h"
#include "Geant/Electron.h"
#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include "Geant/math_wrappers.h"

// from geantV
#include "Geant/Typedefs.h"
#include "Geant/TaskData.h"

using namespace std;
namespace geantphysics {

using geant::Double_v;
using geant::IndexD_v;
using geant::kVecLenD;
using geant::MaskD_v;
using MaskDI_v = vecCore::Mask<IndexD_v>;
using vecCore::AssignMaskLane;
using vecCore::Get;
using vecCore::MaskEmpty;
using vecCore::MaskFull;
using vecCore::Set;

std::vector<double> *SauterGavrilaPhotoElectricModel::fParamHigh[] = {nullptr};
std::vector<double> *SauterGavrilaPhotoElectricModel::fParamLow[]  = {nullptr};
std::vector<double> SauterGavrilaPhotoElectricModel::fBindingEn[];

int SauterGavrilaPhotoElectricModel::fNShells[]          = {0};
int SauterGavrilaPhotoElectricModel::fNShellsUsed[]      = {0};
int SauterGavrilaPhotoElectricModel::fLastSSAliasIndex[] = {0};

Material *SauterGavrilaPhotoElectricModel::fWater         = nullptr;
double SauterGavrilaPhotoElectricModel::fWaterEnergyLimit = 0.0;

XSectionsVector **SauterGavrilaPhotoElectricModel::fShellVectorFull[] = {nullptr};
XSectionsVector **SauterGavrilaPhotoElectricModel::fShellVector[]     = {nullptr};
XSectionsVector *SauterGavrilaPhotoElectricModel::fLECSVector[]       = {nullptr};
XSectionsVector *SauterGavrilaPhotoElectricModel::fCSVector[]         = {nullptr};
bool *SauterGavrilaPhotoElectricModel::fCrossSection                  = nullptr;
bool *SauterGavrilaPhotoElectricModel::fCrossSectionLE                = nullptr;

SauterGavrilaPhotoElectricModel::SauterGavrilaPhotoElectricModel(const std::string &modelname, bool aliasActive)
    : EMModel(modelname)
{

  SetUseSamplingTables(aliasActive);
  fMinPrimEnergy =
      1.e-12 * geant::units::eV; // Minimum of the gamma kinetic energy grid, used to sample the photoelectron direction
  fMaxPrimEnergy = 100 * geant::units::MeV; // Maximum of the gamma kinetic energy grid (after this threshold the e- is
                                            // considered to follow the same direction as the incident gamma)
  fShellMinPrimEnergy    = 1 * geant::units::eV;
  fShellMaxPrimEnergy    = 1 * geant::units::GeV;
  fPrimEnLMin            = 0.;      // will be set in InitSamplingTables if needed
  fPrimEnILDelta         = 0.;      // will be set in InitSamplingTables if needed
  fSamplingPrimEnergies  = nullptr; // will be set in InitSamplingTables if needed
  fLSamplingPrimEnergies = nullptr; // will be set in InitSamplingTables if needed

  fShellSamplingPrimEnergies  = nullptr; // will be set in InitShellSamplingTables if needed
  fShellLSamplingPrimEnergies = nullptr; // will be set in InitShellSamplingTables if needed

  fAliasData         = nullptr; // will be set in InitSamplingTables if needed
  fAliasSampler      = nullptr;
  fShellAliasData    = nullptr; // will be set in InitSamplingTables if needed
  fShellAliasSampler = nullptr;

  fIsBasketizable = true;
}

SauterGavrilaPhotoElectricModel::~SauterGavrilaPhotoElectricModel()
{
  // CLEANING fParamHigh and fParamLow
  for (int i = 0; i < gMaxSizeData; ++i) {
    delete fParamHigh[i];
    delete fParamLow[i];

    fParamHigh[i] = 0;
    fParamLow[i]  = 0;
    fBindingEn[i].clear();
    fSortedBindingEn[i].clear();
    fSortedDoubledBindingEn[i].clear();
    fIndexBaseEn[i].clear();
    fIndexSortedDoubledBindingEn[i].clear();
    fShellSamplingPrimEnergiesNEW[i].clear();
    fShellLSamplingPrimEnergiesNEW[i].clear();
  }

  if (fSamplingPrimEnergies) delete[] fSamplingPrimEnergies;
  if (fLSamplingPrimEnergies) delete[] fLSamplingPrimEnergies;

  if (fShellSamplingPrimEnergies) delete[] fShellSamplingPrimEnergies;
  if (fShellLSamplingPrimEnergies) delete[] fShellLSamplingPrimEnergies;

  if (fAliasData) {
    for (int i = 0; i < fNumSamplingPrimEnergies; ++i) {
      if (fAliasData[i]) {
        delete[] fAliasData[i]->fXdata;
        delete[] fAliasData[i]->fYdata;
        delete[] fAliasData[i]->fAliasW;
        delete[] fAliasData[i]->fAliasIndx;
        delete fAliasData[i];
      }
    }
    delete[] fAliasData;
  }

  if (fShellAliasData) {
    for (int i = 0; i < fNumAliasTables; ++i) {
      if (fShellAliasData[i]) {
        delete[] fShellAliasData[i]->fXdata;
        delete[] fShellAliasData[i]->fYdata;
        delete[] fShellAliasData[i]->fAliasW;
        delete[] fShellAliasData[i]->fAliasIndx;
        delete fShellAliasData[i];
      }
    }
    delete[] fShellAliasData;
  }

  if (fAliasSampler) delete fAliasSampler;
  if (fShellAliasSampler) delete fShellAliasSampler;
}

void SauterGavrilaPhotoElectricModel::Initialize()
{
  EMModel::Initialize();
  fSecondaryInternalCode = Electron::Definition()->GetInternalCode();
  InitializeModel();
}

void SauterGavrilaPhotoElectricModel::InitializeModel()
{
  if (!fCrossSection) {
    fCrossSection = nullptr;
    fCrossSection = new bool[gMaxSizeData];
    for (int i = 0; i < gMaxSizeData; ++i) {
      fCrossSection[i] = false;
    }
  }
  if (!fCrossSectionLE) {
    fCrossSectionLE = nullptr;
    fCrossSectionLE = new bool[gMaxSizeData];
    for (int i = 0; i < gMaxSizeData; ++i) {
      fCrossSectionLE[i] = false;
    }
  }

  for (int i = 0; i < gMaxSizeData; ++i) {
    fSortedBindingEn[i].clear();
    fSortedDoubledBindingEn[i].clear();
    fIndexBaseEn[i].clear();
    fIndexSortedDoubledBindingEn[i].clear();
    fShellSamplingPrimEnergiesNEW[i].clear();
    fShellLSamplingPrimEnergiesNEW[i].clear();
  }

  fVerboseLevel = 1;
  // Uncomment the following 2 lines to run tests:
  //(1)PhysVecSauterGavrilaAliasShellValid
  //(2)PhysVecSauterGavrilaAliasShellBench
  //(3)PhysVecSauterGavrilaRejShellValid
  //(4)PhysVecSauterGavrilaRejShellBench
  //   for(int i=3; i<gMaxSizeData; i++)
  //   ReadData(i);
  LoadData();
  if (GetUseSamplingTables()) {
    InitSamplingTables();
    InitShellSamplingTables();
  }
}

void SauterGavrilaPhotoElectricModel::SetVerboseLevel(int lev)
{
  fVerboseLevel = lev;
}

void SauterGavrilaPhotoElectricModel::LoadData()
{

  int numMatCuts = MaterialCuts::GetTheMaterialCutsTable().size();
  // get list of active region
  std::vector<bool> isActiveInRegion = GetListActiveRegions();
  for (int i = 0; i < numMatCuts; ++i) {
    const MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[i];
    // if this MaterialCuts belongs to a region where this model is active:
    if (isActiveInRegion[matCut->GetRegionIndex()]) {
      // get the list of elements
      const Vector_t<Element *> &theElements = matCut->GetMaterial()->GetElementVector();
      int numElems                           = theElements.size();
      for (int j = 0; j < numElems; ++j) {
        double zet      = theElements[j]->GetZ();
        int elementIndx = std::lrint(zet);
        ReadData(elementIndx);
      }
    }
  }
}

void SauterGavrilaPhotoElectricModel::ReadData(int Z)
{
  using geant::units::barn;
  using geant::units::MeV;
  if (fVerboseLevel > 1) {
    std::cout << "Calling ReadData() of SauterGavrilaPhotoElectricModel" << std::endl;
  }
  if ((fCrossSection[Z]) && ((fCrossSectionLE[Z] && Z > 2) || (!fCrossSectionLE[Z] && Z < 3))) {
    // std::cout<<"Data of "<<Z<<" loaded before!\n";
    return;
  }

  // get the path to the main physics data directory
  char *path = std::getenv("GEANT_PHYSICS_DATA");
  if (!path) {
    std::cerr << "******   ERROR in SauterGavrilaPhotoElectricModel::ReadData() \n"
              << "         GEANT_PHYSICS_DATA is not defined! Set the GEANT_PHYSICS_DATA\n"
              << "         environment variable to the location of Geant data directory!\n"
              << std::endl;
    exit(1);
  }

  std::ostringstream ost;
  ost << path << "/livermore/phot_epics2014/pe-cs-" << Z << ".dat";
  std::ifstream fin(ost.str().c_str());
  if (!fin.is_open()) {
    std::cerr << "SauterGavrilaPhotoElectricModel data file <" << ost.str().c_str() << "> is not opened!" << std::endl;

    return;
  } else {

    fCrossSection[Z] = true;
    if (fVerboseLevel > 2) {
      std::cout << "File " << ost.str().c_str() << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
    }

    fCSVector[Z] = new XSectionsVector;
    fin >> fCSVector[Z]->fEdgeMin >> fCSVector[Z]->fEdgeMax >> fCSVector[Z]->fNumberOfNodes;

    int siz = 0;
    fin >> siz;
    if (fin.fail() || siz <= 0) {
      std::cerr << "SauterGavrilaPhotoElectricModel data file <" << ost.str().c_str() << "> is not opened!"
                << std::endl;
    }

    fCSVector[Z]->fBinVector.reserve(siz);
    fCSVector[Z]->fDataVector.reserve(siz);

    double vBin, vData;

    for (int i = 0; i < siz; i++) {
      vBin  = 0.;
      vData = 0.;
      fin >> vBin >> vData;

      if (fin.fail()) {
        std::cerr << "SauterGavrilaPhotoElectricModel data file <" << ost.str().c_str() << "> is not opened!"
                  << std::endl;
      }

      fCSVector[Z]->fBinVector.push_back(vBin * MeV);
      fCSVector[Z]->fDataVector.push_back(vData * barn);
    }

    // to remove any inconsistency
    fCSVector[Z]->fNumberOfNodes = siz;
    fCSVector[Z]->fEdgeMin       = fCSVector[Z]->fBinVector[0];
    fCSVector[Z]->fEdgeMax       = fCSVector[Z]->fBinVector[fCSVector[Z]->fNumberOfNodes - 1];

    // Use spline interpolator for Cross-sections vector
    fCSVector[Z]->fSplineInt =
        new Spline((fCSVector[Z]->fBinVector.data()), (fCSVector[Z]->fDataVector.data()), fCSVector[Z]->fNumberOfNodes);
    fin.close();
  }

  // read high-energy fit parameters
  fParamHigh[Z] = new std::vector<double>;

  int n1 = 0;
  int n2 = 0;
  double x;
  std::ostringstream ost1;
  ost1 << path << "/livermore/phot_epics2014/pe-high-" << Z << ".dat";
  std::ifstream fin1(ost1.str().c_str());
  if (!fin1.is_open()) {
    std::cerr << "SauterGavrilaPhotoElectricModel data file <" << ost1.str().c_str() << "> is not opened!" << std::endl;
    return;
  } else {
    if (fVerboseLevel > 2) {
      std::cout << "File " << ost1.str().c_str() << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
    }
    fin1 >> n1;
    if (fin1.fail()) {
      return;
    }
    if (0 > n1 || n1 >= INT_MAX) {
      n1 = 0;
    }

    fin1 >> n2;
    if (fin1.fail()) {
      return;
    }
    if (0 > n2 || n2 >= INT_MAX) {
      n2 = 0;
    }

    fin1 >> x;
    if (fin1.fail()) {
      return;
    }

    fNShells[Z] = n1;
    fParamHigh[Z]->reserve(7 * n1 + 1);
    fParamHigh[Z]->push_back(x * MeV);
    for (int i = 0; i < n1; ++i) {
      for (int j = 0; j < 7; ++j) {
        fin1 >> x;
        if (0 == j) {
          x *= MeV;
        } else {
          x *= barn;
        }
        fParamHigh[Z]->push_back(x);
      }
    }
    fin1.close();
  }

  // read low-energy fit parameters
  fParamLow[Z] = new std::vector<double>;
  int n1_low   = 0;
  int n2_low   = 0;
  double x_low;
  std::ostringstream ost1_low;
  ost1_low << path << "/livermore/phot_epics2014/pe-low-" << Z << ".dat";
  std::ifstream fin1_low(ost1_low.str().c_str());
  if (!fin1_low.is_open()) {
    std::cerr << "SauterGavrilaPhotoElectricModel data file <" << ost1_low.str().c_str() << "> is not opened!"
              << std::endl;
    return;
  } else {
    if (fVerboseLevel > 2) {
      std::cout << "File " << ost1_low.str().c_str() << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
    }
    fin1_low >> n1_low;
    if (fin1_low.fail()) {
      return;
    }
    if (0 > n1_low || n1_low >= INT_MAX) {
      n1_low = 0;
    }

    fin1_low >> n2_low;
    if (fin1_low.fail()) {
      return;
    }
    if (0 > n2_low || n2_low >= INT_MAX) {
      n2_low = 0;
    }

    fin1_low >> x_low;
    if (fin1_low.fail()) {
      return;
    }

    fNShells[Z] = n1_low;
    fParamLow[Z]->reserve(7 * n1_low + 1);
    fParamLow[Z]->push_back(x_low * MeV);
    for (int i = 0; i < n1_low; ++i) {
      for (int j = 0; j < 7; ++j) {
        fin1_low >> x_low;
        if (0 == j) {
          x_low *= MeV;
          fBindingEn[Z].push_back(x_low);
        } else {
          x_low *= barn;
        }
        fParamLow[Z]->push_back(x_low);
      }
    }
    fin1_low.close();
  }

  // there is a possibility to use only main shells
  if (gNShellLimit < n2) {
    n2 = gNShellLimit;
  }
  fNShellsUsed[Z]     = n2;
  fShellVector[Z]     = new XSectionsVector *[n2];
  fShellVectorFull[Z] = new XSectionsVector *[n2];
  for (int i = 0; i < n2; i++) {
    fShellVector[Z][i]     = new XSectionsVector;
    fShellVectorFull[Z][i] = new XSectionsVector;
  }
  for (int i = 0; i < n2; i++) {
    fShellVector[Z][i]->fDataVector.clear();
    fShellVector[Z][i]->fBinVector.clear();
    fShellVector[Z][i]->fNumberOfNodes = 0;
    fShellVectorFull[Z][i]->fDataVector.clear();
    fShellVectorFull[Z][i]->fBinVector.clear();
    fShellVectorFull[Z][i]->fNumberOfNodes = 0;
  }

  fNShellsUsed[Z] = n2;

  // If more than one shell is used -> Read sub-shells cross section data
  if (1 < n2) {
    std::ostringstream ost2;
    ost2 << path << "/livermore/phot_epics2014/pe-ss-cs-" << Z << ".dat";
    std::ifstream fin2(ost2.str().c_str());
    if (!fin2.is_open()) {
      std::cerr << "SauterGavrilaPhotoElectricModel data file <" << ost2.str().c_str() << "> is not opened!"
                << std::endl;
      return;
    } else {
      if (fVerboseLevel > 2) {
        std::cout << "File " << ost2.str().c_str() << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
      }

      int n3, n4;
      double y;

      for (int i = 0; i < n2; ++i) {
        fin2 >> x >> y >> n3 >> n4;

        fShellVector[Z][i]->fBinVector.clear();
        fShellVector[Z][i]->fDataVector.clear();
        fShellVector[Z][i]->fBinVector.reserve(n3);
        fShellVector[Z][i]->fDataVector.reserve(n3);
        fShellVector[Z][i]->fEdgeMin = x * MeV;
        fShellVector[Z][i]->fEdgeMax = y * MeV;
        if (fVerboseLevel > 3)
          std::cout << "fShellVector[" << Z << "][" << i << "]->fEdgeMin: " << fShellVector[Z][i]->fEdgeMin
                    << "\t fShellVector[" << Z << "][" << i << "]->fEdgeMax: " << fShellVector[Z][i]->fEdgeMax
                    << std::endl;

        for (int j = 0; j < n3; ++j) {
          fin2 >> x >> y;
          fShellVector[Z][i]->fBinVector.push_back(x * MeV);
          fShellVector[Z][i]->fDataVector.push_back(y * barn);
          fShellVector[Z][i]->fNumberOfNodes++;
        }
        fShellVector[Z][i]->fCompID = n4;
      }

      fin2.close();
    }
    // READ THE DENSER EPICS 2014 subshell cross-sections -> used to Sample the shells
    std::ostringstream ost2_full;
    ost2_full << path << "/livermore/phot_epics2014/pe-denser-ss-cs-" << Z << ".dat";
    std::ifstream fin2_full(ost2_full.str().c_str());

    if (!fin2_full.is_open()) {
      std::cerr << "SauterGavrilaPhotoElectricModel data file <" << ost2_full.str().c_str() << "> is not opened!"
                << std::endl;
      return;
    } else {
      if (fVerboseLevel > 2) {
        std::cout << "File " << ost2_full.str().c_str() << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
      }

      int n3, n4;
      double y;

      for (int i = 0; i < n2; ++i) {
        fin2_full >> x >> y >> n3 >> n4;
        fShellVectorFull[Z][i]->fBinVector.clear();
        fShellVectorFull[Z][i]->fDataVector.clear();
        fShellVectorFull[Z][i]->fBinVector.reserve(n3);
        fShellVectorFull[Z][i]->fDataVector.reserve(n3);
        fShellVectorFull[Z][i]->fEdgeMin = x * MeV;
        fShellVectorFull[Z][i]->fEdgeMax = y * MeV;
        if (fVerboseLevel > 3)
          std::cout << "fShellVectorFull[" << Z << "][" << i << "]->fEdgeMin: " << fShellVectorFull[Z][i]->fEdgeMin
                    << "\t fShellVectorFull[" << Z << "][" << i << "]->fEdgeMax: " << fShellVectorFull[Z][i]->fEdgeMax
                    << std::endl;

        for (int j = 0; j < n3; ++j) {
          fin2_full >> x >> y;
          // std::cout<<j<<"-th element:  \t"<<x<<"\t"<<y<<std::endl;;
          fShellVectorFull[Z][i]->fBinVector.push_back(x * MeV);
          fShellVectorFull[Z][i]->fDataVector.push_back(y * barn);
          fShellVectorFull[Z][i]->fNumberOfNodes++;
        }
        fShellVectorFull[Z][i]->fCompID = n4;
      }

      fin2_full.close();
    }
  }

  // no spline for photoeffect total x-section below K-shell
  if (1 < fNShells[Z]) {

    std::ostringstream ost3;
    ost3 << path << "/livermore/phot_epics2014/pe-le-cs-" << Z << ".dat";
    std::ifstream fin3(ost3.str().c_str());
    if (!fin3.is_open()) {
      std::cerr << "SauterGavrilaPhotoElectricModel data file <" << ost3.str().c_str() << "> is not opened!"
                << std::endl;
      return;
    } else {
      if (fVerboseLevel > 2) {
        std::cout << "File " << ost3.str().c_str() << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
      }

      fCrossSectionLE[Z] = true;
      fLECSVector[Z]     = new XSectionsVector;

      fin3 >> fLECSVector[Z]->fEdgeMin >> fLECSVector[Z]->fEdgeMax >> fLECSVector[Z]->fNumberOfNodes;
      int siz = 0;
      fin3 >> siz;
      if (fin3.fail() || siz <= 0) {
        std::cerr << "SauterGavrilaPhotoElectricModel data file <" << ost3.str().c_str() << "> is not opened!"
                  << std::endl;
      }

      fLECSVector[Z]->fBinVector.reserve(siz);
      fLECSVector[Z]->fDataVector.reserve(siz);

      double vBin, vData;

      for (int i = 0; i < siz; i++) {
        vBin  = 0.;
        vData = 0.;
        fin3 >> vBin >> vData;

        // Scale vector
        vBin *= MeV;
        vData *= barn;

        if (fin3.fail()) {
          std::cerr << "SauterGavrilaPhotoElectricModel data file <" << ost3.str().c_str() << "> is not opened!"
                    << std::endl;
        }

        fLECSVector[Z]->fBinVector.push_back(vBin);
        fLECSVector[Z]->fDataVector.push_back(vData);
      }

      // to remove any inconsistency
      fLECSVector[Z]->fNumberOfNodes = siz;
      fLECSVector[Z]->fEdgeMin       = fLECSVector[Z]->fBinVector[0];
      fLECSVector[Z]->fEdgeMax       = fLECSVector[Z]->fBinVector[fLECSVector[Z]->fNumberOfNodes - 1];

      fin3.close();
    }
    // std::cout<<"pe-le-cs- cross sections for ["<<Z<<"], loaded\n";
  }
}

//____________________
// NB: cosTheta is supposed to contain the dirZ of the incoming photon
void SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirection_Rejection(double gammaEnIn, double &cosTheta,
                                                                             geant::TaskData *td) const
{

  // 1) initialize energy-dependent variables
  // Variable naming according to Eq. (2.24) of Penelope Manual
  // (pag. 44)
  double gamma  = 1.0 + gammaEnIn / geant::units::kElectronMassC2;
  double gamma2 = gamma * gamma;
  double beta   = std::sqrt((gamma2 - 1.0) / gamma2);

  // ac corresponds to "A" of Eq. (2.31)
  //
  double ac = (1.0 / beta) - 1.0;
  double a1 = 0.5 * beta * gamma * (gamma - 1.0) * (gamma - 2.0);
  double a2 = ac + 2.0;
  // gtmax = maximum of the rejection function according to Eq. (2.28), obtained for tsam=0
  double gtmax = 2.0 * (a1 + 1.0 / ac);

  double tsam = 0;
  double gtr  = 0;

  // 2) sampling. Eq. (2.31) of Penelope Manual
  // tsam = 1-Math::Cos(theta)
  // gtr = rejection function according to Eq. (2.28)
  double rndArray[2];
  do {
    td->fRndm->uniform_array(2, rndArray);
    tsam = 2.0 * ac * (2.0 * rndArray[0] + a2 * std::sqrt(rndArray[0])) / (a2 * a2 - 4.0 * rndArray[0]);
    gtr  = (2.0 - tsam) * (a1 + 1.0 / (ac + tsam));

  } while (rndArray[1] * gtmax > gtr);

  cosTheta = 1.0 - tsam;
  return;
}

//____________________
double SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirection_Alias(double primekin, double r1, double r2,
                                                                           double r3) const
{

  // determine primary energy lower grid point
  if (primekin > 100 * geant::units::MeV) {
    return 1.;
  } else {
    // std::cout<<"::::SamplePhotoElectronDirection_Alias::::\n";

    double lGammaEnergy = Math::Log(primekin);
    int gammaEnergyIndx = (int)((lGammaEnergy - fPrimEnLMin) * fPrimEnILDelta);
    //
    if (gammaEnergyIndx >= fNumSamplingPrimEnergies - 1) gammaEnergyIndx = fNumSamplingPrimEnergies - 2;
    //
    double pLowerGammaEner = (fLSamplingPrimEnergies[gammaEnergyIndx + 1] - lGammaEnergy) * fPrimEnILDelta;
    if (r1 > pLowerGammaEner) {
      ++gammaEnergyIndx;
    }

    // sample the outgoing electron cosTheta
    double ecosTheta = fAliasSampler->SampleLinear(
        fAliasData[gammaEnergyIndx]->fXdata, fAliasData[gammaEnergyIndx]->fYdata, fAliasData[gammaEnergyIndx]->fAliasW,
        fAliasData[gammaEnergyIndx]->fAliasIndx, fAliasData[gammaEnergyIndx]->fNumdata, r2, r3);

    // This have to be seen (Transformation)
    // double xsi=fAliasSampler->SampleLinear(fAliasData[gammaEnergyIndx]->fXdata,
    // fAliasData[gammaEnergyIndx]->fYdata,fAliasData[gammaEnergyIndx]->fAliasW,
    // fAliasData[gammaEnergyIndx]->fAliasIndx,fAliasData[gammaEnergyIndx]->fNumdata,r2,r3);

    // double ecosTheta= Math::Exp(xsi)-2;
    return ecosTheta;
  }
}

double SauterGavrilaPhotoElectricModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts *,
                                                               double kinenergy, const Particle *)
{

  double xsec = 0.0;
  if (kinenergy < GetLowEnergyUsageLimit() || kinenergy > GetHighEnergyUsageLimit()) {
    return xsec;
  }
  // compute the parameterized atomic cross section: depends only on target Z and gamma energy.
  xsec = ComputeXSectionPerAtom(elem->GetZ(), kinenergy);
  return xsec;
}

double SauterGavrilaPhotoElectricModel::ComputeXSectionPerAtom(double zeta, double energy) const
{

  int verboseLevel = 0;
  using geant::units::barn;
  using geant::units::keV;
  using geant::units::MeV;

  if (verboseLevel > 3) {
    std::cout << "G4LivermorePhotoElectricModel_new::ComputeCrossSectionPerAtom():"
              << " Z= " << zeta << "  R(keV)= " << energy / keV << std::endl;
  }
  double cs = 0.0;
  int Z     = std::lrint(zeta);
  if (Z < 1 || Z >= gMaxSizeData) {
    return cs;
  }

  // EXTRA CHECK --
  // if( !fCrossSection[Z] && !fCrossSectionLE[Z] )
  //{
  //    ReadData(Z);
  //    if(!fCrossSectionLE[Z] && !fCrossSection[Z]) { return cs; }
  //}

  int idx = fNShells[Z] * 7 - 5; // 7: rows in the parameterization file; 5: number of parameters

  if (energy < (*(fParamHigh[Z]))[idx - 1]) {
    energy = (*(fParamHigh[Z]))[idx - 1];
  }

  double x1 = (MeV) / (energy);
  double x2 = x1 * x1;
  double x3 = x2 * x1;

  // (*) High energy parameterisation
  if (energy >= (*(fParamHigh[Z]))[0]) {
    double x4 = x2 * x2;
    double x5 = x4 * x1;

    cs = x1 * ((*(fParamHigh[Z]))[idx] + x1 * (*(fParamHigh[Z]))[idx + 1] + x2 * (*(fParamHigh[Z]))[idx + 2] +
               x3 * (*(fParamHigh[Z]))[idx + 3] + x4 * (*(fParamHigh[Z]))[idx + 4] + x5 * (*(fParamHigh[Z]))[idx + 5]);

  }
  // (**) Low energy parameterisation
  else if (energy >= (*(fParamLow[Z]))[0]) {
    double x4 = x2 * x2;
    double x5 = x4 * x1; // this variable usage can probably be optimized
    cs        = x1 * ((*(fParamLow[Z]))[idx] + x1 * (*(fParamLow[Z]))[idx + 1] + x2 * (*(fParamLow[Z]))[idx + 2] +
               x3 * (*(fParamLow[Z]))[idx + 3] + x4 * (*(fParamLow[Z]))[idx + 4] + x5 * (*(fParamLow[Z]))[idx + 5]);
  }

  // (***) Tabulated values above k-shell ionization energy
  else if (energy >= (*(fParamHigh[Z]))[1]) {
    // this is an extra-check - if energy<fCSVector[Z]->edgeMin it should go to the next else (****)
    if (energy < fCSVector[Z]->fEdgeMin)
      cs = x3 * fCSVector[Z]->fDataVector[0];
    else {
      size_t idx = 0;
      idx        = fCSVector[Z]->FindCSBinLocation(energy, idx);
      cs         = x3 * fCSVector[Z]->fSplineInt->GetValueAt(energy, idx);
    }

  }

  //(****) Tabulated values below k-shell ionization energy
  else {
    // this check is needed to have a constant cross-section(cs) value for energies below the lowest cs point.
    // in this way the gamma is always absorbed - instead of giving zero cs -
    if (energy < fLECSVector[Z]->fEdgeMin)
      cs = x3 * fLECSVector[Z]->fDataVector[0];
    else {
      size_t idx = 0;
      cs         = x3 * fLECSVector[Z]->GetValue(energy, idx);
    }
  }
  if (verboseLevel > 1) {
    std::cout << "LivermorePhotoElectricModel: E(keV)= " << energy / keV << " Z= " << Z << " cross(barn)= " << cs / barn
              << std::endl;
  }

  return cs;
}

double SauterGavrilaPhotoElectricModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy,
                                                                   const Particle * /*particle*/)
{

  double xsec = 0.0;
  if (kinenergy < GetLowEnergyUsageLimit() || kinenergy > GetHighEnergyUsageLimit()) {
    return xsec;
  }
  // compute the macroscopic cross section as the sum of the atomic cross sections weighted by the number of atoms in
  // in unit volume.
  const Material *mat = matcut->GetMaterial();
  double egamma       = kinenergy;
  // we need the element composition of this material
  const Vector_t<Element *> &theElements  = mat->GetElementVector();
  const double *theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int numElems                            = theElements.size();
  for (int iel = 0; iel < numElems; ++iel) {
    xsec += theAtomicNumDensityVector[iel] * ComputeXSectionPerAtom(theElements[iel]->GetZ(), egamma);
  }
  return xsec;
}

size_t SauterGavrilaPhotoElectricModel::SampleTargetElementIndex(const MaterialCuts *matCut, double gammaekin0,
                                                                 const double prestepmfp, geant::TaskData *td)
{
  size_t index                            = 0;
  const Material *mat                     = matCut->GetMaterial();
  const double *theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  const Vector_t<Element *> &theElements  = mat->GetElementVector();
  size_t num                              = mat->GetNumberOfElements();
  // sample target element index:
  // calculate the cumulative of the partial (per-element) macroscopic cross sections
  // the normalization factor i.e. the sum is known (the inevrse is the mfp)
  double cum = theAtomicNumDensityVector[index] * ComputeXSectionPerAtom(theElements[index]->GetZ(), gammaekin0);
  double rnd = td->fRndm->uniform();
  while (index < num-1 && rnd > cum*prestepmfp) {
    ++index;
    cum += theAtomicNumDensityVector[index] * ComputeXSectionPerAtom(theElements[index]->GetZ(), gammaekin0);
  }
  return index;
}

void SauterGavrilaPhotoElectricModel::TestSampleTargetElementIndex(const MaterialCuts *matcut, double energy,
                                                                   geant::TaskData *td)
{

  std::cout << "testSampleTargetElementIndex\n";
  size_t index = 0;
  double sum   = 0;
  // retrieve the elements vector
  const Vector_t<Element *> &theElements = matcut->GetMaterial()->GetElementVector();
  // retrieve the number of elements in the material
  int num = matcut->GetMaterial()->GetNumberOfElements();
  double xsec[num];
  double xsecSampled[num];

  const Material *mat                     = matcut->GetMaterial();
  const double *theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  for (int i = 0; i < num; i++) {
    xsec[i]        = theAtomicNumDensityVector[i] * ComputeXSectionPerAtom(theElements[i]->GetZ(), energy);
    xsecSampled[i] = 0.;
    sum += xsec[i];
  }
  for (int i = 0; i < 1000000000; i++) {
    index = SampleTargetElementIndex(matcut, energy, 1./sum, td);
    xsecSampled[index]++;
  }
  for (int i = 0; i < num; i++) {
    xsec[i] /= sum;
    xsecSampled[i] /= 1000000000;
  }

  char filename[521];
  sprintf(filename, "SampleTargetElementIndexTest_%s", (matcut->GetMaterial()->GetName()).c_str());
  FILE *f = fopen(filename, "w");

  for (int i = 0; i < num; ++i) {
    fprintf(f, "%d\t%.8g\t%.8g\n", i, xsec[i], xsecSampled[i]);
  }
  fclose(f);
}

int SauterGavrilaPhotoElectricModel::SampleSecondaries(LightTrack &track, geant::TaskData *td)
{
  using geant::units::MeV;
  double gammaekin0 = track.GetKinE();
  // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
  // check if kinetic energy is above fHighEnergyUsageLimit and do nothing if yes;
  if (gammaekin0 < GetLowEnergyUsageLimit() || gammaekin0 > GetHighEnergyUsageLimit()) {
    std::cout << "  --- MUST ADD CHECK1:  " << gammaekin0
              << " and  GetLowEnergyUsageLimit(): " << GetLowEnergyUsageLimit()
              << " - GetHighEnergyUsageLimit(): " << GetHighEnergyUsageLimit() << "\n";
    exit(-1);
    return 0; // numSecondaries is zero since the interaction is not happening
  }

  // interaction is possible so sample target element

  MaterialCuts *matCut                   = MaterialCuts::GetTheMaterialCutsTable()[track.GetMaterialCutCoupleIndex()];
  const Vector_t<Element *> &theElements = matCut->GetMaterial()->GetElementVector();
  size_t targetElemIndx                  = 0;
  if (theElements.size() > 1) {
    // uncomment the following line to test SampleTargetElementIndex
    // testSampleTargetElementIndex(matCut, gammaekin0, td );
    const double preStepMFP = track.GetTotalMFP();
    targetElemIndx = SampleTargetElementIndex(matCut, gammaekin0, preStepMFP, td);
  }
  double zeta = theElements[targetElemIndx]->GetZ();
  int Z       = std::lrint(zeta);
  // if element was not initialised, gamma should be absorbed
  if (!fCrossSectionLE[Z] && !fCrossSection[Z]) {
    track.SetEnergyDeposit(gammaekin0);
    std::cout << "Model not initialized. Element: " << Z << "! Depositing the energy\n";
    std::cout << "GetUseSamplingTables: " << GetUseSamplingTables() << std::endl;
    exit(-1);
    return 0;
  }
  // SAMPLING OF THE SHELL

  double r1       = td->fRndm->uniform();
  size_t shellIdx = 0;
  size_t tmp      = Z;

  /* TEMPORARY COMMENTED OUT UNTILL DEBUGGING THE CRASH IN SampleShellAlias (A.G. Jan 8 2019) */
  //  if (GetUseSamplingTables() && fNShells[tmp] > 1) {
  //    double r2 = td->fRndm->uniform();
  //    SampleShellAlias(gammaekin0, tmp, r1, r2, shellIdx);
  //  } else {
  //    SampleShell(gammaekin0, Z, r1, shellIdx);
  //  }

  if (fNShells[tmp] > 1) SampleShell(gammaekin0, Z, r1, shellIdx);

  // Retrieving ionized shell bindingEnergy
  double bindingEnergy = (*(fParamHigh[Z]))[shellIdx * 7 + 1];
  if (gammaekin0 < bindingEnergy) {
    track.SetEnergyDeposit(gammaekin0);
    std::cout << "Ekin " << gammaekin0 << " lower than bindingEnergy: " << bindingEnergy << std::endl;
    std::cout
        << "SauterGavrilaPhotoElectricModel::SampleSecondaries: must add this check to the vectorized version! \n";
    exit(-1);
    return 0; // numSecondaries
  }

  // since edep is equal to bindingenergy I get rid of it
  // double edep = bindingEnergy;

  //  deexcitation is MISSING for now
  /*
   DEEXCITATION
   */
  //

  // Create the secondary particle: the photoelectron
  double elecKineEnergy = gammaekin0 - bindingEnergy;
  double cosTheta       = track.GetDirZ(); // no need to initialize here cosTheta !
  double sinTheta       = 0.0;
  double phi            = 0.0;
  double eDirX1;
  double eDirY1;
  double eDirZ1;

  if (gammaekin0 <= 100 * geant::units::MeV) {
    if (!GetUseSamplingTables()) {
      SamplePhotoElectronDirection_Rejection(gammaekin0, cosTheta, td);
    } else {
      double *rndArray = td->fDblArray;
      td->fRndm->uniform_array(3, rndArray);
      cosTheta = SamplePhotoElectronDirection_Alias(gammaekin0, rndArray[0], rndArray[1], rndArray[2]);
    }
    sinTheta   = std::sqrt((1 - cosTheta) * (1 + cosTheta));
    double rnd = td->fRndm->uniform();
    phi        = geant::units::kTwoPi * rnd;

    // new photoelectron direction in the scattering frame
    eDirX1 = sinTheta * Math::Cos(phi);
    eDirY1 = sinTheta * Math::Sin(phi);
    eDirZ1 = cosTheta;

    // rotate new photoelectron direction to the lab frame:
    Math::RotateToLabFrame(eDirX1, eDirY1, eDirZ1, track.GetDirX(), track.GetDirY(), track.GetDirZ());

  } else {
    eDirX1 = track.GetDirX();
    eDirY1 = track.GetDirY();
    eDirZ1 = track.GetDirZ();
  }

  // create the secondary particle i.e. the photoelectron
  LightTrack &emTrack = td->fPhysicsData->InsertSecondary();
  emTrack.SetDirX(eDirX1);
  emTrack.SetDirY(eDirY1);
  emTrack.SetDirZ(eDirZ1);
  emTrack.SetKinE(elecKineEnergy);
  emTrack.SetGVcode(fSecondaryInternalCode); // electron GV code
  emTrack.SetMass(geant::units::kElectronMassC2);
  emTrack.SetTrackIndex(track.GetTrackIndex()); // parent Track index

  /*if(fabs(gammaekin0 - elecKineEnergy - esec - edep) > geant::units::eV) {
   std::cout << "### SauterGavrilaPhotoElectricModel dE(eV)= "
   << (gammaekin0 - elecKineEnergy - esec - edep)/geant::units::eV
   << "  shell= " << shellIdx
   << "  E(keV)= " << gammaekin0/geant::units::keV
   << "  Ebind(keV)= " << bindingEnergy/geant::units::keV
   << "  Ee(keV)= " << elecKineEnergy/geant::units::keV
   << "  Esec(keV)= " << esec/geant::units::keV
   << "  Edep(keV)= " << edep/geant::units::keV
   << std::endl;
   }*/

  // always kill primary photon
  track.SetTrackStatus(LTrackStatus::kKill);
  track.SetKinE(0.0);
  // edep is = bindingEnergy
  if (bindingEnergy > 0.0) {
    track.SetEnergyDeposit(bindingEnergy);
  }
  // return with number of secondaries i.e. 1 photoelectron
  return 1;
}

void SauterGavrilaPhotoElectricModel::SampleSecondaries(LightTrack_v &tracks, geant::TaskData *td)
{

  int N                      = tracks.GetNtracks();
  double *kin                = tracks.GetKinEArr();
  auto zed                   = td->fPhysicsData->fPhysicsScratchpad.fIzet;
  auto nshells               = td->fPhysicsData->fPhysicsScratchpad.fNshells;
  auto sampledShells         = td->fPhysicsData->fPhysicsScratchpad.fSampledShells;
  double *cosTheta           = td->fPhysicsData->fPhysicsScratchpad.fDoubleArr;
  double energyDepositionVec = 0.;

  for (int i = 0; i < N; ++i) {
    const MaterialCuts *matCut             = MaterialCuts::GetMaterialCut(tracks.GetMaterialCutCoupleIndex(i));
    const Vector_t<Element *> &theElements = matCut->GetMaterial()->GetElementVector();
    size_t targetElemIndx                  = 0;
    if (theElements.size() > 1) {
      const double preStepMFP = tracks.GetTotalMFP(i);
      targetElemIndx = SampleTargetElementIndex(matCut, kin[i], preStepMFP, td);
    }
    zed[i]     = (int)theElements[targetElemIndx]->GetZ();
    nshells[i] = fNShellsUsed[zed[i]];
  }
  if (GetUseSamplingTables()) {
    /* TEMPORARY INSERTED SCALAR VERSION UNTILL DEBUGGING THE CRASH IN SampleShellAlias (A.G. Jan 8 2019) */
    size_t shellIdx = 0;
    // calling the scalar implementation for the sampling of the shell
    for (int i = 0; i < N; i++) {
      double r1 = td->fRndm->uniform();
      SampleShell(kin[i], zed[i], r1, shellIdx);
      sampledShells[i] = shellIdx;
    }

    for (int i = 0; i < N; i += kVecLenD) {

      Double_v gammaekin;
      // IndexD_v nshells_v, z;
      vecCore::Load(gammaekin, kin + i);
      // vecCore::Load(nshells_v, nshells + i);
      // vecCore::Load(z, zed + i);

      // SAMPLING OF THE SHELL WITH ALIAS
      // IndexD_v sampledShells_v(0);
      // Double_v r1     = td->fRndm->uniformV();
      // Double_v r2     = td->fRndm->uniformV();
      // sampledShells_v = SampleShellAliasVec(gammaekin, z, r1, r2);
      // vecCore::Store(sampledShells_v, sampledShells + i);

      // SAMPLING OF THE ANGLE WITH ALIAS
      Double_v cosTheta_v;
      MaskDI_v activateSamplingAngle(gammaekin <= 100 * geant::units::MeV);
      if (!vecCore::MaskEmpty(activateSamplingAngle)) {
        Double_v r1 = td->fRndm->uniformV();
        Double_v r2 = td->fRndm->uniformV();
        Double_v r3 = td->fRndm->uniformV();
        cosTheta_v  = SamplePhotoElectronDirectionAliasVec(gammaekin, r1, r2, r3);
        vecCore::Store(cosTheta_v, cosTheta + i);
      }
    }

  } else {

    size_t shellIdx = 0;
    // calling the scalar implementation for the sampling of the shell
    for (int i = 0; i < N; i++) {
      double r1 = td->fRndm->uniform();
      SampleShell(kin[i], zed[i], r1, shellIdx);
      sampledShells[i] = shellIdx;
    }
    //            double *rands = td->fDblArray;
    //            //NB: generating random numbers here just for reproducibility issues
    //            td->fRndm->uniform_array(N, rands);
    //            SampleShellVec(kin, zed, sampledShells, N, td, rands);

    // Vectorized Sampling of the Angle
    SamplePhotoElectronDirectionRejVec(kin, cosTheta, N, td);
  }
  int globalCounter = 0;
  for (int i = 0; i < N; i += kVecLenD) {

    Double_v gammaekin_v, cosTheta_v;
    IndexD_v zed_v;
    vecCore::Load(gammaekin_v, kin + i);
    vecCore::Load(zed_v, zed + i);
    vecCore::Load(cosTheta_v, cosTheta + i);
    // Retrieving ionized shell bindingEnergy
    Double_v bindingEnergy_v(0);
    for (int k = 0; k < kVecLenD; ++k) {
      vecCore::Set(bindingEnergy_v, k, (*(fParamHigh[vecCore::Get(zed_v, k)]))[sampledShells[k + i] * 7 + 1]);
    }

    // Create the secondary particle e-
    Double_v eekin = gammaekin_v - bindingEnergy_v;
    MaskDI_v activateSamplingAngle(gammaekin_v <= 100 * geant::units::MeV);

    Double_v eDirX1;
    Double_v eDirY1;
    Double_v eDirZ1;
    if (!vecCore::MaskEmpty(activateSamplingAngle)) {

      Double_v sinTheta = vecCore::math::Sqrt((1.0 - cosTheta_v) * (1.0 + cosTheta_v));
      Double_v phi      = geant::units::kTwoPi * td->fRndm->uniformV();

      // new photoelectron direction in the scattering frame
      eDirX1 = sinTheta * vecCore::math::Cos(phi);
      eDirY1 = sinTheta * vecCore::math::Sin(phi);
      eDirZ1 = cosTheta_v;

      Math::RotateToLabFrame(eDirX1, eDirY1, eDirZ1, tracks.GetDirXVec(i), tracks.GetDirYVec(i), tracks.GetDirZVec(i));
    }
    vecCore::MaskedAssign(eDirX1, !activateSamplingAngle, tracks.GetDirXVec(i));
    vecCore::MaskedAssign(eDirY1, !activateSamplingAngle, tracks.GetDirYVec(i));
    vecCore::MaskedAssign(eDirZ1, !activateSamplingAngle, tracks.GetDirZVec(i));
    LightTrack_v &secondaries = td->fPhysicsData->GetSecondarySOA();
    for (int l = 0; l < kVecLenD; ++l) {
      int idx = secondaries.InsertTrack();
      secondaries.SetKinE(Get(eekin, l), idx);
      secondaries.SetDirX(Get(eDirX1, l), idx);
      secondaries.SetDirY(Get(eDirY1, l), idx);
      secondaries.SetDirZ(Get(eDirZ1, l), idx);
      secondaries.SetGVcode(fSecondaryInternalCode, idx); // electron GV code
      secondaries.SetMass(geant::units::kElectronMassC2, idx);
      secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx); // parent Track index
    }

    // update primary track - always kill primary photon
    for (int l = 0; l < kVecLenD; ++l) {
      bool canDeposit     = kin[i + l] > GetLowEnergyUsageLimit() && kin[i + l] < GetHighEnergyUsageLimit();
      double be           = Get(bindingEnergy_v, l);
      bool notDepositEkin = kin[i + l] > be;
      if (be > 0.0 && canDeposit && notDepositEkin) {
        globalCounter++;
        tracks.SetEnergyDeposit(be, i + l);
        energyDepositionVec += be;
      } else {
        std::cout << "Cannot deposit: kin[" << i + l << "]: " << kin[i + l] << " <  " << GetLowEnergyUsageLimit()
                  << " or >  " << GetHighEnergyUsageLimit() << std::endl;
        tracks.SetEnergyDeposit(kin[i + l], i + l);
        exit(-1);
      }
      tracks.SetKinE(0.0, i + l);
      tracks.SetTrackStatus(LTrackStatus::kKill, i + l);
    }
  }
}
IndexD_v SauterGavrilaPhotoElectricModel::SampleShellAliasVec(Double_v egamma, IndexD_v zed, Double_v r1, Double_v r2)
{

  IndexD_v sampledShells(0);
  MaskDI_v enableSamplingShells = (zed != 1) && (zed != 2);
  if (!vecCore::MaskEmpty(enableSamplingShells)) {
    // this will be directly stored in the track
    Double_v lGammaEnergy_v = vecCore::math::Log(egamma);

    // LOWER bin in the BASE vector
    IndexD_v tableIndexBase_v = (IndexD_v)((lGammaEnergy_v - fShellPrimEnLMin) * Double_v(fShellPrimEnILDelta));

    // These are static informations that can be passed as an argument in Real_v form - TO DO
    Double_v kBindingEn_v;
    for (int k = 0; k < kVecLenD; k++) {
      vecCore::Set(kBindingEn_v, k, fBindingEn[vecCore::Get(zed, k)][0]);
    }

    MaskDI_v lowEn(egamma < kBindingEn_v);
    // IndexD_v tableIndexBinding_v;
    Double_v baseEn_v(999), bindEn_v(999);
    IndexD_v indexBaseEn_v(-1), indexBindingEn_v(-1);

    for (int k = 0; k < kVecLenD; k++) {
      if (vecCore::Get(enableSamplingShells, k)) {
        if (fIndexBaseEn[vecCore::Get(zed, k)].size() > (size_t)(vecCore::Get(tableIndexBase_v, k) + 1)) {
          vecCore::Set(indexBaseEn_v, k, fIndexBaseEn[vecCore::Get(zed, k)][vecCore::Get(tableIndexBase_v, k) + 1] - 1);
        } else {
          vecCore::Set(indexBaseEn_v, k, fIndexBaseEn[vecCore::Get(zed, k)][vecCore::Get(tableIndexBase_v, k)] - 1);
        }
      }
    }

    IndexD_v tableIndex_v(indexBaseEn_v);
    // Only the values of tableIndex_v that need to be changed
    if (!vecCore::MaskEmpty(lowEn)) {
      for (int k = 0; k < kVecLenD; k++) {

        if (vecCore::Get(lowEn, k) && vecCore::Get(enableSamplingShells, k)) {
          vecCore::Set(
              baseEn_v, k,
              fShellSamplingPrimEnergies[vecCore::Get(tableIndexBase_v, k) + 1]); // UPPER VALUE (it could be
                                                                                  // the last meaningful value in
                                                                                  // the vector (for those that
                                                                                  // have the mask set to FALSE)

          // LOWER bin in the BINDING vector
          int tableIndexBinding =
              std::lower_bound(fSortedDoubledBindingEn[vecCore::Get(zed, k)].begin(),
                               fSortedDoubledBindingEn[vecCore::Get(zed, k)].end(), vecCore::Get(egamma, k)) -
              fSortedDoubledBindingEn[vecCore::Get(zed, k)].begin() - 1;
          vecCore::Set(bindEn_v, k,
                       fSortedDoubledBindingEn[vecCore::Get(zed, k)][tableIndexBinding + 1]); // UPPER VALUE
          vecCore::Set(indexBindingEn_v, k,
                       fIndexSortedDoubledBindingEn[vecCore::Get(zed, k)][tableIndexBinding + 1] - 1);
        }
      }
      MaskDI_v checkMinVal(baseEn_v > bindEn_v); // If TRUE, take bindingEnergies, otherwise keep the base energies
      vecCore::MaskedAssign(tableIndex_v, checkMinVal && lowEn, indexBindingEn_v);
    }
    IndexD_v lastSSAliasIndex_v;
    for (int k = 0; k < kVecLenD; k++) {
      if (vecCore::Get(enableSamplingShells, k))
        vecCore::Set(lastSSAliasIndex_v, k, fLastSSAliasIndex[vecCore::Get(zed, k) - 1]);
    }
    Double_v val = (lGammaEnergy_v - fShellPrimEnLMin) *
                   fShellPrimEnILDelta; // To correct - inverse of delta of the log of real en
    // LOWER energy bin index
    IndexD_v indxEgamma = (IndexD_v)val;
    Double_v pIndxHigh  = val - indxEgamma;
    MaskDI_v check(r1 <= pIndxHigh);
    vecCore::MaskedAssign(tableIndex_v, check, tableIndex_v + 1);
    IndexD_v indxTable_v = lastSSAliasIndex_v + tableIndex_v;
    // NB: the SCALAR and the VECTORIZED are almost equivalent
    // SCALAR
    for (int i = 0; i < kVecLenD; i++) {
      if (vecCore::Get(enableSamplingShells, i)) {
        int indxI = vecCore::Get(indxTable_v, i);
        int xsampl =
            fShellAliasSampler->SampleDiscrete(fShellAliasData[indxI]->fAliasW, fShellAliasData[indxI]->fAliasIndx,
                                               fShellAliasData[indxI]->fNumdata, vecCore::Get(r2, i));
        vecCore::Set(sampledShells, i, xsampl);
      }
    }
    // END SCALAR

    //    //VECTORIZED
    //    Real_v aliasW_v, aliasIndx_v, aliasNumdata_v;
    //    //first I need the numData
    //    for (size_t kk=0; kk<kRealSize; kk++)
    //        vecCore::Set(aliasNumdata_v, kk, fShellAliasData[indxTable_v[kk]]->fNumdata);
    //
    //    Real_v rest_v  = r2*aliasNumdata_v;
    //    Real_v indxBin = vecCore::math::Floor(rest_v);
    //    RIndex indxBin_v(indxBin);
    //
    //    for (size_t kk=0; kk<kRealSize; kk++){
    //        vecCore::Set(aliasW_v, kk, fShellAliasData[indxTable_v[kk]]->fAliasW[indxBin_v[kk]]);
    //        vecCore::Set(aliasIndx_v, kk, fShellAliasData[indxTable_v[kk]]->fAliasIndx[indxBin_v[kk]]);
    //
    //    }
    //    RMask check3(aliasW_v<rest_v-indxBin_v);
    //    vecCore::MaskedAssign(indxBin,check3,aliasIndx_v);
    //    //std::cout<<indxBin_v<<std::endl;
    //    RIndex temp(indxBin);
    //    sampledShells=temp;
    //    //END VECTORIZED
  }
  return sampledShells;
}

Double_v SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirectionAliasVec(Double_v egamma, Double_v r1,
                                                                               Double_v r2, Double_v r3)
{
  Double_v ecosT(1.);
  MaskDI_v checkEnergy = egamma > (Double_v)100 * geant::units::MeV;
  // if some energy is below 100 MeV
  if (!vecCore::MaskFull(checkEnergy)) {
    Double_v legamma = vecCore::math::Log(egamma);
    //        //With this implementation is assumed that the model is never built for energies outside the range where
    //        the alias tables where built Double_v val        = (legamma - fPrimEnLMin) * fPrimEnILDelta; IndexD_v
    //        gammaEnergyIndx = (IndexD_v)val; // lower electron energy bin index Double_v pIndxHigh  = val -
    //        gammaEnergyIndx; MaskD_v checkIndex = r1 < pIndxHigh; if (!vecCore::MaskEmpty(checkIndex)) {
    //            vecCore::MaskedAssign(gammaEnergyIndx, checkIndex, gammaEnergyIndx + 1);
    //        }

    IndexD_v gammaEnergyIndx = (IndexD_v)((legamma - fPrimEnLMin) * fPrimEnILDelta);
    MaskDI_v check1          = (gammaEnergyIndx >= (fNumSamplingPrimEnergies - (Double_v)1));
    vecCore::MaskedAssign(gammaEnergyIndx, check1, IndexD_v(fNumSamplingPrimEnergies - 2));
    Double_v fLSamplingPrimEnergies_v;
    for (int i = 0; i < kVecLenD; i++)
      vecCore::Set(fLSamplingPrimEnergies_v, i, fLSamplingPrimEnergies[vecCore::Get(gammaEnergyIndx, i) + 1]);
    Double_v pLowerGammaEner = (fLSamplingPrimEnergies_v - legamma) * fPrimEnILDelta;
    MaskDI_v check2          = r1 > pLowerGammaEner;
    vecCore::MaskedAssign(gammaEnergyIndx, check2, gammaEnergyIndx + 1);

    // going scalar here
    for (int i = 0; i < kVecLenD; ++i) {
      if (vecCore::Get(egamma, i) < 100 * geant::units::MeV) {
        size_t idxI      = vecCore::Get(gammaEnergyIndx, i);
        double ecosTheta = fAliasSampler->SampleLinear(
            fAliasData[idxI]->fXdata, fAliasData[idxI]->fYdata, fAliasData[idxI]->fAliasW, fAliasData[idxI]->fAliasIndx,
            fAliasData[idxI]->fNumdata, vecCore::Get(r2, i), vecCore::Get(r3, i));
        vecCore::Set(ecosT, i, ecosTheta);
      }
    }
  }
  return ecosT;
}

void SauterGavrilaPhotoElectricModel::SampleShellVec(double *egamma, int *zed, int *ss, int N,
                                                     const geant::TaskData * /*td*/, double *randoms)
{

  IndexD_v nshells_v;
  IndexD_v sampledShells_v(0);

  std::vector<double> ehep;
  std::vector<double> elep;
  std::vector<double> etab;

  std::vector<int> zhep;
  std::vector<int> zlep;
  std::vector<int> ztab;

  std::vector<int> indexhep;
  std::vector<int> indexlep;
  std::vector<int> indextab;

  std::vector<int> nshellshep;
  std::vector<int> nshellslep;
  std::vector<int> nshellstab;

  std::vector<double> randhep;
  std::vector<double> randlep;
  std::vector<double> randtab;

  // PREFILTERING
  for (int i = 0; i < N; ++i) {

    if (egamma[i] >= (*(fParamHigh[zed[i]]))[0]) {
      ehep.push_back(egamma[i]);
      zhep.push_back(zed[i]);
      nshellshep.push_back(fNShells[zed[i]]);
      indexhep.push_back(i);
      randhep.push_back(randoms[i]);
    } else if (egamma[i] >= (*(fParamLow[zed[i]]))[0]) {
      elep.push_back(egamma[i]);
      zlep.push_back(zed[i]);
      nshellslep.push_back(fNShells[zed[i]]);
      indexlep.push_back(i);
      randlep.push_back(randoms[i]);
    } else {
      etab.push_back(egamma[i]);
      ztab.push_back(zed[i]);
      indextab.push_back(i);
      randtab.push_back(randoms[i]);
    }
  }

  //**** PROCESS THE LEP
  size_t currlep = 0;
  MaskDI_v lanesDonelep(0); // MaskDI_v::Zero(); // no lanes done
  IndexD_v idxlep;

  for (int l = 0; l < kVecLenD; ++l) {
    vecCore::Set(idxlep, l, currlep++); // indexes initialization
  }
  IndexD_v idxForLoop(0);
  int sampledShellslep[elep.size()];

  while (elep.size() > 3 && (currlep < elep.size() || !vecCore::MaskFull(lanesDonelep))) {
    IndexD_v zeds     = vecCore::Gather<IndexD_v>(zlep.data(), idxlep);
    Double_v egamma_v = vecCore::Gather<Double_v>(elep.data(), idxlep);
    Double_v iegamma  = (geant::units::MeV) / egamma_v;
    Double_v iegamma2 = iegamma * iegamma;
    Double_v iegamma3 = iegamma2 * iegamma;
    Double_v iegamma4 = iegamma2 * iegamma2;
    Double_v iegamma5 = iegamma2 * iegamma3;

    IndexD_v totShells = vecCore::Gather<IndexD_v>(nshellslep.data(), idxlep);
    IndexD_v idx_v     = totShells * 7 - 5;

    Double_v pm1, p0, p1, p2, p3, p4, p5;

    for (int k = 0; k < kVecLenD; k++) {
      if (!vecCore::Get(lanesDonelep, k)) {
        size_t zedk = vecCore::Get(zeds, k);
        size_t idxk = vecCore::Get(idx_v, k);
        vecCore::Set(p0, k, (*(fParamLow[zedk]))[idxk]);
        vecCore::Set(p1, k, (*(fParamLow[zedk]))[idxk + 1]);
        vecCore::Set(p2, k, (*(fParamLow[zedk]))[idxk + 2]);
        vecCore::Set(p3, k, (*(fParamLow[zedk]))[idxk + 3]);
        vecCore::Set(p4, k, (*(fParamLow[zedk]))[idxk + 4]);
        vecCore::Set(p5, k, (*(fParamLow[zedk]))[idxk + 5]);
      }
    }

    Double_v rand_v    = vecCore::Gather<Double_v>(randlep.data(), idxlep);
    Double_v cs0       = rand_v * (p0 + iegamma * p1 + iegamma2 * p2 + iegamma3 * p3 + iegamma4 * p4 + iegamma5 * p5);
    Double_v idxShells = idxForLoop * 7 + 2;
    for (int k = 0; k < kVecLenD; k++) {
      if (!vecCore::Get(lanesDonelep, k)) {
        size_t zedk = vecCore::Get(zeds, k);
        size_t idxk = vecCore::Get(idxShells, k);
        vecCore::Set(pm1, k, (*(fParamLow[zedk]))[idxk] - 1);
        vecCore::Set(p0, k, (*(fParamLow[zedk]))[idxk]);
        vecCore::Set(p1, k, (*(fParamLow[zedk]))[idxk + 1]);
        vecCore::Set(p2, k, (*(fParamLow[zedk]))[idxk + 2]);
        vecCore::Set(p3, k, (*(fParamLow[zedk]))[idxk + 3]);
        vecCore::Set(p4, k, (*(fParamLow[zedk]))[idxk + 4]);
        vecCore::Set(p5, k, (*(fParamLow[zedk]))[idxk + 5]);
      }
    }

    MaskDI_v checkKinE(egamma_v > pm1);
    Double_v cs = p0 + iegamma * p1 + iegamma2 * p2 + iegamma3 * p3 + iegamma4 * p4 + iegamma5 * p5;
    MaskDI_v accepted(cs >= cs0);
    MaskDI_v lastShell(idxForLoop == totShells - 1);
    MaskDI_v checkOut(accepted || lastShell);
    vecCore::MaskedAssign(idxForLoop, !checkOut, idxForLoop + 1);
    // Could scatter directly to the original indexes
    vecCore::Scatter(idxForLoop, sampledShellslep, idxlep);
    vecCore::MaskedAssign(idxForLoop, checkOut, (IndexD_v)0);
    lanesDonelep = lanesDonelep || checkOut;
    for (int l = 0; l < kVecLenD; ++l) {
      auto laneDone = vecCore::Get(checkOut, l);
      if (laneDone) {
        if (currlep < elep.size()) {
          vecCore::Set(idxlep, l, currlep++);
          vecCore::Set(lanesDonelep, l, false);
        } else {
          vecCore::Set(idxlep, l, elep.size());
        }
      }
    }
  }
  for (size_t i = 0; i < elep.size(); i++) {
    ss[indexlep[i]] = sampledShellslep[i];
  }

  //**** PROCESS THE HEP
  size_t currhep = 0;
  MaskDI_v lanesDonehep(0); // MaskDI_v::Zero(); // no lanes done
  IndexD_v idxhep;

  for (int l = 0; l < kVecLenD; ++l) {
    vecCore::Set(idxhep, l, currhep++); // indexes initialization
  }
  idxForLoop = 0;
  int sampledShellshep[ehep.size()];

  while (ehep.size() > 3 && (currhep < ehep.size() || !vecCore::MaskFull(lanesDonehep))) {

    IndexD_v zeds     = vecCore::Gather<IndexD_v>(zhep.data(), idxhep);
    Double_v egamma_v = vecCore::Gather<Double_v>(ehep.data(), idxhep);
    Double_v rand_v   = vecCore::Gather<Double_v>(randhep.data(), idxhep);

    Double_v iegamma   = (geant::units::MeV) / egamma_v;
    Double_v iegamma2  = iegamma * iegamma;
    Double_v iegamma3  = iegamma2 * iegamma;
    Double_v iegamma4  = iegamma2 * iegamma2;
    Double_v iegamma5  = iegamma2 * iegamma3;
    IndexD_v totShells = vecCore::Gather<IndexD_v>(nshellshep.data(), idxhep);
    IndexD_v idx_v     = totShells * 7 - 5;

    Double_v pm1, p0, p1, p2, p3, p4, p5;

    for (int k = 0; k < kVecLenD; k++) {
      if (!vecCore::Get(lanesDonehep, k)) {
        int zedk = vecCore::Get(zeds, k);
        int idxk = vecCore::Get(idx_v, k);
        vecCore::Set(p0, k, (*(fParamHigh[zedk]))[idxk]);
        vecCore::Set(p1, k, (*(fParamHigh[zedk]))[idxk + 1]);
        vecCore::Set(p2, k, (*(fParamHigh[zedk]))[idxk + 2]);
        vecCore::Set(p3, k, (*(fParamHigh[zedk]))[idxk + 3]);
        vecCore::Set(p4, k, (*(fParamHigh[zedk]))[idxk + 4]);
        vecCore::Set(p5, k, (*(fParamHigh[zedk]))[idxk + 5]);
      }
    }
    Double_v cs0 = rand_v * (p0 + iegamma * p1 + iegamma2 * p2 + iegamma3 * p3 + iegamma4 * p4 + iegamma5 * p5);

    Double_v idxShells = idxForLoop * 7 + 2;
    for (int k = 0; k < kVecLenD; k++) {
      if (!vecCore::Get(lanesDonehep, k)) {
        size_t zedk = vecCore::Get(zeds, k);
        size_t idxk = vecCore::Get(idxShells, k);
        vecCore::Set(pm1, k, (*(fParamHigh[zedk]))[idxk - 1]);
        vecCore::Set(p0, k, (*(fParamHigh[zedk]))[idxk]);
        vecCore::Set(p1, k, (*(fParamHigh[zedk]))[idxk + 1]);
        vecCore::Set(p2, k, (*(fParamHigh[zedk]))[idxk + 2]);
        vecCore::Set(p3, k, (*(fParamHigh[zedk]))[idxk + 3]);
        vecCore::Set(p4, k, (*(fParamHigh[zedk]))[idxk + 4]);
        vecCore::Set(p5, k, (*(fParamHigh[zedk]))[idxk + 5]);
      }
    }

    MaskDI_v checkKinE(egamma_v > pm1);
    Double_v cs = p0 + iegamma * p1 + iegamma2 * p2 + iegamma3 * p3 + iegamma4 * p4 + iegamma5 * p5;
    MaskDI_v accepted(cs >= cs0);
    MaskDI_v lastShell(idxForLoop == totShells - 1);
    MaskDI_v checkOut(accepted || lastShell);
    vecCore::MaskedAssign(idxForLoop, !checkOut, idxForLoop + 1);
    // Could scatter directly to the original indexes
    vecCore::Scatter(idxForLoop, sampledShellshep, idxhep);
    vecCore::MaskedAssign(idxForLoop, checkOut, (IndexD_v)0);
    lanesDonehep = lanesDonehep || checkOut;
    for (int l = 0; l < kVecLenD; ++l) {
      auto laneDone = vecCore::Get(checkOut, l);
      if (laneDone) {
        if (currhep < ehep.size()) {
          vecCore::Set(idxhep, l, currhep++);
          vecCore::Set(lanesDonehep, l, false);
        } else {
          vecCore::Set(idxhep, l, ehep.size());
        }
      }
    }
  }
  for (size_t i = 0; i < ehep.size(); i++) {
    ss[indexhep[i]] = sampledShellshep[i];
  }

  //**** PROCESS THE LET && HET in old SCALAR MODE
  for (size_t i = 0; i < etab.size(); i++) {
    size_t nshells  = fNShells[ztab[i]];
    size_t shellIdx = 0;
    if (nshells > 1) {
      double cs  = randtab[i];
      size_t idx = 0;
      // (***) Tabulated values above k-shell ionization energy
      if (etab[i] >= (*(fParamHigh[ztab[i]]))[1]) {
        idx = fCSVector[ztab[i]]->FindCSBinLocation(etab[i], idx);
        cs *= fCSVector[ztab[i]]->fSplineInt->GetValueAt(etab[i], idx);
      }
      //(****) Tabulated values below k-shell ionization energy
      else {
        cs *= fLECSVector[ztab[i]]->GetValue(etab[i], idx);
      }
      for (size_t j = 0; j < nshells; ++j) {
        shellIdx = (size_t)fShellVector[ztab[i]][j]->fCompID;
        if (etab[i] > (*(fParamLow[ztab[i]]))[7 * shellIdx + 1]) {
          size_t idx = 0;
          cs -= fShellVector[ztab[i]][j]->GetValue(etab[i], idx);
        }
        if (cs <= 0.0 || j + 1 == nshells) break;
      }
      ss[indextab[i]] = shellIdx;
    } else
      ss[indextab[i]] = 0;
  }
  //    std::cout<<"*******"<<std::endl;
  //    std::cout<<"N: "<<N<<std::endl;
  //    std::cout<<"lep.size(): "<<elep.size()<<std::endl;
  //    std::cout<<"hep.size(): "<<ehep.size()<<std::endl;
  //    std::cout<<"tab.size(): "<<etab.size()<<std::endl;
}

void SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirectionRejVec(const double *egamma, double *cosTheta, int N,
                                                                         const geant::TaskData *td)
{

  int currN = 0;
  MaskDI_v lanesDone(0); // MaskDI_v::Zero(); // no lanes done
  IndexD_v idx;
  for (int l = 0; l < kVecLenD; ++l) {
    vecCore::Set(idx, l, currN++); // indexes initialization
  }
  while (currN < N || !vecCore::MaskFull(lanesDone)) {

    // 1) initialize energy-dependent variables
    // Variable naming according to Eq. (2.24) of Penelope Manual
    // (pag. 44)
    Double_v gamma  = 1.0 + vecCore::Gather<Double_v>(egamma, idx) / geant::units::kElectronMassC2;
    Double_v gamma2 = gamma * gamma;
    Double_v beta   = std::sqrt((gamma2 - 1.0) / gamma2);

    // ac corresponds to "A" of Eq. (2.31)
    //
    Double_v ac = (1.0 / beta) - 1.0;
    Double_v a1 = 0.5 * beta * gamma * (gamma - 1.0) * (gamma - 2.0);
    Double_v a2 = ac + 2.0;
    // gtmax = maximum of the rejection function according to Eq. (2.28), obtained for tsam=0
    Double_v gtmax = 2.0 * (a1 + 1.0 / ac);
    Double_v tsam  = 0;
    Double_v gtr   = 0;

    // 2) sampling. Eq. (2.31) of Penelope Manual
    // tsam = 1-std::cos(theta)
    // gtr = rejection function according to Eq. (2.28)
    Double_v rnd1 = td->fRndm->uniformV(); // here, wasting random numbers  - to be handled for reproducibility issues
    Double_v rnd2 = td->fRndm->uniformV();

    tsam                = 2.0 * ac * (2.0 * rnd1 + a2 * std::sqrt(rnd1)) / (a2 * a2 - 4.0 * rnd1);
    gtr                 = (2.0 - tsam) * (a1 + 1.0 / (ac + tsam));
    MaskDI_v cond1      = rnd2 * gtmax > gtr;
    Double_v cosTheta_v = 1.0 - tsam;

    // Scatter anyway all the values, but if cond1 is false the lane will be processed again
    if (!vecCore::MaskEmpty(cond1)) {
      vecCore::Scatter(cosTheta_v, cosTheta, idx);
    }
    lanesDone = lanesDone || cond1;
    for (int l = 0; l < kVecLenD; ++l) {
      auto laneDone = vecCore::Get(cond1, l);
      if (laneDone) {
        if (currN < N) {
          vecCore::Set(idx, l, currN++);
          vecCore::Set(lanesDone, l, false);
        } else {
          vecCore::Set(idx, l, N);
        }
      }
    }
  }
}

void SauterGavrilaPhotoElectricModel::SampleShell(double kinE, const int Z, double &rand, size_t &sampledShells)
{
  if (!((fCrossSection[Z]) && ((fCrossSectionLE[Z] && Z > 2) || (!fCrossSectionLE[Z] && Z < 3)))) {
    std::cout << "Data not loaded before!\n";
    ReadData(Z);
  }
  size_t nshells  = fNShells[Z];
  size_t shellIdx = 0;

  if (nshells > 1) {
    // (*) High energy parameterisation
    if (kinE >= (*(fParamHigh[Z]))[0]) {

      double x1 = (geant::units::MeV) / kinE;
      double x2 = x1 * x1;
      double x3 = x2 * x1;
      double x4 = x3 * x1;
      double x5 = x4 * x1;
      int idx   = nshells * 7 - 5;
      // when do sampling common factors are not taken into account
      // so cross section is not real
      double cs0 = rand * ((*(fParamHigh[Z]))[idx] + x1 * (*(fParamHigh[Z]))[idx + 1] +
                           x2 * (*(fParamHigh[Z]))[idx + 2] + x3 * (*(fParamHigh[Z]))[idx + 3] +
                           x4 * (*(fParamHigh[Z]))[idx + 4] + x5 * (*(fParamHigh[Z]))[idx + 5]);
      for (shellIdx = 0; shellIdx < nshells; ++shellIdx)
      // for(shellIdx=0; shellIdx<nn-1; ++shellIdx)//  shellIdx=n-2 --> optimization
      {
        idx = shellIdx * 7 + 2;
        if (kinE > (*(fParamHigh[Z]))[idx - 1]) {
          double cs = (*(fParamHigh[Z]))[idx] + x1 * (*(fParamHigh[Z]))[idx + 1] + x2 * (*(fParamHigh[Z]))[idx + 2] +
                      x3 * (*(fParamHigh[Z]))[idx + 3] + x4 * (*(fParamHigh[Z]))[idx + 4] +
                      x5 * (*(fParamHigh[Z]))[idx + 5];
          if (cs >= cs0) {
            break;
          }
        }
      }
      if (shellIdx >= nshells) {
        shellIdx = nshells - 1;
      } // optimization: this can be taken out.

    }
    // (**) Low energy parameterisation
    else if (kinE >= (*(fParamLow[Z]))[0]) {
      double x1 = (geant::units::MeV) / kinE;
      double x2 = x1 * x1;
      double x3 = x2 * x1;
      double x4 = x3 * x1;
      double x5 = x4 * x1;
      int idx   = nshells * 7 - 5;
      double cs0 =
          rand * ((*(fParamLow[Z]))[idx] + x1 * (*(fParamLow[Z]))[idx + 1] + x2 * (*(fParamLow[Z]))[idx + 2] +
                  x3 * (*(fParamLow[Z]))[idx + 3] + x4 * (*(fParamLow[Z]))[idx + 4] + x5 * (*(fParamLow[Z]))[idx + 5]);
      for (shellIdx = 0; shellIdx < nshells; ++shellIdx) {
        idx = shellIdx * 7 + 2;
        if (kinE > (*(fParamHigh[Z]))[idx - 1]) {
          double cs = (*(fParamLow[Z]))[idx] + x1 * (*(fParamLow[Z]))[idx + 1] + x2 * (*(fParamLow[Z]))[idx + 2] +
                      x3 * (*(fParamLow[Z]))[idx + 3] + x4 * (*(fParamLow[Z]))[idx + 4] +
                      x5 * (*(fParamLow[Z]))[idx + 5];
          if (cs >= cs0) {
            break;
          }
        }
      }
      // this means that if I go out from the for loop without meeting the condition - I select the last shell. Then it
      // is enough to change the condition of the for loop up to  like not needed check
      if (shellIdx >= nshells) {
        shellIdx = nshells - 1;
      }
    } else {
      double cs  = rand;
      size_t idx = 0;
      // (***) Tabulated values above k-shell ionization energy
      if (kinE >= (*(fParamHigh[Z]))[1]) {
        // this is an extra-check - if energy<fCSVector[Z]->edgeMin it should go to the next else (****)
        if (kinE < fCSVector[Z]->fEdgeMin)
          cs *= fCSVector[Z]->fDataVector[0];
        else {
          idx = fCSVector[Z]->FindCSBinLocation(kinE, idx);
          cs *= fCSVector[Z]->fSplineInt->GetValueAt(kinE, idx);
        }
      }
      //(****) Tabulated values below k-shell ionization energy
      else {
        // this check is needed to have a constant cross-section(cs) value for energies below the lowest cs point.
        // in this way the gamma is always absorbed - instead of giving zero cs -
        if (kinE < fLECSVector[Z]->fEdgeMin)
          cs *= fLECSVector[Z]->fDataVector[0];
        else
          cs *= fLECSVector[Z]->GetValue(kinE, idx);
      }
      // size_t j=0;
      for (size_t j = 0; j < nshells; ++j) {
        shellIdx = (size_t)fShellVector[Z][j]->fCompID;
        if (kinE > (*(fParamLow[Z]))[7 * shellIdx + 1]) {
          size_t idx = 0;
          cs -= fShellVector[Z][j]->GetValue(kinE, idx);
        }

        if (cs <= 0.0 || j + 1 == nshells) break;
      }
    }
  } else

  {
    sampledShells = 0;
    return;
  }
  sampledShells = shellIdx;
}

void SauterGavrilaPhotoElectricModel::SampleShellAlias(double kinE, size_t &zed, double &r1, double &r2,
                                                       size_t &sampledShells)
{
  double lGammaEnergy = std::log(kinE);
  int tableIndex;
  int tableIndexBaseEn = (int)((lGammaEnergy - fShellPrimEnLMin) * fShellPrimEnILDelta); // this the lower bin
  if (tableIndexBaseEn >= fNumAliasTables - 1) {
    tableIndexBaseEn = fNumAliasTables - 2;
    std::cout << "SauterGavrilaPhotoElectricModel::SampleShellAlias: CRITICAL, this case is not handled in the "
                 "vectorized implementation.\n";
    exit(-1);
  }

  if (kinE <= fBindingEn[zed][0]) {
    // I can do the search directly with non-log energy values
    int tableIndexBinding =
        std::lower_bound(fSortedDoubledBindingEn[zed].begin(), fSortedDoubledBindingEn[zed].end(), kinE) -
        fSortedDoubledBindingEn[zed].begin(); // this is already the upper bin
    if (((size_t)tableIndexBinding < fSortedDoubledBindingEn[zed].size() - 1) &&
        fSortedDoubledBindingEn[zed][tableIndexBinding] == fSortedDoubledBindingEn[zed][tableIndexBinding - 1])
      tableIndexBinding--;
    if (fShellSamplingPrimEnergies[tableIndexBaseEn + 1] <= fSortedDoubledBindingEn[zed][tableIndexBinding]) {
      // select the Base Energy corresponding index
      tableIndex = fIndexBaseEn[zed][tableIndexBaseEn + 1] - 1; // lower bin in the complete vector
      if (lGammaEnergy > fShellLSamplingPrimEnergiesNEW[zed][tableIndex + 1] ||
          lGammaEnergy < fShellLSamplingPrimEnergiesNEW[zed][tableIndex]) {
        std::cout << "** BASE ** Error\n";
        exit(-1);
      }
    } else {
      // select the Sorted doubles binding energies corresponding index
      tableIndex = fIndexSortedDoubledBindingEn[zed][tableIndexBinding] - 1; // lower bin in the complete vector
      if (lGammaEnergy > fShellLSamplingPrimEnergiesNEW[zed][tableIndex + 1] ||
          lGammaEnergy < fShellLSamplingPrimEnergiesNEW[zed][tableIndex]) {
        std::cout << "** BINDING ** Error\n";
        exit(-1);
      }
    }

  } else {
    tableIndex = fIndexBaseEn[zed][tableIndexBaseEn + 1] - 1;
    if (lGammaEnergy > fShellLSamplingPrimEnergiesNEW[zed][tableIndex + 1] ||
        lGammaEnergy < fShellLSamplingPrimEnergiesNEW[zed][tableIndex]) {
      std::cout << "** HIGH EN **  Error\n";
      exit(-1);
    }
  }

  if ((fShellLSamplingPrimEnergiesNEW[zed][tableIndex] == fShellLSamplingPrimEnergiesNEW[zed][tableIndex + 1]) &&
      lGammaEnergy >= fShellLSamplingPrimEnergiesNEW[zed][tableIndex]) {
    tableIndex++;
    std::cout << "SauterGavrilaPhotoElectricModel::SampleShellAlias::::Attention, this check must be added to the "
                 "vectorized implementation++ "
              << tableIndex << "\n";
    std::cout << lGammaEnergy << std::endl;
    exit(-1);
  }

  double val = (lGammaEnergy - fShellPrimEnLMin) * fShellPrimEnILDelta; // To correct - perché bisogna prendere il
                                                                        // corretto delta (inverso del delta del
                                                                        // logaritmo delle energie//
  int indxEgamma   = (int)val;                                          // lower electron energy bin index
  double pIndxHigh = val - indxEgamma;
  if (r1 <= pIndxHigh) tableIndex++;

  // this has to be tranformed to the localIndex, considering the Z
  int indx      = fLastSSAliasIndex[zed - 1] + tableIndex;
  int xsampl    = fShellAliasSampler->SampleDiscrete(fShellAliasData[indx]->fAliasW, fShellAliasData[indx]->fAliasIndx,
                                                  fShellAliasData[indx]->fNumdata, r2);
  sampledShells = xsampl;
}

int SauterGavrilaPhotoElectricModel::PrepareDiscreteAlias(int Z, double ekin, std::vector<double> &x,
                                                          std::vector<double> &y)
{

  int numShell = fNShells[Z];
  x.resize(numShell);
  y.resize(numShell);
  size_t idx = 0;
  for (int i = 0; i < numShell; i++) {
    x[i] = i;
    if (ekin < fBindingEn[Z][i]) {
      y[i] = 0.;

    } else {
      y[i] = fShellVectorFull[Z][i]->GetValue(ekin, idx); // one could change this with the Exact value..
    }
  }
  return numShell;
}

void SauterGavrilaPhotoElectricModel::BuildOneDiscreteAlias(int Z, int indx, double ekin, bool &flagSpecial,
                                                            int &idSpecialShell)
{

  std::vector<double> x;
  std::vector<double> y;
  int numShell = PrepareDiscreteAlias(Z, ekin, x, y);
  if (flagSpecial) {
    y[idSpecialShell] = 0;
  }
  bool allZeros = true;
  for (int mm = 0; mm < numShell; mm++)
    if (y[mm] != 0) allZeros = false;
  fShellAliasData[indx]             = new LinAlias();
  fShellAliasData[indx]->fNumdata   = numShell;
  fShellAliasData[indx]->fXdata     = new double[numShell];
  fShellAliasData[indx]->fYdata     = new double[numShell];
  fShellAliasData[indx]->fAliasW    = new double[numShell];
  fShellAliasData[indx]->fAliasIndx = new int[numShell];
  // copy data
  for (int i = 0; i < (int)x.size(); i++) {
    fShellAliasData[indx]->fXdata[i] = x[i];
    fShellAliasData[indx]->fYdata[i] = y[i];
    if ((ekin - fBindingEn[Z][i]) < 0 && y[i] != 0) {
      std::cout << "SauterGavrilaPhotoElectricModel::BuildOneDiscreteAlias error  " << (ekin - fBindingEn[Z][i])
                << "\n";
      exit(-1);
    }
  }
  // prepare the alias data for this PDF(x,y)
  if (!allZeros) {
    fShellAliasSampler->PreparDiscreteTable(fShellAliasData[indx]->fYdata, fShellAliasData[indx]->fAliasW,
                                            fShellAliasData[indx]->fAliasIndx, fShellAliasData[indx]->fNumdata);

  } else {
    for (int mm = 0; mm < numShell; mm++) {
      fShellAliasData[indx]->fAliasW[mm]    = 0;
      fShellAliasData[indx]->fAliasIndx[mm] = -1;
    }
  }
}

std::vector<double> merge2Sorted(const std::vector<double> &left, const std::vector<double> &right)
{
  std::vector<double> output;
  std::merge(left.begin(), left.end(), right.begin(), right.end(), std::back_inserter(output));
  return output;
}

int binary_search_find_index(std::vector<double> v, double data)
{
  auto it = std::lower_bound(v.begin(), v.end(), data);
  if (it == v.end() || *it != data) {
    return -1;
  } else {
    std::size_t index = std::distance(v.begin(), it);
    return index;
  }
}

void SauterGavrilaPhotoElectricModel::InitShellSamplingTables()
{

  int Z[gMaxSizeData];
  int nTotShells = 0;
  for (int i = 0; i < gMaxSizeData; i++) {
    Z[i] = 0; // element not present in the simulation

    // Uncomment the following lines to run tests:
    //(1)PhysVecSauterGavrilaAliasShellValid
    //(2)PhysVecSauterGavrilaAliasShellBench
    //(3)PhysVecSauterGavrilaRejShellValid
    //(4)PhysVecSauterGavrilaRejShellBench
    //            Z[i]=i;
    //            nTotShells+= fNShells[i];
  }

  // Comment out the following lines to run tests:
  //(1)PhysVecSauterGavrilaAliasShellValid
  //(2)PhysVecSauterGavrilaAliasShellBench
  //(3)PhysVecSauterGavrilaRejShellValid
  //(4)PhysVecSauterGavrilaRejShellBench
  //*** START
  int numMatCuts = MaterialCuts::GetTheMaterialCutsTable().size();
  // get list of active region
  std::vector<bool> isActiveInRegion = GetListActiveRegions();
  for (int i = 0; i < numMatCuts; ++i) {
    const MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[i];
    // if this MaterialCuts belongs to a region where this model is active:
    if (isActiveInRegion[matCut->GetRegionIndex()]) {
      // get the list of elements
      const Vector_t<Element *> &theElements = matCut->GetMaterial()->GetElementVector();
      int numElems                           = theElements.size();
      for (int j = 0; j < numElems; ++j) {
        double zet      = theElements[j]->GetZ();
        int elementIndx = std::lrint(zet);
        Z[elementIndx]  = elementIndx;
        nTotShells += fNShells[elementIndx];
      }
    }
  }
  //*** END

  int oldNumGridPoints = fNumAliasTables; // fShellNumSamplingPrimEnergies;
  fShellNumSamplingPrimEnergies =
      fShellNumSamplingPrimEnergiesPerDecade * std::lrint(std::log10(fShellMaxPrimEnergy / fShellMinPrimEnergy)) + 1;
  if (fShellNumSamplingPrimEnergies < 2) {
    fShellNumSamplingPrimEnergies = 2;
  }

  fNumAliasTables = gMaxSizeData * fShellNumSamplingPrimEnergies + nTotShells * 2;

  // numElements*fShellNumSamplingPrimEnergies;
  // set up the initial gamma energy grid
  if (fShellSamplingPrimEnergies) {
    delete[] fShellSamplingPrimEnergies;
    delete[] fShellLSamplingPrimEnergies;
    fShellSamplingPrimEnergies  = nullptr;
    fShellLSamplingPrimEnergies = nullptr;
  }
  fShellSamplingPrimEnergies  = new double[fShellNumSamplingPrimEnergies];
  fShellLSamplingPrimEnergies = new double[fShellNumSamplingPrimEnergies];
  fShellPrimEnLMin            = std::log(fShellMinPrimEnergy);
  double delta        = std::log(fShellMaxPrimEnergy / fShellMinPrimEnergy) / (fShellNumSamplingPrimEnergies - 1.0);
  fShellPrimEnILDelta = 1.0 / delta;
  fShellSamplingPrimEnergies[0]                                  = fShellMinPrimEnergy;
  fShellLSamplingPrimEnergies[0]                                 = fShellPrimEnLMin;
  fShellSamplingPrimEnergies[fShellNumSamplingPrimEnergies - 1]  = fShellMaxPrimEnergy;
  fShellLSamplingPrimEnergies[fShellNumSamplingPrimEnergies - 1] = std::log(fShellMaxPrimEnergy);
  // Baseline sampling energies (equal for each element)
  std::vector<double> E;
  E.push_back(fShellMinPrimEnergy);
  for (int i = 1; i < fShellNumSamplingPrimEnergies - 1; ++i) {
    double nextE = fShellPrimEnLMin + i * delta;
    E.push_back(std::exp(nextE));
    fShellLSamplingPrimEnergies[i] = nextE;
    fShellSamplingPrimEnergies[i]  = std::exp(nextE); // TO DO Optimize
  }
  E.push_back(fShellMaxPrimEnergy);

  for (int i = 3; i < gMaxSizeData; ++i) {

    fSortedDoubledBindingEn[i].clear();
    std::vector<double> bindingEDoubled; //= fBindingEn[i];
    for (size_t k = 0; k < fBindingEn[i].size(); k++) {
      bindingEDoubled.push_back(fBindingEn[i][k]);
      fSortedDoubledBindingEn[i].push_back(fBindingEn[i][k]);
    }

    std::vector<double> sortedBindingEnergies = fBindingEn[i];
    std::sort(sortedBindingEnergies.begin(), sortedBindingEnergies.end());
    fSortedBindingEn[i] = sortedBindingEnergies;

    // sorting binding energies for the merge
    std::sort(bindingEDoubled.begin(), bindingEDoubled.end());
    std::sort(fSortedDoubledBindingEn[i].begin(), fSortedDoubledBindingEn[i].end());

    // doubling the binding energies to build different sampling tables
    bindingEDoubled                  = merge2Sorted(bindingEDoubled, sortedBindingEnergies);
    fSortedDoubledBindingEn[i]       = merge2Sorted(fSortedDoubledBindingEn[i], fSortedBindingEn[i]);
    fShellSamplingPrimEnergiesNEW[i] = merge2Sorted(E, fSortedDoubledBindingEn[i]);

    //
    // add store the log of the base energies in fShellLSamplingPrimEnergiesNEW
    for (size_t ii = 0; ii < fShellSamplingPrimEnergiesNEW[i].size(); ii++) {
      fShellLSamplingPrimEnergiesNEW[i].push_back(std::log(fShellSamplingPrimEnergiesNEW[i][ii]));
      // Store the indexes of the two vectors in the new vector fShellSamplingPrimEnergiesNEW
      for (size_t ll = 0; ll < fSortedDoubledBindingEn[i].size(); ll++) {
        if (fShellSamplingPrimEnergiesNEW[i][ii] == fSortedDoubledBindingEn[i][ll]) {
          fIndexSortedDoubledBindingEn[i].push_back(ii);
          if (fSortedDoubledBindingEn[i][ll] == fSortedDoubledBindingEn[i][ll + 1]) {
            ll++;
          }
        }
      }
    }

    for (size_t ii = 0; ii < fShellSamplingPrimEnergiesNEW[i].size(); ii++) {
      for (int ll = 0; ll < fShellNumSamplingPrimEnergies; ll++)
        if (fShellSamplingPrimEnergiesNEW[i][ii] == fShellSamplingPrimEnergies[ll]) {
          fIndexBaseEn[i].push_back(ii);
        }
    }
  }

  //
  // build the sampling tables at each primary gamma energy grid point.
  //
  // prepare the array that stores pointers to sampling data structures
  if (fShellAliasData) {
    for (int i = 0; i < oldNumGridPoints; ++i) {
      if (fShellAliasData[i]) {
        delete[] fShellAliasData[i]->fXdata;
        delete[] fShellAliasData[i]->fYdata;
        delete[] fShellAliasData[i]->fAliasW;
        delete[] fShellAliasData[i]->fAliasIndx;
        delete fShellAliasData[i];
      }
    }
    delete[] fShellAliasData;
  }
  // create new fAliasData array
  fShellAliasData = new LinAlias *[fNumAliasTables];
  for (int i = 0; i < fNumAliasTables; ++i) {
    fShellAliasData[i] = nullptr;
  }
  // create one sampling data structure at each primary gamma energy grid point:
  // -first create an AliasTable object
  if (fShellAliasSampler) {
    delete fShellAliasSampler;
  }
  fLastSSAliasIndex[2] = 0;
  // -the prepare each table one-by-one
  for (int i = 3; i < gMaxSizeData; ++i) {
    int totTablePerElement = fShellNumSamplingPrimEnergies + 2 * fNShells[i];
    fLastSSAliasIndex[i]   = fLastSSAliasIndex[i - 1] + totTablePerElement;
    for (int j = 0; j < totTablePerElement; j++) {
      int localIndex     = fLastSSAliasIndex[i - 1] + j + 1;
      bool flagSpecial   = false;
      int idSpecialShell = 0;
      if (j < totTablePerElement - 2)
        if (fShellSamplingPrimEnergiesNEW[i][j] == fShellSamplingPrimEnergiesNEW[i][j + 1]) {
          flagSpecial = true;
          std::vector<double> sortedBindingEnergies(fBindingEn[i]);
          std::sort(sortedBindingEnergies.begin(), sortedBindingEnergies.end());
          idSpecialShell = binary_search_find_index(sortedBindingEnergies, fShellSamplingPrimEnergiesNEW[i][j]);
          idSpecialShell = fNShells[i] - idSpecialShell - 1;
        }
      if (Z[i] != 0) {
        BuildOneDiscreteAlias(Z[i], localIndex, fShellSamplingPrimEnergiesNEW[i][j], flagSpecial, idSpecialShell);
      }
    }
  }
}

void SauterGavrilaPhotoElectricModel::InitSamplingTables()
{
  // set number of primary gamma energy grid points
  // keep the prev. value of primary energy grid points.
  int oldNumGridPoints = fNumSamplingPrimEnergies;
  fNumSamplingPrimEnergies =
      fNumSamplingPrimEnergiesPerDecade * std::lrint(Math::Log10(fMaxPrimEnergy / fMinPrimEnergy)) + 1;
  if (fNumSamplingPrimEnergies < 2) {
    fNumSamplingPrimEnergies = 2;
  }

  // set up the initial gamma energy grid
  if (fSamplingPrimEnergies) {
    delete[] fSamplingPrimEnergies;
    delete[] fLSamplingPrimEnergies;
    fSamplingPrimEnergies  = nullptr;
    fLSamplingPrimEnergies = nullptr;
  }
  fSamplingPrimEnergies     = new double[fNumSamplingPrimEnergies];
  fLSamplingPrimEnergies    = new double[fNumSamplingPrimEnergies];
  fPrimEnLMin               = Math::Log(fMinPrimEnergy);
  double delta              = Math::Log(fMaxPrimEnergy / fMinPrimEnergy) / (fNumSamplingPrimEnergies - 1.0);
  fPrimEnILDelta            = 1.0 / delta;
  fSamplingPrimEnergies[0]  = fMinPrimEnergy;
  fLSamplingPrimEnergies[0] = fPrimEnLMin;
  fSamplingPrimEnergies[fNumSamplingPrimEnergies - 1]  = fMaxPrimEnergy;
  fLSamplingPrimEnergies[fNumSamplingPrimEnergies - 1] = Math::Log(fMaxPrimEnergy);
  for (int i = 1; i < fNumSamplingPrimEnergies - 1; ++i) {
    fLSamplingPrimEnergies[i] = fPrimEnLMin + i * delta;
    fSamplingPrimEnergies[i]  = Math::Exp(fPrimEnLMin + i * delta);
  }
  //
  // build the sampling tables at each primary gamma energy grid point.
  //
  // prepare the array that stores pointers to sampling data structures
  if (fAliasData) {
    for (int i = 0; i < oldNumGridPoints; ++i) {
      if (fAliasData[i]) {
        delete[] fAliasData[i]->fXdata;
        delete[] fAliasData[i]->fYdata;
        delete[] fAliasData[i]->fAliasW;
        delete[] fAliasData[i]->fAliasIndx;
        delete fAliasData[i];
      }
    }
    delete[] fAliasData;
  }
  // create new fAliasData array
  fAliasData = new LinAlias *[fNumSamplingPrimEnergies];
  for (int i = 0; i < fNumSamplingPrimEnergies; ++i) {
    fAliasData[i] = nullptr;
  }
  // create one sampling data structure at each primary gamma energy grid point:
  // -first create an AliasTable object
  if (fAliasSampler) {
    delete fAliasSampler;
  }
  // -the prepare each table one-by-one
  for (int i = 0; i < fNumSamplingPrimEnergies; ++i) {
    double egamma = fSamplingPrimEnergies[i];
    double tau    = egamma / geant::units::kElectronMassC2;
    BuildOneLinAlias(i, tau);
  }
}

// This method is calculating the differential cross section in the transformed variable xsi
double SauterGavrilaPhotoElectricModel::CalculateDiffCrossSectionLog(double tau, double xsi) const
{

  // Based on Geant4 : G4SauterGavrilaAngularDistribution
  // SauterGavrila approximation for K-shell, correct to the first \alphaZ order
  // input  : energy0  (incoming photon energy)
  // input  : cosTheta (cons(theta) of photo-electron)
  // output : dsigma   (differential cross section, K-shell only)

  // double tau = energy0 / geant::units::kElectronMassC2;

  // std::cout<<"CalculateDiffCrossSectionLog. tau: "<<tau<<" and xsi: "<<xsi<<std::endl;
  double cosTheta = Math::Exp(xsi) - 2;

  // gamma and beta: Lorentz factors of the photoelectron
  double gamma = tau + 1.0;
  double beta  = std::sqrt(tau * (tau + 2.0)) / gamma;

  double z  = 1 - beta * cosTheta;
  double z2 = z * z;
  double z4 = z2 * z2;
  double y  = 1 - cosTheta * cosTheta; // sen^2(theta)

  double dsigmadcostheta = (y / z4) * (1 + 0.5 * gamma * (tau) * (gamma - 2) * z);
  double dsigmadxsi      = dsigmadcostheta * (Math::Exp(xsi));
  // std::cout<<"dsigmadcostheta: "<<dsigmadcostheta<<" and dsigmadxsi: "<<dsigmadxsi<<std::endl;
  return dsigmadxsi;
}

double SauterGavrilaPhotoElectricModel::CalculateDiffCrossSection(double tau, double cosTheta) const
{

  // Based on Geant4 : G4SauterGavrilaAngularDistribution
  // SauterGavrila approximation for K-shell, correct to the first \alphaZ order
  // input  : energy0  (incoming photon energy)
  // input  : cosTheta (cons(theta) of photo-electron)
  // output : dsigma   (differential cross section, K-shell only)

  // double tau = energy0 / geant::units::kElectronMassC2;

  // gamma and beta: Lorentz factors of the photoelectron
  double gamma = tau + 1.0;
  double beta  = std::sqrt(tau * (tau + 2.0)) / gamma;

  double z  = 1 - beta * cosTheta;
  double z2 = z * z;
  double z4 = z2 * z2;
  double y  = 1 - cosTheta * cosTheta; // sen^2(theta)

  double dsigma = (y / z4) * (1 + 0.5 * gamma * (tau) * (gamma - 2) * z);
  return dsigma;
}

int SauterGavrilaPhotoElectricModel::PrepareLinAlias(double tau, std::vector<double> &x, std::vector<double> &y)
{

  int numpoints            = 40;
  int curNumData           = 5;
  double maxErrorThreshold = gsingleTableErrorThreshold;

  x.resize(numpoints);
  y.resize(numpoints);

  // cosTheta variable between [-1, 1]
  // start with 5 points: thetaMin, thetaMin+0.1, 0, thetaMax-0.1, thetaMax
  double thetaMin = -1.;
  double thetaMax = 1.;
  x[0]            = thetaMin;
  x[1]            = thetaMin + 0.1;
  x[2]            = 0.;
  x[3]            = thetaMax - 0.1;
  x[4]            = thetaMax;
  y[0]            = CalculateDiffCrossSection(tau, x[0]);
  y[1]            = CalculateDiffCrossSection(tau, x[1]);
  y[2]            = CalculateDiffCrossSection(tau, x[2]);
  y[3]            = CalculateDiffCrossSection(tau, x[3]);
  y[4]            = CalculateDiffCrossSection(tau, x[4]);
  double maxerr   = 1.0;

  // expand the data up to the required precision level
  while (curNumData < numpoints && maxerr >= maxErrorThreshold) {
    // find the lower index of the bin, where we have the biggest linear interp. error
    maxerr         = 0.0; // value of the current maximum error
    double thexval = 0.0;
    double theyval = 0.0;
    int maxerrindx = 0; // the lower index of the corresponding bin
    for (int i = 0; i < curNumData - 1; ++i) {

      double xx  = 0.5 * (x[i] + x[i + 1]);            // mid x point
      double yy  = 0.5 * (y[i] + y[i + 1]);            // lin. interpolated pdf value at the mid point
      double val = CalculateDiffCrossSection(tau, xx); // real pdf value at the mid point
      double err = std::abs(1. - (yy / val));
      if (err > maxerr) {
        maxerr     = err;
        maxerrindx = i;
        thexval    = xx;
        theyval    = val;
      }
    }
    // extend x,y data by putting a new real value at the mid point of the highest error bin
    // first shift all values to the right
    for (int j = curNumData; j > maxerrindx + 1; --j) {
      x[j] = x[j - 1];
      y[j] = y[j - 1];
    }
    // fill x mid point
    x[maxerrindx + 1] = thexval;
    y[maxerrindx + 1] = theyval;

    // increase number of data
    ++curNumData;
    if (curNumData >= numpoints) {
      numpoints *= 2;
      x.resize(numpoints);
      y.resize(numpoints);
    }
  } // end while

  x.resize(curNumData);
  y.resize(curNumData);
  return curNumData;
}

void SauterGavrilaPhotoElectricModel::BuildOneLinAlias(int indx, double tau)
{

  std::vector<double> x;
  std::vector<double> y;
  fNumSamplingAngles = PrepareLinAlias(tau, x, y);

  fAliasData[indx]             = new LinAlias();
  fAliasData[indx]->fNumdata   = fNumSamplingAngles;
  fAliasData[indx]->fXdata     = new double[fNumSamplingAngles];
  fAliasData[indx]->fYdata     = new double[fNumSamplingAngles];
  fAliasData[indx]->fAliasW    = new double[fNumSamplingAngles];
  fAliasData[indx]->fAliasIndx = new int[fNumSamplingAngles];

  // INTEGRAL CALCULATION : is not needed, the pdf can be not normalized
  // GLIntegral  *gl   = new GLIntegral(pointsForIntegral, -1., 1.0);
  // std::vector<double> glx  = gl->GetAbscissas();
  // std::vector<double> glw  = gl->GetWeights();

  // calculate the integral of the differential cross section
  // double integral=0.;
  // for(int i = 0 ; i < pointsForIntegral ; i++)
  //    integral+= (glw[i]* CalculateDiffCrossSection(tau, glx[i]));

  // for (int i = 0; i < (int)y.size(); ++i) {
  //    y[i]/=integral;
  // }

  // copy data
  for (int i = 0; i < (int)x.size(); i++) {
    fAliasData[indx]->fXdata[i] = x[i];
    fAliasData[indx]->fYdata[i] = y[i];
  }

  // prepare the alias data for this PDF(x,y)
  fAliasSampler->PreparLinearTable(fAliasData[indx]->fXdata, fAliasData[indx]->fYdata, fAliasData[indx]->fAliasW,
                                   fAliasData[indx]->fAliasIndx, fAliasData[indx]->fNumdata);
}

} // namespace geantphysics
