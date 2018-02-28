
#include "Geant/PWATotalXsecTable.h"

#include "Geant/Material.h"
#include "Geant/Element.h"

// from material
#include "Geant/Types.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

namespace geantphysics {

const double PWATotalXsecZ::gPWATotalXsecEnergyGrid[] = {
    // energy bin values for total elastic, first and second transport cross scetions in GeV
    1.00000000e-07, 1.16591440e-07, 1.35935639e-07, 1.58489319e-07, 1.84784980e-07, 2.15443469e-07, 2.51188643e-07,
    2.92864456e-07, 3.41454887e-07, 3.98107171e-07, 4.64158883e-07, 5.41169527e-07, 6.30957344e-07, 7.35642254e-07,
    8.57695899e-07, 1.00000000e-06, 1.16591440e-06, 1.35935639e-06, 1.58489319e-06, 1.84784980e-06, 2.15443469e-06,
    2.51188643e-06, 2.92864456e-06, 3.41454887e-06, 3.98107171e-06, 4.64158883e-06, 5.41169527e-06, 6.30957344e-06,
    7.35642254e-06, 8.57695899e-06, 1.00000000e-05, 1.16591440e-05, 1.35935639e-05, 1.58489319e-05, 1.84784980e-05,
    2.15443469e-05, 2.51188643e-05, 2.92864456e-05, 3.41454887e-05, 3.98107171e-05, 4.64158883e-05, 5.41169527e-05,
    6.30957344e-05, 7.35642254e-05, 8.57695899e-05, 1.00000000e-04, 1.16591440e-04, 1.35935639e-04, 1.58489319e-04,
    1.84784980e-04, 2.15443469e-04, 2.51188643e-04, 2.92864456e-04, 3.41454887e-04, 3.98107171e-04, 4.64158883e-04,
    5.41169527e-04, 6.30957344e-04, 7.35642254e-04, 8.57695899e-04, 1.00000000e-03, 1.16591440e-03, 1.35935639e-03,
    1.58489319e-03, 1.84784980e-03, 2.15443469e-03, 2.51188643e-03, 2.92864456e-03, 3.41454887e-03, 3.98107171e-03,
    4.64158883e-03, 5.41169527e-03, 6.30957344e-03, 7.35642254e-03, 8.57695899e-03, 1.00000000e-02, 1.16591440e-02,
    1.35935639e-02, 1.58489319e-02, 1.84784980e-02, 2.15443469e-02, 2.51188643e-02, 2.92864456e-02, 3.41454887e-02,
    3.98107171e-02, 4.64158883e-02, 5.41169527e-02, 6.30957344e-02, 7.35642254e-02, 8.57695899e-02, 1.00000000e-01,
    1.16591440e-01, 1.35935639e-01, 1.58489319e-01, 1.84784980e-01, 2.15443469e-01, 2.51188643e-01, 2.92864456e-01,
    3.41454887e-01, 3.98107171e-01, 4.64158883e-01, 5.41169527e-01, 6.30957344e-01, 7.35642254e-01, 8.57695899e-01,
    1.00000000e+00};

PWATotalXsecZ::PWATotalXsecZ(int Z)
{
  int nn = gNumTotalXsecBins * 6;
  for (int i = 0; i < nn; ++i) {
    fPWAXsecs[i]     = 0.0;
    fInterpParamA[i] = 0.0;
    fInterpParamB[i] = 0.0;
  }
  LoadPWATotalXsecZ(Z);
}

void PWATotalXsecZ::LoadPWATotalXsecZ(int Z)
{
  double dum;
  char fname[512];
  char *path = getenv("GEANT_PHYSICS_DATA");
  if (!path) {
    std::cerr << "******   ERROR in PWATotalXsecZ::LoadPWATotalXsecZ() \n"
              << "         GEANT_PHYSICS_DATA is not defined! Set the GEANT_PHYSICS_DATA\n"
              << "         environmental variable to the location of Geant data directory!\n"
              << std::endl;
    exit(1);
  }
  std::string pathString(path);
  sprintf(fname, "%s/msc_GS/xsecs/xsecs_%d", path, Z);
  std::ifstream infile(fname, std::ios::in);
  if (!infile.is_open()) {
    std::string strfname(fname);
    std::cerr << "******   ERROR in PWATotalXsecZ::LoadPWATotalXsecZ() \n"
              << "         Cannot open file: " << fname << " \n"
              << std::endl;
    exit(1);
  }
  //
  double dummy;
  for (int i = 0; i < gNumTotalXsecBins; ++i) {
    for (int j = 0; j < 7; ++j) {
      if (j == 0) {
        infile >> dum;
      } else {
        // load pwa xsection that are stored in cm2 units in file and change to
        infile >> dummy;
        fPWAXsecs[(j - 1) * gNumTotalXsecBins + i] = dummy * geant::units::cm2;
      }
    }
  }
  infile.close();
  // compute log-log linear intrp. parameters
  for (int i = 0; i < gNumTotalXsecBins - 1; ++i) {
    for (int k = 0; k < 6; ++k) {
      int j            = k * gNumTotalXsecBins + i;
      double val2      = fPWAXsecs[j + 1];
      double val1      = fPWAXsecs[j];
      fInterpParamA[j] = std::log(val2 / val1) / std::log(gPWATotalXsecEnergyGrid[i + 1] / gPWATotalXsecEnergyGrid[i]);
      fInterpParamB[j] = std::exp(std::log(val1) - fInterpParamA[j] * std::log(gPWATotalXsecEnergyGrid[i]));
    }
  }
}

// Get the index of the lower energy bin edge
int PWATotalXsecZ::GetPWATotalXsecEnergyBinIndex(double energy) const
{
  // log(gPWATotalXsecEnergyGrid[0]);
  constexpr double lne0 = -1.61180956509583e+01;
  // 1./log(gPWATotalXsecEnergyGrid[i+1]/gPWATotalXsecEnergyGrid[i]);
  const double invlnde = 6.51441722854880e+00;
  return (int)((std::log(energy) - lne0) * invlnde);
}

// j-dependent type interploated cross section in Geant4 internal length2 unit
// energy is assumed to be in [GeV] which is the internal geantV energy unit
double PWATotalXsecZ::GetInterpXsec(double energy, int elowindex, int j) const
{
  // protection : out of energy grid range
  if (energy < GetLowestEnergy()) {
    return GetLowestXsecValue(j);
  }
  if (energy >= GetHighestEnergy()) {
    return GetHighestXsecValue(j);
  }
  // normal case log-log linear intrp.
  int k = j * gNumTotalXsecBins + elowindex;
  return std::exp(std::log(energy) * fInterpParamA[k]) * fInterpParamB[k];
}

////////////////////////////////////////////////////////////////////////////////
double PWATotalXsecZ::GetInterpXsec(double energy, int j) const
{
  // protection : out of energy grid range
  if (energy < GetLowestEnergy()) {
    return GetLowestXsecValue(j);
  }
  if (energy >= GetHighestEnergy()) {
    return GetHighestXsecValue(j);
  }
  // normal case log-log linear intrp.
  int elowindex = GetPWATotalXsecEnergyBinIndex(energy);
  int k         = j * gNumTotalXsecBins + elowindex;
  return std::exp(std::log(energy) * fInterpParamA[k]) * fInterpParamB[k];
}

//
//  PWATotalXsecTable
PWATotalXsecZ *PWATotalXsecTable::gPWATotalXsecTable[gNumZet] = {0};

PWATotalXsecTable::~PWATotalXsecTable()
{
  for (int i = 0; i < gNumZet; ++i) {
    if (gPWATotalXsecTable[i]) {
      delete gPWATotalXsecTable[i];
      gPWATotalXsecTable[i] = 0;
    }
  }
}

void PWATotalXsecTable::Initialise()
{
  int isUsedZ[gNumZet] = {0}; // xsec data available up to gNumZet Z-number
  // check used elements
  // get the material table (do not bother if the model is used for that or not)
  const Vector_t<Material *> &theMaterialTable = Material::GetTheMaterialTable();
  // get number of materials in the table
  unsigned int numMaterials = theMaterialTable.size();
  for (unsigned int imat = 0; imat < numMaterials; ++imat) {
    const Vector_t<Element *> &theElemVect = theMaterialTable[imat]->GetElementVector();
    size_t numelems                        = theElemVect.size();
    for (size_t ielem = 0; ielem < numelems; ++ielem) {
      int zet = std::lrint(theElemVect[ielem]->GetZ());
      zet     = zet > gNumZet ? gNumZet : zet;
      if (!isUsedZ[zet - 1]) {
        isUsedZ[zet - 1] = 1;
      }
    }
  }
  //
  for (int i = 0; i < gNumZet; ++i) {
    if (isUsedZ[i] && !gPWATotalXsecTable[i]) { // used but not there yet -> load it
      gPWATotalXsecTable[i] = new PWATotalXsecZ(i + 1);
    } else if (!isUsedZ[i] && gPWATotalXsecTable[i]) { // there but not used now -> delete
      delete gPWATotalXsecTable[i];
      gPWATotalXsecTable[i] = 0;
    }
  }
}

} // namespace geantphysics
