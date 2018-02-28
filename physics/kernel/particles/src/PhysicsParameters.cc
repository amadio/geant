
#include "Geant/PhysicsParameters.h"

#include "Geant/SystemOfUnits.h"

#include <iostream>
#include <iomanip>
#include <cmath>

namespace geantphysics {

std::vector<PhysicsParameters *> PhysicsParameters::gThePhysicsParametersTable;

// set to default values
double PhysicsParameters::fMinAllowedGammaCutEnergy    = 990.0 * geant::units::eV;
double PhysicsParameters::fMaxAllowedGammaCutEnergy    = 10.0 * geant::units::GeV;
double PhysicsParameters::fMinAllowedElectronCutEnergy = 990.0 * geant::units::eV;
double PhysicsParameters::fMaxAllowedElectronCutEnergy = 10.0 * geant::units::GeV;
double PhysicsParameters::fMinAllowedPositronCutEnergy = 990.0 * geant::units::eV;
double PhysicsParameters::fMaxAllowedPositronCutEnergy = 10.0 * geant::units::GeV;

double PhysicsParameters::fDefaultGammaCutInLength    = 1.0 * geant::units::mm;
double PhysicsParameters::fDefaultElectronCutInLength = 1.0 * geant::units::mm;
double PhysicsParameters::fDefaultPositronCutInLength = 1.0 * geant::units::mm;

double PhysicsParameters::fDefaultGammaCutInEnergy    = 100.0 * geant::units::keV;
double PhysicsParameters::fDefaultElectronCutInEnergy = 100.0 * geant::units::keV;
double PhysicsParameters::fDefaultPositronCutInEnergy = 100.0 * geant::units::keV;

PhysicsParameters::PhysicsParameters()
{
  fMinLossTableEnergy        = 100.0 * geant::units::eV;
  fMaxLossTableEnergy        = 100.0 * geant::units::TeV;
  fNumLossTableBins          = 84;
  fNumLossTableBinsPerDecade = 7;

  fIsComputeCSDARange = false;

  fMinLambdaTableEnergy        = fMinLossTableEnergy;
  fMaxLambdaTableEnergy        = fMaxLossTableEnergy;
  fNumLambdaTableBins          = fNumLossTableBins;
  fNumLambdaTableBinsPerDecade = fNumLossTableBinsPerDecade;

  fLowestElectronTrackingEnergy = 1.0 * geant::units::keV;

  fLinearEnergyLossLimit = 0.01;

  fDRoverRange = 0.2;
  fFinalRange  = 1.0 * geant::units::mm;

  gThePhysicsParametersTable.push_back(this);
}

const PhysicsParameters *PhysicsParameters::GetPhysicsParametersForRegion(int regionindx)
{
  PhysicsParameters *physPar = nullptr;
  for (unsigned i = 0; i < gThePhysicsParametersTable.size(); ++i) {
    if (gThePhysicsParametersTable[i]->IsActiveRegion(regionindx)) {
      physPar = gThePhysicsParametersTable[i];
      break;
    }
  }
  return physPar;
}

void PhysicsParameters::Clear()
{
  for (unsigned i = 0; i < gThePhysicsParametersTable.size(); ++i) {
    delete gThePhysicsParametersTable[i];
  }
  gThePhysicsParametersTable.clear();
}

void PhysicsParameters::SetMinLossTableEnergy(double val)
{
  if (val > 1.e-3 * geant::units::eV && val < fMaxLossTableEnergy) {
    fMinLossTableEnergy = val;
    // update number of bins
    fNumLossTableBins = fNumLossTableBinsPerDecade * std::lrint(std::log10(fMaxLossTableEnergy / fMinLossTableEnergy));
  } else {
    std::cerr << "  **** WARNING: PhysicsParameters::SetMinLossTableEnergy() " << std::endl
              << "    Value of fMinLossTableEnergy is out of range: " << val / geant::units::keV
              << " [keV] so it's ignored!" << std::endl;
  }
}

void PhysicsParameters::SetMaxLossTableEnergy(double val)
{
  if (val > fMinLossTableEnergy && val < 1.e+7 * geant::units::TeV) {
    fMaxLossTableEnergy = val;
    // update number of bins
    fNumLossTableBins = fNumLossTableBinsPerDecade * std::lrint(std::log10(fMaxLossTableEnergy / fMinLossTableEnergy));
  } else {
    std::cerr << "  **** WARNING: PhysicsParameters::SetMaxLossTableEnergy() " << std::endl
              << "    Value of fMaxLossTableEnergy is out of range: " << val / geant::units::GeV
              << " [GeV] so it's ignored!" << std::endl;
  }
}

// minimum 5 for the better constrained splie
void PhysicsParameters::SetNumLossTableBins(int val)
{
  if (val > 4 && val < 10000000) {
    fNumLossTableBins          = val;
    fNumLossTableBinsPerDecade = std::lrint(fNumLossTableBins / std::log10(fMaxLossTableEnergy / fMinLossTableEnergy));
  } else {
    std::cerr << "  **** WARNING: PhysicsParameters::SetNumLossTableBins() " << std::endl
              << "    Value of fNumLossTableBins is out of range: " << val << " so it's ignored!" << std::endl;
  }
}

// minimum 5 for the better constrained splie
void PhysicsParameters::SetNumLossTableBinsPerDecade(int val)
{
  if (val > 4 && val < 1000000) {
    fNumLossTableBinsPerDecade = val;
    fNumLossTableBins = fNumLossTableBinsPerDecade * std::lrint(std::log10(fMaxLossTableEnergy / fMinLossTableEnergy));
  } else {
    std::cerr << "  **** WARNING: PhysicsParameters::SetNumLossTableBinsPerDecade() " << std::endl
              << "    Value of fNumLossTableBinsPerDecade is out of range: " << val << " so it's ignored!" << std::endl;
  }
}

////
void PhysicsParameters::SetMinLambdaTableEnergy(double val)
{
  if (val > 1.e-3 * geant::units::eV && val < fMaxLambdaTableEnergy) {
    fMinLambdaTableEnergy = val;
    // update number of bins
    fNumLambdaTableBins =
        fNumLambdaTableBinsPerDecade * std::lrint(std::log10(fMaxLambdaTableEnergy / fMinLambdaTableEnergy));
  } else {
    std::cerr << "  **** WARNING: PhysicsParameters::SetMinLambdaTableEnergy() " << std::endl
              << "    Value of fMinLambdaTableEnergy is out of range: " << val / geant::units::keV
              << " [keV] so it's ignored!" << std::endl;
  }
}

void PhysicsParameters::SetMaxLambdaTableEnergy(double val)
{
  if (val > fMinLambdaTableEnergy && val < 1.e+7 * geant::units::TeV) {
    fMaxLambdaTableEnergy = val;
    // update number of bins
    fNumLambdaTableBins =
        fNumLambdaTableBinsPerDecade * std::lrint(std::log10(fMaxLambdaTableEnergy / fMinLambdaTableEnergy));
  } else {
    std::cerr << "  **** WARNING: PhysicsParameters::SetMaxLambdaTableEnergy() " << std::endl
              << "    Value of fMaxLambdaTableEnergy is out of range: " << val / geant::units::GeV
              << " [GeV] so it's ignored!" << std::endl;
  }
}

// minimum 5 for the better constrained splie
void PhysicsParameters::SetNumLambdaTableBins(int val)
{
  if (val > 4 && val < 10000000) {
    fNumLambdaTableBins = val;
    fNumLambdaTableBinsPerDecade =
        std::lrint(fNumLambdaTableBins / std::log10(fMaxLambdaTableEnergy / fMinLambdaTableEnergy));
  } else {
    std::cerr << "  **** WARNING: PhysicsParameters::SetNumLambdaTableBins() " << std::endl
              << "    Value of fNumLambdaTableBins is out of range: " << val << " so it's ignored!" << std::endl;
  }
}

// minimum 5 for the better constrained splie
void PhysicsParameters::SetNumLambdaTableBinsPerDecade(int val)
{
  if (val > 4 && val < 1000000) {
    fNumLambdaTableBinsPerDecade = val;
    fNumLambdaTableBins =
        fNumLambdaTableBinsPerDecade * std::lrint(std::log10(fMaxLambdaTableEnergy / fMinLambdaTableEnergy));
  } else {
    std::cerr << "  **** WARNING: PhysicsParameters::SetNumLambdaTableBinsPerDecade() " << std::endl
              << "    Value of fNumLambdaTableBinsPerDecade is out of range: " << val << " so it's ignored!"
              << std::endl;
  }
}

void PhysicsParameters::SetLowestElectronTrackingEnergy(double val)
{
  if (val >= 0.0) {
    fLowestElectronTrackingEnergy = val;
  } else {
    std::cerr << "  **** WARNING: PhysicsParameters::SetLowestElectronTrackingEnergy " << std::endl
              << "    Value of fLowestElectronTrackingEnergy is out of range: " << val / geant::units::keV
              << " [keV] so it's ignored!" << std::endl;
  }
}

void PhysicsParameters::SetLinearEnergyLossLimit(double val)
{
  if (val > 0.0 && val < 0.5) {
    fLinearEnergyLossLimit = val;
  } else {
    std::cerr << "  **** WARNING: PhysicsParameters::SetLinearEnergyLossLimit " << std::endl
              << "    Value of fLinearEnergyLossLimit is out of range: " << val << " so it's ignored!" << std::endl;
  }
}

void PhysicsParameters::SetStepFunction(double roverrange, double finalrange)
{
  if (roverrange > 0.0 && roverrange <= 1.0 && finalrange > 0.0) {
    fDRoverRange = roverrange;
    fFinalRange  = finalrange;
  } else {
    std::cerr << "  **** WARNING: PhysicsParameters::SetStepFunction " << std::endl
              << "    Values of step function are out of range: " << roverrange << ", " << finalrange / geant::units::mm
              << " [mm]  so they are ignored!" << std::endl;
  }
}

std::ostream &operator<<(std::ostream &flux, PhysicsParameters *physpar)
{
  flux << *physpar;
  return flux;
}

std::ostream &operator<<(std::ostream &flux, PhysicsParameters &physpar)
{
  using geant::units::GeV;
  using geant::units::cm;
  //  using geant::units::eV;
  int textW = 70;
  int valW  = 12;

  flux << std::setw(100) << std::setfill('=') << std::left << "  =" << std::endl
       << std::setw(40) << std::setfill('-') << std::left << "  -"
       << " Physics parameters " << std::setw(40) << std::setfill('-') << std::left << "-" << std::endl
       << std::setw(100) << std::setfill('=') << std::left << "  =" << std::setfill(' ') << std::endl;

  // parameters that are the same in each region
  flux << "    **** (Static)Parameters that must have the same value in each region: \n"
       << "     - fMinAllowedGammaCutEnergy    =   " << std::setw(10) << physpar.GetMinAllowedGammaCutEnergy() / GeV
       << " [GeV]\n"
       << "     - fMaxAllowedGammaCutEnergy    =   " << std::setw(10) << physpar.GetMaxAllowedGammaCutEnergy() / GeV
       << " [GeV]\n"
       << "     - fMinAllowedElectronCutEnergy =   " << std::setw(10) << physpar.GetMinAllowedElectronCutEnergy() / GeV
       << " [GeV]\n"
       << "     - fMaxAllowedElectronCutEnergy =   " << std::setw(10) << physpar.GetMaxAllowedElectronCutEnergy() / GeV
       << " [GeV]\n"
       << "     - fMinAllowedPositronCutEnergy =   " << std::setw(10) << physpar.GetMinAllowedPositronCutEnergy() / GeV
       << " [GeV]\n"
       << "     - fMaxAllowedPositronCutEnergy =   " << std::setw(10) << physpar.GetMaxAllowedPositronCutEnergy() / GeV
       << " [GeV]\n"
       << "     - fDefaultGammaCutInLength     =   " << std::setw(10) << physpar.GetDefaultGammaCutInLength() / cm
       << "  [cm]\n"
       << "     - fDefaultElectronCutInLength  =   " << std::setw(10) << physpar.GetDefaultElectronCutInLength() / cm
       << "  [cm]\n"
       << "     - fDefaultPositronCutInLength  =   " << std::setw(10) << physpar.GetDefaultPositronCutInLength() / cm
       << "  [cm]\n"
       << "     - fDefaultGammaCutInEnergy     =   " << std::setw(10) << physpar.GetDefaultGammaCutInEnergy() / GeV
       << " [GeV]\n"
       << "     - fDefaultElectronCutInEnergy  =   " << std::setw(10) << physpar.GetDefaultElectronCutInEnergy() / GeV
       << " [GeV]\n"
       << "     - fDefaultPositronCutInEnergy  =   " << std::setw(10) << physpar.GetDefaultPositronCutInEnergy() / GeV
       << " [GeV]\n\n";

  flux << "    **** Active in regions: ";
  for (unsigned i = 0; i < physpar.GetListActiveRegions().size(); ++i) {
    if (physpar.IsActiveRegion(i)) flux << std::setw(10) << i;
  }
  flux << "\n\n";

  flux << "    **** Additional parameter values: \n\n";
  flux << std::setw(90) << std::left << std::setfill('-') << "     --- Energy loss table related parameter values "
       << std::setfill(' ') << "\n";
  flux << std::left << std::setw(textW) << "     - Minimum kinetic energy of the energy loss table grid " << std::right
       << std::setw(valW) << physpar.GetMinLossTableEnergy() / geant::units::keV << "   [keV] \n"
       << std::left << std::setw(textW) << "     - Maximum kinetic energy of the energy loss table grid " << std::right
       << std::setw(valW) << physpar.GetMaxLossTableEnergy() / geant::units::TeV << "   [TeV] \n"
       << std::left << std::setw(textW) << "     - Number of energy bins in the energy loss table " << std::right
       << std::setw(valW) << physpar.GetNumLossTableBins() << "\n"
       << std::left << std::setw(textW) << "     - Number of energy bins per decade in the energy loss table "
       << std::right << std::setw(valW) << physpar.GetNumLossTableBinsPerDecade() << "\n"
       << std::endl;

  flux << std::setw(90) << std::left << std::setfill('-') << "     --- Lambda table related parameter values "
       << std::setfill(' ') << "\n";
  flux << std::left << std::setw(textW) << "     - Minimum kinetic energy of the lambda table grid " << std::right
       << std::setw(valW) << physpar.GetMinLambdaTableEnergy() / geant::units::keV << "   [keV] \n"
       << std::left << std::setw(textW) << "     - Maximum kinetic energy of the lambda table grid " << std::right
       << std::setw(valW) << physpar.GetMaxLambdaTableEnergy() / geant::units::TeV << "   [TeV] \n"
       << std::left << std::setw(textW) << "     - Number of energy bins in the lambda table " << std::right
       << std::setw(valW) << physpar.GetNumLambdaTableBins() << "\n"
       << std::left << std::setw(textW) << "     - Number of energy bins per decade in the lambda table " << std::right
       << std::setw(valW) << physpar.GetNumLambdaTableBinsPerDecade() << "\n"
       << std::endl;

  flux << std::setw(90) << std::left << std::setfill('-') << "     --- Electron/positron related parameters "
       << std::setfill(' ') << "\n";
  flux << std::left << std::setw(textW) << "     - Lowest e-/e+ kinetic energy for tracking  " << std::right
       << std::setw(valW) << physpar.GetLowestElectronTrackingEnergy() / geant::units::keV << "   [keV] \n"
       << std::endl;

  flux << std::setw(90) << std::left << std::setfill('-') << "     --- Other parameters " << std::setfill(' ') << "\n";
  flux << std::left << std::setw(textW) << "     - Compute CSDA range tables  " << std::right << std::setw(valW)
       << physpar.GetIsComputeCSDARange() << "\n"
       << std::left << std::setw(textW) << "     - Linear energy loss approximation limit  " << std::right
       << std::setw(valW) << physpar.GetLinearEnergyLossLimit() << "\n"
       << std::left << std::setw(textW) << "     - Continuous step limit function parameters : fDRoverRange "
       << std::right << std::setw(valW) << physpar.GetDRoverRange() << "\n"
       << std::left << std::setw(textW) << "     - Continuous step limit function parameters : fFInalRange "
       << std::right << std::setw(valW) << physpar.GetFinalRange() / geant::units::mm << "   [mm] " << std::endl;

  flux << std::setw(100) << std::setfill('=') << std::left << "  =" << std::setfill(' ') << std::endl << std::endl;

  return flux;
}

} // namespace geantphysics
