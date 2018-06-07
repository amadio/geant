#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Geant/TNudyInterface.h"

using namespace Nudy;
using NudyPhysics::TNudyInterface;
NudyPhysics::TNudyInterface nudyxs;
geantphysics::NudyProcessType pType;

template <typename T>
std::ostream &operator<<(typename std::enable_if<std::is_enum<T>::value, std::ostream>::type &stream, const T &e)
{
  return stream << static_cast<typename std::underlying_type<T>::type>(e);
}

void DumpProcE2R()
{
  ////////////////////////////////////
  // This is an example to show dumping of ENDF data file
  // to ROOT file after computing resonance parameters, linearization etc.
  // DumpEndf2Root(std::string fIN, std::string fOUT, std::string fSUBName,
  // int tA, int tZ, double temp, std::string isotopeN)
  std::string fEndfIN    = "";
  std::string fRootOUT   = "";
  std::string fEndfSubIN = "";
  double temp            = 293.60608;
  std::cout << " ENDF data file name as INPUT (e.g. 094_Pu_240): ";
  std::getline(std::cin, fEndfIN);
  std::cout << " ROOT data file name as OUTPUT <press ENTER for automatic name>: ";
  std::getline(std::cin, fRootOUT);
  std::cout << " ENDFSUB file name < If does not exist, please ENTER <ZZ> (e.g.094_Pu_240) : ";
  std::cout << " We are selecting " << fEndfIN << std::endl;
  fEndfSubIN = fEndfIN;
  // std::getline(std::cin, fEndfSubIN);
  // std::cout << " ENDFSUB data file name for fission data as INPUT: ";
  // std::cin >> fEndfSubIN;
  std::cout << "Temperature (0 means 293.60608 default ) : ";
  std::cin >> temp;
  temp = (temp <= 0.0) ? 293.60608 : temp;
  nudyxs.DumpEndf2Root(fEndfIN, fEndfSubIN, fRootOUT, temp);
  //////////////////////////////////
}

void CalXS()
{
  int projectileCode  = 2112;
  std::string eleName = "Pu";
  int AtomicNumber    = 94;  // 4; //
  int MassNumber      = 240; //
  double temperature  = 293.60608;
  // double fTempInit, fTempOut;
  double EnergyValue = 1.0e+6; // 14.0 * geant::MeV;  // 2.0e+6;  // in terms of eV    // 1.0 * geant::MeV;

  // @brief Here we provide code for projectile say 2112 for neutron, energy of projectile say 1.0 MeV
  // @brief Then we provide temperature which is required for fission etc.
  // @brief Then we provide element name like Pu , U etc. for which cross section is required
  // @brief Z for element name
  // @brief A for Mass number for the element.

  double Eval = EnergyValue; //  /geant::eV;

  pType      = geantphysics::NudyProcessType::kElastic;
  double xse = nudyxs.GetXS(projectileCode, Eval, temperature, eleName, AtomicNumber, MassNumber, pType);
  pType      = geantphysics::NudyProcessType::kFission;
  double xsf = nudyxs.GetXS(projectileCode, Eval, temperature, eleName, AtomicNumber, MassNumber, pType);
  pType      = geantphysics::NudyProcessType::kCapture;
  double xsc = nudyxs.GetXS(projectileCode, Eval, temperature, eleName, AtomicNumber, MassNumber, pType);
  //  std::cout << "Neutron Energy = " << EnergyValue/geant::MeV << " MeV " << " Isotope "
  std::cout << "Neutron Energy = " << EnergyValue << " eV "
            << " Isotope " << AtomicNumber << "-" << eleName << "-" << MassNumber << " Elastic :: XS = " << xse
            << " Fission :: XS = " << xsf << " Capture :: XS = " << xsc << std::endl;
}

int main(int /*argc*/, char ** /*argv*/)
{
  std::cout << " Select choice:----" << std::endl
            << "\t 1. Dump ENDF data file to ROOT after processing, linearizing etc. " << std::endl
            << "\t 2. Compute cross section example for 94-Pu-240 for 1 Mev neutron." << std::endl
            << "\t 3. Quit." << std::endl
            << "\t \t ENTER CHOICE:--->";
  int cInput;
  do {
    std::cin >> cInput;
  } while (cInput < 1 || cInput > 3);

  std::cout << std::endl;
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  switch (cInput) {
  case 1:
    DumpProcE2R();
    exit(0);
  case 2:
    CalXS();
    exit(0);
  case 3:
    exit(0);
  }

  return 0;
}