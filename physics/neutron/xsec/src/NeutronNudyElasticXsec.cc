#include "Geant/Proton.h"
#include "Geant/Neutron.h"
#include "Geant/SystemOfUnits.h"
#include "Geant/PhysicalConstants.h"
#include "Geant/math_wrappers.h"
#include "Geant/TNudyEndfRecoPoint.h"
#include "Geant/NeutronNudyElasticXsec.h"

#include <cmath>
#include <iostream>

namespace geantphysics {

NeutronNudyElasticXsec::NeutronNudyElasticXsec()
{
  this->SetName("NeutronNudyElasticXsec");

  std::vector<int> projVec;
  projVec.push_back(1);
  projVec.push_back(3);
  projVec.push_back(10);
  projVec.push_back(11);
  projVec.push_back(12);
  projVec.push_back(13);
  projVec.push_back(14);
  projVec.push_back(15);
  projVec.push_back(16);
  projVec.push_back(17);
  char *path = std::getenv("GEANT_PHYSICS_DATA");
  if (!path) {
    std::cerr << "******   ERROR in NeutronNudyElasticXsec() \n"
              << "         GEANT_PHYSICS_DATA is not defined! Set the GEANT_PHYSICS_DATA\n"
              << "         environmental variable to the location of Geant data directory\n"
              << "         It should be .root file processed from ENDF data using EndfToPointRoot\n"
              << "         executable. For more details see EndfToRoot in \n"
              << "         physics/neutron/nudy/EndfToRoot/README\n"
              << "         root file name format is n-Z_A.root!\n"
              << std::endl;
    exit(1);
  }
  std::string tmp = path;
  filename        = tmp + "/neutron/nudy/n-";
  this->SetProjectileCodeVec(projVec);
}

NeutronNudyElasticXsec::~NeutronNudyElasticXsec() {}

double NeutronNudyElasticXsec::GetIsotopeCrossSection(const int /*particleCode*/, const double energyKin,
                                                      const double /*mass*/, const int Z, const int A)
{
  double xsection;
  std::string z   = std::to_string(Z);
  std::string a   = std::to_string(A);
  std::string tmp = filename + z + "_" + a + ".root";
  // std::cout<<"file name "<< tmp << " mass "<< mass << std::endl;
  fRENDF     = tmp.c_str();
  int elemId = Z * 1000 + A;
  recopoint.GetData(elemId, fRENDF);
  xsection = recopoint.GetSigmaPartial(elemId, 2, energyKin / geant::units::eV);
  return xsection * geant::units::barn;
}
} // namespace geantphysics
