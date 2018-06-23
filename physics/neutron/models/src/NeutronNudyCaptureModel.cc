//
#include "Geant/NeutronNudyCaptureModel.h"

#include "Geant/PhysicalConstants.h"
// for Vector_t
#include "Geant/Types.h"

#include "Geant/Isotope.h"

#include "Geant/LightTrack.h"

#include <iostream>
#include "Geant/math_wrappers.h"

// this should be replaced by some VecMath classes
#include "base/Lorentz.h"
#include "base/Vector3D.h"

namespace geantphysics {

NeutronNudyCaptureModel::NeutronNudyCaptureModel(const std::string &modelname) : HadronicFinalStateModel()
{
  this->SetType(HadronicModelType::kNeutronNudy);
  this->SetName(modelname);

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
    std::cerr << "******   ERROR in NeutronNudyXsec() \n"
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

NeutronNudyCaptureModel::~NeutronNudyCaptureModel() {}

void NeutronNudyCaptureModel::Initialize()
{
  HadronicFinalStateModel::Initialize(); //
}

int NeutronNudyCaptureModel::SampleFinalState(LightTrack &track, Isotope *targetisotope, geant::TaskData * /*td*/)
{
  using vecgeom::LorentzRotation;
  using vecgeom::LorentzVector;
  using vecgeom::Vector3D;

  int numSecondaries = 0;
  double Ek          = track.GetKinE();

  // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
  // check if kinetic energy is above fHighEnergyUsageLimit andd o nothing if yes
  if (Ek < GetLowEnergyUsageLimit() || Ek > GetHighEnergyUsageLimit()) {
    return numSecondaries;
  }

  /*
   *    std::cout << "NeutronNudyFissionModel: "
   *    << track.GetGVcode()
   *    << " Ekin(MeV) = " << Ek/geant::units::MeV
   *    << " scattered off Z= " << targetisotope->GetZ()
   *    << " N= " << targetisotope->GetN()
   *    << std::endl;
   */
  int Z           = targetisotope->GetZ();
  int A           = targetisotope->GetN();
  std::string z   = std::to_string(Z);
  std::string a   = std::to_string(A);
  std::string tmp = filename + z + "_" + a + ".root";
  // int particlePDG = Particle::GetParticleByInternalCode(particleCode)->GetPDGCode();
  fRENDF     = tmp.c_str();
  int elemId = Z * 1000 + A;
  recopoint.GetData(elemId, fRENDF);
  // projectile mass
  // double mass1 = track.GetMass();
  // target mass
  // double mass2 = targetisotope->GetIsoMass();
  // momentum in lab frame
  // double plab = std::sqrt(Ek * (Ek + 2 * mass1));
  // killing primary track after Capture
  track.SetTrackStatus(LTrackStatus::kKill);
  track.SetKinE(0.0);
  return numSecondaries;
}
} // namespace geantphysics
