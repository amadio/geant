//
#include "Geant/NeutronNudyFissionModel.h"

#include "Geant/PhysicalConstants.h"
// for Vector_t
#include "Geant/Types.h"
#include "Geant/PhysicsData.h"

#include "Geant/Isotope.h"

#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"

#include <iostream>
#include "Geant/math_wrappers.h"

// this should be replaced by some VecMath classes
#include "base/Lorentz.h"
#include "base/Vector3D.h"

namespace geantphysics {

NeutronNudyFissionModel::NeutronNudyFissionModel(const std::string &modelname) : HadronicFinalStateModel()
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
  filename = tmp +"/neutron/nudy/n-";
  this->SetProjectileCodeVec(projVec);
}

NeutronNudyFissionModel::~NeutronNudyFissionModel()
{
}

void NeutronNudyFissionModel::Initialize()
{
  HadronicFinalStateModel::Initialize(); //
}

int NeutronNudyFissionModel::SampleFinalState(LightTrack &track, Isotope *targetisotope, geant::TaskData *td)
{
  using vecgeom::Vector3D;
  using vecgeom::LorentzVector;
  using vecgeom::LorentzRotation;

  int numSecondaries = 0;
  double Ek          = track.GetKinE();

  // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
  // check if kinetic energy is above fHighEnergyUsageLimit andd o nothing if yes
  if (Ek < GetLowEnergyUsageLimit() || Ek > GetHighEnergyUsageLimit()) {
    return numSecondaries;
  }

  /*
    std::cout << "NeutronNudyFissionModel: "
    << track.GetGVcode()
    << " Ekin(MeV) = " << Ek/geant::units::MeV
    << " scattered off Z= " << targetisotope->GetZ()
    << " N= " << targetisotope->GetN()
    << std::endl;
  */
  int Z = targetisotope->GetZ();
  int A = targetisotope->GetN();
  std::string z = std::to_string(Z);
  std::string a = std::to_string(A);
  std::string tmp = filename + z + "_" + a + ".root";
  // int particlePDG = Particle::GetParticleByInternalCode(particleCode)->GetPDGCode();
  fRENDF = tmp.c_str();
  int elemId = Z * 1000 + A;
  recopoint.GetData(elemId,fRENDF);
  // projectile mass
  double mass1 = track.GetMass();
  // target mass
  // double mass2 = targetisotope->GetIsoMass();
  // momentum in lab frame
  double plab = std::sqrt(Ek * (Ek + 2 * mass1));
  double nut = recopoint.GetNuTotal(elemId,  Ek/geant::units::eV);
  
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::poisson_distribution<int> distribution(nut);
  int number = distribution(generator);
  while (number ==0 ) {
    number = distribution(generator);
  }
  numSecondaries = number;
  for (int i = 0; i != number; ++i) {
    double cosLab = recopoint.GetCos4(elemId, 2, Ek/geant::units::eV);
    if (cosLab == -99) {
      cosLab = 2 * td->fRndm->uniform() - 1;
    }
    double sinLab = 0;
    // problem in sampling 
    if (cosLab > 1.0 || cosLab < -1.0) {
      std::cout << "GVNeutronNudyFission WARNING (1 - cost)= " << 1 - cosLab << " after Fission of " << track.GetGVcode()
      << " p(GeV/c)= " << plab / geant::units::GeV << " on an ion Z= " << targetisotope->GetZ()
      << " N= " << targetisotope->GetN() << std::endl;
      cosLab = 1.0;
      sinLab = 0.0;
      
      // normal situation
    } else {
      sinLab = std::sqrt((1.0 - cosLab) * (1.0 + cosLab));
    }
    double phi = geant::units::kTwoPi * td->fRndm->uniform();
    double sinphi, cosphi;
    Math::SinCos(phi, sinphi, cosphi);
    
    double nDirX = sinLab * cosphi;
    double nDirY = sinLab * sinphi;
    double nDirZ = cosLab;
    // rotate in the origional particle frame : This need to be checked: Harphool
    // energy of the neutron from the ENDF distributions
    double energy = recopoint.GetEnergy5(elemId, 18, Ek/geant::units::eV) * geant::units::eV;
    // killing primary track after fission
    track.SetTrackStatus(LTrackStatus::kKill);
    track.SetKinE(0.0);
    
    auto &secondarySoA = td->fPhysicsData->InsertSecondary();
    // Add neutron
//    int idx = secondarySoA.InsertTrack();
    secondarySoA.SetDirX(nDirX);
    secondarySoA.SetDirY(nDirY);
    secondarySoA.SetDirZ(nDirZ);
    secondarySoA.SetKinE(energy);
    secondarySoA.SetGVcode(track.GetGVcode());
    secondarySoA.SetMass(track.GetMass());
    secondarySoA.SetTrackIndex(track.GetTrackIndex()); // parent Track index
  }
  return numSecondaries;
}

} // namespace geantphysics
