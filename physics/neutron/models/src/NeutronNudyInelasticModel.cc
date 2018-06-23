//
#include "Geant/NeutronNudyInelasticModel.h"

#include "Geant/PhysicalConstants.h"
// for Vector_t
#include "Geant/Types.h"

#include "Geant/Isotope.h"

#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"

#include <iostream>
#include "Geant/math_wrappers.h"

// this should be replaced by some VecMath classes
#include "base/Lorentz.h"
#include "base/Vector3D.h"

namespace geantphysics {

NeutronNudyInelasticModel::NeutronNudyInelasticModel(const std::string &modelname) : HadronicFinalStateModel()
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

NeutronNudyInelasticModel::~NeutronNudyInelasticModel() {}

void NeutronNudyInelasticModel::Initialize()
{
  HadronicFinalStateModel::Initialize(); //
}

int NeutronNudyInelasticModel::SampleFinalState(LightTrack &track, Isotope *targetisotope, geant::TaskData *td)
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
   *    std::cout << "NeutronNudyInelasticModel: "
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
  double mass1 = track.GetMass();
  // target mass
  double mass2 = targetisotope->GetIsoMass();
  double M2m   = mass2 / mass1;
  // momentum in lab frame
  // double plab   = std::sqrt(Ek * (Ek + 2 * mass1));
  int idMt   = recopoint.GetElementId(elemId);
  int mtsize = recopoint.fMtValues[idMt].size();
  std::vector<double> xsec, MtNo;
  double sumXec = 0;
  double sumPar = 0;
  int MT        = 0;
  double eLab   = 0;
  double cosLab = 0;
  int LCT;
  double eCm, cosCm, sinLab;
  for (int i = 0; i != mtsize; ++i) {
    if (recopoint.fMtValues[idMt][i] != 2 && recopoint.fMtValues[idMt][i] != 18 && recopoint.fMtValues[idMt][i] != 19 &&
        recopoint.fMtValues[idMt][i] != 20 && recopoint.fMtValues[idMt][i] != 21 &&
        recopoint.fMtValues[idMt][i] != 38 && recopoint.fMtValues[idMt][i] != 102 &&
        recopoint.fMtValues[idMt][i] != 5) {
      double crs = recopoint.GetSigmaPartial(elemId, recopoint.fMtValues[idMt][i], Ek / geant::units::eV);
      sumXec += crs;
      xsec.push_back(crs);
      MtNo.push_back(recopoint.fMtValues[idMt][i]);
    }
  }
  double rnd = td->fRndm->uniform();
  for (int i = 0, xsecSize = xsec.size(); i != xsecSize; ++i) {
    sumPar += xsec[i] / sumXec;
    if (rnd <= sumPar) {
      MT = MtNo[i];
      break;
    }
  }
  if (MT > 50 && MT < 92) {
    numSecondaries = 0;
    cosCm          = recopoint.GetCos4(elemId, MT, Ek / geant::units::eV);
    LCT            = recopoint.GetCos4Lct(elemId, MT);
    eCm            = recopoint.GetEnergy5(elemId, MT, Ek / geant::units::eV) * geant::units::eV;
    if (eCm < 0) eCm = recopoint.GetEnergy5(elemId, 5, Ek / geant::units::eV) * geant::units::eV;
    switch (LCT) {
    case 2: {
      eLab   = eCm + (Ek + 2 * cosCm * (M2m + 1) * std::sqrt(Ek * eCm)) / ((M2m + 1) * (M2m + 1));
      cosLab = cosCm * std::sqrt(eCm / eLab) + std::sqrt(Ek / eLab) * (1 / (M2m + 1));
    } break;
    case 1: {
      cosLab = cosCm;
      eLab   = eCm;
    } break;
    }
    sinLab = 0;
    // problem in sampling
    if (cosLab > 1.0 || cosLab < -1.0) {
      std::cout << "GVNeutronNudyInelastic WARNING cost= " << cosLab << " after inelastic scattering of "
                << track.GetGVcode() << " E(GeV)= " << eLab << " on an ion Z= " << targetisotope->GetZ()
                << " N= " << targetisotope->GetN() << " for MT = " << MT << std::endl;
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

    Math::RotateToLabFrame(nDirX, nDirY, nDirZ, track.GetDirX(), track.GetDirY(), track.GetDirZ());
    track.SetDirection(nDirX, nDirY, nDirZ);
    track.SetKinE(eLab);
    // std::cout << MT <<"  "<< cosLab <<" "<< eLab << std::endl;
    return numSecondaries;
  }
  switch (MT) {
  case 16: {
    numSecondaries = 2;
  } break;
  case 17: {
    numSecondaries = 3;
  } break;
  case 37: {
    numSecondaries = 4;
  } break;
  default:
    break;
  }
  for (int i = 0; i != numSecondaries; ++i) {
    double cosLab = recopoint.GetCos4(elemId, MT, Ek / geant::units::eV);
    if (cosLab == -99) {
      cosLab = 2 * td->fRndm->uniform() - 1;
    }
    eLab          = recopoint.GetEnergy5(elemId, MT, Ek / geant::units::eV) * geant::units::eV;
    double sinLab = 0;
    // problem in sampling
    if (cosLab > 1.0) {
      std::cout << "GVNeutronNudyInelastic WARNING cost= " << cosLab << " after reaction of " << track.GetGVcode()
                << " E(GeV)= " << eLab << " on an ion Z= " << targetisotope->GetZ() << " N= " << targetisotope->GetN()
                << " for MT = " << MT << std::endl;
      cosLab = 1.0;
      sinLab = 0.0;
    } else if (cosLab < -1.0) {
      std::cout << "GVNeutronNudyInelastic WARNING cost= " << cosLab << " after reaction of " << track.GetGVcode()
                << " E(GeV)= " << eLab << " on an ion Z= " << targetisotope->GetZ() << " N= " << targetisotope->GetN()
                << " for MT = " << MT << std::endl;
      cosLab = -1.0;
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
    Math::RotateToLabFrame(nDirX, nDirY, nDirZ, track.GetDirX(), track.GetDirY(), track.GetDirZ());
    // energy of the neutron from the ENDF distributions
    auto &secondary = td->fPhysicsData->InsertSecondary();
    // Add neutron
    secondary.SetDirX(nDirX);
    secondary.SetDirY(nDirY);
    secondary.SetDirZ(nDirZ);
    secondary.SetKinE(eLab);
    secondary.SetGVcode(track.GetGVcode());
    secondary.SetMass(track.GetMass());
    secondary.SetTrackIndex(track.GetTrackIndex()); // parent Track index
    Ek = Ek - eLab;
    if (Ek < 0) Ek = 0;
    // std::cout << MT <<"  "<< cosLab <<" "<< Ek <<"  "<< eLab << std::endl;
  }
  // killing primary track
  track.SetTrackStatus(LTrackStatus::kKill);
  track.SetKinE(0.0);
  return numSecondaries;
}
} // namespace geantphysics
