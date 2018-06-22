//
#include "Geant/NeutronNudyElasticModel.h"

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

NeutronNudyElasticModel::NeutronNudyElasticModel(const std::string &modelname) : HadronicFinalStateModel()
{
  this->SetType(HadronicModelType::kElastic);
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

NeutronNudyElasticModel::~NeutronNudyElasticModel()
{
}

void NeutronNudyElasticModel::Initialize()
{
  HadronicFinalStateModel::Initialize(); //
}

int NeutronNudyElasticModel::SampleFinalState(LightTrack &track, Isotope *targetisotope, geant::TaskData *td)
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
    std::cout << "NeutronNudyElasticModel: "
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
  double mass2 = targetisotope->GetIsoMass();
  double M2m   = mass2 / mass1;
  // momentum in lab frame
  double plab = std::sqrt(Ek * (Ek + 2 * mass1));
  
  double cosCm = recopoint.GetCos4(elemId, 2, Ek/geant::units::eV);
  double cosLab = (1 + M2m * cosCm) / std::sqrt(1 + M2m * M2m + 2 * M2m * cosCm);
  double sinLab = 0;
  // problem in sampling 
  if (cosLab > 1.0 || cosLab < -1.0) {
    std::cout << "GVNeutronNudyElastic WARNING cost= " << cosCm << " after scattering of " << track.GetGVcode()
    << " p(GeV/c)= " << plab / geant::units::GeV << " on an ion Z= " << targetisotope->GetZ()
    << " N= " << targetisotope->GetN() << std::endl;
    cosLab = 1.0;
    sinLab = 0.0;
    
    // normal situation
  } else {
    sinLab = std::sqrt((1.0 - cosLab) * (1.0 + cosLab));
  }
  LorentzVector<double> lv1(plab * track.GetDirX(), plab * track.GetDirY(), plab * track.GetDirZ(), Ek + mass1);
  LorentzRotation<double> fromZ;
  fromZ.rotateZ(lv1.Phi());
  fromZ.rotateY(lv1.Theta());
  
  LorentzVector<double> lv(0.0, 0.0, 0.0, mass2);
  lv += lv1;
  
  Vector3D<double> bst = lv.BoostVector();
  
  lv1.Boost(-bst);
  
  Vector3D<double> p1 = lv1.vect();
  // momentum in CMS frame
  double momentumCMS = p1.Mag();
  double phi         = td->fRndm->uniform() * 2 * geant::units::kPi;
  double sinCm       = std::sqrt((1.0 - cosCm) * (1.0 + cosCm));
  Vector3D<double> v1(sinCm * Math::Cos(phi), sinCm * Math::Sin(phi), cosCm);
  
  v1 *= momentumCMS;
  
  LorentzVector<double> nlv1(v1.x(), v1.y(), v1.z(), std::sqrt(momentumCMS * momentumCMS + mass1 * mass1));
  
  // rotation back from the Z axis
  nlv1 *= fromZ;
  
  nlv1.Boost(bst);
  Vector3D<double> newdir = nlv1.vect().Unit();
  // double eFinal1 = nlv1.e() - mass1;
  
  lv -= nlv1;
  
  // recoil energy
  double erec = lv.e() - mass2;
  // std::cout<<"lorenz direction "<< newdir.x()<<"  "<< newdir.y() <<" "<< newdir.z() <<std::endl;
  // std::cout<<"lorenz energy "<< eFinal1 <<" recoil energy "<< erec  <<std::endl;

  phi = geant::units::kTwoPi * td->fRndm->uniform();
  double sinphi, cosphi;
  Math::SinCos(phi, sinphi, cosphi);

  double nDirX = sinCm * cosphi;
  double nDirY = sinCm * sinphi;
  double nDirZ = cosCm;
  // std::cout<<"direction cm "<< nDirX <<"  "<< nDirY <<" "<< nDirZ <<std::endl;
  
   nDirX = sinLab * cosphi;
   nDirY = sinLab * sinphi;
   nDirZ = cosLab;
  //std::cout<<"direction Lab "<< nDirX <<"  "<< nDirY <<" "<< nDirZ <<std::endl;
  
  Math::RotateToLabFrame(nDirX, nDirY, nDirZ, track.GetDirX(), track.GetDirY(), track.GetDirZ());

  // std::cout<<"direction rotation "<< nDirX <<"  "<< nDirY <<" "<< nDirZ <<std::endl;
  
  int alpha = (M2m - 1) * (M2m - 1) / ((M2m + 1) * (M2m + 1));
  double ELab = (0.5 * (Ek/geant::units::eV) * ((1 - alpha) * cosCm + 1 + alpha)) * geant::units::eV;
  
  double eFinal = Ek - ELab;

  if (eFinal <= GetLowEnergyUsageLimit()) {
    if (eFinal < 0.0) {
      std::cout << "GVNeutronNudyElastic WARNING Efinal= " << eFinal << " after scattering of " << track.GetGVcode()
                << " p(GeV/c)= " << plab / geant::units::GeV << " on an ion Z= " << targetisotope->GetZ()
                << " N= " << targetisotope->GetN() << std::endl;
    }
    track.SetDirection(nDirX, nDirY, nDirZ);
    track.SetKinE(0);
  } else {
     track.SetDirection(nDirX, nDirY, nDirZ);
     track.SetKinE(ELab);
  }
  //
  if (erec > 0.0) track.SetEnergyDeposit(erec);
  // std::cout<<"energy "<< ELab <<" recoil energy "<< eFinal  <<std::endl;
  
  return numSecondaries;
}

} // namespace geantphysics
