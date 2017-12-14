//
#include "DiffuseElasticModel.h"

#include "PhysicalConstants.h"
// for Vector_t
#include "Types.h"

#include "Isotope.h"

#include "LightTrack.h"

#include <iostream>

// this should be replaced by some VecMath classes
#include "base/Lorentz.h"
#include "base/Vector3D.h"

namespace geantphysics {

  DiffuseElasticModel::DiffuseElasticModel(const std::string &modelname) : HadronicFinalStateModel() {
    this->SetType(HadronicModelType::kElastic);
    this->SetName(modelname);
    
    std::vector< int > projVec;
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
    
    this->SetProjectileCodeVec(projVec);
  }

  DiffuseElasticModel::~DiffuseElasticModel() {
  }

  void DiffuseElasticModel::Initialize() {
    HadronicFinalStateModel::Initialize();  //
  }

  
  int DiffuseElasticModel::SampleFinalState(LightTrack &track, Isotope* targetisotope,
					    Geant::GeantTaskData *td)
  {
    using vecgeom::Vector3D;
    using vecgeom::LorentzVector;
    using vecgeom::LorentzRotation;
    
    int    numSecondaries      = 0;
    double Ek                  = track.GetKinE();
    double eps0                = geant::kElectronMassC2/Ek;

    // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
    // check if kinetic energy is above fHighEnergyUsageLimit andd o nothing if yes
    if (Ek < GetLowEnergyUsageLimit() || Ek > GetHighEnergyUsageLimit() || eps0 > 0.5) {
      return numSecondaries;
    }

    /*
      std::cout << "DiffuseElasticModel: " 
      << track.GetGVcode()
      << " Ekin(MeV) = " << ekin/geant::MeV
      << " scattered off Z= " << targetisotope->GetZ() 
      << " N= " << targetisotope->GetN()
      << std::endl;
    */

    // projectile mass
    double mass1 = track.GetMass();
    // target mass
    double mass2 = targetisotope->GetIsoMass();
    // momentum in lab frame
    double plab = std::sqrt(Ek * (Ek + 2 * mass1));
    
    LorentzVector<double> lv1(plab * track.GetDirX(), plab * track.GetDirY(), plab * track.GetDirZ(), Ek + mass1);

    LorentzRotation<double> fromZ;
    fromZ.rotateZ(lv1.Phi());
    fromZ.rotateY(lv1.Theta());

    LorentzVector<double> lv(0.0,0.0,0.0,mass2);   
    lv += lv1;
     
    Vector3D<double> bst = lv.BoostVector();
    
    lv1.Boost(-bst);

    Vector3D<double> p1 = lv1.vect();
    // momentum in CMS frame
    double momentumCMS = p1.Mag();
    double tmax = 4.0 * momentumCMS * momentumCMS;
  
    // Sampling in CM system
    double t    = SampleInvariantT(mass1, plab, targetisotope, td);
    double phi  = td->fRndm->uniform() * 2 * geant::kPi;
    double cost = 1. - 2.0 * t / tmax;
    double sint;

    // problem in sampling
    if(cost > 1.0 || cost < -1.0) {
      std::cout << "G4HadronElastic WARNING (1 - cost)= " << 1 - cost
                << " after scattering of " 
                << track.GetGVcode()
                << " p(GeV/c)= " << plab/geant::GeV
                << " on an ion Z= " << targetisotope->GetZ() << " N= " << targetisotope->GetN()
                << std::endl;
      cost = 1.0;
      sint = 0.0;

      // normal situation
    } else  {
      sint = std::sqrt((1.0 - cost) * (1.0 + cost));
    }

    
    Vector3D<double> v1(sint * std::cos(phi),sint * std::sin(phi),cost);

    v1 *= momentumCMS;
    
    LorentzVector<double> nlv1(v1.x(),v1.y(),v1.z(),
			 std::sqrt(momentumCMS * momentumCMS + mass1 * mass1));

    // rotation back from the Z axis
    nlv1 *= fromZ;

    nlv1.Boost(bst); 

    double eFinal = nlv1.e() - mass1;
            
    if(eFinal <= GetLowEnergyUsageLimit()) {
      if(eFinal < 0.0) {
        std::cout << "G4HadronElastic WARNING Efinal= " << eFinal
                  << " after scattering of " 
                  << track.GetGVcode()
                  << " p(GeV/c)= " << plab/geant::GeV
                  << " on an ion Z= " << targetisotope->GetZ() << " N= " << targetisotope->GetN()
                  << std::endl;
      }
    
      nlv1 = LorentzVector<double>(0.0,0.0,0.0,mass1);

    } else {

      // rotate back to lab frame

      Vector3D<double> newdir = nlv1.vect().Unit();
      track.SetDirection(newdir.x(), newdir.y(), newdir.z());
      track.SetKinE(eFinal);
      
    }  
      
    lv -= nlv1;
    
    // recoil energy
    double erec =  lv.e() - mass2;

    //
    if(erec > 0.0) track.SetEnergyDeposit(erec);
    
    /*
      std::cout << "Recoil: " <<" m= " << mass2 << " Erec(MeV)= " << erec
      << " 4-mom: " << lv 
      << std::endl;
    */
    return numSecondaries;
  }
  

  double DiffuseElasticModel::SampleInvariantT(double mass, 
					       double plab, Isotope* targetisotope, Geant::GeantTaskData *td)
  {
    static const double GeV2 = geant::GeV * geant::GeV;

    double A = (double)targetisotope->GetN();

    double m1 = mass;
    double m12 = m1 * m1;
    double mass2 = targetisotope->GetIsoMass();
    double momentumCMS = plab * mass2/std::sqrt(m12 + mass2 * mass2 + 2. * mass2 * std::sqrt(m12 + plab * plab));
  
    double tmax = 4.0 * momentumCMS * momentumCMS / GeV2;
    
    double aa, bb, cc;
    double dd = 10.;

    if (A <= 62) {
      bb = 14.5 * std::pow(A, 2.0/3.0);
      aa = std::pow(A, 1.63)/bb;
      cc = 1.4 * std::pow(A, 1.0/3.0)/dd;
    } else {
      bb = 60. * std::pow(A, 1.0/3.0);
      aa = std::pow(A, 1.33)/bb;
      cc = 0.4 * std::pow(A, 0.4)/dd;
    }
    double q1 = 1.0 - std::exp(-bb * tmax);
    double q2 = 1.0 - std::exp(-dd * tmax);
    double s1 = q1 * aa;
    double s2 = q2 * cc;

    if((s1 + s2) * td->fRndm->uniform() < s2) {
      q1 = q2;
      bb = dd;
    }
  
    return -GeV2 * std::log(1.0 - td->fRndm->uniform() * q1)/bb;
  }


} // namespace geantphysics
