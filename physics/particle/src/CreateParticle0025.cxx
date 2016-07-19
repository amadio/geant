#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize off
#endif
#include "Particle.h"
#ifdef GEANT_NVCC
#include "base/Vector.h"
template <typename T>
using vector = vecgeom::Vector<T>;
#else
using std::vector;
#endif
namespace geant {
   inline namespace GEANT_IMPL_NAMESPACE {


//________________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void CreateParticle0025() {

   // Creating N(1900)+
   new Particle("N(1900)+", 42124, 1, "Unknown", 100, 1, 1.9, 0.5, -100, 0, -100, -1, -1);

   // Creating N(1710)+
   new Particle("N(1710)+", 42212, 1, "Unknown", 100, 1, 1.71, 0.1, -100, 0, -100, -1, -1);

   // Creating lambda(1800)
   new Particle("lambda(1800)", 43122, 1, "Unknown", 100, 0, 1.8, 0.3, -100, 0, -100, -1, -1);

   // Creating N(2090)0
   new Particle("N(2090)0", 52114, 1, "Unknown", 100, 0, 2.08, 0.35, -100, 0, -100, -1, -1);

   // Creating N(2090)+
   new Particle("N(2090)+", 52214, 1, "Unknown", 100, 1, 2.08, 0.35, -100, 0, -100, -1, -1);

   // Creating lambda(1810)
   new Particle("lambda(1810)", 53122, 1, "Unknown", 100, 0, 1.81, 0.15, -100, 0, -100, -1, -1);

   // Creating pi(1300)0
   new Particle("pi(1300)0", 100111, 1, "Unknown", 100, 0, 1.3, 0.4, -100, 0, -100, -1, -1);

   // Creating rho(1450)0
   new Particle("rho(1450)0", 100113, 1, "Unknown", 100, 0, 1.465, 0.4, -100, 0, -100, -1, -1);

   // Creating pi(1300)+
   new Particle("pi(1300)+", 100211, 1, "Unknown", 100, 1, 1.3, 0.4, -100, 0, -100, -1, -1);

   // Creating rho(1450)+
   new Particle("rho(1450)+", 100213, 1, "Unknown", 100, 1, 1.465, 0.4, -100, 0, -100, -1, -1);

   // Creating eta(1295)
   new Particle("eta(1295)", 100221, 1, "Unknown", 100, 0, 1.294, 0.055, -100, 0, -100, -1, -1);

   // Creating omega(1420)
   new Particle("omega(1420)", 100223, 1, "Unknown", 100, 0, 1.425, 0.215, -100, 0, -100, -1, -1);

   // Creating k(1460)0
   new Particle("k(1460)0", 100311, 1, "Unknown", 100, 0, 1.46, 0.26, -100, 0, -100, -1, -1);

   // Creating k_star(1410)0
   new Particle("k_star(1410)0", 100313, 1, "Unknown", 100, 0, 1.414, 0.232, -100, 0, -100, -1, -1);

   // Creating k2_star(1980)0
   new Particle("k2_star(1980)0", 100315, 1, "Unknown", 100, 0, 1.973, 0.373, -100, 0, -100, -1, -1);

   // Creating k(1460)+
   new Particle("k(1460)+", 100321, 1, "Unknown", 100, 1, 1.46, 0.26, -100, 0, -100, -1, -1);

   // Creating k_star(1410)+
   new Particle("k_star(1410)+", 100323, 1, "Unknown", 100, 1, 1.414, 0.232, -100, 0, -100, -1, -1);

   // Creating k2_star(1980)+
   new Particle("k2_star(1980)+", 100325, 1, "Unknown", 100, 1, 1.973, 0.373, -100, 0, -100, -1, -1);

   // Creating eta(1475)
   new Particle("eta(1475)", 100331, 1, "Unknown", 100, 0, 1.476, 0.085, -100, 0, -100, -1, -1);

   // Creating phi(1680)
   new Particle("phi(1680)", 100333, 1, "Unknown", 100, 0, 1.68, 0.15, -100, 0, -100, -1, -1);

   // Creating psi'
   new Particle("psi'", 100443, 0, "Unknown", 100, 0, 3.68609, 0, -100, -1, -100, -1, -1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(100443));
   part->AddDecay(Particle::Decay(0, 0.324,  vector<int>{443,211,-211}));
   part->AddDecay(Particle::Decay(12, 0.1866,  vector<int>{82,-82}));
   part->AddDecay(Particle::Decay(0, 0.184,  vector<int>{443,111,111}));
   part->AddDecay(Particle::Decay(0, 0.093,  vector<int>{10441,22}));
   part->AddDecay(Particle::Decay(0, 0.087,  vector<int>{20443,22}));
   part->AddDecay(Particle::Decay(0, 0.078,  vector<int>{445,22}));
   part->AddDecay(Particle::Decay(0, 0.027,  vector<int>{443,221}));
   part->AddDecay(Particle::Decay(0, 0.0083,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(0, 0.0083,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(0, 0.0028,  vector<int>{441,22}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{443,111}));

   // Creating Upsilon'
   new Particle("Upsilon'", 100553, 0, "Unknown", 100, 0, 10.0233, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(100553));
   part->AddDecay(Particle::Decay(4, 0.425,  vector<int>{21,21,21}));
   part->AddDecay(Particle::Decay(0, 0.185,  vector<int>{553,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.088,  vector<int>{553,111,111}));
   part->AddDecay(Particle::Decay(0, 0.067,  vector<int>{20553,22}));
   part->AddDecay(Particle::Decay(0, 0.066,  vector<int>{555,22}));
   part->AddDecay(Particle::Decay(0, 0.043,  vector<int>{10551,22}));
   part->AddDecay(Particle::Decay(32, 0.024,  vector<int>{2,-2}));
   part->AddDecay(Particle::Decay(32, 0.024,  vector<int>{4,-4}));
   part->AddDecay(Particle::Decay(4, 0.02,  vector<int>{22,21,21}));
   part->AddDecay(Particle::Decay(0, 0.014,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(0, 0.014,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(0, 0.014,  vector<int>{15,-15}));
   part->AddDecay(Particle::Decay(32, 0.008,  vector<int>{1,-1}));
   part->AddDecay(Particle::Decay(32, 0.008,  vector<int>{3,-3}));

   // Creating ~d_L
   new Particle("~d_L", 1000001, 1, "Sparticle", 100, -0.333333, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000001));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000002,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000002,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,1}));

   // Creating ~u_L
   new Particle("~u_L", 1000002, 1, "Sparticle", 100, 0.666667, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000002));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000001,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000001,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,2}));

   // Creating ~s_L
   new Particle("~s_L", 1000003, 1, "Sparticle", 100, -0.333333, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000003));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000004,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000004,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,3}));

   // Creating ~c_L
   new Particle("~c_L", 1000004, 1, "Sparticle", 100, 0.666667, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000004));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000003,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000003,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,4}));

   // Creating ~b_1
   new Particle("~b_1", 1000005, 1, "Sparticle", 100, -0.333333, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000005));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,5}));

   // Creating ~t_1
   new Particle("~t_1", 1000006, 1, "Sparticle", 100, 0.666667, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000006));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000005,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000005,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-15,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,16,5}));

   // Creating ~e_L-
   new Particle("~e_L-", 1000011, 1, "Sparticle", 100, -1, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000011));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000012,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000012,-37}));

   // Creating ~nu_eL
   new Particle("~nu_eL", 1000012, 1, "Sparticle", 100, 0, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000012));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000011,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000011,37}));
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
