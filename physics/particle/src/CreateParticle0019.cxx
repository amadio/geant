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
void CreateParticle0019() {

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

   // Creating ~mu_L-
   new Particle("~mu_L-", 1000013, 1, "Sparticle", 100, -1, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000013));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000014,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000014,-37}));

   // Creating ~nu_muL
   new Particle("~nu_muL", 1000014, 1, "Sparticle", 100, 0, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000014));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000013,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000013,37}));

   // Creating ~tau_1-
   new Particle("~tau_1-", 1000015, 1, "Sparticle", 100, -1, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000015));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000016,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000016,-37}));

   // Creating ~nu_tauL
   new Particle("~nu_tauL", 1000016, 1, "Sparticle", 100, 0, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000016));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000015,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000015,37}));

   // Creating ~g
   new Particle("~g", 1000021, 0, "Sparticle", 100, 0, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000021));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,21}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000001,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000001,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000002,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000002,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000003,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000003,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000004,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000004,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000005,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000005,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,4,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,6,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,4,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,6,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,4,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,6,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,4,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,6,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,5,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-5,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,5,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-5,6}));

   // Creating ~chi_10
   new Particle("~chi_10", 1000022, 0, "Sparticle", 100, 0, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000022));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,22}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{4,-1,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1,-3,12}));

   // Creating ~chi_20
   new Particle("~chi_20", 1000023, 0, "Sparticle", 100, 0, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000023));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,22}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,22}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,11,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,13,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,15,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,12,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,14,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,16,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,4,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000001,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000001,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000002,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000002,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000003,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000003,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000004,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000004,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000005,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000005,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000011,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000011,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000012,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000012,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000013,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000013,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000014,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000014,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000014,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000015,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000015,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000016,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000016,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,4,-4}));

   // Creating ~chi_1+
   new Particle("~chi_1+", 1000024, 1, "Sparticle", 100, 1, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000024));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000002,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000001,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000004,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000003,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000005,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000012,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000011,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000014,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000013,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000016,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000015,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,-3,4}));

   // Creating ~chi_30
   new Particle("~chi_30", 1000025, 0, "Sparticle", 100, 0, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000025));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,22}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,22}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,11,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,13,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,15,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,12,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,14,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,16,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,4,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,22}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,11,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,13,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,15,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,12,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,14,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,16,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,4,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000001,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000001,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000002,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000002,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000003,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000003,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000004,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000004,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000005,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000005,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000011,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000011,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000012,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000012,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000013,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000013,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000014,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000014,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000014,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000015,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000015,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000016,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000016,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,4,-4}));

   // Creating ~chi_40
   new Particle("~chi_40", 1000035, 0, "Sparticle", 100, 0, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000035));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,22}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,22}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,11,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,13,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,15,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,12,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,14,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,16,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,4,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,22}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,11,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,13,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,15,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,12,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,14,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,16,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,4,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,22}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,11,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,13,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,15,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,12,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,14,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,16,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,4,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000001,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000001,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000002,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000002,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000003,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000003,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000004,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000004,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000005,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000005,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000011,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000011,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000012,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000012,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000013,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000013,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000014,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000014,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000014,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000015,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000015,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000016,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000016,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,4,-4}));

   // Creating ~chi_2+
   new Particle("~chi_2+", 1000037, 1, "Sparticle", 100, 1, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1000037));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,11,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,13,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,15,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,12,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,14,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,16,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,1,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,3,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,5,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,2,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,4,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,-13,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,-15,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,-3,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000002,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000001,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000004,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000003,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000005,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000012,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000011,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000014,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000013,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000016,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000015,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,-1,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,-3,4}));

   // Creating ~gravitino
   new Particle("~gravitino", 1000039, 0, "Sparticle", 100, 0, 500, 0, -100, -1, -100, -1, -1);

   // Creating ~d_R
   new Particle("~d_R", 2000001, 1, "Sparticle", 100, -0.333333, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(2000001));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000002,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000002,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,1}));

   // Creating ~u_R
   new Particle("~u_R", 2000002, 1, "Sparticle", 100, 0.666667, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(2000002));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000001,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000001,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,2}));

   // Creating ~s_R
   new Particle("~s_R", 2000003, 1, "Sparticle", 100, -0.333333, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(2000003));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000004,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000004,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,3}));

   // Creating ~c_R
   new Particle("~c_R", 2000004, 1, "Sparticle", 100, 0.666667, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(2000004));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000003,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000003,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,4}));

   // Creating ~b_2
   new Particle("~b_2", 2000005, 1, "Sparticle", 100, -0.333333, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(2000005));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,5}));

   // Creating ~t_2
   new Particle("~t_2", 2000006, 1, "Sparticle", 100, 0.666667, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(2000006));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000005,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000005,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,6}));

   // Creating ~e_R-
   new Particle("~e_R-", 2000011, 1, "Sparticle", 100, -1, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(2000011));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000012,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000012,-37}));

   // Creating ~nu_eR
   new Particle("~nu_eR", 2000012, 1, "Sparticle", 100, 0, 500, 0, -100, -1, -100, -1, -1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
