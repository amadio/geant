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
void CreateParticle0000() {

   // Creating Alpha_bar
   new Particle("Alpha_bar", -1000020040, 0, "ion", 100, -2, 3.7284, 1.6916e-33, 100, 100, 0, 100, 1);

   // Creating HE3_bar
   new Particle("HE3_bar", -1000020030, 0, "ion", 100, -2, 2.80941, 0, 100, 100, 0, 100, 1);

   // Creating Triton_bar
   new Particle("Triton_bar", -1000010030, 0, "ion", 100, -1, 2.80941, 1.6916e-33, 100, 100, 0, 100, 1);

   // Creating Deuteron_bar
   new Particle("Deuteron_bar", -1000010020, 0, "ion", 100, -1, 1.87106, 0, 100, 100, 0, 100, 1);

   // Creating N(2250)+_bar
   new Particle("N(2250)+_bar", -100012210, 0, "Unknown", 100, -1, 2.275, 0.5, 100, 100, 0, 100, 1);

   // Creating N(2250)0_bar
   new Particle("N(2250)0_bar", -100012110, 0, "Unknown", 100, 0, 2.275, 0.5, 100, 100, 0, 100, 1);

   // Creating N(2220)+_bar
   new Particle("N(2220)+_bar", -100002210, 0, "Unknown", 100, -1, 2.25, 0.4, 100, 100, 0, 100, 1);

   // Creating N(2220)0_bar
   new Particle("N(2220)0_bar", -100002110, 0, "Unknown", 100, 0, 2.25, 0.4, 100, 100, 0, 100, 1);

   // Creating GenericIon_bar
   new Particle("GenericIon_bar", -50000060, 0, "Unknown", 100, -0.333333, 0.938272, 0, 100, 100, 0, 100, 1);

   // Creating Cherenkov_bar
   new Particle("Cherenkov_bar", -50000050, 0, "Unknown", 100, 0, 0, 0, 100, 100, 0, 100, 1);

   // Creating f2(2010)_bar
   new Particle("f2(2010)_bar", -9060225, 0, "Unknown", 100, 0, 2.01, 0.2, 100, 100, 0, 100, 1);

   // Creating f2(1810)_bar
   new Particle("f2(1810)_bar", -9030225, 0, "Unknown", 100, 0, 1.815, 0.197, 100, 100, 0, 100, 1);

   // Creating f0(1500)_bar
   new Particle("f0(1500)_bar", -9030221, 0, "Unknown", 100, 0, 1.505, 0.109, 100, 100, 0, 100, 1);

   // Creating eta(1405)_bar
   new Particle("eta(1405)_bar", -9020221, 0, "Unknown", 100, 0, 1.4098, 0.0511, 100, 100, 0, 100, 1);

   // Creating f0(980)_bar
   new Particle("f0(980)_bar", -9010221, 0, "Unknown", 100, 0, 0.98, 0.07, 100, 100, 0, 100, 1);

   // Creating f0(600)_bar
   new Particle("f0(600)_bar", -9000221, 0, "Unknown", 100, 0, 0.8, 0.8, 100, 100, 0, 100, 1);

   // Creating a0(980)-_bar
   new Particle("a0(980)-_bar", -9000211, 0, "Unknown", 100, -1, 0.98, 0.06, 100, 100, 0, 100, 1);

   // Creating a0(980)0_bar
   new Particle("a0(980)0_bar", -9000111, 0, "Unknown", 100, 0, 0.98, 0.075, 100, 100, 0, 100, 1);

   // Creating nu*_e0_bar
   new Particle("nu*_e0_bar", -4000012, 0, "Excited", 100, 0, 400, 0.41917, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-4000012));
   part->AddDecay(Particle::Decay(0, 0.610139,  vector<int>{-24,-11}));
   part->AddDecay(Particle::Decay(0, 0.389861,  vector<int>{-23,-12}));

   // Creating e*+
   new Particle("e*+", -4000011, 0, "Excited", 100, 1, 400, 0.42901, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4000011));
   part->AddDecay(Particle::Decay(0, 0.596149,  vector<int>{24,-12}));
   part->AddDecay(Particle::Decay(0, 0.294414,  vector<int>{-22,-11}));
   part->AddDecay(Particle::Decay(0, 0.109437,  vector<int>{-23,-11}));

   // Creating u*_bar
   new Particle("u*_bar", -4000002, 0, "Excited", 100, -0.666667, 400, 2.65499, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4000002));
   part->AddDecay(Particle::Decay(0, 0.853166,  vector<int>{-21,-2}));
   part->AddDecay(Particle::Decay(0, 0.0963291,  vector<int>{-24,-1}));
   part->AddDecay(Particle::Decay(0, 0.029361,  vector<int>{-23,-2}));
   part->AddDecay(Particle::Decay(0, 0.021144,  vector<int>{-22,-2}));

   // Creating d*_bar
   new Particle("d*_bar", -4000001, 0, "Excited", 100, 0.333333, 400, 2.65171, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4000001));
   part->AddDecay(Particle::Decay(53, 0.85422,  vector<int>{-21,-1}));
   part->AddDecay(Particle::Decay(0, 0.096449,  vector<int>{24,-2}));
   part->AddDecay(Particle::Decay(0, 0.044039,  vector<int>{-23,-1}));
   part->AddDecay(Particle::Decay(0, 0.005292,  vector<int>{-22,-1}));

   // Creating ~nu_tauR_bar
   new Particle("~nu_tauR_bar", -2000016, 0, "Sparticle", 100, 0, 500, 0, 100, 100, 1, 100, 1);

   // Creating ~tau_2+
   new Particle("~tau_2+", -2000015, 0, "Sparticle", 100, 1, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2000015));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,-23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,-25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,-35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,-36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000016,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000016,37}));

   // Creating ~nu_muR_bar
   new Particle("~nu_muR_bar", -2000014, 0, "Sparticle", 100, 0, 500, 0, 100, 100, 1, 100, 1);

   // Creating ~mu_R+
   new Particle("~mu_R+", -2000013, 0, "Sparticle", 100, 1, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2000013));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,-23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,-25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,-35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,-36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000014,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000014,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000014,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000014,37}));

   // Creating ~nu_eR_bar
   new Particle("~nu_eR_bar", -2000012, 0, "Sparticle", 100, 0, 500, 0, 100, 100, 1, 100, 1);

   // Creating ~e_R+
   new Particle("~e_R+", -2000011, 0, "Sparticle", 100, 1, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2000011));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,-23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,-25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,-35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,-36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000012,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000012,37}));

   // Creating ~t_2_bar
   new Particle("~t_2_bar", -2000006, 0, "Sparticle", 100, -0.666667, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2000006));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,-23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,-25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,-35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,-36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000005,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000005,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-6}));

   // Creating ~b_2_bar
   new Particle("~b_2_bar", -2000005, 0, "Sparticle", 100, 0.333333, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2000005));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,-23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,-25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,-35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,-36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-5}));
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
