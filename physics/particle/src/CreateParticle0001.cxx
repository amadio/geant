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
void CreateParticle0001() {

   // Creating ~c_R_bar
   new Particle("~c_R_bar", -2000004, 0, "Sparticle", 100, -0.666667, 500, 1, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-2000004));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,-23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,-25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,-35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,-36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000003,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000003,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-4}));

   // Creating ~s_R_bar
   new Particle("~s_R_bar", -2000003, 0, "Sparticle", 100, 0.333333, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2000003));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,-23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,-25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,-35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,-36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000004,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000004,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-3}));

   // Creating ~u_R_bar
   new Particle("~u_R_bar", -2000002, 0, "Sparticle", 100, -0.666667, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2000002));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,-23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,-25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,-35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,-36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000001,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000001,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-2}));

   // Creating ~d_R_bar
   new Particle("~d_R_bar", -2000001, 0, "Sparticle", 100, 0.333333, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2000001));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,-23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,-25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,-35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,-36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000002,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000002,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-1}));

   // Creating ~chi_2-
   new Particle("~chi_2-", -1000037, 0, "Sparticle", 100, -1, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000037));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-11,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-13,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-15,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-12,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-14,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-16,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-1,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-3,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-5,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-2,2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-4,4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000002,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000001,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000004,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000003,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000005,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000012,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000011,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000014,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000014,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000013,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000016,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000015,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,3,-4}));

   // Creating ~chi_1-
   new Particle("~chi_1-", -1000024, 0, "Sparticle", 100, -1, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000024));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,13,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,15,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,3,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000002,1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000001,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000004,3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000003,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000005,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000012,11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000011,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000014,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000014,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000013,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000016,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000015,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,1,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,3,-4}));

   // Creating ~nu_tauL_bar
   new Particle("~nu_tauL_bar", -1000016, 0, "Sparticle", 100, 0, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000016));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000015,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000015,-37}));

   // Creating ~tau_1+
   new Particle("~tau_1+", -1000015, 0, "Sparticle", 100, 1, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000015));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000016,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000016,37}));

   // Creating ~nu_muL_bar
   new Particle("~nu_muL_bar", -1000014, 0, "Sparticle", 100, 0, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000014));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000013,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000013,-37}));

   // Creating ~mu_L+
   new Particle("~mu_L+", -1000013, 0, "Sparticle", 100, 1, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000013));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000014,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000014,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000014,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000014,37}));

   // Creating ~nu_eL_bar
   new Particle("~nu_eL_bar", -1000012, 0, "Sparticle", 100, 0, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000012));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000011,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000011,-37}));

   // Creating ~e_L+
   new Particle("~e_L+", -1000011, 0, "Sparticle", 100, 1, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000011));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000012,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000012,37}));

   // Creating ~t_1_bar
   new Particle("~t_1_bar", -1000006, 0, "Sparticle", 100, -0.666667, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000006));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000005,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000005,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,15,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-16,-5}));

   // Creating ~b_1_bar
   new Particle("~b_1_bar", -1000005, 0, "Sparticle", 100, 0.333333, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000005));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-5}));

   // Creating ~c_L_bar
   new Particle("~c_L_bar", -1000004, 0, "Sparticle", 100, -0.666667, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000004));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000003,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000003,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-4}));

   // Creating ~s_L_bar
   new Particle("~s_L_bar", -1000003, 0, "Sparticle", 100, 0.333333, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000003));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000004,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000004,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-3}));

   // Creating ~u_L_bar
   new Particle("~u_L_bar", -1000002, 0, "Sparticle", 100, -0.666667, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000002));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000001,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000001,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-2}));

   // Creating ~d_L_bar
   new Particle("~d_L_bar", -1000001, 0, "Sparticle", 100, 0.333333, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000001));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000002,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000002,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-1}));

   // Creating phi(1680)_bar
   new Particle("phi(1680)_bar", -100333, 0, "Unknown", 100, 0, 1.68, 0.15, 100, 100, 0, 100, 1);

   // Creating eta(1475)_bar
   new Particle("eta(1475)_bar", -100331, 0, "Unknown", 100, 0, 1.476, 0.085, 100, 100, 0, 100, 1);

   // Creating k2_star(1980)-_bar
   new Particle("k2_star(1980)-_bar", -100325, 0, "Unknown", 100, -1, 1.973, 0.373, 100, 100, 0, 100, 1);

   // Creating k_star(1410)-_bar
   new Particle("k_star(1410)-_bar", -100323, 0, "Unknown", 100, -1, 1.414, 0.232, 100, 100, 0, 100, 1);

   // Creating k(1460)-_bar
   new Particle("k(1460)-_bar", -100321, 0, "Unknown", 100, -1, 1.46, 0.26, 100, 100, 0, 100, 1);

   // Creating k2_star(1980)0_bar
   new Particle("k2_star(1980)0_bar", -100315, 0, "Unknown", 100, 0, 1.973, 0.373, 100, 100, 0, 100, 1);

   // Creating k_star(1410)0_bar
   new Particle("k_star(1410)0_bar", -100313, 0, "Unknown", 100, 0, 1.414, 0.232, 100, 100, 0, 100, 1);

   // Creating k(1460)0_bar
   new Particle("k(1460)0_bar", -100311, 0, "Unknown", 100, 0, 1.46, 0.26, 100, 100, 0, 100, 1);

   // Creating omega(1420)_bar
   new Particle("omega(1420)_bar", -100223, 0, "Unknown", 100, 0, 1.425, 0.215, 100, 100, 0, 100, 1);

   // Creating eta(1295)_bar
   new Particle("eta(1295)_bar", -100221, 0, "Unknown", 100, 0, 1.294, 0.055, 100, 100, 0, 100, 1);

   // Creating rho(1450)-_bar
   new Particle("rho(1450)-_bar", -100213, 0, "Unknown", 100, -1, 1.465, 0.4, 100, 100, 0, 100, 1);

   // Creating pi(1300)-_bar
   new Particle("pi(1300)-_bar", -100211, 0, "Unknown", 100, -1, 1.3, 0.4, 100, 100, 0, 100, 1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
