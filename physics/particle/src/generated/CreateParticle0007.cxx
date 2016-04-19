// This files was autogenerated by geant::Particle::ReadFile

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
void CreateParticle0007() {
   vector<int> daughters;
   Particle *part = nullptr;

   // Creating D_1s-
   new Particle("D_1s-", -10433, 0, "Unknown", 100, -1, 2.5353, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10433));
   daughters.clear();
   daughters.push_back(-423);
   daughters.push_back(-321);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(-413);
   daughters.push_back(-311);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));

   // Creating D*_0s-
   new Particle("D*_0s-", -10431, 0, "Unknown", 100, -1, 2.3178, 0.0046, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10431));
   daughters.clear();
   daughters.push_back(-431);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.8,  daughters));
   daughters.clear();
   daughters.push_back(-431);
   daughters.push_back(-22);
   part->AddDecay(Particle::Decay(0, 0.2,  daughters));

   // Creating D_10_bar
   new Particle("D_10_bar", -10423, 0, "Unknown", 100, 0, 2.4223, 0.02, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10423));
   daughters.clear();
   daughters.push_back(-413);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(-423);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating D*_00_bar
   new Particle("D*_00_bar", -10421, 0, "Unknown", 100, 0, 2.272, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10421));
   daughters.clear();
   daughters.push_back(-411);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(-421);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating D_1-
   new Particle("D_1-", -10413, 0, "Unknown", 100, -1, 2.424, 0.02, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10413));
   daughters.clear();
   daughters.push_back(-423);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(-413);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating D*_0-
   new Particle("D*_0-", -10411, 0, "Unknown", 100, -1, 2.272, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10411));
   daughters.clear();
   daughters.push_back(-421);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(-411);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating eta2(1870)_bar
   new Particle("eta2(1870)_bar", -10335, 0, "Unknown", 100, 0, 1.842, 0.225, 100, 100, 0, 100, 1);

   // Creating k2(1770)-_bar
   new Particle("k2(1770)-_bar", -10325, 0, "Unknown", 100, -1, 1.773, 0.186, 100, 100, 0, 100, 1);

   // Creating K_1-
   new Particle("K_1-", -10323, 0, "Unknown", 100, -1, 1.272, 0.09, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10323));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.313,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-213);
   part->AddDecay(Particle::Decay(0, 0.28,  daughters));
   daughters.clear();
   daughters.push_back(-323);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.157,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(-113);
   part->AddDecay(Particle::Decay(0, 0.14,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(-223);
   part->AddDecay(Particle::Decay(0, 0.11,  daughters));

   // Creating K*_0-
   new Particle("K*_0-", -10321, 0, "Unknown", 100, -1, 1.42, 0.287, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10321));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating k2(1770)0_bar
   new Particle("k2(1770)0_bar", -10315, 0, "Unknown", 100, 0, 1.773, 0.186, 100, 100, 0, 100, 1);

   // Creating K_10_bar
   new Particle("K_10_bar", -10313, 0, "Unknown", 100, 0, 1.272, 0.09, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10313));
   daughters.clear();
   daughters.push_back(-323);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.313,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(213);
   part->AddDecay(Particle::Decay(0, 0.28,  daughters));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.157,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-113);
   part->AddDecay(Particle::Decay(0, 0.14,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-223);
   part->AddDecay(Particle::Decay(0, 0.11,  daughters));

   // Creating K*_00_bar
   new Particle("K*_00_bar", -10311, 0, "Unknown", 100, 0, 1.42, 0.287, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10311));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating eta2(1645)_bar
   new Particle("eta2(1645)_bar", -10225, 0, "Unknown", 100, 0, 1.617, 0.181, 100, 100, 0, 100, 1);

   // Creating pi2(1670)-_bar
   new Particle("pi2(1670)-_bar", -10215, 0, "Unknown", 100, -1, 1.6722, 0.26, 100, 100, 0, 100, 1);

   // Creating b_1-
   new Particle("b_1-", -10213, 0, "Unknown", 100, -1, 1.2295, 0.142, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10213));
   daughters.clear();
   daughters.push_back(-223);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 1,  daughters));

   // Creating a_0-
   new Particle("a_0-", -10211, 0, "Unknown", 100, -1, 0.9835, 0.06, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10211));
   daughters.clear();
   daughters.push_back(-221);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 1,  daughters));

   // Creating pi2(1670)0_bar
   new Particle("pi2(1670)0_bar", -10115, 0, "Unknown", 100, 0, 1.6722, 0.26, 100, 100, 0, 100, 1);

   // Creating Omega*_bbb+
   new Particle("Omega*_bbb+", -5554, 0, "B-Baryon", 100, 1, 15.1106, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5554));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.14,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-13);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-4);
   daughters.push_back(-1);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(16);
   daughters.push_back(-15);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.04,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.015,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-4);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.01,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.005,  daughters));

   // Creating Omega*_bbc0_bar
   new Particle("Omega*_bbc0_bar", -5544, 0, "B-Baryon", 100, 0, 11.7115, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5544));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.14,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-13);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-4);
   daughters.push_back(-1);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(16);
   daughters.push_back(-15);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.04,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.015,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-4);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.01,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.005,  daughters));

   // Creating Omega_bbc0_bar
   new Particle("Omega_bbc0_bar", -5542, 0, "B-Baryon", 100, 0, 11.7077, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5542));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.14,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-13);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.105,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-4);
   daughters.push_back(-1);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(16);
   daughters.push_back(-15);
   daughters.push_back(-4);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.04,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-1);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.015,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-4);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.01,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-3);
   daughters.push_back(-2);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.005,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif