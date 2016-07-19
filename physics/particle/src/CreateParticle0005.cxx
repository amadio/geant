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
void CreateParticle0005() {

   // Creating delta(1930)-_bar
   new Particle("delta(1930)-_bar", -11116, 0, "Unknown", 100, 1, 1.96, 0.36, 100, 100, 0, 100, 1);

   // Creating delta(1700)-_bar
   new Particle("delta(1700)-_bar", -11114, 0, "Unknown", 100, 1, 1.7, 0.3, 100, 100, 0, 100, 1);

   // Creating delta(1900)-_bar
   new Particle("delta(1900)-_bar", -11112, 0, "Unknown", 100, 1, 1.9, 0.2, 100, 100, 0, 100, 1);

   // Creating B_1c-
   new Particle("B_1c-", -10543, 0, "Unknown", 100, -1, 7.3, 0.05, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-10543));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-513,-411}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-523,-421}));

   // Creating B*_0c-
   new Particle("B*_0c-", -10541, 0, "Unknown", 100, -1, 7.25, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10541));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-511,-411}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-521,-421}));

   // Creating B_1s0_bar
   new Particle("B_1s0_bar", -10533, 0, "Unknown", 100, 0, 5.97, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10533));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-523,321}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-513,311}));

   // Creating B*_0s0_bar
   new Particle("B*_0s0_bar", -10531, 0, "Unknown", 100, 0, 5.92, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10531));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-521,321}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-511,311}));

   // Creating B_1-
   new Particle("B_1-", -10523, 0, "Unknown", 100, -1, 5.73, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10523));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-513,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-523,-111}));

   // Creating B*_0-
   new Particle("B*_0-", -10521, 0, "Unknown", 100, -1, 5.68, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10521));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-511,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-521,-111}));

   // Creating B_10_bar
   new Particle("B_10_bar", -10513, 0, "Unknown", 100, 0, 5.73, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10513));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-523,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-513,-111}));

   // Creating B*_00_bar
   new Particle("B*_00_bar", -10511, 0, "Unknown", 100, 0, 5.68, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10511));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-521,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-511,-111}));

   // Creating D_1s-
   new Particle("D_1s-", -10433, 0, "Unknown", 100, -1, 2.5353, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10433));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-423,-321}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-413,-311}));

   // Creating D*_0s-
   new Particle("D*_0s-", -10431, 0, "Unknown", 100, -1, 2.3178, 0.0046, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10431));
   part->AddDecay(Particle::Decay(0, 0.8,  vector<int>{-431,-111}));
   part->AddDecay(Particle::Decay(0, 0.2,  vector<int>{-431,-22}));

   // Creating D_10_bar
   new Particle("D_10_bar", -10423, 0, "Unknown", 100, 0, 2.4223, 0.02, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10423));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-413,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-423,-111}));

   // Creating D*_00_bar
   new Particle("D*_00_bar", -10421, 0, "Unknown", 100, 0, 2.272, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10421));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-411,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-421,-111}));

   // Creating D_1-
   new Particle("D_1-", -10413, 0, "Unknown", 100, -1, 2.424, 0.02, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10413));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-423,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-413,-111}));

   // Creating D*_0-
   new Particle("D*_0-", -10411, 0, "Unknown", 100, -1, 2.272, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10411));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-421,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-411,-111}));

   // Creating eta2(1870)_bar
   new Particle("eta2(1870)_bar", -10335, 0, "Unknown", 100, 0, 1.842, 0.225, 100, 100, 0, 100, 1);

   // Creating k2(1770)-_bar
   new Particle("k2(1770)-_bar", -10325, 0, "Unknown", 100, -1, 1.773, 0.186, 100, 100, 0, 100, 1);

   // Creating K_1-
   new Particle("K_1-", -10323, 0, "Unknown", 100, -1, 1.272, 0.09, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10323));
   part->AddDecay(Particle::Decay(0, 0.313,  vector<int>{-313,-211}));
   part->AddDecay(Particle::Decay(0, 0.28,  vector<int>{-311,-213}));
   part->AddDecay(Particle::Decay(0, 0.157,  vector<int>{-323,-111}));
   part->AddDecay(Particle::Decay(0, 0.14,  vector<int>{-321,-113}));
   part->AddDecay(Particle::Decay(0, 0.11,  vector<int>{-321,-223}));

   // Creating K*_0-
   new Particle("K*_0-", -10321, 0, "Unknown", 100, -1, 1.42, 0.287, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10321));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-311,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-321,-111}));

   // Creating k2(1770)0_bar
   new Particle("k2(1770)0_bar", -10315, 0, "Unknown", 100, 0, 1.773, 0.186, 100, 100, 0, 100, 1);

   // Creating K_10_bar
   new Particle("K_10_bar", -10313, 0, "Unknown", 100, 0, 1.272, 0.09, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10313));
   part->AddDecay(Particle::Decay(0, 0.313,  vector<int>{-323,211}));
   part->AddDecay(Particle::Decay(0, 0.28,  vector<int>{-321,213}));
   part->AddDecay(Particle::Decay(0, 0.157,  vector<int>{-313,-111}));
   part->AddDecay(Particle::Decay(0, 0.14,  vector<int>{-311,-113}));
   part->AddDecay(Particle::Decay(0, 0.11,  vector<int>{-311,-223}));

   // Creating K*_00_bar
   new Particle("K*_00_bar", -10311, 0, "Unknown", 100, 0, 1.42, 0.287, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10311));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-321,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-311,-111}));

   // Creating eta2(1645)_bar
   new Particle("eta2(1645)_bar", -10225, 0, "Unknown", 100, 0, 1.617, 0.181, 100, 100, 0, 100, 1);

   // Creating pi2(1670)-_bar
   new Particle("pi2(1670)-_bar", -10215, 0, "Unknown", 100, -1, 1.6722, 0.26, 100, 100, 0, 100, 1);

   // Creating b_1-
   new Particle("b_1-", -10213, 0, "Unknown", 100, -1, 1.2295, 0.142, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10213));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-223,-211}));

   // Creating a_0-
   new Particle("a_0-", -10211, 0, "Unknown", 100, -1, 0.9835, 0.06, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10211));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-221,-211}));

   // Creating pi2(1670)0_bar
   new Particle("pi2(1670)0_bar", -10115, 0, "Unknown", 100, 0, 1.6722, 0.26, 100, 100, 0, 100, 1);

   // Creating Omega*_bbb+
   new Particle("Omega*_bbb+", -5554, 0, "B-Baryon", 100, 1, 15.1106, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5554));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
